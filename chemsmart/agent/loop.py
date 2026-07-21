from __future__ import annotations

import json
from copy import deepcopy
from dataclasses import dataclass
from typing import Any, Callable

from chemsmart.agent.handles import (
    HandleStore,
    is_handle_id,
    json_safe,
    store_result_handle,
)
from chemsmart.agent.permissions import (
    ApprovalDecision,
    PermissionMode,
    PermissionPolicy,
    ResolvedDecision,
    ResolvedPermission,
    RuntimePermissionMode,
    resolve_xtb_run_local,
)
from chemsmart.agent.provider_adapter import (
    ToolOutcome,
    ToolRequest,
)
from chemsmart.agent.registry import ToolRegistry
from chemsmart.agent.services.tool_loop_runner import ToolLoopRunner

ASK_USER_TOOL_NAME = "ask_user"
ASK_USER_TOOL_DEF: dict[str, Any] = {
    "type": "function",
    "function": {
        "name": ASK_USER_TOOL_NAME,
        "description": (
            "Ask the user a clarifying question when truly ambiguous "
            "(no concrete target / multiple candidates / unclear intent)."
        ),
        "parameters": {
            "type": "object",
            "properties": {
                "question": {
                    "type": "string",
                    "description": "The clarifying question to ask the user.",
                },
                "options": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": (
                        "Optional candidate answers to help the user choose."
                    ),
                },
            },
            "required": ["question"],
            "additionalProperties": False,
        },
    },
}


@dataclass(frozen=True)
class ToolLoopBudgets:
    max_model_steps_per_turn: int = 12
    max_total_tool_calls_per_turn: int = 32
    max_consecutive_tool_errors: int = 4
    max_same_signature_retries: int = 2
    max_provider_errors_per_turn: int = 2
    log_provider_turn_raw: bool = False


class ToolLoop:
    def __init__(
        self,
        provider: Any,
        registry: ToolRegistry,
        handle_store: HandleStore | None,
        decision_log: Any,
        budgets: ToolLoopBudgets | None = None,
        policy: PermissionPolicy | None = None,
        approver: Callable[[ToolRequest], ApprovalDecision] | None = None,
        lifecycle: Any | None = None,
    ) -> None:
        self.provider = provider
        self.registry = registry
        self.handle_store = handle_store
        self.decision_log = decision_log
        self.budgets = budgets or ToolLoopBudgets()
        self.policy = policy or PermissionPolicy(mode=PermissionMode.DRIVING)
        self.approver = approver
        self.lifecycle = lifecycle

    def run_turn(
        self,
        messages: list[dict[str, Any]],
        tool_defs: list[dict[str, Any]] | None = None,
        mode: RuntimePermissionMode | None = None,
        allowed_tool_names: set[str] | None = None,
    ) -> dict[str, Any]:
        return ToolLoopRunner(self).run(
            messages, tool_defs, mode, allowed_tool_names
        )

    def _run_one_request(
        self,
        step: int,
        request: ToolRequest,
    ) -> tuple[ToolOutcome, bool]:
        resolved = self.policy.resolve(request)
        if (
            request.name == "run_local"
            and resolved.reason == "always_require_approval"
        ):
            override = self._xtb_run_policy_override(request)
            if override is not None:
                resolved = override
        if resolved.decision == ResolvedDecision.AUTO_DENY:
            self._runtime_permission(
                request,
                decision="denied",
                reason=resolved.reason,
            )
            return (
                self._deny_request(
                    step,
                    request,
                    reason=resolved.reason,
                ),
                False,
            )

        if resolved.decision == ResolvedDecision.NEEDS_USER:
            self._runtime_permission(
                request,
                decision="needs_user",
                reason=resolved.reason,
            )
            if self.approver is None:
                self._runtime_permission(
                    request,
                    decision="denied",
                    reason="no_approver",
                )
                return (
                    self._deny_request(
                        step,
                        request,
                        reason="no_approver",
                    ),
                    False,
                )

            approval = self.approver(request)
            if approval == ApprovalDecision.DENY:
                self._runtime_permission(
                    request,
                    decision="denied",
                    reason="user_denied",
                )
                return (
                    self._deny_request(
                        step,
                        request,
                        reason="user_denied",
                    ),
                    False,
                )

            self.policy.record(request.name, approval)
            self._runtime_permission(
                request,
                decision="approved",
                reason=(
                    "user_session_approval"
                    if approval == ApprovalDecision.ALLOW_SESSION
                    else "user_once_approval"
                ),
            )
            self._log_approved(
                step,
                request,
                scope=(
                    "session"
                    if approval == ApprovalDecision.ALLOW_SESSION
                    else "once"
                ),
                source="approver",
            )
            return self._execute_request(step, request), True

        self._runtime_permission(
            request,
            decision="approved",
            reason=resolved.reason,
        )
        self._log_approved(
            step,
            request,
            scope="session" if resolved.reason == "session_rule" else "auto",
            source=resolved.reason,
        )
        return self._execute_request(step, request), True

    def _execute_request(
        self,
        step: int,
        request: ToolRequest,
    ) -> ToolOutcome:
        if self.lifecycle is not None:
            try:
                self.lifecycle.before_tool(
                    request_id=request.request_id,
                    tool_name=request.name,
                    arguments=request.arguments,
                )
            except Exception as exc:
                return self._error_outcome(
                    step,
                    request,
                    error_type=exc.__class__.__name__,
                    error_message=str(exc),
                )
        try:
            resolved_args = _resolve_handles(
                request.arguments,
                handle_store=self.handle_store,
            )
        except KeyError as exc:
            return self._error_outcome(
                step,
                request,
                error_type="UnknownHandle",
                error_message=f"Unknown handle {exc.args[0]!r}",
            )

        try:
            result = self.registry.call(request.name, resolved_args)
        except Exception as exc:
            return self._error_outcome(
                step,
                request,
                error_type=exc.__class__.__name__,
                error_message=str(exc),
            )
        if _is_tool_error(result):
            error = result.get("error") or {}
            if not isinstance(error, dict):
                error = {
                    "type": "ToolError",
                    "message": str(error or f"{request.name} failed"),
                }
            return self._error_outcome(
                step,
                request,
                error_type=str(error.get("type") or "ToolError"),
                error_message=str(
                    error.get("message") or f"{request.name} failed"
                ),
                raw_result=result,
            )

        handle_id = self._store_result_handle(request.name, result)
        display_result = _display_result(result, handle_id=handle_id)
        outcome = ToolOutcome(
            request_id=request.request_id,
            provider_call_id=request.provider_call_id,
            name=request.name,
            status="ok",
            result=result,
            display_result=display_result,
            raw_result=result,
            handle_id=handle_id,
        )
        self.decision_log.write(
            "tool_use_result",
            {
                "step": step,
                "provider_call_id": request.provider_call_id,
                "tool": request.name,
                "status": outcome.status,
                "description": self._tool_description(request.name),
                "payload": display_result,
                "handle_id": handle_id,
            },
        )
        if self.lifecycle is not None:
            self.lifecycle.after_tool(
                request_id=request.request_id,
                tool_name=request.name,
                result=result,
            )
        return outcome

    def _xtb_run_policy_override(
        self,
        request: ToolRequest,
    ) -> ResolvedPermission | None:
        """Apply the rule-based xtb_real_runs policy to a run_local request.

        Permission resolution happens before handle resolution, so the loop
        (which owns the handle store) is the earliest point that can tell
        whether the requested job is an xTB job and how large it is. Returns
        None for ask mode, non-xTB jobs, unknown handles, or failed guards —
        all of which fall back to the normal approval flow.
        """

        policy_value = getattr(self.policy, "xtb_real_runs", "ask")
        handle_id = (request.arguments or {}).get("job")
        if not isinstance(handle_id, str) or self.handle_store is None:
            return None
        try:
            job = self.handle_store.get(handle_id)
        except KeyError:
            return None
        molecule = getattr(job, "molecule", None)
        symbols = getattr(molecule, "symbols", None)
        return resolve_xtb_run_local(
            policy_value,
            job_program=getattr(job, "PROGRAM", None),
            atom_count=len(symbols) if symbols is not None else None,
        )

    def _deny_request(
        self,
        step: int,
        request: ToolRequest,
        *,
        reason: str,
    ) -> ToolOutcome:
        self.decision_log.write(
            "tool_use_denied",
            {
                "step": step,
                "provider_call_id": request.provider_call_id,
                "tool": request.name,
                "description": self._tool_description(request.name),
                "mode": self.policy.mode.value,
                "reason": reason,
            },
        )
        outcome = ToolOutcome(
            request_id=request.request_id,
            provider_call_id=request.provider_call_id,
            name=request.name,
            status="denied",
            error_type="PermissionDenied",
            error_message=reason,
        )
        self.decision_log.write(
            "tool_use_result",
            {
                "step": step,
                "provider_call_id": request.provider_call_id,
                "tool": request.name,
                "status": outcome.status,
                "description": self._tool_description(request.name),
                "reason": reason,
                "payload": {
                    "ok": False,
                    "error": {
                        "type": outcome.error_type,
                        "message": reason,
                        "tool": request.name,
                    },
                },
                "handle_id": None,
            },
        )
        return outcome

    def _log_approved(
        self,
        step: int,
        request: ToolRequest,
        *,
        scope: str,
        source: str,
    ) -> None:
        self.decision_log.write(
            "tool_use_approved",
            {
                "step": step,
                "provider_call_id": request.provider_call_id,
                "tool": request.name,
                "description": self._tool_description(request.name),
                "scope": scope,
                "source": source,
            },
        )

    def _error_outcome(
        self,
        step: int,
        request: ToolRequest,
        *,
        error_type: str,
        error_message: str,
        raw_result: Any = None,
    ) -> ToolOutcome:
        if self.lifecycle is not None:
            self.lifecycle.tool_failed(
                request_id=request.request_id,
                tool_name=request.name,
                error_type=error_type,
                error_message=error_message,
                result=raw_result,
            )
        outcome = ToolOutcome(
            request_id=request.request_id,
            provider_call_id=request.provider_call_id,
            name=request.name,
            status="error",
            error_type=error_type,
            error_message=error_message,
            raw_result=raw_result,
        )
        self.decision_log.write(
            "tool_use_result",
            {
                "step": step,
                "provider_call_id": request.provider_call_id,
                "tool": request.name,
                "status": outcome.status,
                "description": self._tool_description(request.name),
                "payload": {
                    "ok": False,
                    "error": {
                        "type": error_type,
                        "message": error_message,
                        "tool": request.name,
                    },
                },
                "handle_id": None,
            },
        )
        return outcome

    def _runtime_permission(
        self,
        request: ToolRequest,
        *,
        decision: str,
        reason: str,
    ) -> None:
        if self.lifecycle is None:
            return
        self.lifecycle.permission(
            request_id=request.request_id,
            tool_name=request.name,
            decision=decision,
            reason=reason,
        )

    def _skipped_outcome(
        self,
        step: int,
        request: ToolRequest,
        *,
        reason: str,
    ) -> ToolOutcome:
        self.decision_log.write(
            "tool_use_skipped",
            {
                "step": step,
                "provider_call_id": request.provider_call_id,
                "tool": request.name,
                "reason": reason,
            },
        )
        outcome = ToolOutcome(
            request_id=request.request_id,
            provider_call_id=request.provider_call_id,
            name=request.name,
            status="skipped",
            error_type="ToolSkipped",
            error_message=reason,
        )
        self.decision_log.write(
            "tool_use_result",
            {
                "step": step,
                "provider_call_id": request.provider_call_id,
                "tool": request.name,
                "status": outcome.status,
                "description": self._tool_description(request.name),
                "reason": reason,
                "payload": {
                    "ok": False,
                    "error": {
                        "type": "ToolSkipped",
                        "message": reason,
                        "tool": request.name,
                    },
                },
                "handle_id": None,
            },
        )
        return outcome

    def _normalized_args(self, request: ToolRequest) -> dict[str, Any]:
        if hasattr(self.registry, "normalize_args"):
            return self.registry.normalize_args(
                request.name, request.arguments
            )
        return dict(request.arguments)

    def _tool_description(self, tool_name: str) -> str:
        if hasattr(self.registry, "describe_tool"):
            return self.registry.describe_tool(tool_name)
        return tool_name

    def _tool_defs_for_mode(
        self,
        provider_name: str,
        mode: RuntimePermissionMode,
    ) -> list[dict[str, Any]]:
        if hasattr(self.registry, "assemble_tool_pool"):
            tool_pool = self.registry.assemble_tool_pool(mode)
            if hasattr(self.registry, "tool_defs_for_provider"):
                try:
                    return with_virtual_tool_defs(
                        provider_name,
                        self.registry.tool_defs_for_provider(
                            provider_name, tool_pool
                        ),
                    )
                except TypeError:
                    pass
            if provider_name == "anthropic":
                return with_virtual_tool_defs(
                    provider_name,
                    [tool.anthropic_tool_def() for tool in tool_pool],
                )
            return with_virtual_tool_defs(
                provider_name,
                [tool.openai_tool_def() for tool in tool_pool],
            )
        return self._tool_defs_for_provider(provider_name)

    def _tool_defs_for_provider(
        self,
        provider_name: str,
    ) -> list[dict[str, Any]]:
        return registry_tool_defs_for_provider(self.registry, provider_name)

    def _filter_tool_defs(
        self,
        provider_name: str,
        tool_defs: list[dict[str, Any]],
        allowed_tool_names: set[str],
    ) -> list[dict[str, Any]]:
        return _filter_tool_defs(provider_name, tool_defs, allowed_tool_names)

    def _store_result_handle(
        self,
        tool_name: str,
        result: Any,
    ) -> str | None:
        return store_result_handle(
            self.handle_store,
            tool_name,
            result,
            summary=json_safe(result),
        )


def _assistant_message(
    provider_name: str,
    response_dict: dict[str, Any],
) -> dict[str, Any]:
    if provider_name == "anthropic":
        return {
            "role": "assistant",
            "content": deepcopy(response_dict.get("content") or []),
        }

    choice = (response_dict.get("choices") or [{}])[0] or {}
    message = deepcopy(choice.get("message") or {})
    message.setdefault("role", "assistant")
    return message


def _canonical_args_json(request: ToolRequest) -> str:
    return json.dumps(request.arguments, sort_keys=True)


def _display_result(result: Any, *, handle_id: str | None) -> Any:
    if handle_id is None:
        return json_safe(result)
    payload: dict[str, Any] = {"handle_id": handle_id}
    summary = json_safe(result)
    if isinstance(summary, dict):
        payload["summary"] = summary
    return payload


def with_virtual_tool_defs(
    provider_name: str,
    tool_defs: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    defs = [deepcopy(tool_def) for tool_def in tool_defs]
    if any(
        _tool_def_name(tool_def) == ASK_USER_TOOL_NAME for tool_def in defs
    ):
        return defs
    defs.append(_ask_user_tool_def_for_provider(provider_name))
    return defs


def registry_tool_defs_for_provider(
    registry: ToolRegistry,
    provider_name: str,
) -> list[dict[str, Any]]:
    """Return provider-native registry definitions plus virtual tools."""

    if hasattr(registry, "tool_defs_for_provider"):
        tool_defs = registry.tool_defs_for_provider(provider_name)
    else:
        tool_defs = registry.openai_tool_defs()
    return with_virtual_tool_defs(provider_name, tool_defs)


def _filter_tool_defs(
    provider_name: str,
    tool_defs: list[dict[str, Any]],
    allowed_tool_names: set[str],
) -> list[dict[str, Any]]:
    filtered = [
        tool_def
        for tool_def in tool_defs
        if _tool_def_name(tool_def) in allowed_tool_names
    ]
    return with_virtual_tool_defs(provider_name, filtered)


def _ask_user_tool_def_for_provider(provider_name: str) -> dict[str, Any]:
    if provider_name != "anthropic":
        return deepcopy(ASK_USER_TOOL_DEF)

    function: dict[str, Any] = ASK_USER_TOOL_DEF.get("function") or {}
    return {
        "name": function.get("name"),
        "description": function.get("description"),
        "input_schema": deepcopy(function.get("parameters") or {}),
    }


def _tool_def_name(tool_def: dict[str, Any]) -> str | None:
    if not isinstance(tool_def, dict):
        return None
    function = tool_def.get("function")
    if isinstance(function, dict):
        name = function.get("name")
        return name if isinstance(name, str) else None
    name = tool_def.get("name")
    return name if isinstance(name, str) else None


def _resolve_handles(
    value: Any,
    *,
    handle_store: HandleStore | None,
) -> Any:
    if (
        isinstance(value, str)
        and handle_store is not None
        and is_handle_id(value)
    ):
        return handle_store.get(value)
    if isinstance(value, list):
        return [
            _resolve_handles(item, handle_store=handle_store) for item in value
        ]
    if isinstance(value, dict):
        return {
            key: _resolve_handles(item, handle_store=handle_store)
            for key, item in value.items()
        }
    return value


def _is_tool_error(result: Any) -> bool:
    return (
        isinstance(result, dict)
        and result.get("ok") is False
        and "error" in result
    )
