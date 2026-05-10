from __future__ import annotations

import json
from collections import Counter
from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

from chemsmart.agent.handles import HandleStore, is_handle_id
from chemsmart.agent.permissions import (
    ApprovalDecision,
    PermissionMode,
    PermissionPolicy,
    ResolvedDecision,
)
from chemsmart.agent.provider_adapter import (
    ToolOutcome,
    ToolRequest,
    build_tool_result_messages,
    normalize_response,
)
from chemsmart.agent.providers import (
    DEFAULT_TIMEOUT_S,
    extract_response_usage,
)
from chemsmart.agent.registry import ToolRegistry


@dataclass(frozen=True)
class ToolLoopBudgets:
    max_model_steps_per_turn: int = 12
    max_total_tool_calls_per_turn: int = 32
    max_consecutive_tool_errors: int = 4
    max_same_signature_retries: int = 2
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
    ) -> None:
        self.provider = provider
        self.registry = registry
        self.handle_store = handle_store
        self.decision_log = decision_log
        self.budgets = budgets or ToolLoopBudgets()
        self.policy = policy or PermissionPolicy(mode=PermissionMode.DRIVING)
        self.approver = approver

    def run_turn(
        self,
        messages: list[dict[str, Any]],
        tool_defs: list[dict[str, Any]],
    ) -> dict[str, Any]:
        provider_name = getattr(self.provider, "name", None) or "openai"
        history = [deepcopy(message) for message in messages]
        assistant_text = ""
        stop_reason = None
        limit_reason = None
        model_steps = 0
        tool_calls = 0
        consecutive_tool_errors = 0
        total_input_tokens = 0
        total_output_tokens = 0
        signature_counts: Counter[tuple[str, str]] = Counter()
        tool_requests: list[ToolRequest] = []
        tool_outcomes: list[ToolOutcome] = []
        approvals_count = 0
        denials_count = 0

        while True:
            if model_steps >= self.budgets.max_model_steps_per_turn:
                limit_reason = "max_model_steps"
                break

            model_steps += 1
            response = self.provider.chat(
                history,
                tools=tool_defs,
                timeout_s=DEFAULT_TIMEOUT_S,
            )
            response_dict = _response_payload(response)
            if self.budgets.log_provider_turn_raw:
                self.decision_log.write(
                    "provider_turn_raw",
                    {
                        "step": model_steps,
                        "response_dict": response_dict,
                    },
                )

            assistant_message = _assistant_message(
                provider_name, response_dict
            )
            assistant_text, requests, stop_reason = normalize_response(
                provider_name,
                response_dict,
            )
            usage = extract_response_usage(response_dict)
            total_input_tokens += int(usage["input_tokens"] or 0)
            total_output_tokens += int(usage["output_tokens"] or 0)
            self.decision_log.write(
                "assistant_turn",
                {
                    "step": model_steps,
                    "assistant_text": assistant_text,
                    "stop_reason": stop_reason,
                    "usage": usage,
                },
            )

            if not requests:
                history.append(assistant_message)
                break

            history.append(assistant_message)
            step_outcomes: list[ToolOutcome] = []
            should_stop = False

            stop_index: int | None = None
            for request_index, request in enumerate(requests):
                tool_requests.append(request)
                self.decision_log.write(
                    "tool_use_request",
                    {
                        "step": model_steps,
                        "provider_call_id": request.provider_call_id,
                        "tool": request.name,
                        "args": request.arguments,
                        "raw": request.raw,
                    },
                )

                signature = (request.name, _canonical_args_json(request))
                signature_counts[signature] += 1
                if (
                    signature_counts[signature]
                    > self.budgets.max_same_signature_retries
                ):
                    limit_reason = "repeat_signature"
                    outcome = self._skipped_outcome(
                        model_steps,
                        request,
                        reason="repeat_signature",
                    )
                    step_outcomes.append(outcome)
                    tool_outcomes.append(outcome)
                    should_stop = True
                    stop_index = request_index
                    break

                if tool_calls >= self.budgets.max_total_tool_calls_per_turn:
                    limit_reason = "max_tool_calls"
                    outcome = self._skipped_outcome(
                        model_steps,
                        request,
                        reason="max_tool_calls",
                    )
                    step_outcomes.append(outcome)
                    tool_outcomes.append(outcome)
                    should_stop = True
                    stop_index = request_index
                    break

                outcome, approved = self._run_one_request(
                    model_steps,
                    request,
                )
                step_outcomes.append(outcome)
                tool_outcomes.append(outcome)
                if approved:
                    approvals_count += 1
                    tool_calls += 1
                elif outcome.status == "denied":
                    denials_count += 1

                if outcome.status == "error":
                    consecutive_tool_errors += 1
                else:
                    consecutive_tool_errors = 0

                if (
                    consecutive_tool_errors
                    >= self.budgets.max_consecutive_tool_errors
                ):
                    limit_reason = "max_consecutive_errors"
                    should_stop = True
                    stop_index = request_index
                    break

            if should_stop and stop_index is not None:
                for request in requests[stop_index + 1 :]:
                    tool_requests.append(request)
                    self.decision_log.write(
                        "tool_use_request",
                        {
                            "step": model_steps,
                            "provider_call_id": request.provider_call_id,
                            "tool": request.name,
                            "args": request.arguments,
                            "raw": request.raw,
                        },
                    )
                    outcome = self._skipped_outcome(
                        model_steps,
                        request,
                        reason=limit_reason or "loop_stopped",
                    )
                    step_outcomes.append(outcome)
                    tool_outcomes.append(outcome)

            if step_outcomes:
                history.extend(
                    build_tool_result_messages(provider_name, step_outcomes)
                )

            if should_stop:
                break

        if limit_reason is not None:
            self.decision_log.write(
                "loop_limit_exceeded",
                {
                    "limit_reason": limit_reason,
                    "model_steps": model_steps,
                    "tool_calls": tool_calls,
                },
            )

        return {
            "assistant_text": assistant_text,
            "tool_requests": tool_requests,
            "tool_outcomes": tool_outcomes,
            "stop_reason": stop_reason,
            "model_steps": model_steps,
            "limit_reason": limit_reason,
            "total_input_tokens": total_input_tokens,
            "total_output_tokens": total_output_tokens,
            "messages": history,
            "approvals_count": approvals_count,
            "denials_count": denials_count,
        }

    def _run_one_request(
        self,
        step: int,
        request: ToolRequest,
    ) -> tuple[ToolOutcome, bool]:
        resolved = self.policy.resolve(request)
        if resolved.decision == ResolvedDecision.AUTO_DENY:
            return (
                self._deny_request(
                    step,
                    request,
                    reason=resolved.reason,
                ),
                False,
            )

        if resolved.decision == ResolvedDecision.NEEDS_USER:
            if self.approver is None:
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
                return (
                    self._deny_request(
                        step,
                        request,
                        reason="user_denied",
                    ),
                    False,
                )

            self.policy.record(request.name, approval)
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
                "payload": display_result,
                "handle_id": handle_id,
            },
        )
        return outcome

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

    def _store_result_handle(
        self,
        tool_name: str,
        result: Any,
    ) -> str | None:
        if self.handle_store is None:
            return None
        kind = _result_handle_kind(tool_name, result)
        if kind is None:
            return None
        return self.handle_store.put(
            kind=kind,
            obj=result,
            summary=_json_safe(result),
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
        return _json_safe(result)
    payload = {"handle_id": handle_id}
    summary = _json_safe(result)
    if isinstance(summary, dict):
        payload["summary"] = summary
    return payload


def _response_payload(response: Any) -> dict[str, Any]:
    if isinstance(response, dict):
        return response
    if hasattr(response, "model_dump"):
        payload = response.model_dump()
        if isinstance(payload, dict):
            return payload
    raise TypeError("Provider response must be a dict-like payload")


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


def _result_handle_kind(tool_name: str, result: Any) -> str | None:
    if tool_name == "build_molecule":
        return "mol"
    if tool_name == "build_gaussian_settings":
        return "gset"
    if tool_name == "build_orca_settings":
        return "oset"
    if tool_name == "build_job":
        return "job"
    if tool_name == "dry_run_input" and isinstance(result, dict):
        return "dryrun"
    if tool_name == "validate_runtime" and isinstance(result, dict):
        return "runtime"
    if tool_name == "run_local" and isinstance(result, dict):
        return "runresult"
    if tool_name == "extract_optimized_geometry":
        return "geom"
    if (
        tool_name == "submit_hpc"
        and isinstance(result, dict)
        and "job_id" in result
    ):
        return "submit"
    if tool_name == "recommend_method" and isinstance(result, dict):
        return "recmethod"
    return None


def _json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, bytes):
        return {"type": "bytes", "length": len(value)}
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return {"type": value.__class__.__name__, "repr": repr(value)}
