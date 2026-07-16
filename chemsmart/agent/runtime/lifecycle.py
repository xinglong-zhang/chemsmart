"""Tool lifecycle hooks that emit compact public runtime evidence."""

from __future__ import annotations

import hashlib
import json
from typing import Any, Protocol

from chemsmart.agent.runtime.contracts import RuntimeV2Mode
from chemsmart.agent.runtime.events import EventKind
from chemsmart.agent.runtime.receipts import collect_artifact_refs
from chemsmart.agent.runtime.tool_catalog import ToolSelection


class EventEmitter(Protocol):
    def emit(
        self,
        kind: EventKind,
        payload: dict[str, Any],
        *,
        idempotency_key: str = "",
    ) -> Any: ...


class ToolExposureViolation(RuntimeError):
    pass


class RuntimeLifecycle:
    def __init__(
        self,
        *,
        emitter: EventEmitter,
        selection: ToolSelection,
        mode: RuntimeV2Mode,
    ) -> None:
        self.emitter = emitter
        self.selection = selection
        self.mode = mode

    def before_tool(
        self,
        *,
        request_id: str,
        tool_name: str,
        arguments: dict[str, Any],
    ) -> None:
        if tool_name not in self.selection.direct:
            payload = {
                "rule_id": "runtime.tool.not_exposed",
                "tool": tool_name,
                "phase": self.selection.phase.value,
            }
            self.emitter.emit(EventKind.SHADOW_VIOLATION, payload)
            if self.mode is RuntimeV2Mode.ACTIVE:
                raise ToolExposureViolation(
                    f"tool {tool_name!r} is not exposed in "
                    f"phase {self.selection.phase.value!r}"
                )
        canonical = json.dumps(arguments, sort_keys=True, default=str)
        self.emitter.emit(
            EventKind.TOOL_STARTED,
            {
                "request_id": request_id,
                "tool": tool_name,
                "arg_keys": sorted(arguments),
                "signature_hash": hashlib.sha256(
                    canonical.encode()
                ).hexdigest(),
            },
            idempotency_key=f"tool-start:{request_id}",
        )

    def permission(
        self,
        *,
        request_id: str,
        tool_name: str,
        decision: str,
        reason: str,
    ) -> None:
        self.emitter.emit(
            EventKind.PERMISSION_RESOLVED,
            {
                "request_id": request_id,
                "tool": tool_name,
                "decision": decision,
                "reason": reason,
            },
            idempotency_key=f"permission:{request_id}:{decision}",
        )

    def after_tool(
        self,
        *,
        request_id: str,
        tool_name: str,
        result: Any,
    ) -> None:
        payload = _success_payload(request_id, tool_name, result)
        self.emitter.emit(
            EventKind.TOOL_SUCCEEDED,
            payload,
            idempotency_key=f"tool-result:{request_id}",
        )
        if isinstance(result, dict):
            _emit_state_delta(self.emitter, tool_name, result)
        for receipt in collect_artifact_refs(result, producer_tool=tool_name):
            self.emitter.emit(
                EventKind.ARTIFACT_RECORDED,
                receipt.model_dump(mode="json"),
                idempotency_key=f"artifact:{receipt.sha256}",
            )

    def tool_failed(
        self,
        *,
        request_id: str,
        tool_name: str,
        error_type: str,
        error_message: str,
        result: Any = None,
    ) -> None:
        self.emitter.emit(
            EventKind.TOOL_FAILED,
            {
                "request_id": request_id,
                "tool": tool_name,
                "error_type": error_type,
                "message": error_message[:500],
                "rule_ids": list(_rule_ids(result)),
            },
            idempotency_key=f"tool-result:{request_id}",
        )


def _success_payload(
    request_id: str,
    tool_name: str,
    result: Any,
) -> dict[str, Any]:
    payload: dict[str, Any] = {
        "request_id": request_id,
        "tool": tool_name,
        "status": "ok",
        "rule_ids": list(_rule_ids(result)),
    }
    if isinstance(result, dict):
        payload["result_keys"] = sorted(str(key) for key in result)
        if isinstance(result.get("state_delta"), dict):
            payload["state_delta"] = result["state_delta"]
    return payload


def _emit_state_delta(
    emitter: EventEmitter,
    tool_name: str,
    result: dict[str, Any],
) -> None:
    state_delta = result.get("state_delta")
    project = (
        state_delta.get("project") if isinstance(state_delta, dict) else None
    )
    if isinstance(project, dict) and project.get("selected"):
        emitter.emit(
            EventKind.PROJECT_SELECTED,
            {
                "name": str(project.get("project") or ""),
                "program": str(project.get("program") or ""),
                "path": str(project.get("path") or ""),
                "sha256": str(project.get("sha256") or ""),
            },
            idempotency_key=(
                f"project:{project.get('sha256') or project.get('path')}"
            ),
        )
    if tool_name in {"synthesize_command", "repair_command"}:
        command = str(result.get("command") or "").strip()
        if command:
            emitter.emit(
                EventKind.COMMAND_SYNTHESIZED,
                {
                    "command": command,
                    "semantic_verdict": _nested_value(
                        result, "semantic", "verdict"
                    ),
                    "intent_verdict": _nested_value(
                        result, "intent", "verdict"
                    ),
                },
                idempotency_key=f"command:{hashlib.sha256(command.encode()).hexdigest()}",
            )
        if result.get("status") == "needs_clarification":
            emitter.emit(
                EventKind.CLARIFICATION_REQUESTED,
                {"slots": list(result.get("missing_info") or [])},
            )


def _rule_ids(value: Any) -> tuple[str, ...]:
    found: list[str] = []

    def visit(node: Any) -> None:
        if isinstance(node, dict):
            for key, child in node.items():
                if key == "rule_id" and isinstance(child, str):
                    found.append(child)
                elif key in {"failed_rule_ids", "rule_ids"} and isinstance(
                    child, (list, tuple)
                ):
                    found.extend(str(item) for item in child)
                elif key in {"semantic", "intent", "validation", "error"}:
                    visit(child)
        elif isinstance(node, (list, tuple)):
            for child in node:
                visit(child)

    visit(value)
    return tuple(dict.fromkeys(found))


def _nested_value(value: dict[str, Any], key: str, child: str) -> Any:
    nested = value.get(key)
    return nested.get(child) if isinstance(nested, dict) else None


__all__ = ["RuntimeLifecycle", "ToolExposureViolation"]
