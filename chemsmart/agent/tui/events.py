"""Decision-log adapters for TUI cells."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from chemsmart.agent.core import (
    CriticVerdict,
    Plan,
    _restore_json_result,
    render_plan,
)
from chemsmart.io.molecules.structure import Molecule


@dataclass(slots=True)
class AgentEvent:
    kind: str


@dataclass(slots=True)
class RequestEvent(AgentEvent):
    request: str


@dataclass(slots=True)
class PlanEvent(AgentEvent):
    plan: Plan
    text: str


@dataclass(slots=True)
class AssistantTurnEvent(AgentEvent):
    text: str
    stop_reason: str | None


@dataclass(slots=True)
class ToolCallEvent(AgentEvent):
    step_index: int
    tool: str
    status: str
    args: dict[str, Any]
    payload: Any | None = None


@dataclass(slots=True)
class ToolPreviewEvent(AgentEvent):
    step_index: int
    tool: str
    status: str
    args: dict[str, Any]
    payload: Any | None = None


@dataclass(slots=True)
class ToolUseEvent(AgentEvent):
    step_index: int
    tool: str
    status: str
    args: dict[str, Any]
    description: str
    provider_call_id: str
    reason: str | None = None
    payload: Any | None = None
    queue_index: int | None = None
    queue_total: int | None = None
    scope: str | None = None


@dataclass(slots=True)
class MethodEvent(AgentEvent):
    recommendation: dict[str, Any]


@dataclass(slots=True)
class DryRunInputEvent(AgentEvent):
    step_index: int
    inputfile: str | None
    content: str


@dataclass(slots=True)
class RuntimeValidationEvent(AgentEvent):
    step_index: int
    validation: dict[str, Any]


@dataclass(slots=True)
class SubmissionPreviewEvent(AgentEvent):
    step_index: int
    preview: dict[str, Any]


@dataclass(slots=True)
class GeometryHandoffEvent(AgentEvent):
    molecule: Molecule
    session_dir: str | None


@dataclass(slots=True)
class CriticVerdictEvent(AgentEvent):
    verdict: CriticVerdict


@dataclass(slots=True)
class ErrorEvent(AgentEvent):
    title: str
    message: str
    details: dict[str, Any]


@dataclass(slots=True)
class SessionSummaryEvent(AgentEvent):
    blocked: bool
    block_reason: str | None
    total_steps_executed: int
    total_steps_planned: int


@dataclass(slots=True)
class IgnoredEvent(AgentEvent):
    payload: dict[str, Any]


_AGENT_EVENT_TYPES = (
    RequestEvent
    | PlanEvent
    | AssistantTurnEvent
    | ToolCallEvent
    | ToolPreviewEvent
    | ToolUseEvent
    | MethodEvent
    | DryRunInputEvent
    | RuntimeValidationEvent
    | SubmissionPreviewEvent
    | GeometryHandoffEvent
    | CriticVerdictEvent
    | ErrorEvent
    | SessionSummaryEvent
    | IgnoredEvent
)


def parse_decision_event(
    entry: dict[str, Any],
    *,
    session_dir: Path | None = None,
) -> _AGENT_EVENT_TYPES:
    kind = str(entry.get("kind") or "")
    payload = entry.get("payload") or {}

    if kind == "request":
        return RequestEvent(
            kind=kind, request=str(payload.get("request") or "")
        )

    if kind == "plan":
        plan = Plan.model_validate(payload)
        return PlanEvent(kind=kind, plan=plan, text=render_plan(plan))

    if kind == "assistant_turn":
        return AssistantTurnEvent(
            kind=kind,
            text=str(payload.get("assistant_text") or ""),
            stop_reason=(
                str(payload.get("stop_reason"))
                if payload.get("stop_reason") is not None
                else None
            ),
        )

    if kind == "tool_call":
        return ToolCallEvent(
            kind=kind,
            step_index=int(payload.get("step_index") or 0),
            tool=str(payload.get("tool") or ""),
            status="running",
            args=dict(payload.get("args") or {}),
        )

    if kind == "tool_preview":
        return ToolPreviewEvent(
            kind=kind,
            step_index=int(payload.get("step_index") or 0),
            tool=str(payload.get("tool") or ""),
            status="running",
            args=dict(payload.get("args") or {}),
        )

    if kind == "tool_use_request":
        return ToolUseEvent(
            kind=kind,
            step_index=int(payload.get("step") or 0),
            tool=str(payload.get("tool") or ""),
            status="pending",
            args=dict(
                payload.get("normalized_args") or payload.get("args") or {}
            ),
            description=str(payload.get("description") or ""),
            provider_call_id=str(payload.get("provider_call_id") or ""),
            queue_index=_coerce_int(payload.get("queue_index")),
            queue_total=_coerce_int(payload.get("queue_total")),
        )

    if kind == "tool_use_approved":
        return ToolUseEvent(
            kind=kind,
            step_index=int(payload.get("step") or 0),
            tool=str(payload.get("tool") or ""),
            status="approved",
            args={},
            description=str(payload.get("description") or ""),
            provider_call_id=str(payload.get("provider_call_id") or ""),
            scope=(
                str(payload.get("scope"))
                if payload.get("scope") is not None
                else None
            ),
        )

    if kind == "tool_use_denied":
        return ToolUseEvent(
            kind=kind,
            step_index=int(payload.get("step") or 0),
            tool=str(payload.get("tool") or ""),
            status="denied",
            args={},
            description=str(payload.get("description") or ""),
            provider_call_id=str(payload.get("provider_call_id") or ""),
            reason=(
                str(payload.get("reason"))
                if payload.get("reason") is not None
                else None
            ),
        )

    if kind == "tool_use_result":
        return ToolUseEvent(
            kind=kind,
            step_index=int(payload.get("step") or 0),
            tool=str(payload.get("tool") or ""),
            status=str(payload.get("status") or "ok"),
            args={},
            description=str(payload.get("description") or ""),
            provider_call_id=str(payload.get("provider_call_id") or ""),
            payload=payload.get("payload"),
            reason=(
                str(payload.get("reason"))
                if payload.get("reason") is not None
                else None
            ),
        )

    if kind in {"tool_result", "tool_preview_result"}:
        tool = payload.get("tool")
        tool_payload = payload.get("payload") or {}
        step_index = int(payload.get("step_index") or 0)
        if kind == "tool_preview_result" and tool != "submit_hpc":
            return ToolPreviewEvent(
                kind=kind,
                step_index=step_index,
                tool=str(tool or ""),
                status="done",
                args={},
                payload=tool_payload,
            )
        if tool in {
            "build_molecule",
            "build_gaussian_settings",
            "build_job",
            "run_local",
        }:
            return ToolCallEvent(
                kind=kind,
                step_index=step_index,
                tool=str(tool or ""),
                status="done",
                args={},
                payload=tool_payload,
            )
        if tool == "recommend_method":
            return MethodEvent(kind=kind, recommendation=tool_payload)
        if tool == "dry_run_input":
            return DryRunInputEvent(
                kind=kind,
                step_index=step_index,
                inputfile=tool_payload.get("inputfile"),
                content=str(tool_payload.get("content") or ""),
            )
        if tool == "validate_runtime":
            return RuntimeValidationEvent(
                kind=kind,
                step_index=step_index,
                validation=tool_payload,
            )
        if tool == "submit_hpc":
            return SubmissionPreviewEvent(
                kind=kind,
                step_index=step_index,
                preview=tool_payload,
            )
        if tool == "extract_optimized_geometry":
            molecule = _restore_json_result(tool_payload)
            if isinstance(molecule, Molecule):
                return GeometryHandoffEvent(
                    kind=kind,
                    molecule=molecule,
                    session_dir=str(session_dir) if session_dir else None,
                )

    if kind == "critic_verdict":
        return CriticVerdictEvent(
            kind=kind,
            verdict=CriticVerdict.model_validate(payload),
        )

    if kind in {"tool_error", "llm_error"}:
        return ErrorEvent(
            kind=kind,
            title=str(payload.get("error_type") or kind),
            message=str(payload.get("message") or "Unknown error"),
            details=payload,
        )

    if kind == "session_summary":
        return SessionSummaryEvent(
            kind=kind,
            blocked=bool(payload.get("blocked")),
            block_reason=payload.get("block_reason"),
            total_steps_executed=int(payload.get("total_steps_executed") or 0),
            total_steps_planned=int(payload.get("total_steps_planned") or 0),
        )

    return IgnoredEvent(kind=kind, payload=payload)


def session_completed(session_dir: Path) -> bool:
    metadata_path = session_dir / "session_metadata.json"
    return metadata_path.exists()


def _coerce_int(value: Any) -> int | None:
    if isinstance(value, int):
        return value
    if isinstance(value, str) and value.isdigit():
        return int(value)
    return None
