"""Decision-log adapters for Phase 1 TUI cells."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from chemsmart.agent.core import CriticVerdict, Plan, render_plan


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
class DryRunInputEvent(AgentEvent):
    inputfile: str | None
    content: str


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
    | DryRunInputEvent
    | CriticVerdictEvent
    | ErrorEvent
    | SessionSummaryEvent
    | IgnoredEvent
)


def parse_decision_event(entry: dict[str, Any]) -> _AGENT_EVENT_TYPES:
    kind = str(entry.get("kind") or "")
    payload = entry.get("payload") or {}

    if kind == "request":
        return RequestEvent(
            kind=kind, request=str(payload.get("request") or "")
        )

    if kind == "plan":
        plan = Plan.model_validate(payload)
        return PlanEvent(kind=kind, plan=plan, text=render_plan(plan))

    if kind == "tool_result" and payload.get("tool") == "dry_run_input":
        tool_payload = payload.get("payload") or {}
        return DryRunInputEvent(
            kind=kind,
            inputfile=tool_payload.get("inputfile"),
            content=str(tool_payload.get("content") or ""),
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
