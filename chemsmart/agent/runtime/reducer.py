"""Pure reducer for reconstructing runtime state from append-only events."""

from __future__ import annotations

from typing import Any, Iterable

from pydantic import BaseModel, ConfigDict, Field

from chemsmart.agent.runtime.contracts import (
    ArtifactRef,
    TaskPhase,
    WorkspaceRef,
)
from chemsmart.agent.runtime.events import EventKind, RuntimeEvent


class RuntimeState(BaseModel):
    model_config = ConfigDict(extra="forbid")

    session_id: str = ""
    turn_id: str = ""
    cwd: str = ""
    phase: TaskPhase = TaskPhase.ROUTE
    request: str = ""
    active_project: WorkspaceRef | None = None
    active_server: WorkspaceRef | None = None
    previous_command: str = ""
    unresolved_slots: list[str] = Field(default_factory=list)
    exposed_tools: list[str] = Field(default_factory=list)
    active_tool_calls: dict[str, str] = Field(default_factory=dict)
    completed_tools: list[str] = Field(default_factory=list)
    completed_tool_receipts: list[dict[str, str]] = Field(default_factory=list)
    artifacts: list[ArtifactRef] = Field(default_factory=list)
    pending_approval: str = ""
    blocked_reason: str = ""
    last_failure_rule_ids: list[str] = Field(default_factory=list)
    shadow_violations: list[str] = Field(default_factory=list)
    latest_sequence: int = 0
    latest_event_hash: str = ""


def reduce_events(
    events: Iterable[RuntimeEvent],
    initial: RuntimeState | None = None,
) -> RuntimeState:
    state = initial.model_copy(deep=True) if initial else RuntimeState()
    for event in events:
        state = apply_event(state, event)
    return state


def apply_event(state: RuntimeState, event: RuntimeEvent) -> RuntimeState:
    updates: dict[str, Any] = {
        "latest_sequence": event.sequence,
        "latest_event_hash": event.event_hash,
    }
    handler = _EVENT_HANDLERS.get(event.kind)
    if handler is not None:
        updates.update(handler(state, event, event.payload))
    return state.model_copy(update=updates, deep=True)


def _on_session_started(
    state: RuntimeState, event: RuntimeEvent, payload: dict[str, Any]
) -> dict[str, Any]:
    return {
        "session_id": event.session_id,
        "cwd": str(payload.get("cwd") or ""),
    }


def _on_turn_started(
    state: RuntimeState, event: RuntimeEvent, payload: dict[str, Any]
) -> dict[str, Any]:
    return {
        "turn_id": event.turn_id,
        "request": str(payload.get("request") or ""),
        "phase": TaskPhase(str(payload.get("phase") or TaskPhase.ROUTE.value)),
        "active_tool_calls": {},
        "completed_tools": [],
        "completed_tool_receipts": [],
        "blocked_reason": "",
        "last_failure_rule_ids": [],
    }


def _on_exposure_planned(
    state: RuntimeState, event: RuntimeEvent, payload: dict[str, Any]
) -> dict[str, Any]:
    updates: dict[str, Any] = {
        "exposed_tools": [str(item) for item in payload.get("tools") or []]
    }
    if payload.get("phase"):
        updates["phase"] = TaskPhase(str(payload["phase"]))
    return updates


def _on_tool_started(
    state: RuntimeState, event: RuntimeEvent, payload: dict[str, Any]
) -> dict[str, Any]:
    calls = dict(state.active_tool_calls)
    calls[str(payload.get("request_id") or event.event_id)] = str(
        payload.get("tool") or ""
    )
    return {"active_tool_calls": calls}


def _on_tool_finished(
    state: RuntimeState, event: RuntimeEvent, payload: dict[str, Any]
) -> dict[str, Any]:
    calls = dict(state.active_tool_calls)
    calls.pop(str(payload.get("request_id") or ""), None)
    updates: dict[str, Any] = {"active_tool_calls": calls}
    if event.kind is EventKind.TOOL_SUCCEEDED:
        tool_name = str(payload.get("tool") or "")
        updates["completed_tools"] = [*state.completed_tools, tool_name]
        updates["completed_tool_receipts"] = [
            *state.completed_tool_receipts,
            {
                "tool": tool_name,
                "verdict": str(payload.get("verdict") or ""),
            },
        ]
    else:
        updates["last_failure_rule_ids"] = [
            str(item) for item in payload.get("rule_ids") or []
        ]
    return updates


def _on_permission_resolved(
    state: RuntimeState, event: RuntimeEvent, payload: dict[str, Any]
) -> dict[str, Any]:
    return {
        "pending_approval": (
            str(payload.get("tool") or "")
            if payload.get("decision") == "needs_user"
            else ""
        )
    }


def _on_project_selected(
    state: RuntimeState, event: RuntimeEvent, payload: dict[str, Any]
) -> dict[str, Any]:
    return {"active_project": WorkspaceRef.model_validate(payload)}


def _on_command_synthesized(
    state: RuntimeState, event: RuntimeEvent, payload: dict[str, Any]
) -> dict[str, Any]:
    return {"previous_command": str(payload.get("command") or "")}


def _on_clarification_requested(
    state: RuntimeState, event: RuntimeEvent, payload: dict[str, Any]
) -> dict[str, Any]:
    return {
        "phase": TaskPhase.WAITING_USER,
        "unresolved_slots": [str(item) for item in payload.get("slots") or []],
    }


def _on_artifact_recorded(
    state: RuntimeState, event: RuntimeEvent, payload: dict[str, Any]
) -> dict[str, Any]:
    artifact = ArtifactRef.model_validate(payload)
    return {"artifacts": [*state.artifacts, artifact]}


def _on_shadow_violation(
    state: RuntimeState, event: RuntimeEvent, payload: dict[str, Any]
) -> dict[str, Any]:
    return {
        "shadow_violations": [
            *state.shadow_violations,
            str(payload.get("rule_id") or "runtime.shadow.tool_exposure"),
        ]
    }


def _on_turn_completed(
    state: RuntimeState, event: RuntimeEvent, payload: dict[str, Any]
) -> dict[str, Any]:
    return {
        "phase": TaskPhase.COMPLETE,
        "unresolved_slots": [],
        "pending_approval": "",
    }


def _on_turn_blocked(
    state: RuntimeState, event: RuntimeEvent, payload: dict[str, Any]
) -> dict[str, Any]:
    return {
        "phase": TaskPhase.BLOCKED,
        "blocked_reason": str(payload.get("reason") or "blocked"),
    }


_EVENT_HANDLERS = {
    EventKind.SESSION_STARTED: _on_session_started,
    EventKind.TURN_STARTED: _on_turn_started,
    EventKind.EXPOSURE_PLANNED: _on_exposure_planned,
    EventKind.TOOL_STARTED: _on_tool_started,
    EventKind.TOOL_SUCCEEDED: _on_tool_finished,
    EventKind.TOOL_FAILED: _on_tool_finished,
    EventKind.PERMISSION_RESOLVED: _on_permission_resolved,
    EventKind.PROJECT_SELECTED: _on_project_selected,
    EventKind.COMMAND_SYNTHESIZED: _on_command_synthesized,
    EventKind.CLARIFICATION_REQUESTED: _on_clarification_requested,
    EventKind.ARTIFACT_RECORDED: _on_artifact_recorded,
    EventKind.SHADOW_VIOLATION: _on_shadow_violation,
    EventKind.TURN_COMPLETED: _on_turn_completed,
    EventKind.TURN_BLOCKED: _on_turn_blocked,
}


__all__ = ["RuntimeState", "apply_event", "reduce_events"]
