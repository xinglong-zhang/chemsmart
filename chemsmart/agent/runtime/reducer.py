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
    payload = event.payload
    updates: dict[str, Any] = {
        "latest_sequence": event.sequence,
        "latest_event_hash": event.event_hash,
    }
    if event.kind is EventKind.SESSION_STARTED:
        updates.update(
            session_id=event.session_id,
            cwd=str(payload.get("cwd") or ""),
        )
    elif event.kind is EventKind.TURN_STARTED:
        updates.update(
            turn_id=event.turn_id,
            request=str(payload.get("request") or ""),
            phase=TaskPhase(
                str(payload.get("phase") or TaskPhase.ROUTE.value)
            ),
            blocked_reason="",
            last_failure_rule_ids=[],
        )
    elif event.kind is EventKind.EXPOSURE_PLANNED:
        updates["exposed_tools"] = [
            str(item) for item in payload.get("tools") or []
        ]
        if payload.get("phase"):
            updates["phase"] = TaskPhase(str(payload["phase"]))
    elif event.kind is EventKind.TOOL_STARTED:
        calls = dict(state.active_tool_calls)
        calls[str(payload.get("request_id") or event.event_id)] = str(
            payload.get("tool") or ""
        )
        updates["active_tool_calls"] = calls
    elif event.kind in {EventKind.TOOL_SUCCEEDED, EventKind.TOOL_FAILED}:
        calls = dict(state.active_tool_calls)
        calls.pop(str(payload.get("request_id") or ""), None)
        updates["active_tool_calls"] = calls
        if event.kind is EventKind.TOOL_SUCCEEDED:
            updates["completed_tools"] = [
                *state.completed_tools,
                str(payload.get("tool") or ""),
            ]
        else:
            updates["last_failure_rule_ids"] = [
                str(item) for item in payload.get("rule_ids") or []
            ]
    elif event.kind is EventKind.PERMISSION_RESOLVED:
        updates["pending_approval"] = (
            str(payload.get("tool") or "")
            if payload.get("decision") == "needs_user"
            else ""
        )
    elif event.kind is EventKind.PROJECT_SELECTED:
        updates["active_project"] = WorkspaceRef.model_validate(payload)
    elif event.kind is EventKind.COMMAND_SYNTHESIZED:
        updates["previous_command"] = str(payload.get("command") or "")
    elif event.kind is EventKind.CLARIFICATION_REQUESTED:
        updates["phase"] = TaskPhase.WAITING_USER
        updates["unresolved_slots"] = [
            str(item) for item in payload.get("slots") or []
        ]
    elif event.kind is EventKind.ARTIFACT_RECORDED:
        artifact = ArtifactRef.model_validate(payload)
        updates["artifacts"] = [*state.artifacts, artifact]
    elif event.kind is EventKind.SHADOW_VIOLATION:
        updates["shadow_violations"] = [
            *state.shadow_violations,
            str(payload.get("rule_id") or "runtime.shadow.tool_exposure"),
        ]
    elif event.kind is EventKind.TURN_COMPLETED:
        updates.update(
            phase=TaskPhase.COMPLETE,
            unresolved_slots=[],
            pending_approval="",
        )
    elif event.kind is EventKind.TURN_BLOCKED:
        updates.update(
            phase=TaskPhase.BLOCKED,
            blocked_reason=str(payload.get("reason") or "blocked"),
        )
    return state.model_copy(update=updates, deep=True)


__all__ = ["RuntimeState", "apply_event", "reduce_events"]
