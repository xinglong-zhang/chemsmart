"""Pure replay reducer for legacy and additive Runtime V2 events."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable

from chemsmart.agent._contracts import ContractError, canonical_sha256
from chemsmart.agent.runtime.events import (
    CAPABILITY_QUERIED,
    COMMAND_COMPILED,
    COMMAND_INSPECTED,
    ENGINE_BOUND,
    ENVIRONMENT_QUERIED,
    OPTIMIZED_GEOMETRY_HANDED_OFF,
    PERMISSION_RESOLVED,
    PROGRAM_BOUND,
    PROGRAM_EXECUTED,
    PROGRAM_PREFLIGHTED,
    PROJECT_VALIDATED,
    PROJECT_PROMOTED,
    PROVIDER_TURN_OBSERVED,
    RESULT_VERIFIED,
    SCIENTIFIC_DECISION_RECORDED,
    RUNTIME_TERMINATED,
    SAFE_PREVIEWED,
    SUBSTITUTION_ASSESSED,
    VALIDATOR_OBSERVED,
    WORKFLOW_PLANNED,
    RuntimeEvent,
)


@dataclass
class RuntimeState:
    """Replay state; new fields default empty for historical event streams."""

    session_id: str = ""
    turn_id: str = ""
    phase: str = "route"
    terminal_state: str = ""
    blocked_reason: str = ""
    latest_sequence: int = 0
    latest_event_hash: str = ""
    capability_receipts: list[str] = field(default_factory=list)
    environment_receipts: list[str] = field(default_factory=list)
    program_bindings: list[str] = field(default_factory=list)
    engine_bindings: list[str] = field(default_factory=list)
    project_receipts: list[str] = field(default_factory=list)
    project_promotion_receipts: list[str] = field(default_factory=list)
    scientific_decision_receipts: list[str] = field(default_factory=list)
    workflow_receipts: list[str] = field(default_factory=list)
    command_invocations: list[str] = field(default_factory=list)
    command_inspections: list[str] = field(default_factory=list)
    safe_preview_receipts: list[str] = field(default_factory=list)
    validator_receipts: list[str] = field(default_factory=list)
    preflight_receipts: list[str] = field(default_factory=list)
    substitution_receipts: list[str] = field(default_factory=list)
    result_verification_receipts: list[str] = field(default_factory=list)
    execution_receipts: list[str] = field(default_factory=list)
    optimized_geometry_handoff_receipts: list[str] = field(default_factory=list)
    permission_receipts: list[str] = field(default_factory=list)
    legacy_permission_events: int = 0
    provider_turn_receipts: list[str] = field(default_factory=list)
    api_attempt_receipts: list[str] = field(default_factory=list)
    unknown_event_kinds: list[str] = field(default_factory=list)
    seen_idempotency_keys: dict[str, str] = field(default_factory=dict)
    cwd: str = ""
    request: str = ""
    active_project: dict = field(default_factory=dict)
    active_server: dict = field(default_factory=dict)
    previous_command: str = ""
    unresolved_slots: list[str] = field(default_factory=list)
    exposed_tools: list[str] = field(default_factory=list)
    active_tool_calls: dict[str, str] = field(default_factory=dict)
    completed_tools: list[str] = field(default_factory=list)
    completed_tool_receipts: list[dict[str, str]] = field(default_factory=list)
    artifacts: list[dict] = field(default_factory=list)
    pending_approval: str = ""
    last_failure_rule_ids: list[str] = field(default_factory=list)
    shadow_violations: list[str] = field(default_factory=list)


def reduce_event(state: RuntimeState, event: RuntimeEvent) -> RuntimeState:
    if state.terminal_state:
        raise ContractError("runtime terminal state is absorbing")
    if state.session_id and event.session_id != state.session_id:
        raise ContractError("runtime event session changed during replay")
    if state.latest_sequence and event.sequence != state.latest_sequence + 1:
        raise ContractError("runtime event sequence is not contiguous")
    if not state.latest_sequence and event.sequence != 1:
        raise ContractError("runtime replay must begin at sequence 1")
    if state.latest_event_hash and event.previous_hash != state.latest_event_hash:
        raise ContractError("runtime event hash chain is broken")
    idempotency_identity = canonical_sha256(
        {"kind": event.kind, "payload": event.payload}
    )
    if event.idempotency_key and event.idempotency_key in state.seen_idempotency_keys:
        if state.seen_idempotency_keys[event.idempotency_key] != idempotency_identity:
            raise ContractError(
                "idempotency key was reused with a different action payload"
            )
        # The duplicate domain action is not applied, but the verified event
        # envelope still advances so the next hash-chain link remains valid.
        state.latest_sequence = event.sequence
        state.latest_event_hash = event.event_hash
        return state

    state.session_id = event.session_id
    state.turn_id = event.turn_id
    digest = str(
        event.payload.get("receipt_sha256")
        or event.payload.get("binding_sha256")
        or ""
    )
    targets = {
        CAPABILITY_QUERIED: state.capability_receipts,
        ENVIRONMENT_QUERIED: state.environment_receipts,
        PROGRAM_BOUND: state.program_bindings,
        ENGINE_BOUND: state.engine_bindings,
        PROJECT_VALIDATED: state.project_receipts,
        PROJECT_PROMOTED: state.project_promotion_receipts,
        SCIENTIFIC_DECISION_RECORDED: state.scientific_decision_receipts,
        WORKFLOW_PLANNED: state.workflow_receipts,
        COMMAND_COMPILED: state.command_invocations,
        COMMAND_INSPECTED: state.command_inspections,
        SAFE_PREVIEWED: state.safe_preview_receipts,
        VALIDATOR_OBSERVED: state.validator_receipts,
        PROGRAM_PREFLIGHTED: state.preflight_receipts,
        SUBSTITUTION_ASSESSED: state.substitution_receipts,
        RESULT_VERIFIED: state.result_verification_receipts,
        PROGRAM_EXECUTED: state.execution_receipts,
        OPTIMIZED_GEOMETRY_HANDED_OFF: (
            state.optimized_geometry_handoff_receipts
        ),
        PERMISSION_RESOLVED: state.permission_receipts,
        PROVIDER_TURN_OBSERVED: state.provider_turn_receipts,
        "api_attempt_observed": state.api_attempt_receipts,
    }
    if event.kind == PERMISSION_RESOLVED and event.schema_version == 1:
        state.legacy_permission_events += 1
    elif event.kind in targets:
        targets[event.kind].append(digest)
    elif event.kind == RUNTIME_TERMINATED:
        state.terminal_state = str(event.payload["terminal_state"])
        state.blocked_reason = str(event.payload.get("reason", ""))
        state.phase = "terminal"
    elif event.kind == "session_started":
        state.phase = str(event.payload.get("phase", state.phase))
        state.cwd = str(event.payload.get("cwd") or "")
    elif event.kind == "turn_started":
        state.request = str(event.payload.get("request") or "")
        state.phase = str(event.payload.get("phase") or "route")
        state.active_tool_calls = {}
        state.completed_tools = []
        state.completed_tool_receipts = []
        state.blocked_reason = ""
        state.last_failure_rule_ids = []
    elif event.kind == "exposure_planned":
        state.exposed_tools = [
            str(item) for item in event.payload.get("tools") or []
        ]
        if event.payload.get("phase"):
            state.phase = str(event.payload["phase"])
    elif event.kind == "tool_started":
        request_id = str(event.payload.get("request_id") or event.event_id)
        state.active_tool_calls[request_id] = str(
            event.payload.get("tool") or ""
        )
    elif event.kind in {"tool_succeeded", "tool_failed", "tool_request_rejected"}:
        request_id = str(event.payload.get("request_id") or "")
        state.active_tool_calls.pop(request_id, None)
        if event.kind == "tool_succeeded":
            tool = str(event.payload.get("tool") or "")
            state.completed_tools.append(tool)
            state.completed_tool_receipts.append(
                {
                    "tool": tool,
                    "verdict": str(event.payload.get("verdict") or ""),
                    "typed_command_status": str(
                        event.payload.get("typed_command_status") or ""
                    ),
                    "typed_receipt_status": str(
                        event.payload.get("typed_receipt_status") or ""
                    ),
                }
            )
        else:
            state.last_failure_rule_ids = [
                str(item) for item in event.payload.get("rule_ids") or []
            ]
    elif event.kind == "project_selected":
        state.active_project = dict(event.payload)
    elif event.kind == "command_synthesized":
        state.previous_command = str(event.payload.get("command") or "")
    elif event.kind == "clarification_requested":
        state.phase = "waiting_user"
        state.unresolved_slots = [
            str(item) for item in event.payload.get("slots") or []
        ]
    elif event.kind == "artifact_recorded":
        state.artifacts.append(dict(event.payload))
    elif event.kind == "shadow_violation":
        state.shadow_violations.append(
            str(event.payload.get("rule_id") or "runtime.shadow.tool_exposure")
        )
    elif event.kind in {"turn_completed", "turn_blocked"}:
        state.phase = "complete" if event.kind == "turn_completed" else "blocked"
        if event.kind == "turn_completed":
            state.unresolved_slots = []
            state.pending_approval = ""
        else:
            state.blocked_reason = str(event.payload.get("reason") or "blocked")
    else:
        state.unknown_event_kinds.append(event.kind)

    if event.idempotency_key:
        state.seen_idempotency_keys[event.idempotency_key] = idempotency_identity
    state.latest_sequence = event.sequence
    state.latest_event_hash = event.event_hash
    return state


def replay_events(
    events: Iterable[RuntimeEvent], *, verify_hashes: bool = True
) -> RuntimeState:
    state = RuntimeState()
    for event in events:
        if verify_hashes and not event.verify_hash():
            raise ContractError("runtime event hash verification failed")
        reduce_event(state, event)
    return state


__all__ = ["RuntimeState", "reduce_event", "replay_events"]
