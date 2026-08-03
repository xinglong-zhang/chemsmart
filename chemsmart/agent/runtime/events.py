"""Additive hash-chained Runtime V2 events with v1 replay compatibility."""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from datetime import datetime, timezone
from enum import Enum
from typing import Any, Mapping
from uuid import uuid4

from chemsmart.agent._contracts import (
    ContractError,
    canonical_data,
    require_sha256,
)


class EventKind(str, Enum):
    SESSION_STARTED = "session_started"
    TURN_STARTED = "turn_started"
    EXPOSURE_PLANNED = "exposure_planned"
    TOOL_REQUEST_REJECTED = "tool_request_rejected"
    TOOL_STARTED = "tool_started"
    TOOL_SUCCEEDED = "tool_succeeded"
    TOOL_FAILED = "tool_failed"
    PERMISSION_RESOLVED = "permission_resolved"
    PROJECT_SELECTED = "project_selected"
    COMMAND_SYNTHESIZED = "command_synthesized"
    CLARIFICATION_REQUESTED = "clarification_requested"
    ARTIFACT_RECORDED = "artifact_recorded"
    SHADOW_VIOLATION = "shadow_violation"
    TURN_COMPLETED = "turn_completed"
    TURN_BLOCKED = "turn_blocked"
    CAPABILITY_QUERIED = "program_capability_queried"
    ENVIRONMENT_QUERIED = "program_environment_queried"
    PROGRAM_BOUND = "program_binding_resolved"
    ENGINE_BOUND = "engine_binding_resolved"
    PROJECT_VALIDATED = "project_validated"
    PROJECT_PROMOTED = "project_artifact_promoted"
    SCIENTIFIC_DECISION_RECORDED = "scientific_decision_recorded"
    WORKFLOW_PLANNED = "command_workflow_planned"
    COMMAND_COMPILED = "command_compiled"
    COMMAND_INSPECTED = "command_inspected"
    SAFE_PREVIEWED = "safe_preview_observed"
    VALIDATOR_OBSERVED = "program_validator_observed"
    PROGRAM_PREFLIGHTED = "program_node_preflighted"
    SUBSTITUTION_ASSESSED = "program_substitution_assessed"
    RESULT_VERIFIED = "program_result_verified"
    PROGRAM_EXECUTED = "program_execution_observed"
    OPTIMIZED_GEOMETRY_HANDED_OFF = "optimized_geometry_handed_off"
    PROVIDER_TURN_OBSERVED = "provider_turn_observed"
    API_ATTEMPT_OBSERVED = "api_attempt_observed"
    RUNTIME_TERMINATED = "runtime_terminated"


CAPABILITY_QUERIED = EventKind.CAPABILITY_QUERIED.value
ENVIRONMENT_QUERIED = EventKind.ENVIRONMENT_QUERIED.value
PROGRAM_BOUND = EventKind.PROGRAM_BOUND.value
ENGINE_BOUND = EventKind.ENGINE_BOUND.value
PROJECT_VALIDATED = EventKind.PROJECT_VALIDATED.value
PROJECT_PROMOTED = EventKind.PROJECT_PROMOTED.value
SCIENTIFIC_DECISION_RECORDED = EventKind.SCIENTIFIC_DECISION_RECORDED.value
WORKFLOW_PLANNED = EventKind.WORKFLOW_PLANNED.value
COMMAND_COMPILED = EventKind.COMMAND_COMPILED.value
COMMAND_INSPECTED = EventKind.COMMAND_INSPECTED.value
SAFE_PREVIEWED = EventKind.SAFE_PREVIEWED.value
VALIDATOR_OBSERVED = EventKind.VALIDATOR_OBSERVED.value
PROGRAM_PREFLIGHTED = EventKind.PROGRAM_PREFLIGHTED.value
SUBSTITUTION_ASSESSED = EventKind.SUBSTITUTION_ASSESSED.value
RESULT_VERIFIED = EventKind.RESULT_VERIFIED.value
PROGRAM_EXECUTED = EventKind.PROGRAM_EXECUTED.value
OPTIMIZED_GEOMETRY_HANDED_OFF = EventKind.OPTIMIZED_GEOMETRY_HANDED_OFF.value
PERMISSION_RESOLVED = EventKind.PERMISSION_RESOLVED.value
PROVIDER_TURN_OBSERVED = EventKind.PROVIDER_TURN_OBSERVED.value
API_ATTEMPT_OBSERVED = EventKind.API_ATTEMPT_OBSERVED.value
RUNTIME_TERMINATED = EventKind.RUNTIME_TERMINATED.value

_RECEIPT_EVENTS = frozenset(
    {
        CAPABILITY_QUERIED,
        ENVIRONMENT_QUERIED,
        PROGRAM_BOUND,
        ENGINE_BOUND,
        PROJECT_VALIDATED,
        PROJECT_PROMOTED,
        SCIENTIFIC_DECISION_RECORDED,
        WORKFLOW_PLANNED,
        COMMAND_COMPILED,
        COMMAND_INSPECTED,
        SAFE_PREVIEWED,
        VALIDATOR_OBSERVED,
        PROGRAM_PREFLIGHTED,
        SUBSTITUTION_ASSESSED,
        RESULT_VERIFIED,
        PROGRAM_EXECUTED,
        OPTIMIZED_GEOMETRY_HANDED_OFF,
        PERMISSION_RESOLVED,
        PROVIDER_TURN_OBSERVED,
        API_ATTEMPT_OBSERVED,
    }
)


@dataclass(frozen=True)
class RuntimeEvent:
    schema_version: int
    sequence: int
    event_id: str
    session_id: str
    turn_id: str
    kind: str
    timestamp: str
    payload: dict[str, Any] = field(default_factory=dict)
    idempotency_key: str = ""
    previous_hash: str = ""
    event_hash: str = ""

    def __post_init__(self) -> None:
        if self.schema_version not in {1, 2}:
            raise ContractError("unsupported runtime event schema")
        if self.sequence < 1:
            raise ContractError("runtime event sequence must be positive")
        if not self.event_id or not self.session_id or not self.kind:
            raise ContractError("runtime event identity must not be empty")
        _validate_payload(self.schema_version, self.kind, self.payload)

    @classmethod
    def create(
        cls,
        *,
        sequence: int,
        session_id: str,
        turn_id: str,
        kind: str,
        payload: Mapping[str, Any] | None = None,
        previous_hash: str = "",
        idempotency_key: str = "",
        timestamp: str | None = None,
        event_id: str | None = None,
    ) -> "RuntimeEvent":
        body = {
            "schema_version": 2,
            "sequence": sequence,
            "event_id": event_id or uuid4().hex,
            "session_id": session_id,
            "turn_id": turn_id,
            "kind": str(getattr(kind, "value", kind)),
            "timestamp": timestamp or datetime.now(timezone.utc).isoformat(),
            "payload": canonical_data(dict(payload or {})),
            "idempotency_key": idempotency_key,
            "previous_hash": previous_hash,
        }
        return cls(**body, event_hash=_event_hash(body, schema_version=2))

    @classmethod
    def from_dict(cls, raw: Mapping[str, Any]) -> "RuntimeEvent":
        """Load historical events while defaulting only absent envelope fields."""

        return cls(
            schema_version=int(raw.get("schema_version", 1)),
            sequence=int(raw["sequence"]),
            event_id=str(raw["event_id"]),
            session_id=str(raw["session_id"]),
            turn_id=str(raw.get("turn_id", "")),
            kind=str(getattr(raw.get("kind"), "value", raw.get("kind"))),
            timestamp=str(raw.get("timestamp", "")),
            payload=dict(raw.get("payload") or {}),
            idempotency_key=str(raw.get("idempotency_key", "")),
            previous_hash=str(raw.get("previous_hash", "")),
            event_hash=str(raw.get("event_hash", "")),
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "sequence": self.sequence,
            "event_id": self.event_id,
            "session_id": self.session_id,
            "turn_id": self.turn_id,
            "kind": self.kind,
            "timestamp": self.timestamp,
            "payload": canonical_data(self.payload),
            "idempotency_key": self.idempotency_key,
            "previous_hash": self.previous_hash,
            "event_hash": self.event_hash,
        }

    def verify_hash(self) -> bool:
        body = self.to_dict()
        body.pop("event_hash")
        return self.event_hash == _event_hash(
            body, schema_version=self.schema_version
        )


def _validate_payload(
    schema_version: int, kind: str, payload: Mapping[str, Any]
) -> None:
    if schema_version == 2 and kind in _RECEIPT_EVENTS:
        digest = payload.get("receipt_sha256") or payload.get("binding_sha256")
        try:
            require_sha256(str(digest or ""), f"{kind} receipt")
        except ContractError as exc:
            raise ContractError(
                f"{kind} requires a receipt or binding digest"
            ) from exc
    if schema_version == 2:
        _validate_typed_receipt_payload(kind, payload)
    if kind == RUNTIME_TERMINATED and payload.get("terminal_state") not in {
        "complete",
        "planned",
        "failed",
        "blocked",
        "waiting_for_approval",
    }:
        raise ContractError("terminal event requires an explicit terminal state")
    if kind == RUNTIME_TERMINATED and payload.get("terminal_state") == "complete":
        required = tuple(payload.get("required_receipt_sha256s") or ())
        green = tuple(payload.get("green_receipt_sha256s") or ())
        gate = str(payload.get("completion_gate_sha256") or "")
        if not required or required != tuple(sorted(set(required))):
            raise ContractError(
                "complete terminal event requires sorted receipt gates"
            )
        if green != required:
            raise ContractError(
                "complete terminal event requires every gate to be green"
            )
        if len(gate) != 64:
            raise ContractError("complete terminal event requires host gate digest")
    if kind == RUNTIME_TERMINATED and payload.get("terminal_state") == "planned":
        require_sha256(
            str(payload.get("plan_receipt_sha256") or ""),
            "plan_receipt_sha256",
        )


def _validate_typed_receipt_payload(
    kind: str, payload: Mapping[str, Any]
) -> None:
    enum_fields = {
        CAPABILITY_QUERIED: (
            "status",
            {
                "supported",
                "preview_only",
                "reference_only",
                "disabled",
                "unknown_program",
                "unsupported_jobtype",
                "unsupported_engine",
                "cli_schema_mismatch",
            },
        ),
        ENVIRONMENT_QUERIED: (
            "status",
            {"available", "missing", "not_declared", "not_queried"},
        ),
        PROGRAM_BOUND: ("state", {"resolved", "preview_only", "blocked"}),
        ENGINE_BOUND: ("state", {"resolved", "preview_only", "blocked"}),
        PROJECT_VALIDATED: (
            "status",
            {"valid", "invalid", "loader_unavailable"},
        ),
        PROJECT_PROMOTED: (
            "status",
            {"materialized", "validated", "rejected"},
        ),
        SCIENTIFIC_DECISION_RECORDED: ("status", {"recorded"}),
        WORKFLOW_PLANNED: ("status", {"planned"}),
        COMMAND_COMPILED: ("status", {"compiled"}),
        COMMAND_INSPECTED: ("status", {"valid", "invalid"}),
        SAFE_PREVIEWED: ("status", {"previewed", "failed"}),
        VALIDATOR_OBSERVED: ("status", {"valid", "invalid"}),
        PROGRAM_PREFLIGHTED: (
            "plan_state",
            {"blocked", "ready_for_safe_preview", "previewed"},
        ),
        RESULT_VERIFIED: (
            "status",
            {"valid", "invalid", "verifier_unavailable"},
        ),
        PROGRAM_EXECUTED: (
            "execution_state",
            {"not_started", "running", "engine_complete", "failed", "ambiguous"},
        ),
        OPTIMIZED_GEOMETRY_HANDED_OFF: (
            "status",
            {"validated_handoff"},
        ),
        PERMISSION_RESOLVED: (
            "decision",
            {"auto_allow", "allow_once", "deny", "needs_user"},
        ),
    }
    field_and_values = enum_fields.get(kind)
    if field_and_values is not None:
        field, allowed = field_and_values
        if payload.get(field) not in allowed:
            raise ContractError(f"{kind} requires a typed {field}")
    if kind == SAFE_PREVIEWED:
        if payload.get("status") == "previewed" and (
            payload.get("program_validation_status") != "valid"
            or _nonnegative_int(payload, "critical_finding_count") != 0
        ):
            raise ContractError("previewed event requires green program validation")
    if kind == VALIDATOR_OBSERVED:
        if payload.get("status") == "valid" and _nonnegative_int(
            payload, "critical_finding_count"
        ) != 0:
            raise ContractError("valid validator event cannot carry findings")
    if kind == PROGRAM_PREFLIGHTED:
        critical_count = _nonnegative_int(payload, "critical_finding_count")
        if payload.get("plan_state") == "previewed":
            require_sha256(
                str(payload.get("safe_preview_receipt_sha256") or ""),
                "safe_preview_receipt_sha256",
            )
            if critical_count:
                raise ContractError("previewed preflight cannot carry findings")


def _nonnegative_int(payload: Mapping[str, Any], field: str) -> int:
    value = payload.get(field)
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ContractError(f"event requires non-negative integer {field}")
    return value


def _event_hash(body: Mapping[str, Any], *, schema_version: int) -> str:
    if schema_version == 1:
        encoded = json.dumps(
            body, sort_keys=True, separators=(",", ":")
        ).encode("utf-8")
    else:
        encoded = json.dumps(
            canonical_data(body),
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=False,
        ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


__all__ = [
    "CAPABILITY_QUERIED",
    "API_ATTEMPT_OBSERVED",
    "COMMAND_COMPILED",
    "COMMAND_INSPECTED",
    "ENGINE_BOUND",
    "ENVIRONMENT_QUERIED",
    "PERMISSION_RESOLVED",
    "PROGRAM_EXECUTED",
    "PROGRAM_BOUND",
    "PROGRAM_PREFLIGHTED",
    "PROJECT_VALIDATED",
    "PROJECT_PROMOTED",
    "SCIENTIFIC_DECISION_RECORDED",
    "PROVIDER_TURN_OBSERVED",
    "RESULT_VERIFIED",
    "OPTIMIZED_GEOMETRY_HANDED_OFF",
    "RUNTIME_TERMINATED",
    "SAFE_PREVIEWED",
    "VALIDATOR_OBSERVED",
    "WORKFLOW_PLANNED",
    "RuntimeEvent",
    "SUBSTITUTION_ASSESSED",
    "EventKind",
]
