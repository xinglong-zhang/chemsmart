"""Versioned append-only runtime event definitions."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from enum import Enum
from typing import Any
from uuid import uuid4

from pydantic import BaseModel, ConfigDict, Field


class EventKind(str, Enum):
    SESSION_STARTED = "session_started"
    TURN_STARTED = "turn_started"
    EXPOSURE_PLANNED = "exposure_planned"
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


class RuntimeEvent(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    schema_version: int = 1
    sequence: int = Field(ge=1)
    event_id: str
    session_id: str
    turn_id: str
    kind: EventKind
    timestamp: str
    payload: dict[str, Any] = Field(default_factory=dict)
    idempotency_key: str = ""
    previous_hash: str = ""
    event_hash: str

    @classmethod
    def create(
        cls,
        *,
        sequence: int,
        session_id: str,
        turn_id: str,
        kind: EventKind,
        payload: dict[str, Any],
        previous_hash: str,
        idempotency_key: str = "",
    ) -> "RuntimeEvent":
        body = {
            "schema_version": 1,
            "sequence": sequence,
            "event_id": uuid4().hex,
            "session_id": session_id,
            "turn_id": turn_id,
            "kind": kind.value,
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "payload": payload,
            "idempotency_key": idempotency_key,
            "previous_hash": previous_hash,
        }
        event_hash = hashlib.sha256(
            json.dumps(body, sort_keys=True, separators=(",", ":")).encode()
        ).hexdigest()
        return cls(**body, event_hash=event_hash)

    def verify_hash(self) -> bool:
        body = self.model_dump(exclude={"event_hash"}, mode="json")
        digest = hashlib.sha256(
            json.dumps(body, sort_keys=True, separators=(",", ":")).encode()
        ).hexdigest()
        return digest == self.event_hash


__all__ = ["EventKind", "RuntimeEvent"]
