from __future__ import annotations

import hashlib
import json

from chemsmart.agent.runtime.events import (
    CAPABILITY_QUERIED,
    RUNTIME_TERMINATED,
    RuntimeEvent,
)
from chemsmart.agent.runtime.reducer import replay_events


def _legacy_event(*, sequence, kind, payload, previous_hash=""):
    body = {
        "schema_version": 1,
        "sequence": sequence,
        "event_id": f"legacy-{sequence}",
        "session_id": "legacy-session",
        "turn_id": "legacy-turn",
        "kind": kind,
        "timestamp": f"2026-01-01T00:00:0{sequence}+00:00",
        "payload": payload,
        "idempotency_key": "",
        "previous_hash": previous_hash,
    }
    digest = hashlib.sha256(
        json.dumps(body, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    return RuntimeEvent.from_dict({**body, "event_hash": digest})


def test_legacy_permission_event_replays_without_new_receipt_shape():
    started = _legacy_event(
        sequence=1, kind="session_started", payload={"phase": "route"}
    )
    permission = _legacy_event(
        sequence=2,
        kind="permission_resolved",
        payload={"decision": "approved", "tool": "read_project"},
        previous_hash=started.event_hash,
    )
    blocked = _legacy_event(
        sequence=3,
        kind="turn_blocked",
        payload={"reason": "missing geometry"},
        previous_hash=permission.event_hash,
    )

    state = replay_events((started, permission, blocked))

    assert state.legacy_permission_events == 1
    assert state.permission_receipts == []
    assert state.phase == "blocked"
    assert state.latest_event_hash == blocked.event_hash


def test_duplicate_idempotency_advances_chain_without_reapplying_state():
    started = RuntimeEvent.create(
        sequence=1,
        session_id="session",
        turn_id="turn",
        kind="session_started",
    )
    first = RuntimeEvent.create(
        sequence=2,
        session_id="session",
        turn_id="turn",
        kind=CAPABILITY_QUERIED,
        payload={"receipt_sha256": "a" * 64, "status": "supported"},
        previous_hash=started.event_hash,
        idempotency_key="capability-query",
    )
    duplicate = RuntimeEvent.create(
        sequence=3,
        session_id="session",
        turn_id="turn",
        kind=CAPABILITY_QUERIED,
        payload={"receipt_sha256": "a" * 64, "status": "supported"},
        previous_hash=first.event_hash,
        idempotency_key="capability-query",
    )
    terminal = RuntimeEvent.create(
        sequence=4,
        session_id="session",
        turn_id="turn",
        kind=RUNTIME_TERMINATED,
        payload={"terminal_state": "blocked", "reason": "fixture"},
        previous_hash=duplicate.event_hash,
    )

    state = replay_events((started, first, duplicate, terminal))

    assert state.capability_receipts == ["a" * 64]
    assert state.latest_sequence == 4
    assert state.terminal_state == "blocked"
