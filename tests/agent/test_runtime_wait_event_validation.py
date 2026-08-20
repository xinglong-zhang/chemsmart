from __future__ import annotations

import json
from typing import Any

import pytest

from chemsmart.agent._contracts import ContractError, canonical_sha256
from chemsmart.agent.runtime.event_store import RuntimeEventStore
from chemsmart.agent.runtime.events import (
    TOOL_WAITING,
    TOOL_WOKE,
    RuntimeEvent,
)

_MISSING = object()
_WAITING_PAYLOAD = {
    "request_id": "execute-one",
    "tool": "execute_approved_program_node",
    "wait_kind": "approved_program_tool_dispatch",
    "timeout_seconds": 120.0,
    "provider_calls_while_waiting": 0,
}
_WOKE_PAYLOAD = {
    "request_id": "execute-one",
    "tool": "execute_approved_program_node",
    "wake_reason": "result",
    "waited_seconds": 0.0,
    "provider_calls_while_waiting": 0,
}


def _payload(
    baseline: dict[str, Any], field_name: str, value: Any
) -> dict[str, Any]:
    payload = dict(baseline)
    if value is _MISSING:
        payload.pop(field_name, None)
    else:
        payload[field_name] = value
    return payload


def _event(*, kind: str, payload: dict[str, Any]) -> RuntimeEvent:
    return RuntimeEvent(
        schema_version=2,
        sequence=1,
        event_id="event-one",
        session_id="session-one",
        turn_id="turn-one",
        kind=kind,
        timestamp="2026-08-13T00:00:00+00:00",
        payload=payload,
    )


@pytest.mark.parametrize(
    ("kind", "baseline", "field_name", "value"),
    (
        (TOOL_WAITING, _WAITING_PAYLOAD, "request_id", _MISSING),
        (TOOL_WAITING, _WAITING_PAYLOAD, "request_id", "  "),
        (TOOL_WAITING, _WAITING_PAYLOAD, "tool", ""),
        (TOOL_WAITING, _WAITING_PAYLOAD, "tool", 7),
        (TOOL_WOKE, _WOKE_PAYLOAD, "request_id", ""),
        (TOOL_WOKE, _WOKE_PAYLOAD, "tool", _MISSING),
        (TOOL_WOKE, _WOKE_PAYLOAD, "tool", "\t"),
    ),
)
def test_wait_events_reject_empty_or_non_string_identity(
    kind, baseline, field_name, value
):
    with pytest.raises(ContractError, match="requires a non-empty"):
        _event(kind=kind, payload=_payload(baseline, field_name, value))


@pytest.mark.parametrize(
    "kind,baseline",
    ((TOOL_WAITING, _WAITING_PAYLOAD), (TOOL_WOKE, _WOKE_PAYLOAD)),
)
@pytest.mark.parametrize(
    "value", (_MISSING, False, True, "0", -1, 1, float("nan"), float("inf"))
)
def test_wait_events_require_native_numeric_zero_provider_calls(
    kind, baseline, value
):
    payload = _payload(baseline, "provider_calls_while_waiting", value)
    with pytest.raises(ContractError, match="zero provider calls"):
        _event(kind=kind, payload=payload)


@pytest.mark.parametrize(
    "value",
    (_MISSING, False, True, "1", 0, -1, float("nan"), float("inf")),
)
def test_tool_waiting_requires_finite_positive_native_number(value):
    payload = _payload(_WAITING_PAYLOAD, "timeout_seconds", value)
    with pytest.raises(ContractError, match="positive timeout"):
        _event(kind=TOOL_WAITING, payload=payload)


@pytest.mark.parametrize(
    "value", (_MISSING, False, True, "0", -1, float("nan"), float("inf"))
)
def test_tool_woke_requires_finite_non_negative_native_number(value):
    payload = _payload(_WOKE_PAYLOAD, "waited_seconds", value)
    with pytest.raises(ContractError, match="non-negative wait time"):
        _event(kind=TOOL_WOKE, payload=payload)


@pytest.mark.parametrize(
    ("kind", "payload"),
    (
        (TOOL_WAITING, {**_WAITING_PAYLOAD, "timeout_seconds": 1}),
        (
            TOOL_WAITING,
            {
                **_WAITING_PAYLOAD,
                "timeout_seconds": 1.5,
                "provider_calls_while_waiting": 0.0,
            },
        ),
        (TOOL_WOKE, {**_WOKE_PAYLOAD, "waited_seconds": 0}),
        (
            TOOL_WOKE,
            {
                **_WOKE_PAYLOAD,
                "waited_seconds": 0.25,
                "provider_calls_while_waiting": 0.0,
            },
        ),
    ),
)
def test_wait_events_accept_only_valid_native_numeric_boundaries(
    kind, payload
):
    assert _event(kind=kind, payload=payload).payload == payload


def test_persisted_wait_pair_replays_after_runtime_restart(tmp_path):
    path = tmp_path / "events" / "runtime.jsonl"
    store = RuntimeEventStore(path, session_id="session-one")
    store.append(turn_id="turn-one", kind="session_started")
    store.append(
        turn_id="turn-one", kind=TOOL_WAITING, payload=_WAITING_PAYLOAD
    )
    store.append(turn_id="turn-one", kind=TOOL_WOKE, payload=_WOKE_PAYLOAD)

    state = RuntimeEventStore(path, session_id="session-one").state()

    assert state.active_waits == {}
    assert state.completed_waits == [_WOKE_PAYLOAD]


@pytest.mark.parametrize(
    ("kind", "payload"),
    (
        (TOOL_WAITING, {**_WAITING_PAYLOAD, "timeout_seconds": "120"}),
        (
            TOOL_WOKE,
            {**_WOKE_PAYLOAD, "provider_calls_while_waiting": False},
        ),
    ),
)
def test_persisted_malformed_wait_event_is_rejected_during_replay(
    tmp_path, kind, payload
):
    started = RuntimeEvent.create(
        sequence=1,
        session_id="session-one",
        turn_id="turn-one",
        kind="session_started",
    )
    raw = {
        "schema_version": 2,
        "sequence": 2,
        "event_id": "malformed-wait",
        "session_id": "session-one",
        "turn_id": "turn-one",
        "kind": kind,
        "timestamp": "2026-08-13T00:00:01+00:00",
        "payload": payload,
        "idempotency_key": "",
        "previous_hash": started.event_hash,
    }
    raw["event_hash"] = canonical_sha256(raw)
    path = tmp_path / "events" / "runtime.jsonl"
    path.parent.mkdir()
    path.write_text(
        "\n".join(
            (
                json.dumps(started.to_dict(), sort_keys=True),
                json.dumps(raw, sort_keys=True),
            )
        )
        + "\n",
        encoding="utf-8",
    )

    store = RuntimeEventStore(path, session_id="session-one")
    with pytest.raises(ContractError, match="JSONL line 2"):
        store.state()
