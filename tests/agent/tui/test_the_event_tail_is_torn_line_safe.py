"""The live feed reads the same evidence the audit reads, without racing it.

The writer appends whole JSON lines under its own lock; the reader must
never parse a torn tail, never crash on a malformed line, and never re-read
consumed bytes.
"""

from __future__ import annotations

import json
from pathlib import Path

from chemsmart.agent.tui.monitor import (
    EventTailerV1,
    execution_signal,
    planning_feed_update,
)


def _line(kind: str, payload: dict) -> bytes:
    return (
        json.dumps({"kind": kind, "payload": payload}).encode("utf-8") + b"\n"
    )


def test_a_torn_tail_waits_and_is_delivered_whole(tmp_path: Path):
    stream = tmp_path / "events.jsonl"
    tailer = EventTailerV1(stream)

    assert tailer.poll() == ()  # not yet created

    whole = _line("tool_started", {"tool": "preview_command", "request_id": "c1"})
    stream.write_bytes(whole[:20])  # writer mid-line
    assert tailer.poll() == ()

    stream.write_bytes(whole)  # line completed
    events = tailer.poll()
    assert len(events) == 1
    assert events[0]["kind"] == "tool_started"

    # Nothing is re-read once consumed.
    assert tailer.poll() == ()


def test_a_malformed_line_becomes_a_gap_not_a_crash(tmp_path: Path):
    stream = tmp_path / "events.jsonl"
    stream.write_bytes(
        _line("tool_started", {"tool": "compile_command", "request_id": "c1"})
        + b"{this is not json}\n"
        + _line("tool_succeeded", {"tool": "compile_command", "request_id": "c1"})
    )
    tailer = EventTailerV1(stream)
    events = tailer.poll()

    kinds = [event["kind"] for event in events]
    assert kinds == ["tool_started", "__display_gap__", "tool_succeeded"]
    gap_row = planning_feed_update(events[1])
    assert gap_row is not None and gap_row.state == "note"


def test_planning_rows_read_as_sentences(tmp_path: Path):
    started = planning_feed_update(
        {
            "kind": "tool_started",
            "payload": {"tool": "preview_command", "request_id": "c9"},
        }
    )
    assert started is not None
    assert started.text == "~ running the safe preview..."
    assert started.state == "running"

    settled = planning_feed_update(
        {
            "kind": "tool_succeeded",
            "payload": {
                "tool": "extract_result_quantities",
                "request_id": "c9",
                "canonical_result": {
                    "result": {"quantities": [{"quantity_id": "e_scf"}]}
                },
            },
        }
    )
    assert settled is not None
    assert settled.text == "extracted e_scf"
    assert settled.state == "finished"

    refused = planning_feed_update(
        {
            "kind": "tool_failed",
            "payload": {
                "tool": "plan_scientific_workflow",
                "error_class": "ContractError",
                "canonical_result": {"message": "unknown workflow ID"},
            },
        }
    )
    assert refused is not None
    assert refused.state == "failed"
    assert "unknown workflow ID" in refused.text
    # The Python exception class is bookkeeping; the host message suffices.
    assert "ContractError" not in refused.text

    silent = planning_feed_update({"kind": "artifact_recorded", "payload": {}})
    assert silent is None
