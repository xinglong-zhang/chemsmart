"""A goal's word reads every stream its ledger names, however it names it.

The settlement names a pre-registered expectation the physics left from
every stream of the goal (R11 truth, Repair A), and the goal's streams were
the ones its ledger names by a session row or a run row. A ledger names a
session's stream a third way: the analysis evidence a cycle recorded when
it read results and decided. Every goal written before the session row
existed (2026-09-17) names its sessions only that way, and so does a
session that raises today. R11 truth-4's final-pin census found ax41
goal-ino3-r17 settling ``unreachable_from_evidence`` without naming the
three expectations its cycle-3 session scored diverged -- a ferrocene-
referenced potential of -0.494 V against -0.3..0.9, a quartet-doublet gap
of 103 kJ/mol against 5..90, a sulfur spin of 0.002 against 0.05..0.45 --
because that session was named only as analysis evidence. Here a woken
goal is settled over both ledger shapes through the goal loop.
"""

from __future__ import annotations

import pytest

from chemsmart.agent.runtime.event_store import RuntimeEventStore
from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1

from .test_a_goal_word_names_every_expectation_the_physics_left import (
    _POINTS,
    _TASK,
    _claims_and_completion,
    _named,
    _rows,
)
from .test_the_goal_loop_recovers_or_returns import _loop, _planning_session

pytestmark = pytest.mark.capability("tool:declare_requested_observable")

_BARRIER = "barrier-forward"


def _declaring_and_scoring_session(tmp_path):
    """Cycle 1: its host declares a band and a count, claims the barrier
    outside the band, decides, and passes a completion that scores the
    barrier diverged; the count stays undelivered, so the goal is woken
    once more."""

    build = tmp_path / "session-live-1"
    host = CommandCompiledToolHostV1(
        event_store=RuntimeEventStore(
            build / "events.jsonl", session_id="live-1"
        ),
        artifacts={},
        task_spec_sha256s=(_TASK,),
        approved_workspace=build / "ws",
    )
    reply = host.dispatch(
        turn_id="t1",
        tool_name="declare_requested_observable",
        arguments={
            "observables": [
                {
                    "observable_id": _BARRIER,
                    "unit": "kcal/mol",
                    "meaning": "E(TS) - E(reactant)",
                    "expected_low": 2.0,
                    "expected_high": 8.0,
                    "expectation_basis": "a torsional barrier of a few "
                    "kcal/mol",
                },
                {
                    "observable_id": _POINTS,
                    "unit": "1",
                    "meaning": "points on the forward IRC branch",
                },
            ]
        },
    )
    assert reply["status"] == "ok", reply
    declared = tuple(host.requested_observable_declarations.values())
    _claims_and_completion(
        host, [(_BARRIER, 11.2, "kcal/mol")], status="passed"
    )
    decided = host.dispatch(
        turn_id="t1",
        tool_name="record_scientific_decision",
        arguments={
            "decision_id": "d-barrier",
            "task_spec_sha256": _TASK,
            "assumptions": ["the barrier is the claimed number"],
            "method_rationale": "claim the barrier, then compute the count",
            "alternatives": [],
            "uncertainties": ["the barrier lies above the declared band"],
            "diagnostics": [],
            "stage_order": ["claim"],
            "evidence_refs": [],
            "findings": [],
        },
    )
    assert decided["status"] == "ok", decided
    step = _planning_session(
        "live-1", terminal="complete", wake_rows=_rows(build / "events.jsonl")
    )
    completion = [
        row["payload"]["receipt_sha256"]
        for row in _rows(build / "events.jsonl")
        if row["kind"] == "analysis_completion_evaluated"
    ][-1]
    return step, declared, completion


def _delivering_session(tmp_path, declared):
    """Cycle 2: its host claims the count and passes; the goal settles."""

    build = tmp_path / "session-live-2"
    host = CommandCompiledToolHostV1(
        event_store=RuntimeEventStore(
            build / "events.jsonl", session_id="live-2"
        ),
        artifacts={},
        task_spec_sha256s=(_TASK,),
        approved_workspace=build / "ws",
        approved_requested_observable_declarations=declared,
    )
    _claims_and_completion(host, [(_POINTS, 42, "1")], status="passed")
    return _planning_session(
        "live-2", terminal="complete", wake_rows=_rows(build / "events.jsonl")
    )


@pytest.mark.parametrize("shape", ["archived", "live"])
def test_a_goal_word_names_what_a_session_named_as_evidence_scored(
    tmp_path, shape
):
    """Cycle 1's session scores the barrier diverged at 11.2 kcal/mol and
    cycle 2's settles without claiming it again. Archived shape: no session
    result carries its id, so the ledger names cycle 1's stream only as the
    analysis evidence it recorded. Live shape: each result carries its id,
    so a session row names it too. The barrier is delivered from cycle 1 in
    both, and the word names the expectation the physics left, citing the
    completion that scored it."""

    first, declared, completion = _declaring_and_scoring_session(tmp_path)
    second = _delivering_session(tmp_path, declared)
    if shape == "live":
        first, second = _named(first, "live-1"), _named(second, "live-2")

    result = _loop(tmp_path, sessions=[first, second], executes=[])

    ledger = _rows(
        tmp_path
        / "ws"
        / ".chemsmart-agent"
        / "goals"
        / "goal-t1"
        / "ledger.jsonl"
    )
    kinds = [row["kind"] for row in ledger]
    evidence_rows = [
        row["payload"]
        for row in ledger
        if row["kind"] == "analysis_evidence_recorded"
    ]
    assert {"cycle": 1, "evidence": "runs/live-1"} in evidence_rows, kinds
    assert "rewake_opened" in kinds, kinds
    if shape == "archived":
        assert "session_stream_recorded" not in kinds, kinds
    else:
        assert "session_stream_recorded" in kinds, kinds
    (settled,) = [row for row in ledger if row["kind"] == "goal_settled"]
    reasons = " ".join(settled["payload"]["reasons"])
    assert result.settlement == "achieved_with_observations", reasons
    assert f"falsified_expectation:{_BARRIER}" in reasons
    assert completion in (
        settled["payload"]["evidence"].get("receipt_sha256s") or ()
    )
