"""An executed goal used to reach `achieved` with its precision open.

ROUND 10 wired the sufficiency judgement into the analysis-only
settlement and the re-wake, and into nothing on the executed path:
`_settle` built `open_delivery` from failed verdicts, stale and
unclaimed quantities and undelivered ids, and never consulted
`unresolved_requirement_ids`. So the one path the whole mechanism was
built for -- a goal that actually runs engines -- certified a delivery
whose number did not answer the precision its declaration asked for.

It shipped green because the witness that was supposed to cover it
imported `judge_sufficiency` and built an `_AnalysisDelivery` by hand.
These tests drive the goal loop instead: a real execute step writes the
run stream, the driver reads it, and the settlement is observed.
"""

from __future__ import annotations

import json

import pytest

from .test_the_goal_loop_recovers_or_returns import (
    _engine_stream,
    _loop,
    _planning_session,
    _review_payload,
)

pytestmark = pytest.mark.capability("rule:wake.restate_observable")


def _declared(tolerance=2.0):
    row = {
        "observable_id": "dg_solv",
        "unit": "kJ/mol",
        "dimension": (1, 0, 0, 0, 0, 0),
        "meaning": "solvation free energy",
    }
    if tolerance is not None:
        row["required_tolerance"] = tolerance
        row["tolerance_basis"] = "the author's stated design threshold"
    return [
        {
            "kind": "requested_observable_declared",
            "payload": {"observables": [row]},
        }
    ]


def _sufficiency(state, uncertainty):
    return {
        "observable_id": "dg_solv",
        "unit": "kJ/mol",
        "required_tolerance": 2.0,
        "uncertainty": uncertainty,
        "uncertainty_basis": "asserted",
        "state": state,
    }


def _record_claims(run_directory, sufficiency):
    """Append the claims the executor's analysis walk would have written.

    Through the real RuntimeEventStore, because the run stream is
    hash-chained and its event contract checks the record digest and the
    envelope: a hand-written row is not a run.
    """

    from chemsmart.agent._contracts import canonical_sha256
    from chemsmart.agent.runtime.event_store import RuntimeEventStore

    record = {
        "task_spec_sha256": "a" * 64,
        "status": "recorded",
        "claims": [
            {
                "claim_id": "dg_solv",
                "quantity_id": "dg_solv",
                "display_value": 1.0,
                "display_unit": "kJ/mol",
                "dimension": [1, 0, 0, 0, 0, 0],
                "source_receipt_sha256": "4" * 64,
            }
        ],
    }
    store = RuntimeEventStore(
        run_directory / "events.jsonl", session_id="water-session"
    )
    store.append(
        turn_id="turn-1",
        kind="analysis_claims_recorded",
        payload={
            "status": "recorded",
            "task_spec_sha256": "a" * 64,
            "receipt_sha256": canonical_sha256(record),
            "record": record,
            "source_receipt_sha256s": ["4" * 64],
            "claim_ids": ["dg_solv"],
            "critical_finding_count": 0,
            "sufficiency": list(sufficiency),
        },
        idempotency_key="analysis-claims:test",
    )


def _executed_run(tmp_path, *, sufficiency):
    """A run that delivered its number, with the assessment the
    executor's own analysis walk recorded beside it."""

    def step(run_directory):
        _engine_stream(tmp_path, run_directory, failed=False)
        _record_claims(run_directory, sufficiency)
        from types import SimpleNamespace

        return SimpleNamespace(status="completed", analysis_status="completed")

    return step


def _ledger(tmp_path):
    path = (
        tmp_path / "ws" / ".chemsmart-agent" / "goals" / "goal-t1"
    ) / "ledger.jsonl"
    return [
        json.loads(line)
        for line in path.read_text(encoding="utf-8").splitlines()
    ]


def test_a_short_requirement_reopens_an_executed_goal(tmp_path):
    """Budget in hand: the number stays delivered and a cycle opens."""

    result = _loop(
        tmp_path,
        sessions=[
            _planning_session(
                "live-1",
                review=_review_payload(),
                wake_rows=_declared(),
            ),
            _planning_session("live-2", terminal="planned"),
            _planning_session("live-3", terminal="planned"),
        ],
        executes=[
            _executed_run(tmp_path, sufficiency=[_sufficiency("short", 6.0)])
        ],
        max_revisions=3,
    )
    kinds = [entry["kind"] for entry in _ledger(tmp_path)]
    assert "recovery_opened" in kinds, kinds
    opened = next(
        entry
        for entry in _ledger(tmp_path)
        if entry["kind"] == "recovery_opened"
    )
    assert opened["payload"]["unresolved_requirement_ids"] == ["dg_solv"]
    assert result.settlement != "achieved"


def test_a_met_requirement_leaves_an_executed_goal_alone(tmp_path):
    """The control: the same run, resolved, still settles achieved."""

    result = _loop(
        tmp_path,
        sessions=[
            _planning_session(
                "live-1",
                review=_review_payload(),
                wake_rows=_declared(),
            ),
        ],
        executes=[
            _executed_run(tmp_path, sufficiency=[_sufficiency("met", 1.0)])
        ],
        max_revisions=3,
    )
    assert result.settlement == "achieved", result.reasons
    kinds = [entry["kind"] for entry in _ledger(tmp_path)]
    assert "recovery_opened" not in kinds


def test_a_goal_with_no_tolerance_is_untouched(tmp_path):
    """Nothing here narrows a goal that asked for no precision."""

    result = _loop(
        tmp_path,
        sessions=[
            _planning_session(
                "live-1",
                review=_review_payload(),
                wake_rows=_declared(tolerance=None),
            ),
        ],
        executes=[_executed_run(tmp_path, sufficiency=[])],
        max_revisions=3,
    )
    assert result.settlement == "achieved", result.reasons


def test_an_executed_delivery_is_asked_to_assess_what_came_back(tmp_path):
    """The route must fit the state it answers.

    An executed delivery lands `unstated` by construction, because a
    planned claim node is written before its numbers exist. Offering
    that session "plan the calculation that narrows the term you named"
    is a route to nowhere -- it has named no term. The wake asks it to
    read its own result and say what the number is worth, which costs
    no engine call.
    """

    from chemsmart.agent.driver import (
        SUFFICIENCY_SHORT_ROUTE,
        SUFFICIENCY_UNSTATED_ROUTE,
    )

    contexts = []

    def capture(inner):
        def step(workspace, kwargs):
            contexts.append(kwargs["goal_context"])
            return inner(workspace, kwargs)

        return step

    _loop(
        tmp_path,
        sessions=[
            capture(
                _planning_session(
                    "live-1",
                    review=_review_payload(),
                    wake_rows=_declared(),
                )
            ),
            capture(_planning_session("live-2", terminal="planned")),
            capture(_planning_session("live-3", terminal="planned")),
        ],
        executes=[
            _executed_run(
                tmp_path,
                sufficiency=[_sufficiency("unstated", None)],
            )
        ],
        max_revisions=3,
    )
    menu = contexts[1]["repair_menu"]
    assert menu.get("requirement_unstated") == SUFFICIENCY_UNSTATED_ROUTE
    assert "requirement_short" not in menu
    assert "no engine call" in menu["requirement_unstated"]
    # And the wake names the requirement that is open.
    deliverables = contexts[1]["deliverables"]
    assert deliverables["unresolved_requirement_ids"] == ("dg_solv",)
    assert SUFFICIENCY_SHORT_ROUTE not in str(menu)


_REFUSAL = [
    {
        "kind": "result_quantities_extracted",
        "payload": {"receipt_sha256": "e" * 64},
    },
    {
        "kind": "scientific_decision_recorded",
        "payload": {
            "receipt_sha256": "d" * 64,
            "record": {
                "evidence_refs": ["receipt:" + "e" * 64],
                "uncertainties": [],
            },
            "unreachable_observables": [
                {
                    "observable_id": "dg_solv",
                    "statement": "no method in this envelope reaches "
                    "2 kJ/mol",
                    "selector": "",
                    "jobtype": "",
                    "blocked_node_id": "",
                    "receipt_sha256s": ["e" * 64],
                    "verified": True,
                    "basis": "the requirement stands open",
                }
            ],
        },
    },
]


def test_a_sessions_refusal_closes_the_requirement_the_run_recorded(
    tmp_path,
):
    """The refusal and the assessment it answers live in two streams.

    The provider-free walker records extraction, thermochemistry,
    expressions, verdicts and claims and never a scientific decision,
    so a delivery that subtracts only its own stream's refusals can
    never see the session's. A requirement the host itself certified as
    refused re-opened the goal and spent a recovery cycle on it.
    """

    result = _loop(
        tmp_path,
        sessions=[
            _planning_session(
                "live-1",
                review=_review_payload(),
                wake_rows=_declared() + _REFUSAL,
            ),
            _planning_session(
                "live-2", terminal="planned", wake_rows=_REFUSAL
            ),
        ],
        executes=[
            _executed_run(tmp_path, sufficiency=[_sufficiency("short", 6.0)])
        ],
        max_revisions=3,
    )
    assert result.settlement == "unreachable_from_evidence"
    ledger = (
        tmp_path
        / "ws"
        / ".chemsmart-agent"
        / "goals"
        / "goal-t1"
        / "ledger.jsonl"
    )
    kinds = [
        json.loads(line)["kind"]
        for line in ledger.read_text(encoding="utf-8").splitlines()
    ]
    assert "recovery_opened" not in kinds
