"""A woken cycle that delivers nothing while budget remains is woken once
more, carrying a failure report; a second silence settles.

NOVEL-2's po2 ended its woken cycle with a plan nothing could run and no
claim, and one silent session ended a goal with sixteen engine calls,
eight revisions and three and a half hours unspent (2026-09-04). The
owner ruled: one re-wake with the failure report, charged one revision.
"""

from __future__ import annotations

import json

import pytest

from .test_the_goal_loop_recovers_or_returns import (
    _delivery_rows,
    _execute,
    _loop,
    _planning_session,
    _review_payload,
)

pytestmark = pytest.mark.capability("rule:wake.refusal_is_a_deliverable")


def _declared(*ids):
    return [
        {
            "kind": "requested_observable_declared",
            "payload": {
                "observables": [
                    {
                        "observable_id": observable_id,
                        "unit": "kJ/mol",
                        "dimension": (1, 0, 0, 0, 0, 0),
                        "meaning": observable_id,
                    }
                    for observable_id in ids
                ]
            },
        }
    ]


def _claims(*ids):
    return [
        {
            "kind": "analysis_claims_recorded",
            "payload": {
                "receipt_sha256": "3" * 64,
                "record": {
                    "claims": [
                        {
                            "claim_id": claim_id,
                            "quantity_id": claim_id,
                            "display_value": 1.0,
                            "display_unit": "kcal/mol",
                            "source_receipt_sha256": "4" * 64,
                        }
                        for claim_id in ids
                    ]
                },
            },
        }
    ]


def _ledger(tmp_path):
    path = (
        tmp_path / "ws" / ".chemsmart-agent" / "goals" / "goal-t1"
    ) / "ledger.jsonl"
    return [
        json.loads(line)
        for line in path.read_text(encoding="utf-8").splitlines()
    ]


def test_a_silent_woken_cycle_is_woken_once_more_with_the_report(tmp_path):
    contexts = []

    def capture(inner):
        def step(workspace, kwargs):
            contexts.append(kwargs["goal_context"])
            return inner(workspace, kwargs)

        return step

    result = _loop(
        tmp_path,
        sessions=[
            capture(
                _planning_session(
                    "live-1",
                    review=_review_payload(),
                    wake_rows=_declared("dg_solv"),
                )
            ),
            capture(_planning_session("live-2", terminal="planned")),
            capture(
                _planning_session(
                    "live-3", terminal="planned", wake_rows=_delivery_rows()
                )
            ),
        ],
        executes=[
            _execute(
                tmp_path, failed=True, status="partial", analysis="partial"
            ),
        ],
        max_revisions=3,
    )
    assert result.settlement == "achieved"
    assert result.cycles == 3
    assert "failure_report" not in contexts[1]
    report = contexts[2]["failure_report"]
    assert report["gate"] == "goal.cycle_delivers_or_returns"
    assert "ended 'planned'" in report["diagnosis"]
    assert "record_analysis_claims" in report["route"]
    assert "charged one revision" in report["cost"]
    kinds = [entry["kind"] for entry in _ledger(tmp_path)]
    assert kinds.count("rewake_opened") == 1
    # The re-wake is charged: three revisions granted, one spent here.
    assert contexts[2]["budgets"]["revisions_remaining"] == 2


def test_a_second_silence_settles(tmp_path):
    result = _loop(
        tmp_path,
        sessions=[
            _planning_session(
                "live-1",
                review=_review_payload(),
                wake_rows=_declared("dg_solv"),
            ),
            _planning_session("live-2", terminal="planned"),
            _planning_session("live-3", terminal="planned"),
        ],
        executes=[
            _execute(
                tmp_path, failed=True, status="partial", analysis="partial"
            ),
        ],
        max_revisions=3,
    )
    assert result.settlement == "returned_to_human"
    assert result.cycles == 3
    kinds = [entry["kind"] for entry in _ledger(tmp_path)]
    assert kinds.count("rewake_opened") == 1


def test_no_revision_left_means_no_rewake(tmp_path):
    result = _loop(
        tmp_path,
        sessions=[
            _planning_session(
                "live-1",
                review=_review_payload(),
                wake_rows=_declared("dg_solv"),
            ),
            _planning_session("live-2", terminal="planned"),
        ],
        executes=[
            _execute(
                tmp_path, failed=True, status="partial", analysis="partial"
            ),
        ],
        max_revisions=0,
    )
    assert result.cycles <= 2
    kinds = [entry["kind"] for entry in _ledger(tmp_path)]
    assert "rewake_opened" not in kinds


def test_claims_under_other_names_do_not_count_as_delivery(tmp_path):
    """po3 (NOVEL-3) claimed four wall observations and stopped with the
    declared observable undelivered, 28 calls and 7 revisions in hand;
    the first re-wake keyed on 'any claim' and did not fire."""

    contexts = []

    def capture(inner):
        def step(workspace, kwargs):
            contexts.append(kwargs["goal_context"])
            return inner(workspace, kwargs)

        return step

    result = _loop(
        tmp_path,
        sessions=[
            capture(
                _planning_session(
                    "live-1",
                    review=_review_payload(),
                    wake_rows=_declared("ddg-activation-353k"),
                )
            ),
            capture(
                _planning_session(
                    "live-2",
                    terminal="blocked",
                    wake_rows=_claims("barrier-pos-a", "ddg-elec-scan-edge"),
                )
            ),
            capture(
                _planning_session(
                    "live-3", terminal="planned", wake_rows=_delivery_rows()
                )
            ),
        ],
        executes=[
            _execute(
                tmp_path, failed=True, status="partial", analysis="partial"
            ),
        ],
        max_revisions=3,
    )
    report = contexts[2]["failure_report"]
    assert "ddg-activation-353k" in report["diagnosis"]
    assert "barrier-pos-a" in report["diagnosis"]
    assert "under other names" in report["diagnosis"]
    kinds = [entry["kind"] for entry in _ledger(tmp_path)]
    assert kinds.count("rewake_opened") == 1
    assert result.cycles == 3


def test_an_id_delivered_in_an_earlier_cycle_needs_no_rewake(tmp_path):
    """A declared id claimed in cycle 2 under its id is delivered for the
    goal; a later silent cycle is not re-woken for it."""

    result = _loop(
        tmp_path,
        sessions=[
            _planning_session(
                "live-1",
                review=_review_payload(),
                wake_rows=_declared("gap"),
            ),
            _planning_session(
                "live-2", terminal="planned", wake_rows=_claims("gap")
            ),
            _planning_session("live-3", terminal="planned"),
        ],
        executes=[
            _execute(
                tmp_path, failed=True, status="partial", analysis="partial"
            ),
        ],
        max_revisions=3,
    )
    kinds = [entry["kind"] for entry in _ledger(tmp_path)]
    assert "rewake_opened" not in kinds
    assert result.cycles == 2


def test_a_first_cycle_that_planned_nothing_is_woken_once(tmp_path):
    """REACH-1 ino3's shape: cycle 1 ended 'complete' with every declared
    observable undelivered and the whole grant in hand. A cycle-1
    session has no previous outcome, which used to veto the re-wake."""

    contexts = []

    def capture(inner):
        def step(workspace, kwargs):
            contexts.append(kwargs["goal_context"])
            return inner(workspace, kwargs)

        return step

    result = _loop(
        tmp_path,
        sessions=[
            capture(
                _planning_session(
                    "live-1",
                    terminal="complete",
                    wake_rows=_declared("dg_solv"),
                )
            ),
            capture(
                _planning_session(
                    "live-2", terminal="planned", wake_rows=_delivery_rows()
                )
            ),
        ],
        executes=[],
        max_revisions=3,
    )
    kinds = [entry["kind"] for entry in _ledger(tmp_path)]
    assert kinds.count("rewake_opened") == 1
    report = contexts[1]["failure_report"]
    assert "no workflow was planned" in report["diagnosis"]
    assert "ended 'complete'" in report["diagnosis"]
    assert result.settlement == "achieved"
    assert contexts[0]["budgets"]["excursion_calls_remaining"] == 0
