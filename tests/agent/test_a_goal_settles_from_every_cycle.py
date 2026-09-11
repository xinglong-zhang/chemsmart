"""A settlement reads the goal, not one stream.

NOVEL-3 ino3 claimed two declared observables under their exact ids in
cycle 2; cycle 4's settlement read one stream and named all seven
undelivered, while the workspace record held both rows (2026-09-05). And
a woken session that cited a receipt from the previous cycle's executed
chain was refused: the session host held only its own receipts (po1,
ino2).
"""

from __future__ import annotations

import json

import pytest

from chemsmart.agent._contracts import RoutedContractError
from chemsmart.agent.driver import _analysis_delivery, _goal_delivered_ids
from chemsmart.agent.goal import GoalLedger
from chemsmart.agent.workspace_record import record_run

from .test_a_guide_opens_when_something_asks import _host

pytestmark = pytest.mark.capability(
    "rule:wake.claim_by_id_costs_no_engine_call"
)


def _rows(*, claims=(), completion_misses=()):
    rows = []
    if claims:
        rows.append(
            {
                "kind": "analysis_claims_recorded",
                "payload": {
                    "receipt_sha256": "5" * 64,
                    "record": {
                        "claims": [
                            {
                                "claim_id": claim_id,
                                "quantity_id": claim_id,
                                "display_value": 0.35,
                                "display_unit": "eV",
                                "source_receipt_sha256": "6" * 64,
                            }
                            for claim_id in claims
                        ]
                    },
                },
            }
        )
    if completion_misses:
        rows.append(
            {
                "kind": "analysis_completion_evaluated",
                "payload": {
                    "receipt_sha256": "7" * 64,
                    "status": "passed",
                    "limitation_output_ids": [
                        f"declared_observable:{observable_id}"
                        for observable_id in completion_misses
                    ],
                    "declared_observable_misses": [
                        f"declared observable '{observable_id}' (eV) has no "
                        "delivered claim named "
                        f"'{observable_id}' of matching dimension"
                        for observable_id in completion_misses
                    ],
                },
            }
        )
    return rows


def _write(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join(json.dumps(row) for row in rows) + "\n", encoding="utf-8"
    )
    return path


def test_an_id_delivered_in_an_earlier_cycle_is_not_undelivered(tmp_path):
    workspace = tmp_path / "ws"
    cycle_two = _write(
        tmp_path / "c2.jsonl", _rows(claims=("quartet-doublet-gap",))
    )
    record_run(
        workspace,
        goal_id="goal-ino3",
        cycle=2,
        run_events_path=cycle_two,
        run="goals/goal-ino3/runs/cycle-2",
    )
    delivered = _goal_delivered_ids(workspace, "goal-ino3")
    assert delivered["quartet-doublet-gap"]["cycle"] == 2
    assert _goal_delivered_ids(workspace, "goal-other") == {}

    cycle_four = _write(
        tmp_path / "c4.jsonl",
        _rows(
            claims=("charge-delta-ni",),
            completion_misses=("quartet-doublet-gap", "spin-nickel"),
        ),
    )
    alone = _analysis_delivery(cycle_four)
    assert alone.undelivered_declared_ids == (
        "quartet-doublet-gap",
        "spin-nickel",
    )
    joined = _analysis_delivery(cycle_four, goal_delivered_ids=delivered)
    assert joined.undelivered_declared_ids == ("spin-nickel",)
    assert joined.delivered_in_earlier_cycles == (
        "quartet-doublet-gap (delivered in cycle 2)",
    )
    (miss,) = joined.open_declared_misses
    assert "'spin-nickel'" in miss and "of matching dimension" in miss


def test_the_settlement_names_the_cycle_and_the_full_miss(tmp_path):
    from chemsmart.agent.driver import _settle_from_delivery
    from chemsmart.agent.goal import GoalRecordV1

    workspace = tmp_path / "ws"
    record_run(
        workspace,
        goal_id="goal-ino3",
        cycle=2,
        run_events_path=_write(
            tmp_path / "c2.jsonl", _rows(claims=("quartet-doublet-gap",))
        ),
        run="goals/goal-ino3/runs/cycle-2",
    )
    goal_dir = workspace / ".chemsmart-agent" / "goals" / "goal-ino3"
    ledger = GoalLedger(goal_dir)
    ledger.create(
        GoalRecordV1(
            schema_version="chemsmart.goal.v1",
            goal_id="goal-ino3",
            task_spec_sha256="a" * 64,
            scientific_identity_sha256="",
            conditions={"solvents": (), "thermochemistry": ()},
            envelope={
                "allowed_program_engines": (),
                "max_engine_calls": 30,
                "episode_wall_time_seconds": 21600.0,
                "max_excursion_calls": 0,
            },
            max_revisions=8,
            granted_by="claude-owner-delegated-reviewer",
            initial_review_sha256="",
            created_at="2026-09-05T00:00:00+00:00",
        )
    )
    stream = _write(
        tmp_path / "c4.jsonl",
        _rows(
            claims=("charge-delta-ni",),
            completion_misses=("quartet-doublet-gap", "spin-nickel"),
        ),
    )
    result = _settle_from_delivery(
        ledger,
        goal_id="goal-ino3",
        cycles=4,
        revisions_admitted=2,
        events_path=stream,
        terminal="complete",
        workspace=workspace,
    )
    assert result.settlement == "returned_to_human"
    text = " ".join(result.reasons)
    assert "in any cycle: spin-nickel" in text
    assert "quartet-doublet-gap (delivered in cycle 2)" in text
    assert "of matching dimension" in text


def test_a_decision_may_cite_a_receipt_a_recorded_run_minted(tmp_path):
    workspace = tmp_path / "ws"
    _write(
        workspace
        / ".chemsmart-agent"
        / "goals"
        / "goal-x"
        / "runs"
        / "cycle-1"
        / "events.jsonl",
        [
            {
                "kind": "thermochemistry_derived",
                "payload": {"receipt_sha256": "c" * 64},
            }
        ],
    )
    host = _host(tmp_path, run_evidence_root=workspace)
    assert host._recorded_run_receipt("c" * 64)
    assert not host._recorded_run_receipt("d" * 64)
    base = {
        "decision_id": "d1",
        "assumptions": ["a"],
        "method_rationale": "r",
        "alternatives": ["b"],
        "uncertainties": ["u"],
        "diagnostics": ["g"],
        "stage_order": ["s"],
        "evidence_refs": [],
    }
    accepted = host.dispatch(
        turn_id="t1",
        tool_name="record_scientific_decision",
        arguments={**base, "postprocessing_receipt_sha256s": ["c" * 64]},
    )
    assert accepted["status"] == "ok"
    with pytest.raises(RoutedContractError) as refused:
        host.dispatch(
            turn_id="t2",
            tool_name="record_scientific_decision",
            arguments={
                **base,
                "decision_id": "d2",
                "postprocessing_receipt_sha256s": ["d" * 64],
            },
        )
    report = refused.value.failure_report
    assert report["gate"] == "decision.receipt_is_one_the_host_minted"
    assert "inspect_run" in report["route"]
