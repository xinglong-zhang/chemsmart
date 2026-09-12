"""An achieved goal qualifies every validated node of every cycle it ran.

The qualification rows were written from the settling cycle's outcome
alone, so a goal that executed in cycle one and settled ``achieved`` in
an analysis-only cycle two wrote no ``qualified`` row for the nodes it
ran (PySCF round, g2, 2026-09-12: two validated PySCF nodes, zero rows,
while g1, settled in its first cycle, wrote two). Whether a capability is
qualified cannot depend on which cycle a goal happened to settle in.
"""

from __future__ import annotations

import pytest

from chemsmart.agent.driver import _qualification_entries_for_goal
from chemsmart.agent.terminal_states import (
    NodeTerminalStateV1,
    RunOutcomeV1,
)

pytestmark = pytest.mark.capability("program_jobtype:*")


def _outcome(*nodes):
    return RunOutcomeV1(
        run_id="r",
        workflow_id="w",
        plan_sha256="1" * 64,
        approval_sha256="2" * 64,
        workflow_state="validated" if nodes else "completed",
        nodes=tuple(
            NodeTerminalStateV1(
                node_id=node_id, program="pyscf", jobtype=jobtype, state=state
            )
            for node_id, jobtype, state in nodes
        ),
    )


def test_rows_come_from_every_recorded_cycle():
    ran = _outcome(("opt", "opt", "validated"), ("hess", "hess", "validated"))
    analysis_only = _outcome()
    ledger = (
        {
            "kind": "run_recorded",
            "payload": {"cycle": 1, "run": "goals/g/runs/cycle-1"},
        },
        {
            "kind": "run_recorded",
            "payload": {"cycle": 2, "run": "goals/g/runs/cycle-2"},
        },
    )
    entries = _qualification_entries_for_goal(
        ledger_entries=ledger,
        goal_id="g",
        current_run="goals/g/runs/cycle-2",
        current_outcome=analysis_only,
        read_outcome={"goals/g/runs/cycle-1": ran}.__getitem__,
    )
    assert sorted((e["id"], e["run"]) for e in entries) == [
        ("pyscf:cpu:hess", "goals/g/runs/cycle-1"),
        ("pyscf:cpu:opt", "goals/g/runs/cycle-1"),
    ]


def test_an_analysis_only_settlement_qualifies_the_cycle_that_ran(tmp_path):
    """The rows were written only by the executed run's settlement.

    A goal that ran in cycle one and settled ``achieved`` from an
    analysis-only cycle two -- the exact shape the reader above was
    repaired for -- went through the other settlement and wrote nothing
    (PySCF round g5, 2026-09-12). Driven through the analysis-only
    settlement over a stream the real receipt writer built.
    """

    import json
    import shutil

    from chemsmart.agent._contracts import canonical_sha256
    from chemsmart.agent.driver import _settle_from_delivery
    from chemsmart.agent.execution import build_program_execution_receipt
    from chemsmart.agent.goal import GoalLedger
    from chemsmart.agent.runtime.event_store import RuntimeEventStore
    from chemsmart.agent.workspace_record import record_run

    from .test_runtime_v2_launch_fence import _reserve

    # A validated node, through the real launch fence and receipt writer.
    built = tmp_path / "build" / "events" / "runtime.jsonl"
    store = RuntimeEventStore(built, session_id="water-session")
    _, plan, _materialized, _approval, invocation = _reserve(
        store, tmp_path / "build"
    )
    receipt = build_program_execution_receipt(
        invocation,
        execution_state="validated",
        exit_status=0,
        child_exit_status=0,
        engine_complete=True,
        validated=True,
        validator_receipt_sha256s=("e" * 64,),
        result_validation_receipt_sha256="e" * 64,
        started_at="2026-08-04T00:00:00+00:00",
        finished_at="2026-08-04T00:00:05+00:00",
    )
    store.record_program_execution_receipt(
        turn_id="turn-1",
        workflow_id=plan.workflow_id,
        run_id="run.water-approval",
        receipt=receipt,
    )
    # The validator's own record, which is where a node's program and
    # job type come from on read.
    verification = {
        "node_id": "sp-initial",
        "program": "orca",
        "jobtype": "sp",
        "state": "valid",
        "observations": {"jobtype": "sp", "program": "orca"},
        "output_artifacts": [],
    }
    digest = canonical_sha256(verification)
    store.append(
        turn_id="turn-1",
        kind="program_result_verified",
        payload={
            "node_id": "sp-initial",
            "status": "valid",
            "critical_finding_count": 0,
            "receipt_sha256": digest,
            "record": {**verification, "receipt_sha256": digest},
        },
    )
    workspace = tmp_path / "ws"
    run_dir = (
        workspace / ".chemsmart-agent" / "goals" / "g" / "runs" / "cycle-1"
    )
    run_dir.mkdir(parents=True)
    shutil.copy2(built, run_dir / "events.jsonl")
    record_run(
        workspace,
        goal_id="g",
        cycle=1,
        run_events_path=run_dir / "events.jsonl",
        run="goals/g/runs/cycle-1",
    )
    ledger = GoalLedger(workspace / ".chemsmart-agent" / "goals" / "g")
    ledger.append(
        "run_recorded",
        {
            "cycle": 1,
            "run": "goals/g/runs/cycle-1",
            "engine_calls_consumed": 1,
        },
    )
    session = tmp_path / "cycle-2.jsonl"
    session.write_text(
        "\n".join(
            json.dumps(row)
            for row in (
                {
                    "kind": "analysis_claims_recorded",
                    "payload": {
                        "receipt_sha256": "5" * 64,
                        "record": {
                            "claims": [
                                {
                                    "claim_id": "e",
                                    "quantity_id": "e",
                                    "display_value": -76.0,
                                    "display_unit": "hartree",
                                    "dimension": [1, 0, 0, 0, 0, 0],
                                    "source_receipt_sha256": "6" * 64,
                                }
                            ]
                        },
                    },
                },
                {
                    "kind": "analysis_completion_evaluated",
                    "payload": {
                        "receipt_sha256": "7" * 64,
                        "status": "passed",
                        "limitation_output_ids": [],
                        "declared_observable_misses": [],
                    },
                },
            )
        )
        + "\n",
        encoding="utf-8",
    )
    result = _settle_from_delivery(
        ledger,
        goal_id="g",
        cycles=2,
        revisions_admitted=1,
        events_path=session,
        terminal="complete",
        workspace=workspace,
    )
    assert result.settlement in {
        "achieved",
        "achieved_with_observations",
    }, result.reasons
    qualified = [
        entry["payload"]
        for entry in ledger.entries()
        if entry["kind"] == "qualified"
    ]
    assert [row["run"] for row in qualified] == [
        "goals/g/runs/cycle-1"
    ], qualified
