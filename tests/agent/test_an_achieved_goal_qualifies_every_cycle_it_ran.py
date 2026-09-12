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
