"""History is evidence; it is not state.

Two mirrored defects from one missing idea. Within a stream, every
sufficiency row was kept and an id was open if *any* row was open, so a
requirement assessed `unstated`, re-claimed with a proper uncertainty
and assessed `met` stayed open for ever -- and re-claiming is the very
route the re-wake offers, so the repair could not discharge its own
condition. Across streams the opposite: the workspace record dropped
the assessment entirely, so a later cycle reading an earlier claim
found none, and a missing assessment read as no open requirement.

One idea fixes both: the *current* assessment of each requirement,
latest wins, with earlier ones kept as evidence.
"""

from __future__ import annotations

import pytest

from chemsmart.agent.delivery import (
    current_assessments,
    unresolved_requirement_ids,
)
from chemsmart.agent.workspace_record import (
    record_run,
    render_workspace_record,
)

pytestmark = pytest.mark.capability("rule:wake.restate_observable")


def _row(observable_id, state, uncertainty=None):
    return {
        "observable_id": observable_id,
        "unit": "kJ/mol",
        "required_tolerance": 2.0,
        "uncertainty": uncertainty,
        "uncertainty_basis": "asserted",
        "state": state,
    }


def test_a_later_assessment_closes_an_earlier_shortfall():
    rows = (
        _row("dg", "unstated"),
        _row("other", "met", 1.0),
        _row("dg", "met", 1.0),
    )
    assert unresolved_requirement_ids(rows) == ()
    assert current_assessments(rows)["dg"]["state"] == "met"


def test_a_later_assessment_can_also_re_open_one():
    """Latest wins in both directions; nothing is sticky either way."""

    rows = (_row("dg", "met", 1.0), _row("dg", "short", 6.0))
    assert unresolved_requirement_ids(rows) == ("dg",)


def test_an_assessment_survives_into_the_workspace_record(tmp_path):
    """The cross-cycle half: the record carries the assessment beside
    the number it judges, so a later cycle can see it at all."""

    from chemsmart.agent._contracts import canonical_sha256
    from chemsmart.agent.runtime.event_store import RuntimeEventStore

    run = tmp_path / "run"
    run.mkdir()
    record = {
        "task_spec_sha256": "a" * 64,
        "status": "recorded",
        "claims": [
            {
                "claim_id": "dg",
                "quantity_id": "dg",
                "display_value": -10.0,
                "display_unit": "kJ/mol",
                "dimension": [1, 0, 0, 0, 0, 0],
                "source_receipt_sha256": "4" * 64,
            }
        ],
    }
    store = RuntimeEventStore(run / "events.jsonl", session_id="s1")
    store.append(
        turn_id="turn-1",
        kind="analysis_claims_recorded",
        payload={
            "status": "recorded",
            "task_spec_sha256": "a" * 64,
            "receipt_sha256": canonical_sha256(record),
            "record": record,
            "source_receipt_sha256s": ["4" * 64],
            "claim_ids": ["dg"],
            "critical_finding_count": 0,
            "sufficiency": [_row("dg", "short", 6.0)],
        },
        idempotency_key="analysis-claims:one",
    )

    workspace = tmp_path / "ws"
    workspace.mkdir()
    record_run(
        workspace,
        goal_id="goal-t1",
        cycle=1,
        run="goals/goal-t1/runs/cycle-1",
        run_events_path=run / "events.jsonl",
    )
    rendered = render_workspace_record(workspace)
    claim = next(row for row in rendered["claims"] if row["claim_id"] == "dg")
    assert claim["sufficiency"]["state"] == "short"
    assert claim["sufficiency"]["required_tolerance"] == 2.0
