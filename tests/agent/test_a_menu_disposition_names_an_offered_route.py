"""What a session did with the repair menu is a typed field, verified.

REACH-1 po3 cycle 3 (2026-09-06) rejected all four routes the wake's
repair menu offered, each with a mechanism, in prose the next cycle
never saw: the wake carries the menu and not what was done with it, so
cycle 4 inherited the list and not the argument. The owner ruled (R3,
2026-09-06) the disposition is a typed field of the decision record.
The host verifies only what it can know -- the route is one the menu
offered this cycle, the receipts are ones it minted -- and the next
wake shows the dispositions beside the menu it re-offers.
"""

from __future__ import annotations

import pytest

from chemsmart.agent.driver import (
    _recorded_dispositions,
    _session_dispositions,
    _wake_context,
)
from chemsmart.agent.goal import GoalLedger
from chemsmart.agent.tool_runtime import RoutedContractError

from .test_a_guide_opens_when_something_asks import _host
from .test_an_excursion_is_a_grant_line_never_the_deliverable import _goal

_BASE = {
    "decision_id": "d1",
    "assumptions": ["a"],
    "method_rationale": "r",
    "alternatives": ["b"],
    "uncertainties": ["u"],
    "diagnostics": ["g"],
    "stage_order": ["s"],
    "evidence_refs": [],
}


def _decide(host, dispositions, decision_id="d1"):
    return host.dispatch(
        turn_id=f"t-{decision_id}",
        tool_name="record_scientific_decision",
        arguments={
            **_BASE,
            "decision_id": decision_id,
            "menu_route_dispositions": dispositions,
        },
    )


@pytest.mark.capability("rule:wake.menu_dispositions_are_recorded")
def test_a_disposition_names_a_route_the_menu_offered(tmp_path):
    host = _host(
        tmp_path, offered_repair_routes=("failed_wrong_stationary_point",)
    )
    with pytest.raises(RoutedContractError) as refused:
        _decide(
            host,
            [
                {
                    "route": "scf_nonconverged",
                    "disposition": "rejected",
                    "reason": "no SCF failed",
                }
            ],
        )
    report = refused.value.failure_report
    assert report["gate"] == "decision.route_is_one_the_menu_offered"
    assert "failed_wrong_stationary_point" in report["diagnosis"]

    accepted = _decide(
        host,
        [
            {
                "route": "failed_wrong_stationary_point",
                "disposition": "rejected",
                "reason": "the saddle is on the C4 channel; a mode step "
                "along the printed mode stays there",
            }
        ],
    )
    assert accepted["status"] == "ok"
    (recorded,) = _session_dispositions(tmp_path / "events.jsonl")
    assert recorded["route"] == "failed_wrong_stationary_point"
    assert recorded["disposition"] == "rejected"
    assert recorded["receipt_sha256s"] == []


@pytest.mark.capability("rule:wake.menu_dispositions_are_recorded")
def test_a_wake_with_no_menu_takes_no_disposition(tmp_path):
    host = _host(tmp_path)
    with pytest.raises(RoutedContractError) as refused:
        _decide(
            host,
            [{"route": "timeout", "disposition": "taken", "reason": "x"}],
        )
    assert (
        "carried no repair menu" in refused.value.failure_report["diagnosis"]
    )


@pytest.mark.capability("rule:wake.menu_dispositions_are_recorded")
def test_the_next_wake_shows_what_the_previous_cycle_decided(tmp_path):
    ledger = GoalLedger(tmp_path / "goal")
    ledger.create(_goal())
    ledger.append(
        "repair_menu_dispositions",
        {
            "cycle": 2,
            "dispositions": [
                {
                    "route": "failed_wrong_stationary_point",
                    "disposition": "rejected",
                    "reason": "wrong channel",
                    "receipt_sha256s": [],
                }
            ],
        },
    )
    ledger.append(
        "repair_menu_dispositions",
        {
            "cycle": 3,
            "dispositions": [
                {
                    "route": "failed_wrong_stationary_point",
                    "disposition": "taken",
                    "reason": "seeded from the scan maximum",
                    "receipt_sha256s": [],
                }
            ],
        },
    )
    (latest,) = _recorded_dispositions(ledger)
    assert latest["cycle"] == 3
    assert latest["disposition"] == "taken"
    wake = _wake_context(ledger.load(), ledger, None)
    assert wake["repair_menu_dispositions"] == (latest,)
