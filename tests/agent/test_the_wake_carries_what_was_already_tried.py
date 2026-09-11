"""A goal's next cycle sees what its earlier cycles already attempted.

REACH-1 po3 (2026-09-06) diagnosed at cycle 3 that its dual-contact
seed had drifted to a separated pair, rejected four repair routes with
mechanisms, and then at cycle 4 seeded a structurally identical dual
contact. The wake carried budgets, anomalies and a repair menu; it did
not carry what had been tried. Nothing new is asked of the model here:
the register is the stream's own words, kept where the next cycle
reads them.
"""

from __future__ import annotations

import json

import pytest

from chemsmart.agent.driver import (
    _recorded_approaches,
    _run_approaches,
    _session_approaches,
    _wake_context,
)
from chemsmart.agent.goal import GoalLedger

from .test_an_excursion_is_a_grant_line_never_the_deliverable import _goal


@pytest.mark.capability("rule:wake.approaches_already_tried")
def test_a_rejected_alternative_is_kept_in_the_sessions_own_words(tmp_path):
    events = tmp_path / "events.jsonl"
    events.write_text(
        json.dumps(
            {
                "kind": "scientific_decision_recorded",
                "payload": {
                    "record": {
                        "alternatives": [
                            "Restart from the reached geometry: rejected "
                            "because the structure is a separated pair",
                            "Loosen the convergence criterion: rejected "
                            "because it buys the word, not the structure",
                        ]
                    }
                },
            }
        )
        + "\n",
        encoding="utf-8",
    )
    approaches = _session_approaches(events)
    assert len(approaches) == 2
    assert approaches[0]["outcome"] == "rejected_by_the_session"
    assert "separated pair" in approaches[0]["approach"]


@pytest.mark.capability("rule:wake.approaches_already_tried")
def test_a_failed_node_is_kept_with_the_numbers_that_typed_it(tmp_path):
    run = tmp_path / "run.jsonl"
    run.write_text(
        "\n".join(
            json.dumps(item)
            for item in (
                {
                    "kind": "program_result_verified",
                    "payload": {"node_id": "ts-c4", "status": "invalid"},
                },
                {
                    "kind": "program_result_verified",
                    "payload": {"node_id": "scan-c5", "status": "valid"},
                },
                {
                    "kind": "anomaly_observed",
                    "payload": {
                        "node_id": "ts-c4",
                        "record": {
                            "signal_id": "stationary_point.unexpected_order",
                            "values": {
                                "expected_imaginary_modes": 1,
                                "observed_imaginary_modes": 0,
                            },
                        },
                    },
                },
            )
        )
        + "\n",
        encoding="utf-8",
    )
    (approach,) = _run_approaches(run)
    assert approach["approach"] == "ts-c4"
    assert approach["outcome"] == "invalid"
    assert "observed_imaginary_modes=0" in approach["mechanism"]


@pytest.mark.capability("rule:wake.approaches_already_tried")
def test_the_next_wake_shows_them_most_recent_first(tmp_path):
    ledger = GoalLedger(tmp_path / "goal")
    ledger.create(_goal())
    for cycle, node in ((1, "ts-c4b"), (2, "ts-c4-r4")):
        ledger.append(
            "approaches_recorded",
            {
                "cycle": cycle,
                "approaches": [
                    {
                        "approach": node,
                        "outcome": "invalid",
                        "mechanism": "the fragments drifted apart",
                    }
                ],
            },
        )
    recorded = _recorded_approaches(ledger)
    assert [item["approach"] for item in recorded] == ["ts-c4-r4", "ts-c4b"]
    assert recorded[0]["cycle"] == 2

    wake = _wake_context(ledger.load(), ledger, None)
    assert wake["approaches_tried"] == recorded
    assert "approaches_tried is what this goal has already attempted" in (
        wake["authority"]
    )
