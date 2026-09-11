"""The workspace record is what a later cycle and a human actually read.

Two defects earned these pins. The projection's idempotency guard lived
on a driver instance, so any second process re-appended a stream it had
already recorded. And the record kept a bare
``uncertainty_evidence_backed: true`` while dropping the citation the
host had resolved to grant it -- so the one word that closes a precision
contract could not be audited by the reader it exists for.
"""

from __future__ import annotations

import json

import pytest

from chemsmart.agent.workspace_record import (
    read_workspace_record,
    record_run,
)

pytestmark = pytest.mark.capability("tool:record_analysis_claims")


def test_the_record_is_appended_once_per_stream(tmp_path):
    """The idempotency guard was a set on a driver instance.

    A second process -- a restart, a wake, any new driver -- started
    empty and appended the same stream's rows again. Every row already
    carries what identifies it, so the guard belongs where the write
    is.
    """

    workspace = tmp_path / "ws"
    events = tmp_path / "events.jsonl"
    events.write_text(
        json.dumps(
            {
                "kind": "analysis_claims_recorded",
                "payload": {
                    "receipt_sha256": "3" * 64,
                    "record": {
                        "claims": [
                            {
                                "claim_id": "dg",
                                "quantity_id": "dg",
                                "display_value": 1.0,
                                "display_unit": "kJ/mol",
                                "dimension": [1, 0, 0, 0, 0, 0],
                                "source_receipt_sha256": "4" * 64,
                            }
                        ]
                    },
                },
            }
        )
        + "\n",
        encoding="utf-8",
    )
    first = record_run(
        workspace, goal_id="g", cycle=1, run_events_path=events, run="r"
    )
    second = record_run(
        workspace, goal_id="g", cycle=1, run_events_path=events, run="r"
    )
    assert (first, second) == (1, 0)
    rows = [
        row
        for row in read_workspace_record(workspace)
        if row.get("kind") == "claim"
    ]
    assert len(rows) == 1


def test_the_assessment_carries_what_the_host_resolved(tmp_path):
    """`met` is auditable, not protected.

    No host check can decide whether a stated magnitude is really the
    uncertainty of a claim -- that is chemistry, and the host does not
    judge it. What the host owes is that its own word can be audited,
    and the record kept a bare `uncertainty_evidence_backed: true` while
    dropping the citation the host had resolved to grant it.
    """

    workspace = tmp_path / "ws"
    events = tmp_path / "events.jsonl"
    events.write_text(
        json.dumps(
            {
                "kind": "analysis_claims_recorded",
                "payload": {
                    "receipt_sha256": "3" * 64,
                    "record": {
                        "claims": [
                            {
                                "claim_id": "dg",
                                "quantity_id": "dg",
                                "display_value": 1.0,
                                "display_unit": "kJ/mol",
                                "dimension": [1, 0, 0, 0, 0, 0],
                                "source_receipt_sha256": "4" * 64,
                                "uncertainty": 0.4,
                                "uncertainty_basis": "measured",
                                "uncertainty_reference": "5" * 64 + ":spread",
                                "uncertainty_components": [
                                    {
                                        "meaning": "the treatment spread",
                                        "magnitude": 0.4,
                                        "basis": "measured",
                                        "reference": "5" * 64,
                                    }
                                ],
                            }
                        ]
                    },
                    "sufficiency": [
                        {
                            "observable_id": "dg",
                            "unit": "kJ/mol",
                            "required_tolerance": 2.0,
                            "uncertainty": 0.4,
                            "uncertainty_basis": "measured",
                            "uncertainty_evidence_backed": True,
                            "meets_tolerance": True,
                            "state": "met",
                        }
                    ],
                },
            }
        )
        + "\n",
        encoding="utf-8",
    )
    record_run(
        workspace, goal_id="g", cycle=1, run_events_path=events, run="r"
    )
    (row,) = [
        item
        for item in read_workspace_record(workspace)
        if item.get("kind") == "claim"
    ]
    assessment = row["sufficiency"]
    assert assessment["state"] == "met"
    assert assessment["uncertainty_reference"].endswith(":spread")
    assert len(assessment["uncertainty_components"]) == 1
