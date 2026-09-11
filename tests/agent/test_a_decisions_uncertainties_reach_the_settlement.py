"""What the recorded decision said it was unsure of is carried where the
human reads first.

A woken session wrote that the author's 2 kJ/mol design threshold sits
inside the method's 1-2 kJ/mol error; it reached no receipt, no claim
and no settlement (NOVEL-2 po2, 2026-09-04). The decision record has
always carried uncertainties; now the settlement does too.
"""

from __future__ import annotations

import json

import pytest

from chemsmart.agent.driver import (
    _achieved_word,
    _analysis_delivery,
    _settlement_evidence,
)

pytestmark = pytest.mark.capability("tool:record_scientific_decision")


def test_uncertainties_ride_the_evidence_and_the_word(tmp_path):
    rows = [
        {
            "kind": "analysis_claims_recorded",
            "payload": {
                "receipt_sha256": "1" * 64,
                "record": {
                    "claims": [
                        {
                            "claim_id": "dg",
                            "quantity_id": "dg",
                            "source_receipt_sha256": "2" * 64,
                        }
                    ]
                },
            },
        },
        {
            "kind": "scientific_decision_recorded",
            "payload": {
                "receipt_sha256": "3" * 64,
                "record": {
                    "decision_id": "d1",
                    "uncertainties": [
                        "method error 1-2 kJ/mol, the size of the "
                        "2 kJ/mol design threshold"
                    ],
                    "evidence_refs": [],
                },
            },
        },
    ]
    stream = tmp_path / "events.jsonl"
    stream.write_text(
        "\n".join(json.dumps(row) for row in rows) + "\n", encoding="utf-8"
    )
    delivery = _analysis_delivery(stream)
    assert delivery.decision_uncertainties == (
        "method error 1-2 kJ/mol, the size of the 2 kJ/mol design threshold",
    )
    evidence = _settlement_evidence(delivery)
    assert (
        evidence["decision_uncertainties"] == delivery.decision_uncertainties
    )
    _word, reasons = _achieved_word(delivery, ())
    assert any("states its uncertainties" in reason for reason in reasons)
