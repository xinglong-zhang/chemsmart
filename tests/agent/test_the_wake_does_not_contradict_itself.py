"""The wake's deliverables and its rules use the key the gate uses.

NOVEL-2 po2's cycle-2 wake listed six ids as delivered (built from
quantity_id) and the same six as undelivered (built from claim_id), and
its rule text said the gate matched "by kind and unit, never value" while
the gate matched by id. One record, one key.
"""

from __future__ import annotations

import json

import pytest

from chemsmart.agent.driver import _analysis_delivery
from chemsmart.agent.rules import rules_by_id

pytestmark = pytest.mark.capability("rule:wake.restate_observable")


def test_the_restate_rule_says_what_the_gate_does():
    text = rules_by_id()["wake.restate_observable"].text
    assert "by id" in text
    assert "claim_id" in text
    assert "kind and unit" not in text


def test_a_claim_is_delivered_under_both_ids_it_carries(tmp_path):
    stream = tmp_path / "events.jsonl"
    stream.write_text(
        json.dumps(
            {
                "kind": "analysis_claims_recorded",
                "payload": {
                    "receipt_sha256": "1" * 64,
                    "record": {
                        "claims": [
                            {
                                "claim_id": "dg_gauche_anti_sulfone_mecn",
                                "quantity_id": "expr-out",
                                "source_receipt_sha256": "2" * 64,
                            }
                        ]
                    },
                },
            }
        )
        + "\n",
        encoding="utf-8",
    )
    delivery = _analysis_delivery(stream)
    assert set(delivery.delivered_quantity_ids) == {
        "dg_gauche_anti_sulfone_mecn",
        "expr-out",
    }
