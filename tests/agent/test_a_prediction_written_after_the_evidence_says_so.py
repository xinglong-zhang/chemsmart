"""A prediction written after the numbers exist is marked as such.

OPEN-1 ino3 (2026-09-07) extracted its spin populations, then declared
six observables with bands drawn around the values it had just read,
and the completion printed twelve expectation rows saying `agreed` with
nothing to separate them from the two written before any number
existed. The verdict does not move -- being right after the fact is
still being right, and a diverged post-hoc row is still a result -- but
a reader must not mistake a restatement for a pre-registration.
"""

from __future__ import annotations

from types import SimpleNamespace

import pytest

from .test_a_guide_opens_when_something_asks import _host


def _declare(host, observable_id, turn):
    return host.dispatch(
        turn_id=turn,
        tool_name="declare_requested_observable",
        arguments={
            "observables": [
                {
                    "observable_id": observable_id,
                    "unit": "kcal/mol",
                    "meaning": "a barrier difference",
                    "expected_sign": "positive",
                    "expectation_basis": "the ester controls orientation",
                }
            ]
        },
    )["result"]


@pytest.mark.capability("rule:declare.diagnostic_has_standing")
def test_a_declaration_before_any_number_is_a_pre_registration(tmp_path):
    host = _host(tmp_path)
    (record,) = _declare(host, "before", "t1")["declared"]
    assert "declared_after_evidence" not in record


@pytest.mark.capability("rule:declare.diagnostic_has_standing")
def test_a_declaration_after_a_receipt_says_so_on_its_row(tmp_path):
    host = _host(tmp_path)
    _declare(host, "before", "t1")
    host.quantity_extractions["r1"] = SimpleNamespace(receipt_sha256="a" * 64)
    (record,) = _declare(host, "after", "t2")["declared"]
    assert record["declared_after_evidence"] is True

    task = next(iter(host.task_spec_sha256s))
    host.analysis_claim_records["c1"] = SimpleNamespace(
        task_spec_sha256=task,
        claims=(
            SimpleNamespace(
                claim_id="after",
                quantity_id="after",
                dimension=(1, 0, 0, 0, 0, 0),
                display_value=6.3,
                display_unit="kcal/mol",
            ),
        ),
    )
    rows = {
        row["observable_id"]: row
        for row in host._declared_observable_predictions(task_spec_sha256=task)
    }
    assert rows["after"]["agreement"] == "agreed"
    assert rows["after"]["declared_after_evidence"] is True
    assert "declared_after_evidence" not in rows["before"]
