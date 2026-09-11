"""A declared sign that the declared band already implies is recorded as
implied and never printed as agreed.

NOVEL-2 ino1 declared expected_sign positive on a gap defined as the
next state minus the ground state, band 0..60: that sign cannot fail,
so it predicted nothing (2026-09-04). Sibling of the shipped rule that
a sign test on zero can never agree.
"""

from __future__ import annotations

from types import SimpleNamespace

import pytest

from .test_a_guide_opens_when_something_asks import _host

pytestmark = pytest.mark.capability("tool:declare_requested_observable")


def test_an_implied_sign_is_marked_and_the_band_decides(tmp_path):
    host = _host(tmp_path)
    result = host.dispatch(
        turn_id="t1",
        tool_name="declare_requested_observable",
        arguments={
            "observables": [
                {
                    "observable_id": "trans_gap_kjmol",
                    "unit": "kJ/mol",
                    "meaning": "next state minus ground state",
                    "expectation_basis": "near-crossover d6",
                    "expected_sign": "positive",
                    "expected_low": 0.0,
                    "expected_high": 60.0,
                }
            ]
        },
    )["result"]
    (declared,) = result["declared"]
    assert declared["sign_implied_by_band"] is True

    claim = SimpleNamespace(
        claim_id="trans_gap_kjmol",
        quantity_id="trans_gap_kjmol",
        dimension=(1, 0, 0, 0, 0, 0),
        display_value=46.2,
        display_unit="kJ/mol",
    )
    host.analysis_claim_records["r1"] = SimpleNamespace(
        task_spec_sha256=next(iter(host.task_spec_sha256s)), claims=(claim,)
    )
    (row,) = host._declared_observable_predictions(
        task_spec_sha256=next(iter(host.task_spec_sha256s))
    )
    assert row["sign_implied_by_band"] is True
    assert row["agreement"] == "agreed"


def test_a_sign_the_band_does_not_imply_is_a_prediction(tmp_path):
    host = _host(tmp_path)
    result = host.dispatch(
        turn_id="t1",
        tool_name="declare_requested_observable",
        arguments={
            "observables": [
                {
                    "observable_id": "dg",
                    "unit": "kJ/mol",
                    "meaning": "gauche minus anti",
                    "expectation_basis": "the gauche effect",
                    "expected_sign": "negative",
                    "expected_low": -4.0,
                    "expected_high": 4.0,
                }
            ]
        },
    )["result"]
    (declared,) = result["declared"]
    assert "sign_implied_by_band" not in declared
