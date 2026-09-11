"""A diagnostic is the session's own prediction, scored and never owed.

REACH-1 po3 (2026-09-06) diagnosed, from anomaly receipts, that its
saddle search would land on the wrong channel, rejected every route the
repair menu offered with a mechanism, and had nowhere to commit that
prediction: a declaration was a deliverable, so declaring the diagnosis
would have held the goal open on it. The owner ruled (R3, 2026-09-06)
that the agent's own predictions get standing: declared with
``role: diagnostic``, joined by id, scored like any expectation, and
never a limitation -- an undelivered diagnostic settles nothing, a
diverged one is an observation the settlement word carries, and one
inside the method's own resolution prints indeterminate.
"""

from __future__ import annotations

from types import SimpleNamespace

import pytest

from chemsmart.agent.driver import _required_declared_ids
from chemsmart.agent.goal import GoalLedger
from chemsmart.agent.tool_runtime import ContractError

from .test_a_guide_opens_when_something_asks import _host
from .test_an_excursion_is_a_grant_line_never_the_deliverable import _goal

_TASK = "a" * 64

_REQUESTED = {
    "observable_id": "dg_c4_minus_c5",
    "unit": "kJ/mol",
    "dimension": (1, 0, 0, 0, 0, 0),
    "meaning": "G(C4 ester TS) - G(C5 ester TS)",
}

_DIAGNOSTIC = {
    "observable_id": "ts_c5_barrier_minus_c4",
    "unit": "kJ/mol",
    "dimension": (1, 0, 0, 0, 0, 0),
    "meaning": "the C5-side saddle lies above the C4-side saddle",
    "role": "diagnostic",
    "expectation_basis": "the anomaly receipts show the C4 channel",
    "expected_sign": "positive",
    "failure_update_rule": "seed the C5 saddle from the scan maximum",
}


def _claim(observable_id, value):
    return SimpleNamespace(
        claim_id=observable_id,
        quantity_id=observable_id,
        dimension=(1, 0, 0, 0, 0, 0),
        display_value=value,
        display_unit="kJ/mol",
    )


def _seed(host, *claims):
    host.analysis_claim_records["r1"] = SimpleNamespace(
        task_spec_sha256=_TASK, claims=tuple(claims)
    )


@pytest.mark.capability("rule:declare.diagnostic_has_standing")
def test_an_undelivered_diagnostic_is_no_limitation(tmp_path):
    host = _host(
        tmp_path,
        approved_requested_observable_declarations=[_REQUESTED, _DIAGNOSTIC],
    )
    _seed(host, _claim("dg_c4_minus_c5", 6.3))
    assert host._declared_observable_completion(task_spec_sha256=_TASK) == (
        (),
        (),
    )
    (row,) = host._declared_observable_predictions(task_spec_sha256=_TASK)
    assert row["observable_id"] == "ts_c5_barrier_minus_c4"
    assert row["role"] == "diagnostic"
    assert row["agreement"] == "not_comparable"
    assert row["failure_update_rule"] == _DIAGNOSTIC["failure_update_rule"]


@pytest.mark.capability("rule:declare.diagnostic_has_standing")
def test_a_delivered_diagnostic_is_scored_and_a_diverged_one_stands(
    tmp_path,
):
    host = _host(
        tmp_path,
        approved_requested_observable_declarations=[_REQUESTED, _DIAGNOSTIC],
    )
    _seed(
        host,
        _claim("dg_c4_minus_c5", 6.3),
        _claim("ts_c5_barrier_minus_c4", -9.5),
    )
    (row,) = host._declared_observable_predictions(task_spec_sha256=_TASK)
    assert row["agreement"] == "diverged"
    assert host._declared_observable_completion(task_spec_sha256=_TASK) == (
        (),
        (),
    )


@pytest.mark.capability("rule:declare.diagnostic_has_standing")
def test_a_value_inside_the_methods_resolution_is_indeterminate(tmp_path):
    resolved = {**_DIAGNOSTIC, "method_resolution": 4.0}
    host = _host(
        tmp_path,
        approved_requested_observable_declarations=[resolved],
    )
    _seed(host, _claim("ts_c5_barrier_minus_c4", -2.1))
    (row,) = host._declared_observable_predictions(task_spec_sha256=_TASK)
    assert row["agreement"] == "indeterminate"
    assert row["within_method_resolution"] is True
    assert row["method_resolution"] == 4.0


@pytest.mark.capability("rule:declare.diagnostic_has_standing")
def test_a_diagnostic_needs_a_prediction_and_an_update_rule(tmp_path):
    host = _host(tmp_path)

    def declare(item):
        return host.dispatch(
            turn_id="t1",
            tool_name="declare_requested_observable",
            arguments={"observables": [item]},
        )

    bare = {
        "observable_id": "which_channel",
        "unit": "kJ/mol",
        "meaning": "C5 minus C4 saddle",
        "role": "diagnostic",
    }
    with pytest.raises(ContractError, match="predicts nothing"):
        declare(bare)
    predicted = {
        **bare,
        "expected_sign": "positive",
        "expectation_basis": "the receipts show the C4 channel",
    }
    with pytest.raises(ContractError, match="failure_update_rule"):
        declare(predicted)
    result = declare(
        {**predicted, "failure_update_rule": "seed from the scan maximum"}
    )["result"]
    (declared,) = result["declared"]
    assert declared["role"] == "diagnostic"
    assert declared["failure_update_rule"] == "seed from the scan maximum"


@pytest.mark.capability("rule:declare.diagnostic_has_standing")
def test_a_diagnostic_never_holds_the_goal_open(tmp_path):
    ledger = GoalLedger(tmp_path / "goal")
    ledger.create(_goal())
    ledger.append(
        "observables_declared",
        {"cycle": 1, "observables": [_REQUESTED, _DIAGNOSTIC]},
    )
    assert _required_declared_ids(ledger) == ("dg_c4_minus_c5",)
