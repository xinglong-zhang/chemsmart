"""A declared observable is answered by the id the claim carries.

Two live losses earned these pins. A chain rendered six correct numbers
whose receipts carried the declared id in ``quantity_id`` and the plan's
short input label in ``claim_id``; the completion gate read one field and
called all six undelivered (NOVEL-2 po2, 2026-09-04). The plan behind it
had declared the six ids as its claim node's outputs and bound them under
six other input ids, and the plan-time gate accepted that: an output no
input carries is a claim nothing renders.
"""

from __future__ import annotations

from types import SimpleNamespace

import pytest

from chemsmart.agent.scientific_toolchain import (
    AnalysisInputIntentV1,
    AnalysisNodeIntentV1,
    AnalysisOutputIntentV1,
    ScientificToolchainContractError,
)


def _host(tmp_path, declarations):
    from tests.agent.test_a_guide_opens_when_something_asks import _host

    return _host(
        tmp_path, approved_requested_observable_declarations=declarations
    )


_DG = {
    "observable_id": "dg_gauche_anti_sulfone_mecn",
    "unit": "kJ/mol",
    "dimension": (1, 0, 0, 0, 0, 0),
    "meaning": "G(gauche) - G(anti) in acetonitrile",
    "expectation_basis": "opposing hyperconjugation and dipole terms",
    "expected_low": -4.0,
    "expected_high": 4.0,
}


@pytest.mark.capability("rule:plan.claim_carries_declared_id")
def test_a_claim_standing_on_the_declared_quantity_answers_it(tmp_path):
    host = _host(tmp_path, [_DG])
    claim = SimpleNamespace(
        claim_id="dg_sulfone_mecn",
        quantity_id="dg_gauche_anti_sulfone_mecn",
        dimension=(1, 0, 0, 0, 0, 0),
        display_value=-0.906,
        display_unit="kJ/mol",
    )
    host.analysis_claim_records["r1"] = SimpleNamespace(
        task_spec_sha256="a" * 64, claims=(claim,)
    )
    assert host._declared_observable_completion(task_spec_sha256="a" * 64) == (
        (),
        (),
    )
    assert host._declared_observable_join_fields == {
        "dg_gauche_anti_sulfone_mecn": "quantity_id"
    }
    (row,) = host._declared_observable_predictions(task_spec_sha256="a" * 64)
    assert row["delivered_claim_id"] == "dg_sulfone_mecn"
    assert row["agreement"] == "agreed"


@pytest.mark.capability("rule:plan.claim_carries_declared_id")
def test_the_claim_id_still_joins_first_and_a_miss_names_both_fields(
    tmp_path,
):
    host = _host(tmp_path, [_DG])
    named = SimpleNamespace(
        claim_id="dg_gauche_anti_sulfone_mecn",
        quantity_id="expr-out",
        dimension=(1, 0, 0, 0, 0, 0),
        display_value=-0.906,
        display_unit="kJ/mol",
    )
    host.analysis_claim_records["r1"] = SimpleNamespace(
        task_spec_sha256="a" * 64, claims=(named,)
    )
    assert host._declared_observable_completion(task_spec_sha256="a" * 64) == (
        (),
        (),
    )
    assert host._declared_observable_join_fields == {
        "dg_gauche_anti_sulfone_mecn": "claim_id"
    }
    host.analysis_claim_records.clear()
    host.analysis_claim_records["r2"] = SimpleNamespace(
        task_spec_sha256="a" * 64,
        claims=(
            SimpleNamespace(
                claim_id="other",
                quantity_id="also-other",
                dimension=(1, 0, 0, 0, 0, 0),
                display_value=1.0,
                display_unit="kJ/mol",
            ),
        ),
    )
    misses, limitations = host._declared_observable_completion(
        task_spec_sha256="a" * 64
    )
    assert limitations == ("declared_observable:dg_gauche_anti_sulfone_mecn",)
    assert "as its claim_id" in misses[0]
    assert "receipt quantity" in misses[0]


def _claim_node(input_ids, output_ids):
    return AnalysisNodeIntentV1(
        node_id="claim_dg",
        analysis_kind="claim_rendering",
        dependencies=("expr_dg",),
        selectors=(),
        expression_nodes=(),
        expression_output_node_ids=(),
        temperature_k=None,
        pressure_atm=None,
        support_state="planned",
        blocked_reason="",
        validation_rules=(),
        inputs=tuple(
            AnalysisInputIntentV1(
                input_id=input_id,
                source_kind="analysis_output",
                producer_node_id="expr_dg",
                producer_output_id=output_id,
            )
            for input_id, output_id in zip(input_ids, output_ids)
        ),
        outputs=tuple(
            AnalysisOutputIntentV1(
                output_id=output_id, quantity_kind="energy", unit="kJ/mol"
            )
            for output_id in output_ids
        ),
    )


@pytest.mark.capability("rule:plan.claim_carries_declared_id")
def test_a_claim_node_output_no_input_carries_is_refused_when_planned():
    with pytest.raises(ScientificToolchainContractError) as refused:
        _claim_node(("dg_sulfone_mecn",), ("dg_gauche_anti_sulfone_mecn",))
    message = str(refused.value)
    assert "dg_gauche_anti_sulfone_mecn" in message
    assert "rendered under its input_id" in message
    assert "observable_id verbatim" in message


@pytest.mark.capability("rule:plan.claim_carries_declared_id")
def test_a_claim_node_whose_outputs_are_its_inputs_plans():
    node = _claim_node(
        ("dg_gauche_anti_sulfone_mecn",), ("dg_gauche_anti_sulfone_mecn",)
    )
    assert node.outputs[0].output_id == node.inputs[0].input_id
