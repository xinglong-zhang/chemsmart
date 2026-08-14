from __future__ import annotations

import pytest

from chemsmart.agent._contracts import canonical_data
from chemsmart.agent.scientific_toolchain import (
    AnalysisInputIntentV1,
    AnalysisNodeIntentV1,
    AnalysisOutputIntentV1,
    RegisteredResultInputIntentV1,
    ScientificToolchainContractError,
)
from chemsmart.agent.tool_specs import build_command_compiled_tool_surface
from chemsmart.analysis.result_quantities import (
    ThermochemistryReceiptV1,
    canonical_quantity_sha256,
    thermochemistry_receipt_from_record,
)


def _thermochemistry_node(**overrides):
    values = {
        "node_id": "derive-free-energy",
        "analysis_kind": "thermochemistry",
        "dependencies": ("frequency-calculation",),
        "inputs": (
            AnalysisInputIntentV1(
                input_id="frequency-result",
                source_kind="program_output",
                producer_node_id="frequency-calculation",
                producer_output_id="structured-result",
            ),
        ),
        "selectors": (),
        "outputs": (
            AnalysisOutputIntentV1(
                output_id="quasi-harmonic-free-energy",
                quantity_kind="gibbs_free_energy",
                unit="hartree",
            ),
        ),
        "expression_nodes": (),
        "expression_output_node_ids": (),
        "temperature_k": 298.15,
        "pressure_atm": 1.0,
        "support_state": "planned",
        "blocked_reason": "",
        "concentration_mol_l": 1.0,
        "entropy_method": "grimme",
        "entropy_cutoff_cm1": 100.0,
        "enthalpy_cutoff_cm1": 100.0,
        "alpha": 4,
        "use_weighted_mass": False,
        "frequency_scale_factor": 0.985,
    }
    values.update(overrides)
    return AnalysisNodeIntentV1(**values)


def test_quasi_harmonic_controls_are_hashable_typed_dag_state():
    record = canonical_data(_thermochemistry_node())

    assert record["entropy_method"] == "grimme"
    assert record["entropy_cutoff_cm1"] == 100.0
    assert record["enthalpy_cutoff_cm1"] == 100.0
    assert record["concentration_mol_l"] == 1.0
    assert record["alpha"] == 4
    assert record["frequency_scale_factor"] == 0.985


def test_quasi_harmonic_entropy_requires_its_typed_cutoff():
    with pytest.raises(
        ScientificToolchainContractError,
        match="grimme entropy requires entropy_cutoff_cm1",
    ):
        _thermochemistry_node(entropy_cutoff_cm1=None)


def test_planned_thermochemistry_requires_exactly_one_result_input():
    with pytest.raises(
        ScientificToolchainContractError,
        match="exactly one result input",
    ):
        _thermochemistry_node(inputs=(), dependencies=())

    registered = RegisteredResultInputIntentV1(
        input_id="registered-result",
        artifact_id="frequency-result",
    )
    node = _thermochemistry_node(inputs=(registered,), dependencies=())
    assert node.inputs == (registered,)

    with pytest.raises(
        ScientificToolchainContractError,
        match="exactly one result input",
    ):
        _thermochemistry_node(
            inputs=(
                registered,
                RegisteredResultInputIntentV1(
                    input_id="second-result",
                    artifact_id="second-frequency-result",
                ),
            ),
            dependencies=(),
        )


def test_blocked_thermochemistry_can_preserve_intent_without_a_result():
    node = _thermochemistry_node(
        dependencies=(),
        inputs=(),
        support_state="blocked_unsupported",
        blocked_reason="required result reader is unavailable",
    )

    assert node.inputs == ()
    assert node.support_state == "blocked_unsupported"


def test_plan_tool_schema_exposes_thermochemistry_controls():
    tool = next(
        item
        for item in build_command_compiled_tool_surface().tool_definitions
        if item["function"]["name"] == "plan_scientific_workflow"
    )
    properties = tool["function"]["parameters"]["properties"]
    analysis = properties["analysis_nodes"]["items"]["properties"]

    for name in (
        "concentration_mol_l",
        "entropy_method",
        "entropy_cutoff_cm1",
        "enthalpy_cutoff_cm1",
        "alpha",
        "use_weighted_mass",
        "frequency_scale_factor",
    ):
        assert name in analysis

    assert "thermochemistry root" in analysis["artifact_id"]["description"]
    source_kind = analysis["inputs"]["items"]["properties"]["source_kind"]
    assert "immediate producer" in source_kind["description"]


def test_extended_thermochemistry_receipt_replays_with_all_controls():
    body = {
        "schema_version": "chemsmart.thermochemistry-receipt.v1",
        "artifact_id": "frequency-result",
        "artifact_sha256": "b" * 64,
        "program": "orca",
        "engine_id": "orca-6.1.1",
        "temperature_k": 298.15,
        "pressure_atm": 1.0,
        "quantities": (),
        "assumptions": ("typed quasi-harmonic controls",),
        "status": "derived",
        "concentration_mol_l": 1.0,
        "entropy_method": "grimme",
        "entropy_cutoff_cm1": 100.0,
        "enthalpy_cutoff_cm1": 100.0,
        "alpha": 4,
        "use_weighted_mass": False,
        "frequency_scale_factor": 0.985,
    }
    digest = canonical_quantity_sha256(body)
    receipt = ThermochemistryReceiptV1(**body, receipt_sha256=digest)

    replayed = thermochemistry_receipt_from_record(body, receipt_sha256=digest)

    assert replayed == receipt
    assert replayed.frequency_scale_factor == 0.985
    assert replayed.concentration_mol_l == 1.0


def test_legacy_pyscf_rrho_receipt_digest_remains_replayable():
    legacy_body = {
        "schema_version": "chemsmart.thermochemistry-receipt.v1",
        "artifact_id": "legacy-frequency-result",
        "artifact_sha256": "c" * 64,
        "program": "pyscf",
        "engine_id": "pyscf-cpu",
        "temperature_k": 298.15,
        "pressure_atm": 1.0,
        "quantities": (),
        "assumptions": ("legacy harmonic RRHO",),
        "status": "derived",
    }
    legacy_digest = canonical_quantity_sha256(legacy_body)

    replayed = thermochemistry_receipt_from_record(
        legacy_body,
        receipt_sha256=legacy_digest,
    )

    assert replayed.receipt_sha256 == legacy_digest
    assert replayed.entropy_method == "rrho"
    assert replayed.frequency_scale_factor == 1.0
    assert replayed.concentration_mol_l is None
