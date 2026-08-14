"""Typed analysis must consume a registered result, not a native handoff."""

from dataclasses import replace

import pytest

from chemsmart.agent.scientific_toolchain import (
    AnalysisInputIntentV1,
    AnalysisNodeIntentV1,
    AnalysisOutputIntentV1,
    RegisteredResultInputIntentV1,
    ScientificToolchainContractError,
    build_scientific_toolchain_plan,
)
from chemsmart.agent.workflows import (
    ArtifactInputIntentV1,
    ArtifactOutputIntentV1,
    CommandNodeIntentV1,
)


def _calculation(artifact_class: str) -> CommandNodeIntentV1:
    return CommandNodeIntentV1(
        node_id="orca-opt",
        program="orca",
        jobtype="opt",
        project_role="minimum",
        dependencies=(),
        inputs=(
            ArtifactInputIntentV1(
                binding_id="filename",
                artifact_class="geometry_xyz",
                artifact_id="geometry-start",
                producer_node_id="",
                producer_output_id="",
            ),
        ),
        expected_outputs=(
            ArtifactOutputIntentV1(
                output_id="analysis-result",
                artifact_class=artifact_class,
            ),
        ),
        unresolved_fields=(),
    )


def _thermochemistry() -> AnalysisNodeIntentV1:
    return AnalysisNodeIntentV1(
        node_id="thermochemistry",
        analysis_kind="thermochemistry",
        dependencies=("orca-opt",),
        inputs=(
            AnalysisInputIntentV1(
                input_id="result",
                source_kind="program_output",
                producer_node_id="orca-opt",
                producer_output_id="analysis-result",
            ),
        ),
        selectors=(),
        outputs=(
            AnalysisOutputIntentV1(
                output_id="gibbs",
                quantity_kind="gibbs_energy",
                unit="hartree",
            ),
        ),
        expression_nodes=(),
        expression_output_node_ids=(),
        temperature_k=298.15,
        pressure_atm=1.0,
        support_state="planned",
        blocked_reason="",
    )


def _plan(artifact_class: str):
    return build_scientific_toolchain_plan(
        plan_id="result-artifact-probe",
        workflow_id="result-artifact-probe",
        command_workflow_draft_sha256="a" * 64,
        calculation_nodes=(_calculation(artifact_class),),
        calculation_observables={"orca-opt": ("vibrational_frequencies",)},
        analysis_nodes=(_thermochemistry(),),
        required_output_ids=("gibbs",),
    )


def test_thermochemistry_requires_the_registered_program_result_artifact():
    with pytest.raises(
        ScientificToolchainContractError,
        match=r"orca_hessian.*orca_output",
    ):
        _plan("orca_hessian")

    plan = _plan("orca_output")
    assert plan.analysis_nodes[0].inputs[0].producer_output_id == (
        "analysis-result"
    )


def test_registered_result_can_be_the_thermochemistry_root():
    thermochemistry = replace(
        _thermochemistry(),
        dependencies=(),
        inputs=(
            RegisteredResultInputIntentV1(
                input_id="result",
                artifact_id="registered-frequency-result",
            ),
        ),
    )

    plan = build_scientific_toolchain_plan(
        plan_id="registered-thermochemistry",
        workflow_id="registered-thermochemistry",
        command_workflow_draft_sha256="b" * 64,
        calculation_nodes=(),
        calculation_observables={},
        analysis_nodes=(thermochemistry,),
        required_output_ids=("gibbs",),
    )

    assert plan.analysis_nodes[0].inputs == thermochemistry.inputs


def test_analysis_quantity_cannot_replace_thermochemistry_result_artifact():
    quantity_source = AnalysisNodeIntentV1(
        node_id="derived-quantity",
        analysis_kind="claim_rendering",
        dependencies=(),
        inputs=(),
        selectors=(),
        outputs=(
            AnalysisOutputIntentV1(
                output_id="derived-energy",
                quantity_kind="energy",
                unit="hartree",
            ),
        ),
        expression_nodes=(),
        expression_output_node_ids=(),
        temperature_k=None,
        pressure_atm=None,
        support_state="planned",
        blocked_reason="",
    )
    thermochemistry = replace(
        _thermochemistry(),
        dependencies=("derived-quantity",),
        inputs=(
            AnalysisInputIntentV1(
                input_id="result",
                source_kind="analysis_output",
                producer_node_id="derived-quantity",
                producer_output_id="derived-energy",
            ),
        ),
    )

    with pytest.raises(
        ScientificToolchainContractError,
        match=r"analysis_output.*cannot replace",
    ):
        build_scientific_toolchain_plan(
            plan_id="derived-thermochemistry",
            workflow_id="derived-thermochemistry",
            command_workflow_draft_sha256="c" * 64,
            calculation_nodes=(),
            calculation_observables={},
            analysis_nodes=(quantity_source, thermochemistry),
            required_output_ids=("gibbs",),
        )
