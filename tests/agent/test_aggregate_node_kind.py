"""A DAG must be able to declare the stage that carries the answer.

When every node is "invoke a program", an observable that needs host
arithmetic over finished results -- a complete-basis-set limit, a Boltzmann
average, a state difference -- has no expressible producer. Repairing such a
plan then removes the finding *and* the answer.
"""

import pytest

from chemsmart.agent._contracts import ContractError
from chemsmart.agent.workflows import (
    AGGREGATE_NODE_PROGRAM,
    AGGREGATE_NODE_STAGE,
    ArtifactInputIntentV1,
    ArtifactOutputIntentV1,
    CommandNodeIntentV1,
    build_command_workflow_draft,
)


def _energy_input(binding_id="qz", producer="sp_ccpvqz", output="e_qz"):
    return ArtifactInputIntentV1(
        binding_id=binding_id,
        artifact_class="energy",
        artifact_id="",
        producer_node_id=producer,
        producer_output_id=output,
    )


def _single_point():
    return CommandNodeIntentV1(
        node_id="sp_ccpvqz",
        program="pyscf",
        jobtype="sp",
        project_role="sp",
        dependencies=(),
        inputs=(
            ArtifactInputIntentV1(
                binding_id="geom",
                artifact_class="geometry_xyz",
                artifact_id="start.xyz",
                producer_node_id="",
                producer_output_id="",
            ),
        ),
        expected_outputs=(
            ArtifactOutputIntentV1(output_id="e_qz", artifact_class="energy"),
        ),
        unresolved_fields=(),
    )


def _aggregate(**overrides):
    values = dict(
        node_id="cbs_limit",
        program=AGGREGATE_NODE_PROGRAM,
        jobtype=AGGREGATE_NODE_STAGE,
        project_role="aggregate",
        dependencies=("sp_ccpvqz",),
        inputs=(_energy_input(),),
        expected_outputs=(
            ArtifactOutputIntentV1(output_id="e_cbs", artifact_class="energy"),
        ),
        unresolved_fields=(),
        node_kind="aggregate",
    )
    values.update(overrides)
    return CommandNodeIntentV1(**values)


def test_an_aggregate_stage_plans_alongside_program_calls():
    draft = build_command_workflow_draft(
        workflow_id="wf",
        task_spec_id="t",
        nodes=(_single_point(), _aggregate()),
    )
    assert draft.nodes[1].node_kind == "aggregate"


def test_a_program_call_workflow_keeps_its_existing_identity():
    """A defaulted field must not change the digest of recorded workflows."""

    node = _single_point()
    assert node.node_kind == "program_call"
    first = build_command_workflow_draft(
        workflow_id="wf", task_spec_id="t", nodes=(node,)
    )
    second = build_command_workflow_draft(
        workflow_id="wf", task_spec_id="t", nodes=(node,)
    )
    assert first.draft_sha256 == second.draft_sha256


def test_declaring_an_aggregate_stage_changes_the_workflow_identity():
    plain = build_command_workflow_draft(
        workflow_id="wf", task_spec_id="t", nodes=(_single_point(),)
    )
    with_stage = build_command_workflow_draft(
        workflow_id="wf",
        task_spec_id="t",
        nodes=(_single_point(), _aggregate()),
    )
    assert plain.draft_sha256 != with_stage.draft_sha256


def test_an_aggregate_node_may_not_name_a_program():
    with pytest.raises(ContractError, match="runs no program"):
        _aggregate(program="orca")


def test_the_arithmetic_belongs_to_the_expression_not_the_node():
    """One vocabulary: operations live in the quantity expression."""

    with pytest.raises(ContractError, match="quantity_expression"):
        _aggregate(jobtype="cbs_extrapolation")


def test_an_aggregate_node_must_consume_the_results_it_combines():
    with pytest.raises(ContractError, match="consumes nothing"):
        _aggregate(inputs=(), dependencies=())


def test_a_program_call_cannot_masquerade_as_host_arithmetic():
    with pytest.raises(ContractError, match="not an executable program"):
        _aggregate(node_kind="program_call", jobtype="sp")


def test_an_unknown_node_kind_is_refused_by_name():
    with pytest.raises(ContractError, match="node_kind"):
        _aggregate(node_kind="freestyle")


def test_the_aggregate_program_is_absent_from_the_executable_registry():
    """capabilities.py declares executable programs; aggregation is not one."""

    from chemsmart.settings.capabilities import EXECUTABLE_PROGRAMS

    assert AGGREGATE_NODE_PROGRAM not in EXECUTABLE_PROGRAMS


def test_the_command_tool_schema_keeps_aggregate_in_the_analysis_plane():
    from chemsmart.agent.tool_specs import build_command_compiled_tool_surface

    surface = build_command_compiled_tool_surface()
    plan = next(
        item
        for item in surface.tool_definitions
        if item["function"]["name"] == "plan_command_workflow"
    )
    node = plan["function"]["parameters"]["properties"]["nodes"]["items"]
    assert "aggregate" not in node["properties"]["node_kind"]["enum"]
    assert node["properties"]["node_kind"]["enum"] == ["program_call"]
