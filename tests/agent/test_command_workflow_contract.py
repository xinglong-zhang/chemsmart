import pytest

from chemsmart.agent._contracts import ContractError
from chemsmart.agent.workflows import (
    ArtifactBindingV1,
    ArtifactOutputV1,
    CommandNodeV1,
    build_command_workflow_spec,
)


def _node(*, node_id, dependencies=(), inputs=(), outputs=()):
    return CommandNodeV1(
        node_id=node_id,
        program="pyscf",
        jobtype="sp",
        capability_receipt_sha256="1" * 64,
        program_binding_sha256="2" * 64,
        engine_binding_sha256="3" * 64,
        scientific_identity_sha256="4" * 64,
        project_artifact_id="",
        project_sha256="",
        dependencies=dependencies,
        inputs=inputs,
        expected_outputs=outputs,
        invocation_sha256="",
        preflight_receipt_sha256="",
        state="planned",
    )


def test_workflow_requires_exact_topological_producer_output_binding():
    first = _node(
        node_id="node.one",
        inputs=(
            ArtifactBindingV1(
                "input.geometry",
                "geometry.initial",
                "geometry_xyz",
                "a" * 64,
                "",
                "",
                "",
            ),
        ),
        outputs=(
            ArtifactOutputV1(
                "optimized_geometry", "geometry.optimized", "geometry_xyz"
            ),
        ),
    )
    second = _node(
        node_id="node.two",
        dependencies=("node.one",),
        inputs=(
            ArtifactBindingV1(
                "input.optimized",
                "geometry.optimized",
                "geometry_xyz",
                "",
                "node.one",
                "optimized_geometry",
                "",
            ),
        ),
        outputs=(
            ArtifactOutputV1("energy", "result.energy", "result_record"),
        ),
    )
    workflow = build_command_workflow_spec(
        workflow_id="workflow.one",
        task_spec_sha256="b" * 64,
        live_cli_schema_sha256="c" * 64,
        nodes=(first, second),
    )

    assert workflow.nodes[1].inputs[0].resolved is False

    wrong = _node(
        node_id="node.bad",
        dependencies=("node.one",),
        inputs=(
            ArtifactBindingV1(
                "input.wrong",
                "geometry.other",
                "geometry_xyz",
                "",
                "node.one",
                "optimized_geometry",
                "",
            ),
        ),
        outputs=(ArtifactOutputV1("energy", "result.other", "result_record"),),
    )
    with pytest.raises(ContractError, match="exact producer output"):
        build_command_workflow_spec(
            workflow_id="workflow.bad",
            task_spec_sha256="b" * 64,
            live_cli_schema_sha256="c" * 64,
            nodes=(first, wrong),
        )
