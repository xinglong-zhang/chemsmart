from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from chemsmart.agent._contracts import (
    ContractError,
    TrustedArtifactRefV1,
    canonical_data,
    canonical_sha256,
    file_sha256,
)
from chemsmart.agent.commands import build_scientific_identity_binding
from chemsmart.agent.execution import (
    ApprovedNodeBindingV1,
    build_execution_resource_spec,
    build_producer_edge_rule,
    build_program_execution_invocation,
    build_workflow_execution_approval,
    handoff_optimized_xtb_geometry,
)
from chemsmart.agent.runtime.event_store import RuntimeEventStore
from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1
from chemsmart.agent.tool_runtime import _CommandContext
from chemsmart.agent.workflows import (
    ArtifactInputIntentV1,
    ArtifactOutputIntentV1,
    CommandNodeIntentV1,
    ScientificWorkflowEdgeV2,
    ScientificWorkflowNodeV2,
    build_command_workflow_draft,
    build_scientific_workflow_plan,
)
from tests.agent.neutral_workflow_fixture import (
    build_neutral_workflow_fixture,
)


def _artifact(
    path: Path, *, artifact_id: str, kind: str
) -> TrustedArtifactRefV1:
    return TrustedArtifactRefV1(
        artifact_id=artifact_id,
        kind=kind,
        sha256=file_sha256(path),
        size_bytes=path.stat().st_size,
        path=str(path.resolve()),
        cli_value=str(path.resolve()),
    )


def _resources():
    return build_execution_resource_spec(
        execution_target="run",
        cores=4,
        memory_gb=8,
        gpu_count=0,
        scratch_policy="server",
        node_timeout_seconds=600,
    )


@pytest.mark.parametrize(
    ("producer_state", "consumer_state"),
    (
        ((0, 1), (-1, 2)),
        ((-1, 2), (0, 1)),
    ),
    ids=("neutral-geometry-to-anion", "anion-geometry-to-neutral"),
)
def test_cross_state_handoff_preserves_geometry_and_admits_vertical_energy(
    tmp_path, producer_state, consumer_state
):
    initial_path = tmp_path / "pbq-initial.xyz"
    initial_path.write_text(
        "4\nstate-specific PBQ fragment\n"
        "C -1.0 0.0 0.0\nO -2.2 0.0 0.0\n"
        "C 1.0 0.0 0.0\nO 2.2 0.0 0.0\n",
        encoding="utf-8",
    )
    optimized_path = tmp_path / "xtbopt.xyz"
    optimized_path.write_text(
        "4\noptimized exact frame\n"
        "C -1.01 0.02 0.0\nO -2.24 -0.01 0.0\n"
        "C 1.01 -0.02 0.0\nO 2.24 0.01 0.0\n",
        encoding="utf-8",
    )
    project_path = tmp_path / "vertical.yaml"
    project_path.write_text("sp: {}\n", encoding="utf-8")
    initial = _artifact(
        initial_path,
        artifact_id="geometry.pbq.initial",
        kind="geometry_xyz",
    )
    optimized = _artifact(
        optimized_path,
        artifact_id="result.pbq.optimized",
        kind="geometry_xyz",
    )
    project = _artifact(
        project_path,
        artifact_id="project.pbq.vertical",
        kind="project_yaml",
    )
    edge = build_producer_edge_rule(
        producer_node_id="opt-source-state",
        consumer_node_id="sp-vertical-target-state",
        artifact_kind="geometry_xyz",
        selection_rule="validated_optimized_geometry",
    )
    receipt = SimpleNamespace(
        validated=True,
        node_id="opt-source-state",
        receipt_sha256="a" * 64,
        output_artifacts=(optimized,),
    )
    geometry, handoff = handoff_optimized_xtb_geometry(
        producer_receipt=receipt,
        result_artifact=optimized,
        input_artifact=initial,
        producer_edge=edge,
        approved_workspace=tmp_path,
        geometry_artifact_id="geometry.pbq.for-vertical-energy",
        expected_charge=producer_state[0],
        expected_multiplicity=producer_state[1],
        consumer_charge=consumer_state[0],
        consumer_multiplicity=consumer_state[1],
    )

    assert (handoff.charge, handoff.multiplicity) == producer_state
    assert handoff.consumer_state == consumer_state
    assert handoff.symbols == ("C", "O", "C", "O")
    assert geometry.sha256 == file_sha256(geometry.path)
    assert (
        "consumer_charge"
        not in Path(geometry.path).read_text(encoding="utf-8").splitlines()[1]
    )
    observed_rows = [
        line.split()
        for line in Path(geometry.path)
        .read_text(encoding="utf-8")
        .splitlines()[2:]
    ]
    source_rows = [
        line.split()
        for line in optimized_path.read_text(encoding="utf-8").splitlines()[2:]
    ]
    assert [row[0] for row in observed_rows] == [row[0] for row in source_rows]
    assert [
        tuple(float(value) for value in row[1:]) for row in observed_rows
    ] == [tuple(float(value) for value in row[1:]) for row in source_rows]

    target_identity = build_scientific_identity_binding(
        task_spec_sha256="b" * 64,
        geometry_artifact=geometry,
        charge=consumer_state[0],
        multiplicity=consumer_state[1],
    )
    approved_consumer = ApprovedNodeBindingV1(
        node_id="sp-vertical-target-state",
        program="xtb",
        engine="cpu",
        jobtype="sp",
        project_artifact_sha256=project.sha256,
        settings_sha256="c" * 64,
        charge=consumer_state[0],
        multiplicity=consumer_state[1],
        input_mode="producer",
        initial_artifact_id="",
        initial_artifact_sha256="",
        scientific_identity_sha256="",
        producer_edge_sha256=edge.edge_sha256,
    )
    approval = build_workflow_execution_approval(
        approval_id="pbq-four-point",
        workflow_id="pbq-marcus",
        workflow_sha256="d" * 64,
        task_spec_sha256="b" * 64,
        approved_workspace=tmp_path,
        resources=_resources(),
        node_bindings=(
            ApprovedNodeBindingV1(
                node_id="opt-source-state",
                program="xtb",
                engine="cpu",
                jobtype="opt",
                project_artifact_sha256=project.sha256,
                settings_sha256="e" * 64,
                charge=producer_state[0],
                multiplicity=producer_state[1],
                input_mode="initial",
                initial_artifact_id=initial.artifact_id,
                initial_artifact_sha256=initial.sha256,
                scientific_identity_sha256="f" * 64,
                producer_edge_sha256="",
            ),
            approved_consumer,
        ),
        producer_edges=(edge,),
    )
    invocation = build_program_execution_invocation(
        node_id=approved_consumer.node_id,
        approval=approval,
        project_artifact=project,
        input_artifact=geometry,
        scientific_identity_sha256=target_identity.binding_sha256,
        environment_receipt_sha256="1" * 64,
        resources=_resources(),
        argv=("chemsmart", "run", "xtb", "sp"),
        handoff=handoff,
    )

    assert invocation.node_id == "sp-vertical-target-state"
    assert invocation.input_sha256 == geometry.sha256
    assert (
        invocation.scientific_identity_sha256 == target_identity.binding_sha256
    )
    source_state_identity = build_scientific_identity_binding(
        task_spec_sha256="b" * 64,
        geometry_artifact=geometry,
        charge=producer_state[0],
        multiplicity=producer_state[1],
    )
    with pytest.raises(ContractError, match="target state and geometry"):
        build_program_execution_invocation(
            node_id=approved_consumer.node_id,
            approval=approval,
            project_artifact=project,
            input_artifact=geometry,
            scientific_identity_sha256=source_state_identity.binding_sha256,
            environment_receipt_sha256="1" * 64,
            resources=_resources(),
            argv=("chemsmart", "run", "xtb", "sp"),
            handoff=handoff,
        )


@pytest.mark.parametrize(
    "contract",
    (CommandNodeIntentV1, ScientificWorkflowNodeV2),
)
def test_node_state_requires_complete_valid_charge_multiplicity_pair(contract):
    common = (
        {
            "node_id": "vertical-energy",
            "program": "pyscf",
            "jobtype": "sp",
            "project_role": "project.vertical",
            "dependencies": (),
            "inputs": (),
            "expected_outputs": (
                ArtifactOutputIntentV1("energy", "program_result"),
            ),
            "unresolved_fields": (),
        }
        if contract is CommandNodeIntentV1
        else {
            "node_id": "vertical-energy",
            "stage": "sp",
            "requested_program": "pyscf",
            "program": "pyscf",
            "engine": "cpu",
            "project_role": "project.vertical",
            "unresolved_fields": (),
        }
    )
    with pytest.raises(
        ContractError, match="charge and multiplicity together"
    ):
        contract(**common, charge=-1)
    with pytest.raises(ContractError, match="multiplicity must be positive"):
        contract(**common, charge=-1, multiplicity=0)


def test_absent_state_preserves_legacy_command_and_scientific_plan_digests():
    command_node = CommandNodeIntentV1(
        node_id="neutral-sp",
        program="pyscf",
        jobtype="sp",
        project_role="project.neutral",
        dependencies=(),
        inputs=(
            ArtifactInputIntentV1(
                binding_id="geometry",
                artifact_class="geometry_xyz",
                artifact_id="geometry.neutral",
                producer_node_id="",
                producer_output_id="",
            ),
        ),
        expected_outputs=(ArtifactOutputIntentV1("result", "program_result"),),
        unresolved_fields=(),
    )
    draft = build_command_workflow_draft(
        workflow_id="legacy-neutral",
        task_spec_id="2" * 64,
        nodes=(command_node,),
    )
    legacy_command_node = canonical_data(command_node)
    for field in ("node_kind", "charge", "multiplicity"):
        legacy_command_node.pop(field)
    assert draft.draft_sha256 == canonical_sha256(
        {
            "schema_version": "chemsmart.command-workflow-draft.v1",
            "workflow_id": "legacy-neutral",
            "task_spec_id": "2" * 64,
            "nodes": (legacy_command_node,),
            "status": "planned",
        }
    )

    scientific_node = ScientificWorkflowNodeV2(
        node_id="neutral-sp",
        stage="sp",
        requested_program="pyscf",
        program="pyscf",
        engine="cpu",
        project_role="project.neutral",
        unresolved_fields=(),
    )
    plan = build_scientific_workflow_plan(
        workflow_id="legacy-neutral",
        task_spec_sha256="2" * 64,
        scientific_identity_sha256="3" * 64,
        nodes=(scientific_node,),
    )
    legacy_scientific_node = canonical_data(scientific_node)
    for field in ("charge", "multiplicity"):
        legacy_scientific_node.pop(field)
    assert plan.plan_sha256 == canonical_sha256(
        {
            "schema_version": "chemsmart.scientific-workflow-plan.v2",
            "workflow_id": "legacy-neutral",
            "task_spec_sha256": "2" * 64,
            "scientific_identity_sha256": "3" * 64,
            "nodes": (legacy_scientific_node,),
            "edges": (),
            "complexity_factors": (),
            "status": "planned",
        }
    )


def test_model_workflow_state_reaches_command_and_scientific_dags(tmp_path):
    fixture = build_neutral_workflow_fixture(tmp_path / "fixture")
    store = RuntimeEventStore(
        tmp_path / "events" / "runtime.jsonl", session_id="session"
    )
    host = CommandCompiledToolHostV1(
        event_store=store,
        task_spec_sha256s=(fixture.public_context.task_spec_sha256,),
        approved_workspace=tmp_path / "workspace",
        **fixture.host_inputs,
    )
    action = next(
        item
        for item in fixture.public_context.next_actions
        if item.tool_name == "plan_command_workflow"
    )
    fields = dict(action.fields)
    raw_nodes = []
    for raw_node in fields["nodes"]:
        state = (
            {"charge": -1, "multiplicity": 2}
            if raw_node["node_id"] == "node.hess"
            else {}
        )
        raw_nodes.append({**raw_node, **state})

    host.dispatch(
        turn_id="turn-1",
        tool_name="plan_command_workflow",
        arguments={**fields, "nodes": raw_nodes},
    )

    draft = next(iter(host.workflow_drafts.values()))
    plan = next(iter(host.scientific_workflow_plans.values()))
    draft_node = next(
        node for node in draft.nodes if node.node_id == "node.hess"
    )
    scientific_node = next(
        node for node in plan.nodes if node.node_id == "node.hess"
    )
    assert (draft_node.charge, draft_node.multiplicity) == (-1, 2)
    assert (scientific_node.charge, scientific_node.multiplicity) == (-1, 2)


def test_future_bounded_node_context_uses_explicit_target_state(tmp_path):
    geometry_path = tmp_path / "pbq-neutral-frame.xyz"
    geometry_path.write_text(
        "2\nPBQ fragment\nC 0 0 0\nO 1.2 0 0\n", encoding="utf-8"
    )
    geometry = _artifact(
        geometry_path,
        artifact_id="geometry.pbq.neutral-start",
        kind="geometry_xyz",
    )
    producer_identity = build_scientific_identity_binding(
        task_spec_sha256="a" * 64,
        geometry_artifact=geometry,
        charge=0,
        multiplicity=1,
    )
    plan = build_scientific_workflow_plan(
        workflow_id="pbq-marcus",
        task_spec_sha256="a" * 64,
        scientific_identity_sha256=producer_identity.binding_sha256,
        nodes=(
            ScientificWorkflowNodeV2(
                node_id="opt-neutral",
                stage="opt",
                requested_program="pyscf",
                program="pyscf",
                engine="cpu",
                project_role="project.opt",
                unresolved_fields=(),
            ),
            ScientificWorkflowNodeV2(
                node_id="sp-anion-at-neutral-geometry",
                stage="sp",
                requested_program="pyscf",
                program="pyscf",
                engine="cpu",
                project_role="project.sp",
                unresolved_fields=(),
                support_state="unresolved_future",
                charge=-1,
                multiplicity=2,
            ),
        ),
        edges=(
            ScientificWorkflowEdgeV2(
                edge_id="opt-neutral-to-anion-sp",
                source_node_id="opt-neutral",
                target_node_id="sp-anion-at-neutral-geometry",
                edge_kind="data",
                artifact_class="geometry_xyz",
                producer_output_id="optimized-geometry",
                consumer_input_id="filename",
            ),
        ),
    )
    project = SimpleNamespace(artifact_id="project.sp", sha256="1" * 64)
    capability = SimpleNamespace(
        receipt_sha256="2" * 64,
        status=SimpleNamespace(value="supported"),
        query=SimpleNamespace(program="pyscf", jobtype="sp", engine="cpu"),
    )
    validation = SimpleNamespace(
        project_artifact_id=project.artifact_id,
        project_sha256=project.sha256,
        capability_receipt_sha256=capability.receipt_sha256,
        program="pyscf",
        jobtype="sp",
        status="valid",
    )
    engine = SimpleNamespace(
        program="pyscf",
        selected_engine="cpu",
        state="resolved",
        execution_ready=True,
        capability_receipt_sha256=capability.receipt_sha256,
        program_binding_sha256="3" * 64,
        environment_receipt_sha256="4" * 64,
    )
    producer_context = _CommandContext(
        proposal=SimpleNamespace(),
        capability=SimpleNamespace(),
        program_binding=SimpleNamespace(),
        engine_binding=engine,
        project_artifact=SimpleNamespace(),
        project_validation=SimpleNamespace(),
        input_artifact=geometry,
        scientific_identity=producer_identity,
    )
    host = object.__new__(CommandCompiledToolHostV1)
    host.workflow_drafts = {
        "draft": SimpleNamespace(
            workflow_id=plan.workflow_id,
            nodes=(
                SimpleNamespace(
                    node_id="sp-anion-at-neutral-geometry",
                    project_role=project.artifact_id,
                    inputs=(
                        SimpleNamespace(
                            binding_id="filename",
                            artifact_class="geometry_xyz",
                            producer_node_id="opt-neutral",
                            producer_output_id="optimized-geometry",
                        ),
                    ),
                ),
            ),
        )
    }
    host.artifacts = {project.artifact_id: project}
    host.project_validations = {"validation": validation}
    host.engine_bindings = {"engine": engine}
    host.capabilities = {capability.receipt_sha256: capability}
    host.program_bindings = {engine.program_binding_sha256: SimpleNamespace()}
    host._latest_invocation_for_node = lambda _node_id, **_kwargs: (
        SimpleNamespace(),
        producer_context,
    )

    # The future context still carries the producer geometry until execution;
    # admission reads the target state from the planned node and creates the
    # consumer identity only after the validated optimized frame exists.
    context = host._bounded_node_context(
        plan=plan,
        planned_node=plan.nodes[1],
        data_target_ids={"sp-anion-at-neutral-geometry"},
    )
    assert context.input_artifact is geometry
    assert (
        context.scientific_identity.charge,
        context.scientific_identity.multiplicity,
    ) == (0, 1)
    assert (plan.nodes[1].charge, plan.nodes[1].multiplicity) == (-1, 2)
