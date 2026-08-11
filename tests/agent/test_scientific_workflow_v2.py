from __future__ import annotations

import pytest

from chemsmart.agent._contracts import (
    ContractError,
    canonical_data,
    canonical_sha256,
)
from chemsmart.agent.execution import (
    ProgramResultValidationReceiptV1,
    build_execution_resource_spec,
    build_frozen_workflow_approval,
    derive_ready_node_ids,
)
from chemsmart.agent.tool_runtime import (
    _validate_stationary_point_policy_binding,
)
from chemsmart.agent.runtime.event_store import RuntimeEventStore
from chemsmart.agent.runtime.events import EventKind
from chemsmart.agent.runtime.reducer import RuntimeState
from chemsmart.agent.workflows import (
    ArtifactOutputIntentV1,
    CommandNodeIntentV1,
    ContextManifestV1,
    ScientificWorkflowEdgeV2,
    ScientificWorkflowNodeV2,
    SpecialistTaskPacketV1,
    StationaryPointValidationPolicyV1,
    MaterializedNodeV1,
    build_materialized_workflow,
    build_scientific_workflow_plan,
    eligible_specialist_roles,
)


def _node(node_id: str, stage: str):
    return ScientificWorkflowNodeV2(
        node_id=node_id,
        stage=stage,
        requested_program="pyscf",
        program="pyscf",
        engine="cpu",
        project_role="water-project",
        unresolved_fields=(),
    )


def test_command_node_rejects_invalid_unresolved_field_with_local_counterexample():
    with pytest.raises(
        ContractError,
        match=(
            "node 'candidate-hess' has invalid unresolved field "
            "'geometry \(from OPT\)'"
        ),
    ):
        CommandNodeIntentV1(
            node_id="candidate-hess",
            program="pyscf",
            jobtype="hess",
            project_role="candidate.pyscf.hess",
            dependencies=(),
            inputs=(),
            expected_outputs=(
                ArtifactOutputIntentV1(
                    output_id="hessian",
                    artifact_class="hdf5",
                ),
            ),
            unresolved_fields=("geometry (from OPT)",),
        )


def _water_plan():
    return build_scientific_workflow_plan(
        workflow_id="water-workflow",
        task_spec_sha256="a" * 64,
        scientific_identity_sha256="b" * 64,
        nodes=(
            _node("sp-initial", "sp"),
            _node("opt-initial", "opt"),
            _node("hess-optimized", "hess"),
        ),
        edges=(
            ScientificWorkflowEdgeV2(
                edge_id="opt-to-hess",
                source_node_id="opt-initial",
                target_node_id="hess-optimized",
                edge_kind="data",
                artifact_class="geometry_xyz",
                producer_output_id="optimized-geometry",
                consumer_input_id="geometry",
            ),
        ),
    )


def _sp_result_validation() -> ProgramResultValidationReceiptV1:
    body = {
        "schema_version": "chemsmart.program-result-validation-receipt.v1",
        "validator_id": "pyscf-result-validator",
        "validator_schema_version": "chemsmart.pyscf-result-validation.v1",
        "validator_version": "1",
        "invocation_sha256": "7" * 64,
        "node_id": "sp-initial",
        "program": "pyscf",
        "engine": "cpu",
        "jobtype": "sp",
        "input_artifact_sha256": "1" * 64,
        "project_artifact_sha256": "2" * 64,
        "capability_environment_receipt_sha256": "3" * 64,
        "run_environment_receipt_sha256": "",
        "environment_validation_sha256": "",
        "stationary_point_policy_sha256": "",
        "output_artifacts": (),
        "observations": {"state": "validated"},
        "findings": (),
        "state": "valid",
    }
    return ProgramResultValidationReceiptV1(
        **body, receipt_sha256=canonical_sha256(body)
    )


def _materialized(plan, *, resource_sha256="6" * 64):
    materialized = tuple(
        MaterializedNodeV1(
            node_id=node_id,
            input_artifact_sha256="c" * 64,
            project_artifact_sha256="d" * 64,
            project_validation_receipt_sha256="e" * 64,
            environment_receipt_sha256="f" * 64,
            invocation_sha256=invocation,
            preflight_receipt_sha256=preflight,
            state=state,
        )
        for node_id, invocation, preflight, state in (
            ("sp-initial", "1" * 64, "2" * 64, "previewed"),
            ("opt-initial", "3" * 64, "4" * 64, "previewed"),
        )
    )
    return build_materialized_workflow(
        plan=plan,
        live_cli_schema_sha256="5" * 64,
        resource_sha256=resource_sha256,
        nodes=materialized,
        unresolved_node_ids=("hess-optimized",),
        status="partial",
    )


def _approval(plan):
    resources = build_execution_resource_spec(
        execution_target="run",
        cores=4,
        memory_gb=4,
        gpu_count=0,
        scratch_policy="none",
        node_timeout_seconds=600,
    )
    materialized = _materialized(
        plan, resource_sha256=resources.resource_sha256
    )
    return resources, materialized, build_frozen_workflow_approval(
        approval_id="water-approval",
        plan=plan,
        materialized_workflow=materialized,
        resources=resources,
        environment_identity_sha256s=("f" * 64,),
    )


def test_scientific_dag_keeps_sp_and_opt_parallel_but_hess_data_bound():
    plan = _water_plan()

    assert plan.complexity_factors == (
        "multiple_stages",
        "producer_artifact_edge",
    )
    assert eligible_specialist_roles(plan) == ("dag_specialist",)
    resources, partial, approval = _approval(plan)
    assert resources.cores == 4
    materialized = build_materialized_workflow(
        plan=plan,
        live_cli_schema_sha256=partial.live_cli_schema_sha256,
        resource_sha256=resources.resource_sha256,
        nodes=(
            *partial.nodes,
            MaterializedNodeV1(
                node_id="hess-optimized",
                input_artifact_sha256="7" * 64,
                project_artifact_sha256="d" * 64,
                project_validation_receipt_sha256="e" * 64,
                environment_receipt_sha256="f" * 64,
                invocation_sha256="8" * 64,
                preflight_receipt_sha256="",
                state="compiled",
            ),
        ),
        unresolved_node_ids=(),
        status="materialized",
    )
    assert materialized.status == "materialized"
    from chemsmart.agent.execution import build_workflow_run_state

    run = build_workflow_run_state(
        run_id="water-run",
        plan=plan,
        approval=approval,
        approval_consumed=True,
    )
    assert derive_ready_node_ids(plan, run) == (
        "sp-initial",
        "opt-initial",
    )


def test_context_manifest_and_specialist_packet_bind_narrow_tools():
    body = {
        "schema_version": "chemsmart.context-manifest.v1",
        "manifest_id": "context-water",
        "workflow_id": "water-workflow",
        "task_spec_sha256": "a" * 64,
        "scientific_identity_sha256": "b" * 64,
        "source_sha256s": ("c" * 64,),
        "artifact_sha256s": ("d" * 64,),
        "tool_schema_sha256": "e" * 64,
        "allowed_tools": ("inspect_program_capability",),
        "token_budget": 4096,
        "tool_call_budget": 4,
        "wall_time_seconds": 120,
    }
    manifest = ContextManifestV1(
        **body, manifest_sha256=canonical_sha256(body)
    )
    packet_body = {
        "schema_version": "chemsmart.specialist-task-packet.v1",
        "packet_id": "packet-water",
        "workflow_id": "water-workflow",
        "role": "pyscf-specialist",
        "context_manifest_sha256": manifest.manifest_sha256,
        "input_record_sha256s": ("f" * 64,),
        "expected_output_schema": "typed-settings-v1",
        "owner": "coordinator",
        "merge_key": "water-settings",
    }
    packet = SpecialistTaskPacketV1(
        **packet_body, packet_sha256=canonical_sha256(packet_body)
    )

    assert packet.owner == "coordinator"
    assert manifest.allowed_tools == ("inspect_program_capability",)


def test_stationary_point_policy_freezes_expected_mode_count():
    body = {
        "schema_version": "chemsmart.stationary-point-policy.v1",
        "policy_id": "water-minimum",
        "task_spec_sha256": "a" * 64,
        "hessian_node_id": "hess-optimized",
        "stationary_point_kind": "minimum",
        "expected_imaginary_mode_count": 1,
        "imaginary_mode_cutoff_cm1": 20.0,
        "require_finite_modes": True,
        "require_symmetric_hessian": True,
    }
    with pytest.raises(ContractError, match="point kind"):
        StationaryPointValidationPolicyV1(
            **body, policy_sha256=canonical_sha256(body)
        )


def test_stationary_point_policy_binds_exact_plan_task_and_hessian_node():
    plan = _water_plan()
    resources = build_execution_resource_spec(
        execution_target="run",
        cores=4,
        memory_gb=4,
        gpu_count=0,
        scratch_policy="none",
        node_timeout_seconds=600,
    )
    materialized = _materialized(
        plan, resource_sha256=resources.resource_sha256
    )
    policy_body = {
        "schema_version": "chemsmart.stationary-point-policy.v1",
        "policy_id": "water-minimum",
        "task_spec_sha256": plan.task_spec_sha256,
        "hessian_node_id": "hess-optimized",
        "stationary_point_kind": "minimum",
        "expected_imaginary_mode_count": 0,
        "imaginary_mode_cutoff_cm1": 20.0,
        "require_finite_modes": True,
        "require_symmetric_hessian": True,
    }
    policy = StationaryPointValidationPolicyV1(
        **policy_body, policy_sha256=canonical_sha256(policy_body)
    )
    approval = build_frozen_workflow_approval(
        approval_id="water-policy-approval",
        plan=plan,
        materialized_workflow=materialized,
        resources=resources,
        environment_identity_sha256s=("f" * 64,),
        stationary_point_policy=policy,
    )

    _validate_stationary_point_policy_binding(
        approval,
        policy,
        plan=plan,
        hessian_node_id="hess-optimized",
        require_for_hessian=True,
    )
    with pytest.raises(ContractError, match="another Hessian node"):
        _validate_stationary_point_policy_binding(
            approval,
            policy,
            plan=plan,
            hessian_node_id="different-hess",
            require_for_hessian=True,
        )


def test_event_store_consumes_approval_once_and_replays_node_frontier(tmp_path):
    plan = _water_plan()
    _, materialized, approval = _approval(plan)
    store = RuntimeEventStore(
        tmp_path / "events" / "runtime.jsonl", session_id="water-session"
    )
    store.record_materialized_workflow(
        turn_id="turn-1", workflow=materialized
    )
    _, started_event, started = store.consume_and_start_workflow(
        turn_id="turn-1",
        plan=plan,
        approval=approval,
        run_id="water-run",
        node_id="sp-initial",
        invocation_sha256="7" * 64,
        timestamp="2026-08-04T00:00:00+00:00",
    )

    assert started_event.kind == "workflow_execution_started"
    assert derive_ready_node_ids(plan, started) == ("opt-initial",)
    with pytest.raises(ContractError, match="already been consumed"):
        store.consume_and_start_workflow(
            turn_id="turn-2",
            plan=plan,
            approval=approval,
            run_id="water-run-2",
            node_id="opt-initial",
            invocation_sha256="8" * 64,
            timestamp="2026-08-04T00:00:01+00:00",
        )

    _, engine_complete = store.transition_workflow_run_node(
        turn_id="turn-1",
        run_id="water-run",
        node_id="sp-initial",
        new_state="engine_complete",
        execution_receipt_sha256="8" * 64,
        output_artifact_sha256s=("9" * 64,),
        timestamp="2026-08-04T00:00:02+00:00",
    )
    assert engine_complete.state == "running"
    validation = _sp_result_validation()
    with pytest.raises(ContractError, match="typed result validation receipt"):
        store.transition_workflow_run_node(
            turn_id="turn-1",
            run_id="water-run",
            node_id="sp-initial",
            new_state="validated",
            validator_receipt_sha256s=(validation.receipt_sha256,),
            timestamp="2026-08-04T00:00:03+00:00",
        )
    store.append(
        turn_id="turn-1",
        kind=EventKind.RESULT_VERIFIED.value,
        payload={
            "receipt_sha256": validation.receipt_sha256,
            "status": "valid",
            "critical_finding_count": 0,
            "record": canonical_data(validation),
        },
        idempotency_key="result-validation:" + validation.receipt_sha256,
    )
    _, validated = store.transition_workflow_run_node(
        turn_id="turn-1",
        run_id="water-run",
        node_id="sp-initial",
        new_state="validated",
        validator_receipt_sha256s=(validation.receipt_sha256,),
        result_validation_receipt=validation,
        timestamp="2026-08-04T00:00:03+00:00",
    )
    assert validated.nodes[2].output_artifact_sha256s == ("9" * 64,)

    state = store.state()
    assert state.consumed_workflow_approval_ids == ["water-approval"]
    assert state.workflow_run_records["water-run"]["state"] == "running"
    assert RuntimeState().workflow_run_records == {}
