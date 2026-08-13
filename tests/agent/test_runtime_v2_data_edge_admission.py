from __future__ import annotations

import pytest

from chemsmart.agent._contracts import ContractError, canonical_data, canonical_sha256
from chemsmart.agent.execution import (
    FrozenWorkflowApprovalV1,
    OptimizedGeometryHandoffV1,
    ProgramExecutionInvocationV1,
    ProgramExecutionReceiptV1,
    ProgramResultValidationReceiptV1,
    ValidatedDataEdgeBindingV1,
    build_execution_resource_spec,
    build_frozen_workflow_approval,
    build_producer_edge_rule,
    build_validated_data_edge_binding,
    build_workflow_run_state,
    derive_ready_node_ids,
    transition_workflow_node,
)
from chemsmart.agent.runtime.records import frozen_workflow_approval_from_record
from chemsmart.agent.workflows import (
    MaterializedNodeV1,
    ScientificWorkflowEdgeV2,
    ScientificWorkflowNodeV2,
    build_materialized_workflow,
    build_scientific_workflow_plan,
)


def _plan(*, producer_stage: str = "opt"):
    return build_scientific_workflow_plan(
        workflow_id="water-opt-hess",
        task_spec_sha256="a" * 64,
        scientific_identity_sha256="b" * 64,
        nodes=(
            ScientificWorkflowNodeV2(
                node_id="opt-initial",
                stage=producer_stage,
                requested_program="pyscf",
                program="pyscf",
                engine="cpu",
                project_role="water-opt",
                unresolved_fields=(),
            ),
            ScientificWorkflowNodeV2(
                node_id="hess-optimized",
                stage="hess",
                requested_program="pyscf",
                program="pyscf",
                engine="cpu",
                project_role="water-hess",
                unresolved_fields=(),
            ),
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


def _frontier(*, state: str = "previewed", producer_stage: str = "opt"):
    plan = _plan(producer_stage=producer_stage)
    resources = build_execution_resource_spec(
        execution_target="run",
        cores=4,
        memory_gb=4,
        gpu_count=0,
        scratch_policy="none",
        node_timeout_seconds=600,
    )
    node = MaterializedNodeV1(
        node_id="opt-initial",
        input_artifact_sha256="c" * 64,
        project_artifact_sha256="d" * 64,
        project_validation_receipt_sha256="e" * 64,
        environment_receipt_sha256="f" * 64,
        invocation_sha256="1" * 64,
        preflight_receipt_sha256="2" * 64 if state == "previewed" else "",
        state=state,
    )
    materialized = build_materialized_workflow(
        plan=plan,
        live_cli_schema_sha256="3" * 64,
        resource_sha256=resources.resource_sha256,
        nodes=(node,),
        unresolved_node_ids=("hess-optimized",),
        status="partial",
    )
    return plan, resources, materialized


def _approval(*, producer_stage: str = "opt"):
    plan, resources, materialized = _frontier(
        producer_stage=producer_stage
    )
    approval = build_frozen_workflow_approval(
        approval_id="water-approval",
        plan=plan,
        materialized_workflow=materialized,
        resources=resources,
        environment_identity_sha256s=("f" * 64,),
    )
    return plan, materialized, approval


def _result_validation_receipt(
    *, state: str = "valid",
) -> ProgramResultValidationReceiptV1:
    findings = (
        ()
        if state == "valid"
        else (
            "result_validation.state_not_validated",
            "seeded.validation.failed",
        )
    )
    body = {
        "schema_version": "chemsmart.program-result-validation-receipt.v1",
        "validator_id": "pyscf-result-validator",
        "validator_schema_version": "chemsmart.pyscf-result-validation.v1",
        "validator_version": "1",
        "invocation_sha256": "1" * 64,
        "node_id": "opt-initial",
        "program": "pyscf",
        "engine": "cpu",
        "jobtype": "opt",
        "input_artifact_sha256": "c" * 64,
        "project_artifact_sha256": "d" * 64,
        "capability_environment_receipt_sha256": "f" * 64,
        "run_environment_receipt_sha256": "",
        "environment_validation_sha256": "",
        "stationary_point_policy_sha256": "",
        "output_artifacts": (),
        "observations": {"state": "validated" if state == "valid" else "failed"},
        "findings": findings,
        "state": state,
    }
    return ProgramResultValidationReceiptV1(
        **body, receipt_sha256=canonical_sha256(body)
    )


def _validated_opt_run(plan, approval):
    run = build_workflow_run_state(
        run_id="run.water-approval",
        plan=plan,
        approval=approval,
        approval_consumed=True,
    )
    run = transition_workflow_node(
        run,
        node_id="opt-initial",
        new_state="running",
        invocation_sha256="1" * 64,
        timestamp="2026-08-04T00:00:00+00:00",
    )
    run = transition_workflow_node(
        run,
        node_id="opt-initial",
        new_state="engine_complete",
        execution_receipt_sha256="4" * 64,
        output_artifact_sha256s=("5" * 64, "6" * 64),
        timestamp="2026-08-04T00:00:01+00:00",
    )
    validation = _result_validation_receipt()
    return transition_workflow_node(
        run,
        node_id="opt-initial",
        new_state="validated",
        validator_receipt_sha256s=(validation.receipt_sha256,),
        result_validation_receipt=validation,
        timestamp="2026-08-04T00:00:02+00:00",
    )


def _edge_binding(plan, approval):
    rule = approval.producer_edge_rules[0]
    body = {
        "schema_version": "chemsmart.validated-data-edge-binding.v1",
        "run_id": "run.water-approval",
        "workflow_id": plan.workflow_id,
        "plan_sha256": plan.plan_sha256,
        "approval_sha256": approval.approval_sha256,
        "scientific_edge_sha256": canonical_sha256(plan.edges[0]),
        "producer_rule_sha256": rule.rule_sha256,
        "source_node_id": "opt-initial",
        "target_node_id": "hess-optimized",
        "artifact_class": "geometry_xyz",
        "producer_output_id": "optimized-geometry",
        "consumer_input_id": "geometry",
        "selection_rule": "validated_optimized_geometry",
        "producer_execution_receipt_sha256": "4" * 64,
        "producer_validator_receipt_sha256s": (
            _result_validation_receipt().receipt_sha256,
        ),
        "source_artifact_id": "opt-result",
        "source_artifact_sha256": "5" * 64,
        "selected_artifact_id": "optimized-geometry",
        "selected_artifact_sha256": "6" * 64,
        "producer_scientific_identity_sha256": "b" * 64,
        "consumer_scientific_identity_sha256": "8" * 64,
        "atom_order_sha256": "9" * 64,
        "positions_sha256": "a" * 64,
        "charge": 0,
        "multiplicity": 1,
        "handoff_receipt_sha256": "b" * 64,
        "status": "validated",
    }
    return ValidatedDataEdgeBindingV1(
        **body, receipt_sha256=canonical_sha256(body)
    )


def test_frozen_approval_rejects_compiled_materialized_node():
    plan, resources, materialized = _frontier(state="compiled")

    with pytest.raises(ContractError, match="requires exact preview"):
        build_frozen_workflow_approval(
            approval_id="water-approval",
            plan=plan,
            materialized_workflow=materialized,
            resources=resources,
            environment_identity_sha256s=("f" * 64,),
        )


def test_frozen_approval_records_preview_and_exact_future_selection_rule():
    _, _, approval = _approval()

    assert approval.execution_admissible is True
    assert approval.materialized_preview_bindings[0].preflight_receipt_sha256 == (
        "2" * 64
    )
    rule = approval.producer_edge_rules[0]
    assert rule.selection_rule == "validated_optimized_geometry"
    assert rule.environment_receipt_sha256 == "f" * 64
    assert rule.preserve_atom_order is True
    assert rule.preserve_electronic_state is True


def test_ts_producer_is_deferred_under_the_same_exact_geometry_rule():
    plan, _, approval = _approval(producer_stage="ts")

    assert plan.nodes[0].stage == "ts"
    rule = approval.producer_edge_rules[0]
    assert rule.source_node_id == "opt-initial"
    assert rule.producer_output_id == "optimized-geometry"
    assert rule.selection_rule == "validated_optimized_geometry"


def test_future_node_requires_an_exact_environment_when_multiple_are_approved():
    plan, resources, materialized = _frontier()

    with pytest.raises(
        ContractError,
        match="node-specific environment evidence",
    ):
        build_frozen_workflow_approval(
            approval_id="water-approval",
            plan=plan,
            materialized_workflow=materialized,
            resources=resources,
            environment_identity_sha256s=("e" * 64, "f" * 64),
        )

    approval = build_frozen_workflow_approval(
        approval_id="water-approval",
        plan=plan,
        materialized_workflow=materialized,
        resources=resources,
        environment_identity_sha256s=("e" * 64, "f" * 64),
        future_node_environment_identity_sha256s={
            "hess-optimized": "f" * 64,
        },
    )
    assert approval.producer_edge_rules[0].environment_receipt_sha256 == (
        "f" * 64
    )


def test_validated_transition_rejects_bare_or_invalid_validation_receipt():
    plan, _, approval = _approval()
    run = build_workflow_run_state(
        run_id="run.water-approval",
        plan=plan,
        approval=approval,
        approval_consumed=True,
    )
    run = transition_workflow_node(
        run,
        node_id="opt-initial",
        new_state="running",
        invocation_sha256="1" * 64,
        timestamp="2026-08-04T00:00:00+00:00",
    )
    run = transition_workflow_node(
        run,
        node_id="opt-initial",
        new_state="engine_complete",
        execution_receipt_sha256="4" * 64,
        timestamp="2026-08-04T00:00:01+00:00",
    )
    with pytest.raises(ContractError, match="typed result validation receipt"):
        transition_workflow_node(
            run,
            node_id="opt-initial",
            new_state="validated",
            validator_receipt_sha256s=("7" * 64,),
            timestamp="2026-08-04T00:00:02+00:00",
        )
    invalid = _result_validation_receipt(state="invalid")
    with pytest.raises(ContractError, match="requires a valid result"):
        transition_workflow_node(
            run,
            node_id="opt-initial",
            new_state="validated",
            validator_receipt_sha256s=(invalid.receipt_sha256,),
            result_validation_receipt=invalid,
            timestamp="2026-08-04T00:00:02+00:00",
        )


def test_arbitrary_producer_digest_does_not_ready_data_consumer():
    plan, _, approval = _approval()
    run = _validated_opt_run(plan, approval)

    assert derive_ready_node_ids(plan, run) == ()
    assert derive_ready_node_ids(
        plan, run, (_edge_binding(plan, approval),)
    ) == ("hess-optimized",)
    mismatched_body = canonical_data(_edge_binding(plan, approval))
    mismatched_body.pop("receipt_sha256")
    mismatched_body["selected_artifact_sha256"] = "c" * 64
    mismatched_body["producer_validator_receipt_sha256s"] = tuple(
        mismatched_body["producer_validator_receipt_sha256s"]
    )
    mismatched = ValidatedDataEdgeBindingV1(
        **mismatched_body,
        receipt_sha256=canonical_sha256(mismatched_body),
    )
    assert derive_ready_node_ids(plan, run, (mismatched,)) == ()


def test_data_edge_binds_node_identity_not_multi_geometry_plan_aggregate():
    plan, _, approval = _approval()
    producer_identity = "7" * 64
    invocation_body = {
        "schema_version": "chemsmart.program-execution-invocation.v1",
        "node_id": "opt-initial",
        "approval_sha256": "a" * 64,
        "program": "pyscf",
        "engine": "cpu",
        "jobtype": "opt",
        "project_sha256": "d" * 64,
        "input_artifact_id": "geometry.initial",
        "input_sha256": "c" * 64,
        "scientific_identity_sha256": producer_identity,
        "environment_receipt_sha256": "f" * 64,
        "resource_sha256": "e" * 64,
        "workspace": "/tmp/chemsmart-data-edge-test",
        "argv": ("chemsmart", "run", "pyscf", "opt"),
        "idempotency_key": "1" * 64,
        "status": "ready",
    }
    invocation = ProgramExecutionInvocationV1(
        **invocation_body,
        invocation_sha256=canonical_sha256(invocation_body),
    )
    result = {
        "artifact_id": "result.opt",
        "kind": "pyscf_hdf5",
        "sha256": "5" * 64,
        "size_bytes": 1,
        "path": "/tmp/chemsmart-data-edge-test/result.h5",
        "cli_value": "/tmp/chemsmart-data-edge-test/result.h5",
    }
    from chemsmart.agent._contracts import TrustedArtifactRefV1

    result_artifact = TrustedArtifactRefV1(**result)
    receipt_body = {
        "schema_version": "chemsmart.program-execution-receipt.v1",
        "invocation_sha256": invocation.invocation_sha256,
        "node_id": "opt-initial",
        "idempotency_key": invocation.idempotency_key,
        "execution_state": "validated",
        "wrapper_exit_status": 0,
        "child_exit_status": 0,
        "engine_complete": True,
        "validated": True,
        "engine_receipt_sha256": "6" * 64,
        "environment_receipt_sha256": invocation.environment_receipt_sha256,
        "result_validation_receipt_sha256": "4" * 64,
        "output_artifacts": (result_artifact,),
        "validator_receipt_sha256s": ("4" * 64,),
        "findings": (),
        "started_at": "2026-08-10T00:00:00+00:00",
        "finished_at": "2026-08-10T00:00:01+00:00",
    }
    receipt = ProgramExecutionReceiptV1(
        **receipt_body, receipt_sha256=canonical_sha256(receipt_body)
    )
    producer_edge = build_producer_edge_rule(
        producer_node_id="opt-initial",
        consumer_node_id="hess-optimized",
        artifact_kind="geometry_xyz",
        selection_rule="validated_optimized_geometry",
    )
    handoff_body = {
        "schema_version": "chemsmart.optimized-geometry-handoff.v1",
        "producer_node_id": "opt-initial",
        "consumer_node_id": "hess-optimized",
        "producer_edge_sha256": producer_edge.edge_sha256,
        "producer_execution_receipt_sha256": receipt.receipt_sha256,
        "result_artifact_id": result_artifact.artifact_id,
        "result_artifact_sha256": result_artifact.sha256,
        "geometry_artifact_id": "geometry.optimized",
        "geometry_artifact_sha256": "6" * 64,
        "atom_count": 3,
        "symbols": ("O", "H", "H"),
        "positions_sha256": "a" * 64,
        "charge": 0,
        "multiplicity": 1,
        "status": "validated_handoff",
    }
    handoff = OptimizedGeometryHandoffV1(
        **handoff_body, receipt_sha256=canonical_sha256(handoff_body)
    )

    binding = build_validated_data_edge_binding(
        run_id="run.water-approval",
        plan=plan,
        approval=approval,
        scientific_edge=plan.edges[0],
        producer_edge=producer_edge,
        producer_invocation=invocation,
        producer_receipt=receipt,
        handoff=handoff,
        producer_scientific_identity_sha256=producer_identity,
        consumer_scientific_identity_sha256="8" * 64,
    )

    assert producer_identity != plan.scientific_identity_sha256
    assert binding.producer_scientific_identity_sha256 == producer_identity
    with pytest.raises(ContractError, match="execution invocation"):
        build_validated_data_edge_binding(
            run_id="run.water-approval",
            plan=plan,
            approval=approval,
            scientific_edge=plan.edges[0],
            producer_edge=producer_edge,
            producer_invocation=invocation,
            producer_receipt=receipt,
            handoff=handoff,
            producer_scientific_identity_sha256="9" * 64,
            consumer_scientific_identity_sha256="8" * 64,
        )


def test_legacy_frozen_approval_remains_readable_but_not_admissible():
    _, _, approval = _approval()
    record = canonical_data(approval)
    for field in (
        "materialized_preview_bindings",
        "producer_edge_rules",
        "admission_sha256",
    ):
        record.pop(field)
    record.pop("approval_sha256")
    record["approval_sha256"] = canonical_sha256(record)

    observed = frozen_workflow_approval_from_record(record)

    assert isinstance(observed, FrozenWorkflowApprovalV1)
    assert observed.execution_admissible is False
