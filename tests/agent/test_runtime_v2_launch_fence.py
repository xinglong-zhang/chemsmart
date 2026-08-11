from __future__ import annotations

import pytest

from chemsmart.agent._contracts import ContractError, canonical_data, canonical_sha256
from chemsmart.agent.execution import (
    ProgramExecutionInvocationV1,
    build_execution_resource_spec,
    build_frozen_workflow_approval,
    build_program_execution_receipt,
)
from chemsmart.agent.runtime.event_store import RuntimeEventStore
from chemsmart.agent.runtime.records import (
    scientific_workflow_plan_from_record,
)
from chemsmart.agent.workflows import (
    MaterializedNodeV1,
    ScientificWorkflowNodeV2,
    build_materialized_workflow,
    build_scientific_workflow_plan,
)


def _frontier(tmp_path):
    plan = build_scientific_workflow_plan(
        workflow_id="water-workflow",
        task_spec_sha256="a" * 64,
        scientific_identity_sha256="b" * 64,
        nodes=(
            ScientificWorkflowNodeV2(
                node_id="sp-initial",
                stage="sp",
                requested_program="pyscf",
                program="pyscf",
                engine="cpu",
                project_role="water-project",
                unresolved_fields=(),
            ),
        ),
    )
    resources = build_execution_resource_spec(
        execution_target="run",
        cores=4,
        memory_gb=4,
        gpu_count=0,
        scratch_policy="none",
        node_timeout_seconds=600,
    )
    invocation_body = {
        "schema_version": "chemsmart.program-execution-invocation.v1",
        "node_id": "sp-initial",
        "approval_sha256": "4" * 64,
        "program": "pyscf",
        "engine": "cpu",
        "jobtype": "sp",
        "project_sha256": "e" * 64,
        "input_artifact_id": "water-xyz",
        "input_sha256": "d" * 64,
        "scientific_identity_sha256": plan.scientific_identity_sha256,
        "environment_receipt_sha256": "1" * 64,
        "resource_sha256": resources.resource_sha256,
        "workspace": str(tmp_path.resolve()),
        "argv": ("chemsmart", "run", "pyscf", "sp"),
        "idempotency_key": "5" * 64,
        "status": "ready",
    }
    invocation = ProgramExecutionInvocationV1(
        **invocation_body,
        invocation_sha256=canonical_sha256(invocation_body),
    )
    materialized = build_materialized_workflow(
        plan=plan,
        live_cli_schema_sha256="c" * 64,
        resource_sha256=resources.resource_sha256,
        nodes=(
            MaterializedNodeV1(
                node_id="sp-initial",
                input_artifact_sha256="d" * 64,
                project_artifact_sha256="e" * 64,
                project_validation_receipt_sha256="f" * 64,
                environment_receipt_sha256="1" * 64,
                invocation_sha256=invocation.invocation_sha256,
                preflight_receipt_sha256="3" * 64,
                state="previewed",
            ),
        ),
        unresolved_node_ids=(),
        status="ready_for_approval",
    )
    approval = build_frozen_workflow_approval(
        approval_id="water-approval",
        plan=plan,
        materialized_workflow=materialized,
        resources=resources,
        environment_identity_sha256s=("1" * 64,),
    )
    return plan, materialized, approval, invocation


def _reserve(store, tmp_path):
    plan, materialized, approval, invocation = _frontier(tmp_path)
    result = store.reserve_workflow_node_launch(
        turn_id="turn-1",
        plan=plan,
        materialized_workflow=materialized,
        approval=approval,
        invocation=invocation,
        run_id="run.water-approval",
        timestamp="2026-08-04T00:00:00+00:00",
    )
    return result, plan, materialized, approval, invocation


def test_launch_reservation_is_one_atomic_full_frontier_event(tmp_path):
    store = RuntimeEventStore(
        tmp_path / "events" / "runtime.jsonl", session_id="water-session"
    )

    result, plan, materialized, approval, invocation = _reserve(store, tmp_path)

    assert result.status == "reserved"
    events = store.read_events()
    assert len(events) == 1
    assert events[0].kind == "workflow_node_launch_reserved"
    frontier = store.workflow_frontier(
        workflow_id=plan.workflow_id, run_id="run.water-approval"
    )
    assert frontier.plan == plan
    assert frontier.materialized_workflow == materialized
    assert frontier.approval == approval
    assert frontier.run_state == result.run_state
    assert frontier.invocations == (invocation,)
    assert frontier.legacy_incomplete is False
    assert store.state().consumed_workflow_approval_ids == ["water-approval"]


def test_unresolved_reservation_blocks_relaunch(tmp_path):
    store = RuntimeEventStore(
        tmp_path / "events" / "runtime.jsonl", session_id="water-session"
    )
    _, plan, materialized, approval, invocation = _reserve(store, tmp_path)

    with pytest.raises(ContractError, match="remains unresolved"):
        store.reserve_workflow_node_launch(
            turn_id="turn-2",
            plan=plan,
            materialized_workflow=materialized,
            approval=approval,
            invocation=invocation,
            run_id="run.water-approval",
            timestamp="2026-08-04T00:00:01+00:00",
        )


def test_terminal_receipt_replays_exactly_without_second_launch(tmp_path):
    path = tmp_path / "events" / "runtime.jsonl"
    store = RuntimeEventStore(path, session_id="water-session")
    _, plan, materialized, approval, invocation = _reserve(store, tmp_path)
    receipt = build_program_execution_receipt(
        invocation,
        execution_state="failed",
        exit_status=1,
        child_exit_status=1,
        engine_complete=False,
        validated=False,
        findings=("execution.engine.failed",),
        started_at="2026-08-04T00:00:00+00:00",
        finished_at="2026-08-04T00:00:05+00:00",
    )
    store.record_program_execution_receipt(
        turn_id="turn-1",
        workflow_id=plan.workflow_id,
        run_id="run.water-approval",
        receipt=receipt,
    )

    restarted = RuntimeEventStore(path, session_id="water-session")
    replayed = restarted.replayed_execution_receipt(
        workflow_id=plan.workflow_id,
        run_id="run.water-approval",
        node_id="sp-initial",
    )
    assert replayed == receipt
    result = restarted.reserve_workflow_node_launch(
        turn_id="turn-2",
        plan=plan,
        materialized_workflow=materialized,
        approval=approval,
        invocation=invocation,
        run_id="run.water-approval",
        timestamp="2026-08-04T00:00:06+00:00",
    )
    assert result.status == "terminal_replay"
    assert result.execution_receipt == receipt
    assert len(restarted.read_events()) == 2


def test_legacy_two_event_start_is_inspectable_but_cannot_launch(tmp_path):
    store = RuntimeEventStore(
        tmp_path / "events" / "runtime.jsonl", session_id="water-session"
    )
    plan, materialized, approval, invocation = _frontier(tmp_path)
    store.record_materialized_workflow(
        turn_id="turn-1", workflow=materialized
    )
    store.consume_and_start_workflow(
        turn_id="turn-1",
        plan=plan,
        approval=approval,
        run_id="run.water-approval",
        node_id="sp-initial",
        invocation_sha256=invocation.invocation_sha256,
        timestamp="2026-08-04T00:00:00+00:00",
    )

    assert store.state().workflow_run_records["run.water-approval"]
    with pytest.raises(ContractError, match="legacy-incomplete"):
        store.reserve_workflow_node_launch(
            turn_id="turn-2",
            plan=plan,
            materialized_workflow=materialized,
            approval=approval,
            invocation=invocation,
            run_id="run.water-approval",
            timestamp="2026-08-04T00:00:01+00:00",
        )


def test_strict_plan_loader_rejects_extra_and_corrupt_fields(tmp_path):
    plan, _, _, _ = _frontier(tmp_path)
    extra = canonical_data(plan)
    extra["model_note"] = "not authoritative"
    with pytest.raises(ContractError, match="fields differ"):
        scientific_workflow_plan_from_record(extra)

    corrupted = canonical_data(plan)
    corrupted["status"] = "previewed"
    with pytest.raises(ContractError, match="can only be planned"):
        scientific_workflow_plan_from_record(corrupted)
