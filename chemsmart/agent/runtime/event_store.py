"""Secure, locked, append-only Runtime V2 JSONL event store."""

from __future__ import annotations

import hashlib
import json
import os
import stat
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Iterable, Iterator, Mapping, TextIO

try:  # POSIX advisory locking.
    import fcntl as _fcntl
except ImportError:  # pragma: no cover - Windows import boundary
    _fcntl = None

try:  # Windows byte-range locking.
    import msvcrt as _msvcrt
except ImportError:  # pragma: no cover - POSIX import boundary
    _msvcrt = None

from chemsmart.agent._contracts import (
    ContractError,
    canonical_data,
    canonical_sha256,
)
from chemsmart.agent.execution import (
    FrozenWorkflowApprovalV1,
    ProgramExecutionInvocationV1,
    ProgramExecutionReceiptV1,
    ProgramResultValidationReceiptV1,
    ValidatedDataEdgeBindingV1,
    WorkflowExecutionApprovalBundleV1,
    WorkflowReviewResolutionV1,
    WorkflowRunStateV1,
    build_workflow_run_state,
    derive_ready_node_ids,
    transition_workflow_node,
)
from chemsmart.agent.permissions import (
    ApprovalResolutionV1,
    PermissionReceiptV1,
    PermissionRequestV1,
    evaluate_permission,
)
from chemsmart.agent.runtime.events import (
    ANALYSIS_CLAIMS_RECORDED,
    ANALYSIS_COMPLETION_EVALUATED,
    CAPABILITY_QUERIED,
    COMMAND_COMPILED,
    COMMAND_INSPECTED,
    ENGINE_BOUND,
    ENVIRONMENT_QUERIED,
    EXECUTION_BUNDLE_CONSUMED,
    OPTIMIZED_GEOMETRY_HANDED_OFF,
    PERMISSION_RESOLVED,
    PROGRAM_BOUND,
    PROGRAM_EXECUTED,
    PROGRAM_PREFLIGHTED,
    PROJECT_PROMOTED,
    PROJECT_VALIDATED,
    QUANTITY_EXPRESSION_EVALUATED,
    RESULT_QUANTITIES_EXTRACTED,
    RESULT_VERIFIED,
    RUNTIME_TERMINATED,
    SAFE_PREVIEWED,
    SCIENTIFIC_DECISION_RECORDED,
    SCIENTIFIC_VALIDATION_EVALUATED,
    SCIENTIFIC_WORKFLOW_MATERIALIZED,
    SUBSTITUTION_ASSESSED,
    THERMOCHEMISTRY_DERIVED,
    VALIDATOR_OBSERVED,
    WORKFLOW_APPROVAL_CONSUMED,
    WORKFLOW_DATA_EDGE_BOUND,
    WORKFLOW_EXECUTION_STARTED,
    WORKFLOW_LAUNCH_RESERVED,
    WORKFLOW_NODE_STATE_CHANGED,
    WORKFLOW_REVIEW_RESOLVED,
    RuntimeEvent,
)
from chemsmart.agent.runtime.records import (
    LaunchFenceResultV1,
    ReconstructedWorkflowFrontierV1,
    build_workflow_node_launch_reservation,
    canonical_record,
    reconstruct_workflow_frontier,
    workflow_run_state_from_record,
)
from chemsmart.agent.runtime.reducer import RuntimeState, replay_events
from chemsmart.agent.workflows import (
    MaterializedWorkflowV1,
    ScientificWorkflowPlanV2,
)


class RuntimeEventStore:
    """One crash-stable stream with atomic one-shot approval consumption."""

    def __init__(self, path: str | Path, *, session_id: str) -> None:
        # abspath deliberately does not follow the final symlink; secure open
        # rejects links before and during open.
        self.path = Path(os.path.abspath(os.fspath(path)))
        self.session_id = str(session_id)
        if not self.session_id:
            raise ContractError("event store session_id must not be empty")
        self._prepare_private_parent()
        self._lock_path = self.path.with_name("." + self.path.name + ".lock")

    def read_events(self) -> tuple[RuntimeEvent, ...]:
        with self._locked_handle(exclusive=False) as handle:
            return self._read_locked(handle)

    def state(self) -> RuntimeState:
        return replay_events(self.read_events())

    def workflow_frontier(
        self, *, workflow_id: str, run_id: str
    ) -> ReconstructedWorkflowFrontierV1:
        """Return the strict typed frontier reconstructed from this stream."""

        return reconstruct_workflow_frontier(
            self.state(), workflow_id=workflow_id, run_id=run_id
        )

    def append(
        self,
        *,
        turn_id: str,
        kind: str,
        payload: Mapping[str, Any] | None = None,
        idempotency_key: str = "",
    ) -> RuntimeEvent:
        normalized_kind = str(getattr(kind, "value", kind))
        if normalized_kind == RUNTIME_TERMINATED and (payload or {}).get(
            "terminal_state"
        ) in {"complete", "planned"}:
            raise ContractError(
                "host-derived termination must pass RuntimeEventStore.terminate"
            )
        with self._locked_handle(exclusive=True) as handle:
            return self._append_locked(
                handle,
                turn_id=turn_id,
                kind=normalized_kind,
                payload=dict(payload or {}),
                idempotency_key=idempotency_key,
            )

    def consume_approval(
        self,
        *,
        turn_id: str,
        request: PermissionRequestV1,
        approval: ApprovalResolutionV1,
    ) -> tuple[PermissionReceiptV1, RuntimeEvent]:
        """Evaluate and consume an exact one-shot approval under one lock."""

        with self._locked_handle(exclusive=True) as handle:
            events = self._read_locked(handle)
            consumed = any(
                event.kind == PERMISSION_RESOLVED
                and event.payload.get("approval_id") == approval.approval_id
                and event.payload.get("decision") == "allow_once"
                for event in events
            )
            if consumed:
                raise ContractError("approval has already been consumed")
            receipt = evaluate_permission(request, approval=approval)
            payload = {
                "receipt_sha256": receipt.receipt_sha256,
                "permission_request_sha256": request.request_sha256,
                "approval_id": approval.approval_id,
                "approval_resolution_sha256": approval.resolution_sha256,
                "decision": receipt.decision.value,
                "reason": receipt.reason,
            }
            event = self._append_locked(
                handle,
                turn_id=turn_id,
                kind=PERMISSION_RESOLVED,
                payload=payload,
                idempotency_key=(
                    "approval:"
                    + approval.approval_id
                    + ":"
                    + request.request_sha256
                ),
                existing_events=events,
            )
            return receipt, event

    def record_materialized_workflow(
        self,
        *,
        turn_id: str,
        workflow: MaterializedWorkflowV1,
    ) -> RuntimeEvent:
        """Persist the canonical grounded workflow, not only its digest."""

        record = canonical_data(workflow)
        receipt_sha256 = canonical_sha256(record)
        return self.append(
            turn_id=turn_id,
            kind=SCIENTIFIC_WORKFLOW_MATERIALIZED,
            payload={
                "receipt_sha256": receipt_sha256,
                "workflow_id": workflow.workflow_id,
                "plan_sha256": workflow.plan_sha256,
                "status": workflow.status,
                "record": record,
            },
            idempotency_key=(
                "materialized-workflow:" + workflow.materialized_sha256
            ),
        )

    def record_workflow_review_resolution(
        self,
        *,
        turn_id: str,
        resolution: WorkflowReviewResolutionV1,
    ) -> RuntimeEvent:
        """Persist approve/deny/revise/quit against the exact reviewed digest."""

        record = canonical_data(resolution)
        receipt_sha256 = canonical_sha256(record)
        return self.append(
            turn_id=turn_id,
            kind=WORKFLOW_REVIEW_RESOLVED,
            payload={
                "receipt_sha256": receipt_sha256,
                "review_sha256": resolution.review_sha256,
                "resolution_id": resolution.resolution_id,
                "approval_id": resolution.approval_id,
                "decision": resolution.decision,
                "record": record,
            },
            idempotency_key=(
                "workflow-review-resolution:" + resolution.resolution_id
            ),
        )

    def claim_execution_bundle(
        self,
        *,
        turn_id: str,
        bundle: WorkflowExecutionApprovalBundleV1,
    ) -> RuntimeEvent:
        """Atomically consume one compound approval before execution starts."""

        with self._locked_handle(exclusive=True) as handle:
            events = self._read_locked(handle)
            state = replay_events(events)
            if (
                bundle.bundle_sha256 in state.consumed_execution_bundle_sha256s
                or bundle.workflow_approval.approval_id
                in state.consumed_execution_approval_ids
            ):
                raise ContractError(
                    "execution approval bundle has already been consumed"
                )
            record = {
                "schema_version": "chemsmart.execution-bundle-consumption.v1",
                "bundle_sha256": bundle.bundle_sha256,
                "review_sha256": bundle.review_sha256,
                "approval_id": bundle.workflow_approval.approval_id,
                "one_shot": True,
                "status": "consumed",
            }
            return self._append_locked(
                handle,
                turn_id=turn_id,
                kind=EXECUTION_BUNDLE_CONSUMED,
                payload={
                    "receipt_sha256": canonical_sha256(record),
                    "bundle_sha256": bundle.bundle_sha256,
                    "review_sha256": bundle.review_sha256,
                    "approval_id": bundle.workflow_approval.approval_id,
                    "status": "consumed",
                    "record": record,
                },
                idempotency_key=(
                    "execution-bundle-consumed:" + bundle.bundle_sha256
                ),
                existing_events=events,
            )

    def consume_and_start_workflow(
        self,
        *,
        turn_id: str,
        plan: ScientificWorkflowPlanV2,
        approval: FrozenWorkflowApprovalV1,
        run_id: str,
        node_id: str,
        invocation_sha256: str,
        timestamp: str,
    ) -> tuple[RuntimeEvent, RuntimeEvent, WorkflowRunStateV1]:
        """Consume one frozen approval and durably start one ready node.

        Both events are appended while holding the same exclusive stream lock.
        A second caller therefore cannot consume the approval or start the same
        run concurrently.  The host still records this transition before it
        invokes an engine.
        """

        with self._locked_handle(exclusive=True) as handle:
            events = self._read_locked(handle)
            state = replay_events(events)
            if approval.approval_id in state.consumed_workflow_approval_ids:
                raise ContractError(
                    "workflow approval has already been consumed"
                )
            normalized_run_id = str(run_id).strip().lower()
            if not normalized_run_id:
                raise ContractError("run_id must not be empty")
            if normalized_run_id in state.workflow_run_records:
                raise ContractError("workflow run has already started")
            initial = build_workflow_run_state(
                run_id=normalized_run_id,
                plan=plan,
                approval=approval,
                approval_consumed=True,
            )
            if node_id not in derive_ready_node_ids(plan, initial):
                raise ContractError("node is not ready in the frozen workflow")
            started = transition_workflow_node(
                initial,
                node_id=node_id,
                new_state="running",
                plan=plan,
                invocation_sha256=invocation_sha256,
                timestamp=timestamp,
            )
            approval_record = canonical_data(approval)
            approval_receipt = canonical_sha256(approval_record)
            approval_event = self._append_locked(
                handle,
                turn_id=turn_id,
                kind=WORKFLOW_APPROVAL_CONSUMED,
                payload={
                    "receipt_sha256": approval_receipt,
                    "approval_id": approval.approval_id,
                    "approval_sha256": approval.approval_sha256,
                    "workflow_id": approval.workflow_id,
                    "status": "consumed",
                    "record": approval_record,
                },
                idempotency_key=("workflow-approval:" + approval.approval_id),
                existing_events=events,
            )
            run_record = canonical_data(started)
            start_event = self._append_locked(
                handle,
                turn_id=turn_id,
                kind=WORKFLOW_EXECUTION_STARTED,
                payload={
                    "receipt_sha256": canonical_sha256(run_record),
                    "run_id": started.run_id,
                    "workflow_id": started.workflow_id,
                    "approval_id": started.approval_id,
                    "node_id": node_id,
                    "state": "running",
                    "record": run_record,
                },
                idempotency_key="workflow-run-started:" + started.run_id,
                existing_events=events + (approval_event,),
            )
            return approval_event, start_event, started

    def replayed_execution_receipt(
        self,
        *,
        workflow_id: str,
        run_id: str,
        node_id: str,
    ) -> ProgramExecutionReceiptV1 | None:
        """Return a durable terminal receipt or reject an unsafe replay.

        A running reservation is deliberately not interpreted as a retryable
        failure.  The caller must reconcile the original process/output state
        rather than launching the same node again.
        """

        frontier = reconstruct_workflow_frontier(
            self.state(), workflow_id=workflow_id, run_id=run_id
        )
        receipt = frontier.receipt_for_node(node_id)
        reservation = frontier.reservation_for_node(node_id)
        if receipt is not None:
            if reservation is None:
                raise ContractError(
                    "execution receipt lacks a durable launch reservation"
                )
            if receipt.invocation_sha256 != reservation.invocation_sha256:
                raise ContractError(
                    "execution receipt differs from its launch reservation"
                )
            return receipt
        if reservation is not None:
            raise ContractError(
                "launch reservation remains unresolved; relaunch is forbidden"
            )
        if frontier.legacy_incomplete:
            raise ContractError(
                "legacy-incomplete workflow stream is inspectable but cannot launch"
            )
        return None

    def reserve_workflow_node_launch(
        self,
        *,
        turn_id: str,
        plan: ScientificWorkflowPlanV2,
        materialized_workflow: MaterializedWorkflowV1,
        approval: FrozenWorkflowApprovalV1,
        invocation: ProgramExecutionInvocationV1,
        run_id: str,
        timestamp: str,
    ) -> LaunchFenceResultV1:
        """Atomically consume approval and reserve one node before launch."""

        with self._locked_handle(exclusive=True) as handle:
            events = self._read_locked(handle)
            state = replay_events(events)
            frontier = reconstruct_workflow_frontier(
                state, workflow_id=plan.workflow_id, run_id=run_id
            )
            receipt = frontier.receipt_for_node(invocation.node_id)
            if receipt is not None:
                reservation = frontier.reservation_for_node(invocation.node_id)
                if reservation is None or (
                    receipt.invocation_sha256 != invocation.invocation_sha256
                    or reservation.invocation_sha256
                    != invocation.invocation_sha256
                ):
                    raise ContractError(
                        "replayed execution differs from requested invocation"
                    )
                return LaunchFenceResultV1(
                    status="terminal_replay",
                    reservation=None,
                    execution_receipt=receipt,
                    run_state=frontier.run_state,
                )
            prior_reservation = frontier.reservation_for_node(
                invocation.node_id
            )
            if prior_reservation is not None:
                raise ContractError(
                    "launch reservation remains unresolved; relaunch is forbidden"
                )
            if frontier.legacy_incomplete:
                raise ContractError(
                    "legacy-incomplete workflow stream is inspectable but cannot launch"
                )
            for persisted, supplied, label in (
                (frontier.plan, plan, "scientific plan"),
                (
                    frontier.materialized_workflow,
                    materialized_workflow,
                    "materialized workflow",
                ),
                (frontier.approval, approval, "frozen approval"),
            ):
                if persisted is not None and canonical_data(
                    persisted
                ) != canonical_data(supplied):
                    raise ContractError(
                        f"persisted {label} differs from launch"
                    )

            consumes_approval = frontier.run_state is None
            if consumes_approval:
                if (
                    approval.approval_id
                    in state.consumed_workflow_approval_ids
                ):
                    raise ContractError(
                        "workflow approval has already been consumed"
                    )
                current = build_workflow_run_state(
                    run_id=run_id,
                    plan=plan,
                    approval=approval,
                    approval_consumed=True,
                )
            else:
                current = frontier.run_state
                if current is None:  # pragma: no cover - narrowed above
                    raise ContractError("workflow run reconstruction failed")
                if (
                    current.plan_sha256 != plan.plan_sha256
                    or current.approval_sha256 != approval.approval_sha256
                ):
                    raise ContractError(
                        "persisted workflow run differs from launch approval"
                    )
            if invocation.node_id not in derive_ready_node_ids(
                plan, current, frontier.data_edge_bindings
            ):
                raise ContractError(
                    "node is not ready in the replayed workflow"
                )
            started = transition_workflow_node(
                current,
                node_id=invocation.node_id,
                new_state="running",
                plan=plan,
                invocation_sha256=invocation.invocation_sha256,
                timestamp=timestamp,
            )
            reservation = build_workflow_node_launch_reservation(
                run_id=run_id,
                plan=plan,
                materialized_workflow=materialized_workflow,
                approval=approval,
                invocation=invocation,
                data_edge_bindings=frontier.data_edge_bindings,
                consumes_approval=consumes_approval,
                reserved_at=timestamp,
            )
            reservation_record = canonical_record(reservation)
            self._append_locked(
                handle,
                turn_id=turn_id,
                kind=WORKFLOW_LAUNCH_RESERVED,
                payload={
                    "receipt_sha256": canonical_sha256(reservation_record),
                    "workflow_id": plan.workflow_id,
                    "run_id": started.run_id,
                    "node_id": invocation.node_id,
                    "approval_id": approval.approval_id,
                    "state": "running",
                    "record": reservation_record,
                    "scientific_plan_record": canonical_record(plan),
                    "materialized_workflow_record": canonical_record(
                        materialized_workflow
                    ),
                    "frozen_approval_record": canonical_record(approval),
                    "invocation_record": canonical_record(invocation),
                    "run_state_record": canonical_record(started),
                },
                idempotency_key=(
                    "workflow-launch-reservation:"
                    + run_id
                    + ":"
                    + invocation.node_id
                ),
                existing_events=events,
            )
            return LaunchFenceResultV1(
                status="reserved",
                reservation=reservation,
                execution_receipt=None,
                run_state=started,
            )

    def record_validated_data_edge_binding(
        self,
        *,
        turn_id: str,
        binding: ValidatedDataEdgeBindingV1,
    ) -> RuntimeEvent:
        """Persist one exact producer-output selection after validation."""

        with self._locked_handle(exclusive=True) as handle:
            events = self._read_locked(handle)
            state = replay_events(events)
            record = state.workflow_run_records.get(binding.run_id)
            if record is None:
                raise ContractError(
                    "data edge binding requires a workflow run"
                )
            run_state = workflow_run_state_from_record(record)
            source = next(
                (
                    item
                    for item in run_state.nodes
                    if item.node_id == binding.source_node_id
                ),
                None,
            )
            if source is None or source.state != "validated":
                raise ContractError(
                    "data edge producer is not durably validated"
                )
            if (
                source.execution_receipt_sha256
                != binding.producer_execution_receipt_sha256
                or source.validator_receipt_sha256s
                != binding.producer_validator_receipt_sha256s
                or binding.source_artifact_sha256
                not in source.output_artifact_sha256s
                or binding.selected_artifact_sha256
                not in source.output_artifact_sha256s
            ):
                raise ContractError(
                    "data edge binding differs from durable producer evidence"
                )
            frontier = reconstruct_workflow_frontier(
                state,
                workflow_id=binding.workflow_id,
                run_id=binding.run_id,
            )
            if frontier.plan is None or frontier.approval is None:
                raise ContractError(
                    "data edge binding lacks a replayed frontier"
                )
            edge = next(
                (
                    item
                    for item in frontier.plan.edges
                    if canonical_sha256(item) == binding.scientific_edge_sha256
                ),
                None,
            )
            if edge is None:
                raise ContractError(
                    "data edge binding is absent from the plan"
                )
            rules = tuple(
                item
                for item in frontier.approval.producer_edge_rules
                if item.rule_sha256 == binding.producer_rule_sha256
            )
            if len(rules) != 1:
                raise ContractError("data edge binding lacks a frozen rule")
            existing = tuple(
                item
                for item in frontier.data_edge_bindings
                if item.scientific_edge_sha256
                == binding.scientific_edge_sha256
            )
            if existing:
                if existing != (binding,):
                    raise ContractError(
                        "data edge binding conflicts with replay"
                    )
                return next(
                    event
                    for event in events
                    if event.kind == WORKFLOW_DATA_EDGE_BOUND
                    and event.payload.get("receipt_sha256")
                    == binding.receipt_sha256
                )
            return self._append_locked(
                handle,
                turn_id=turn_id,
                kind=WORKFLOW_DATA_EDGE_BOUND,
                payload={
                    "receipt_sha256": binding.receipt_sha256,
                    "run_id": binding.run_id,
                    "workflow_id": binding.workflow_id,
                    "source_node_id": binding.source_node_id,
                    "target_node_id": binding.target_node_id,
                    "status": binding.status,
                    "record": canonical_record(binding),
                },
                idempotency_key=(
                    "workflow-data-edge:"
                    + binding.run_id
                    + ":"
                    + binding.scientific_edge_sha256
                ),
                existing_events=events,
            )

    def record_program_execution_receipt(
        self,
        *,
        turn_id: str,
        workflow_id: str,
        run_id: str,
        receipt: ProgramExecutionReceiptV1,
    ) -> RuntimeEvent:
        """Persist the full terminal receipt against its exact reservation."""

        with self._locked_handle(exclusive=True) as handle:
            events = self._read_locked(handle)
            frontier = reconstruct_workflow_frontier(
                replay_events(events), workflow_id=workflow_id, run_id=run_id
            )
            reservation = frontier.reservation_for_node(receipt.node_id)
            if reservation is None:
                raise ContractError(
                    "execution receipt requires a durable launch reservation"
                )
            if reservation.invocation_sha256 != receipt.invocation_sha256:
                raise ContractError(
                    "execution receipt differs from launch reservation"
                )
            existing = frontier.receipt_for_node(receipt.node_id)
            if existing is not None and existing != receipt:
                raise ContractError("execution receipt conflicts with replay")
            return self._append_locked(
                handle,
                turn_id=turn_id,
                kind=PROGRAM_EXECUTED,
                payload={
                    "receipt_sha256": receipt.receipt_sha256,
                    "execution_state": receipt.execution_state,
                    "wrapper_exit_status": receipt.wrapper_exit_status,
                    "child_exit_status": receipt.child_exit_status,
                    "engine_receipt_sha256": receipt.engine_receipt_sha256,
                    "engine_complete": receipt.engine_complete,
                    "validated": receipt.validated,
                    "result_validation_receipt_sha256": (
                        receipt.result_validation_receipt_sha256
                    ),
                    "workflow_id": workflow_id,
                    "run_id": run_id,
                    "node_id": receipt.node_id,
                    "record": canonical_record(receipt),
                },
                idempotency_key=(
                    "program-execution-receipt:" + receipt.idempotency_key
                ),
                existing_events=events,
            )

    def transition_workflow_run_node(
        self,
        *,
        turn_id: str,
        run_id: str,
        node_id: str,
        new_state: str,
        timestamp: str,
        plan: ScientificWorkflowPlanV2 | None = None,
        invocation_sha256: str = "",
        execution_receipt_sha256: str = "",
        validator_receipt_sha256s: Iterable[str] = (),
        result_validation_receipt: (
            ProgramResultValidationReceiptV1 | None
        ) = None,
        output_artifact_sha256s: Iterable[str] = (),
        failure_rule_ids: Iterable[str] = (),
    ) -> tuple[RuntimeEvent, WorkflowRunStateV1]:
        """Apply and persist one host-validated node-state transition."""

        with self._locked_handle(exclusive=True) as handle:
            events = self._read_locked(handle)
            state = replay_events(events)
            record = state.workflow_run_records.get(str(run_id))
            if record is None:
                raise ContractError("workflow run has not started")
            current = _workflow_run_state_from_record(record)
            frontier = reconstruct_workflow_frontier(
                state,
                workflow_id=current.workflow_id,
                run_id=run_id,
            )
            transition_plan = plan or frontier.plan
            if new_state == "running":
                if transition_plan is None:
                    raise ContractError(
                        "starting another node requires the scientific plan"
                    )
                if node_id not in derive_ready_node_ids(
                    transition_plan, current, frontier.data_edge_bindings
                ):
                    raise ContractError(
                        "node is not ready in the scientific DAG"
                    )
            if new_state == "validated":
                if result_validation_receipt is None:
                    raise ContractError(
                        "validated transition requires a typed result validation receipt"
                    )
                expected_record = canonical_data(result_validation_receipt)
                durable_records = tuple(
                    event.payload.get("record")
                    for event in events
                    if event.kind == RESULT_VERIFIED
                    and event.payload.get("receipt_sha256")
                    == result_validation_receipt.receipt_sha256
                )
                if durable_records != (expected_record,):
                    raise ContractError(
                        "result validation receipt is not durably resolved"
                    )
            updated = transition_workflow_node(
                current,
                node_id=node_id,
                new_state=new_state,
                plan=transition_plan,
                invocation_sha256=invocation_sha256,
                execution_receipt_sha256=execution_receipt_sha256,
                validator_receipt_sha256s=tuple(validator_receipt_sha256s),
                result_validation_receipt=result_validation_receipt,
                output_artifact_sha256s=tuple(output_artifact_sha256s),
                failure_rule_ids=tuple(failure_rule_ids),
                timestamp=timestamp,
            )
            updated_record = canonical_data(updated)
            event = self._append_locked(
                handle,
                turn_id=turn_id,
                kind=WORKFLOW_NODE_STATE_CHANGED,
                payload={
                    "receipt_sha256": canonical_sha256(updated_record),
                    "run_id": updated.run_id,
                    "workflow_id": updated.workflow_id,
                    "node_id": node_id,
                    "node_state": new_state,
                    "workflow_state": updated.state,
                    "record": updated_record,
                },
                idempotency_key=(
                    "workflow-node-state:"
                    + updated.run_id
                    + ":"
                    + node_id
                    + ":"
                    + new_state
                    + ":"
                    + updated.run_state_sha256
                ),
                existing_events=events,
            )
            return event, updated

    def terminate(
        self,
        *,
        turn_id: str,
        terminal_state: str,
        reason: str,
        required_receipt_sha256s: Iterable[str] = (),
    ) -> RuntimeEvent:
        """Create the sole terminal event; green status is host-derived."""

        with self._locked_handle(exclusive=True) as handle:
            events = self._read_locked(handle)
            state = replay_events(events)
            if state.terminal_state:
                raise ContractError("runtime is already terminal")
            payload: dict[str, Any] = {
                "terminal_state": str(terminal_state),
                "reason": str(reason),
            }
            if terminal_state == "complete":
                required = tuple(sorted(set(required_receipt_sha256s)))
                observed = _observed_receipts(state)
                if not required or not set(required).issubset(observed):
                    raise ContractError(
                        "completion requires receipts already observed in stream"
                    )
                green = tuple(
                    digest
                    for digest in required
                    if _receipt_is_green(events, digest)
                )
                if green != required:
                    raise ContractError("a required completion gate is red")
                latest_preflight_is_required = bool(
                    state.preflight_receipts
                    and state.preflight_receipts[-1] in required
                )
                latest_analysis_is_required = bool(
                    state.analysis_completion_receipts
                    and state.analysis_completion_receipts[-1] in required
                )
                if not (
                    latest_preflight_is_required or latest_analysis_is_required
                ):
                    raise ContractError(
                        "completion requires the latest command preflight or "
                        "analysis completion receipt"
                    )
                gate = canonical_sha256(
                    {
                        "session_id": self.session_id,
                        "turn_id": turn_id,
                        "required_receipt_sha256s": required,
                        "green_receipt_sha256s": green,
                        "previous_hash": (
                            events[-1].event_hash if events else ""
                        ),
                    }
                )
                payload.update(
                    {
                        "required_receipt_sha256s": required,
                        "green_receipt_sha256s": green,
                        "completion_gate_sha256": gate,
                    }
                )
            elif terminal_state == "planned":
                required = tuple(sorted(set(required_receipt_sha256s)))
                observed = _observed_receipts(state)
                if (
                    len(required) != 1
                    or required[0] not in observed
                    or not state.workflow_receipts
                    or state.workflow_receipts[-1] != required[0]
                ):
                    raise ContractError(
                        "planned termination requires the latest workflow draft"
                    )
                payload["plan_receipt_sha256"] = required[0]
            return self._append_locked(
                handle,
                turn_id=turn_id,
                kind=RUNTIME_TERMINATED,
                payload=payload,
                idempotency_key="runtime-terminal",
                existing_events=events,
            )

    def persist_public_transcript(
        self, *, turn_id: str, transcript: Iterable[Mapping[str, Any]]
    ) -> dict[str, Any]:
        """Atomically persist the exact sanitized visible transcript."""

        normalized = tuple(canonical_data(dict(item)) for item in transcript)
        transcript_sha256 = canonical_sha256(normalized)
        body = {
            "schema_version": "chemsmart.public-transcript.v1",
            "session_id": self.session_id,
            "turn_id": str(turn_id),
            "transcript_sha256": transcript_sha256,
            "transcript": normalized,
        }
        encoded = (
            json.dumps(
                body,
                sort_keys=True,
                indent=2,
                ensure_ascii=False,
            )
            + "\n"
        ).encode("utf-8")
        artifact_sha256 = hashlib.sha256(encoded).hexdigest()
        suffix = canonical_sha256({"turn_id": str(turn_id)})[:16]
        destination = self.path.with_name(f"public-transcript-{suffix}.json")
        temporary = destination.with_name(
            "." + destination.name + f".tmp-{os.getpid()}"
        )
        with self._locked_handle(exclusive=True):
            if os.path.lexists(destination):
                existing = _secure_open_text(destination)
                try:
                    existing.seek(0)
                    observed = existing.read().encode("utf-8")
                finally:
                    existing.close()
                if observed != encoded:
                    raise ContractError("public transcript artifact conflicts")
            else:
                flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
                if hasattr(os, "O_NOFOLLOW"):
                    flags |= os.O_NOFOLLOW
                descriptor = os.open(temporary, flags, 0o600)
                try:
                    with os.fdopen(descriptor, "wb") as handle:
                        descriptor = -1
                        handle.write(encoded)
                        handle.flush()
                        os.fsync(handle.fileno())
                    os.replace(temporary, destination)
                    os.chmod(destination, 0o600)
                finally:
                    if descriptor >= 0:
                        os.close(descriptor)
                    if temporary.exists():
                        temporary.unlink()
        return {
            "schema_version": "chemsmart.public-transcript-artifact.v1",
            "artifact_id": "public_transcript." + suffix,
            "path": str(destination),
            "artifact_sha256": artifact_sha256,
            "transcript_sha256": transcript_sha256,
        }

    def _prepare_private_parent(self) -> None:
        self.path.parent.mkdir(parents=True, mode=0o700, exist_ok=True)
        if self.path.parent.is_symlink() or not self.path.parent.is_dir():
            raise ContractError("event-store parent must be a real directory")
        try:
            os.chmod(self.path.parent, 0o700)
        except OSError as exc:
            raise ContractError(
                "event-store parent cannot be made private"
            ) from exc

    @contextmanager
    def _locked_handle(self, *, exclusive: bool) -> Iterator[TextIO]:
        lock_handle = _secure_open_text(self._lock_path)
        try:
            _acquire_lock(lock_handle, exclusive=exclusive)
            data_handle = _secure_open_text(self.path)
            try:
                yield data_handle
            finally:
                data_handle.close()
        finally:
            _release_lock(lock_handle)
            lock_handle.close()

    def _read_locked(self, handle: TextIO) -> tuple[RuntimeEvent, ...]:
        handle.seek(0)
        events = []
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            try:
                raw = json.loads(line)
                events.append(RuntimeEvent.from_dict(raw))
            except (
                KeyError,
                TypeError,
                ValueError,
                json.JSONDecodeError,
            ) as exc:
                raise ContractError(
                    f"invalid runtime event at JSONL line {line_number}"
                ) from exc
        state = replay_events(events)
        if state.session_id and state.session_id != self.session_id:
            raise ContractError("event store contains another session")
        return tuple(events)

    def _append_locked(
        self,
        handle: TextIO,
        *,
        turn_id: str,
        kind: str,
        payload: dict[str, Any],
        idempotency_key: str,
        existing_events: tuple[RuntimeEvent, ...] | None = None,
    ) -> RuntimeEvent:
        events = (
            existing_events
            if existing_events is not None
            else self._read_locked(handle)
        )
        if events and replay_events(events).terminal_state:
            raise ContractError("runtime terminal state is absorbing")
        if idempotency_key:
            identity = canonical_sha256({"kind": kind, "payload": payload})
            for existing in events:
                if existing.idempotency_key != idempotency_key:
                    continue
                existing_identity = canonical_sha256(
                    {"kind": existing.kind, "payload": existing.payload}
                )
                if existing_identity != identity:
                    raise ContractError(
                        "idempotency key conflicts with persisted action"
                    )
                return existing
        event = RuntimeEvent.create(
            sequence=len(events) + 1,
            session_id=self.session_id,
            turn_id=str(turn_id),
            kind=kind,
            payload=payload,
            previous_hash=events[-1].event_hash if events else "",
            idempotency_key=idempotency_key,
        )
        handle.seek(0, os.SEEK_END)
        handle.write(
            json.dumps(
                event.to_dict(),
                sort_keys=True,
                separators=(",", ":"),
                ensure_ascii=False,
            )
            + "\n"
        )
        handle.flush()
        os.fsync(handle.fileno())
        return event


def _secure_open_text(path: Path) -> TextIO:
    if os.path.lexists(path) and path.is_symlink():
        raise ContractError("runtime event paths must not be symbolic links")
    flags = os.O_RDWR | os.O_CREAT | os.O_APPEND
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    if hasattr(os, "O_BINARY"):
        flags |= os.O_BINARY
    try:
        descriptor = os.open(path, flags, 0o600)
    except OSError as exc:
        raise ContractError(
            "runtime event file cannot be opened securely"
        ) from exc
    try:
        metadata = os.fstat(descriptor)
        if not stat.S_ISREG(metadata.st_mode):
            raise ContractError("runtime event path must be a regular file")
        os.fchmod(descriptor, 0o600)
        return os.fdopen(descriptor, "a+", encoding="utf-8", newline="\n")
    except Exception:
        os.close(descriptor)
        raise


def _acquire_lock(handle: TextIO, *, exclusive: bool) -> None:
    if _fcntl is not None:
        operation = _fcntl.LOCK_EX if exclusive else _fcntl.LOCK_SH
        _fcntl.flock(handle.fileno(), operation)
        return
    if _msvcrt is None:  # pragma: no cover - unsupported platform
        raise ContractError(
            "no portable file-lock implementation is available"
        )
    handle.seek(0, os.SEEK_END)
    if handle.tell() == 0:
        handle.write("\0")
        handle.flush()
        os.fsync(handle.fileno())
    handle.seek(0)
    operation = _msvcrt.LK_LOCK if exclusive else _msvcrt.LK_RLCK
    _msvcrt.locking(handle.fileno(), operation, 1)


def _release_lock(handle: TextIO) -> None:
    if _fcntl is not None:
        _fcntl.flock(handle.fileno(), _fcntl.LOCK_UN)
        return
    if _msvcrt is not None:  # pragma: no branch - platform boundary
        handle.seek(0)
        _msvcrt.locking(handle.fileno(), _msvcrt.LK_UNLCK, 1)


def _receipt_is_green(
    events: tuple[RuntimeEvent, ...],
    digest: str,
    _seen: frozenset[str] = frozenset(),
) -> bool:
    if digest in _seen:
        return False
    event = next(
        (
            item
            for item in reversed(events)
            if str(
                item.payload.get("receipt_sha256")
                or item.payload.get("binding_sha256")
                or ""
            )
            == digest
        ),
        None,
    )
    if event is None:
        return False
    value = event.payload
    if event.kind == CAPABILITY_QUERIED:
        return value.get("status") in {"supported", "preview_only"}
    if event.kind == ENVIRONMENT_QUERIED:
        return value.get("status") == "available"
    if event.kind in {PROGRAM_BOUND, ENGINE_BOUND}:
        return value.get("state") in {"resolved", "preview_only"}
    if event.kind == PROJECT_VALIDATED:
        return value.get("status") == "valid"
    if event.kind == PROJECT_PROMOTED:
        return value.get("status") == "validated"
    if event.kind == SCIENTIFIC_DECISION_RECORDED:
        return value.get("status") == "recorded"
    if event.kind == COMMAND_COMPILED:
        return value.get("status") == "compiled"
    if event.kind == COMMAND_INSPECTED:
        return value.get("status") == "valid"
    if event.kind == SAFE_PREVIEWED:
        return (
            value.get("status") == "previewed"
            and value.get("program_validation_status") == "valid"
            and _finding_count(value) == 0
        )
    if event.kind == VALIDATOR_OBSERVED:
        return value.get("status") == "valid" and _finding_count(value) == 0
    if event.kind == PROGRAM_PREFLIGHTED:
        preview_digest = str(value.get("safe_preview_receipt_sha256") or "")
        return (
            value.get("plan_state") == "previewed"
            and _finding_count(value) == 0
            and _receipt_is_green(
                events, preview_digest, _seen | frozenset({digest})
            )
        )
    if event.kind == RESULT_VERIFIED:
        return value.get("status") == "valid"
    if event.kind == RESULT_QUANTITIES_EXTRACTED:
        return value.get("status") == "extracted"
    if event.kind == THERMOCHEMISTRY_DERIVED:
        return value.get("status") == "derived"
    if event.kind == QUANTITY_EXPRESSION_EVALUATED:
        return value.get("status") == "derived"
    if event.kind == SCIENTIFIC_VALIDATION_EVALUATED:
        source_receipts = tuple(value.get("source_receipt_sha256s") or ())
        return bool(
            value.get("status") == "evaluated"
            and source_receipts
            and all(
                _receipt_is_green(events, source, _seen | frozenset({digest}))
                for source in source_receipts
            )
        )
    if event.kind == ANALYSIS_CLAIMS_RECORDED:
        source_receipts = tuple(value.get("source_receipt_sha256s") or ())
        return bool(
            value.get("status") == "recorded"
            and _finding_count(value) == 0
            and source_receipts
            and all(
                _receipt_is_green(events, source, _seen | frozenset({digest}))
                for source in source_receipts
            )
        )
    if event.kind == ANALYSIS_COMPLETION_EVALUATED:
        source_receipts = tuple(value.get("source_receipt_sha256s") or ())
        return bool(
            value.get("status") == "passed"
            and _finding_count(value) == 0
            and source_receipts
            and all(
                _receipt_is_green(
                    events,
                    source,
                    _seen | frozenset({digest}),
                )
                for source in source_receipts
            )
        )
    if event.kind == PROGRAM_EXECUTED:
        return (
            value.get("execution_state") == "validated"
            and value.get("engine_complete") is True
            and value.get("validated") is True
        )
    if event.kind == OPTIMIZED_GEOMETRY_HANDED_OFF:
        return value.get("status") == "validated_handoff"
    if event.kind == WORKFLOW_DATA_EDGE_BOUND:
        return value.get("status") == "validated"
    if event.kind == SUBSTITUTION_ASSESSED:
        return value.get("decision") in {"exact", "approved"}
    if event.kind == PERMISSION_RESOLVED:
        return value.get("decision") in {"auto_allow", "allow_once"}
    if event.kind == SCIENTIFIC_WORKFLOW_MATERIALIZED:
        return value.get("status") in {"previewed", "ready_for_approval"}
    if event.kind == WORKFLOW_APPROVAL_CONSUMED:
        return value.get("status") == "consumed"
    if event.kind == WORKFLOW_NODE_STATE_CHANGED:
        return value.get("node_state") == "validated"
    return False


def _finding_count(payload: Mapping[str, Any]) -> int:
    value = payload.get("critical_finding_count")
    return (
        value if isinstance(value, int) and not isinstance(value, bool) else -1
    )


def _observed_receipts(state: RuntimeState) -> set[str]:
    collections = (
        state.capability_receipts,
        state.environment_receipts,
        state.program_bindings,
        state.engine_bindings,
        state.project_receipts,
        state.project_promotion_receipts,
        state.scientific_decision_receipts,
        state.workflow_receipts,
        state.command_invocations,
        state.command_inspections,
        state.safe_preview_receipts,
        state.validator_receipts,
        state.preflight_receipts,
        state.substitution_receipts,
        state.result_verification_receipts,
        state.result_quantity_receipts,
        state.thermochemistry_receipts,
        state.quantity_expression_receipts,
        state.scientific_validation_receipts,
        state.analysis_claim_receipts,
        state.analysis_completion_receipts,
        state.execution_receipts,
        state.optimized_geometry_handoff_receipts,
        state.permission_receipts,
        state.provider_turn_receipts,
        state.api_attempt_receipts,
        state.materialized_workflow_receipts,
        state.workflow_approval_receipts,
        state.workflow_execution_start_receipts,
        state.workflow_launch_reservation_receipts,
        state.workflow_data_edge_binding_receipts,
        state.workflow_node_state_receipts,
    )
    return {
        digest for collection in collections for digest in collection if digest
    }


def _workflow_run_state_from_record(
    record: Mapping[str, Any],
) -> WorkflowRunStateV1:
    return workflow_run_state_from_record(record)


__all__ = ["RuntimeEventStore"]
