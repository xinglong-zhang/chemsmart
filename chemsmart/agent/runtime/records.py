"""Strict canonical-record loading for the Runtime V2 execution frontier.

The event stream stores JSON-compatible records, not live Python objects.  A
restart must therefore reconstruct the exact typed contracts and re-run every
contract invariant before those records can influence execution.  These
loaders reject missing and extra fields instead of silently defaulting them.
"""

from __future__ import annotations

from dataclasses import MISSING, dataclass, fields
from typing import Any, Mapping

from chemsmart.agent._contracts import (
    AuxiliaryArtifactBindingV1,
    ContractError,
    TrustedArtifactRefV1,
    canonical_data,
    canonical_sha256,
    require_identifier,
    require_sha256,
)
from chemsmart.agent.execution import (
    FrozenMaterializedNodePreviewV1,
    FrozenProducerEdgeRuleV1,
    FrozenWorkflowApprovalV1,
    ProgramExecutionInvocationV1,
    ProgramExecutionReceiptV1,
    ValidatedDataEdgeBindingV1,
    WorkflowNodeRunStateV1,
    WorkflowRunStateV1,
)
from chemsmart.agent.workflows import (
    MaterializedNodeV1,
    MaterializedWorkflowV1,
    ScientificWorkflowEdgeV2,
    ScientificWorkflowNodeV2,
    ScientificWorkflowPlanV2,
)


def _strict_record(
    value: Mapping[str, Any], contract: type, label: str
) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise ContractError(f"{label} must be a canonical mapping")
    record = dict(value)
    expected = {item.name for item in fields(contract)}
    # A field carrying a default is optional by construction, so an older
    # record that predates an additive schema extension still loads.  Unknown
    # fields stay refused: forward compatibility is not permission to invent.
    optional = {
        item.name
        for item in fields(contract)
        if item.default is not MISSING or item.default_factory is not MISSING
    }
    observed = set(record)
    if observed != expected:
        missing = sorted(expected.difference(observed).difference(optional))
        extra = sorted(observed.difference(expected))
        detail = []
        if missing:
            detail.append("missing=" + ",".join(missing))
        if extra:
            detail.append("extra=" + ",".join(extra))
        if detail:
            raise ContractError(f"{label} fields differ ({'; '.join(detail)})")
    return record


def _tuple(value: Any, label: str) -> tuple[Any, ...]:
    if isinstance(value, (str, bytes, bytearray)) or not isinstance(
        value, (list, tuple)
    ):
        raise ContractError(f"{label} must be a canonical array")
    return tuple(value)


def _construct(contract: type, record: dict[str, Any], label: str) -> Any:
    try:
        return contract(**record)
    except ContractError:
        raise
    except (TypeError, ValueError) as exc:
        raise ContractError(f"invalid {label}") from exc


def _auxiliary_artifact_bindings(
    value: Any, label: str
) -> tuple[AuxiliaryArtifactBindingV1, ...]:
    return tuple(
        _construct(
            AuxiliaryArtifactBindingV1,
            _strict_record(
                item, AuxiliaryArtifactBindingV1, "auxiliary artifact binding"
            ),
            "auxiliary artifact binding",
        )
        for item in _tuple(value, label)
    )


def scientific_workflow_plan_from_record(
    value: Mapping[str, Any],
) -> ScientificWorkflowPlanV2:
    record = _strict_record(value, ScientificWorkflowPlanV2, "scientific plan")
    nodes = []
    for item in _tuple(record["nodes"], "scientific plan nodes"):
        node = _strict_record(
            item, ScientificWorkflowNodeV2, "scientific node"
        )
        node["unresolved_fields"] = tuple(
            str(value)
            for value in _tuple(
                node["unresolved_fields"], "scientific unresolved fields"
            )
        )
        nodes.append(
            _construct(ScientificWorkflowNodeV2, node, "scientific node")
        )
    record["nodes"] = tuple(nodes)
    record["edges"] = tuple(
        _construct(
            ScientificWorkflowEdgeV2,
            _strict_record(item, ScientificWorkflowEdgeV2, "scientific edge"),
            "scientific edge",
        )
        for item in _tuple(record["edges"], "scientific plan edges")
    )
    record["complexity_factors"] = tuple(
        str(item)
        for item in _tuple(
            record["complexity_factors"], "scientific complexity factors"
        )
    )
    return _construct(
        ScientificWorkflowPlanV2, record, "scientific workflow plan"
    )


def materialized_workflow_from_record(
    value: Mapping[str, Any],
) -> MaterializedWorkflowV1:
    record = _strict_record(
        value, MaterializedWorkflowV1, "materialized workflow"
    )
    nodes = []
    for item in _tuple(record["nodes"], "materialized workflow nodes"):
        node = _strict_record(item, MaterializedNodeV1, "materialized node")
        if "auxiliary_input_bindings" in node:
            node["auxiliary_input_bindings"] = _auxiliary_artifact_bindings(
                node["auxiliary_input_bindings"],
                "materialized node auxiliary inputs",
            )
        nodes.append(_construct(MaterializedNodeV1, node, "materialized node"))
    record["nodes"] = tuple(nodes)
    record["unresolved_node_ids"] = tuple(
        str(item)
        for item in _tuple(
            record["unresolved_node_ids"], "unresolved workflow nodes"
        )
    )
    return _construct(MaterializedWorkflowV1, record, "materialized workflow")


def frozen_workflow_approval_from_record(
    value: Mapping[str, Any],
) -> FrozenWorkflowApprovalV1:
    if not isinstance(value, Mapping):
        raise ContractError(
            "frozen workflow approval must be a canonical mapping"
        )
    normalized = dict(value)
    legacy_fields = {
        item.name for item in fields(FrozenWorkflowApprovalV1)
    }.difference(
        {
            "materialized_preview_bindings",
            "producer_edge_rules",
            "admission_sha256",
        }
    )
    if set(normalized) == legacy_fields:
        normalized.update(
            materialized_preview_bindings=(),
            producer_edge_rules=(),
            admission_sha256="",
        )
    record = _strict_record(
        normalized, FrozenWorkflowApprovalV1, "frozen workflow approval"
    )
    for name in (
        "environment_identity_sha256s",
        "approved_node_ids",
        "producer_edge_sha256s",
    ):
        record[name] = tuple(
            str(item) for item in _tuple(record[name], f"approval {name}")
        )
    previews = []
    for item in _tuple(
        record["materialized_preview_bindings"],
        "approval materialized_preview_bindings",
    ):
        preview = _strict_record(
            item,
            FrozenMaterializedNodePreviewV1,
            "frozen materialized preview",
        )
        if "auxiliary_input_bindings" in preview:
            preview["auxiliary_input_bindings"] = (
                _auxiliary_artifact_bindings(
                    preview["auxiliary_input_bindings"],
                    "frozen preview auxiliary inputs",
                )
            )
        previews.append(
            _construct(
                FrozenMaterializedNodePreviewV1,
                preview,
                "frozen materialized preview",
            )
        )
    record["materialized_preview_bindings"] = tuple(previews)
    record["producer_edge_rules"] = tuple(
        _construct(
            FrozenProducerEdgeRuleV1,
            _strict_record(
                item, FrozenProducerEdgeRuleV1, "frozen producer edge rule"
            ),
            "frozen producer edge rule",
        )
        for item in _tuple(
            record["producer_edge_rules"], "approval producer_edge_rules"
        )
    )
    return _construct(
        FrozenWorkflowApprovalV1, record, "frozen workflow approval"
    )


def validated_data_edge_binding_from_record(
    value: Mapping[str, Any],
) -> ValidatedDataEdgeBindingV1:
    record = _strict_record(
        value, ValidatedDataEdgeBindingV1, "validated data edge binding"
    )
    record["producer_validator_receipt_sha256s"] = tuple(
        str(item)
        for item in _tuple(
            record["producer_validator_receipt_sha256s"],
            "data edge producer validator receipts",
        )
    )
    return _construct(
        ValidatedDataEdgeBindingV1, record, "validated data edge binding"
    )


def workflow_run_state_from_record(
    value: Mapping[str, Any],
) -> WorkflowRunStateV1:
    record = _strict_record(value, WorkflowRunStateV1, "workflow run state")
    nodes = []
    for value_node in _tuple(record["nodes"], "workflow run nodes"):
        node = _strict_record(
            value_node, WorkflowNodeRunStateV1, "workflow node run state"
        )
        for name in (
            "validator_receipt_sha256s",
            "output_artifact_sha256s",
            "failure_rule_ids",
        ):
            node[name] = tuple(
                str(item)
                for item in _tuple(node[name], f"workflow node {name}")
            )
        nodes.append(
            _construct(WorkflowNodeRunStateV1, node, "workflow node run state")
        )
    record["nodes"] = tuple(nodes)
    return _construct(WorkflowRunStateV1, record, "workflow run state")


def program_execution_invocation_from_record(
    value: Mapping[str, Any],
) -> ProgramExecutionInvocationV1:
    record = _strict_record(
        value, ProgramExecutionInvocationV1, "program execution invocation"
    )
    record["argv"] = tuple(
        str(item) for item in _tuple(record["argv"], "execution argv")
    )
    if "auxiliary_input_bindings" in record:
        record["auxiliary_input_bindings"] = _auxiliary_artifact_bindings(
            record["auxiliary_input_bindings"],
            "execution auxiliary inputs",
        )
    return _construct(
        ProgramExecutionInvocationV1, record, "program execution invocation"
    )


def program_execution_receipt_from_record(
    value: Mapping[str, Any],
) -> ProgramExecutionReceiptV1:
    record = _strict_record(
        value, ProgramExecutionReceiptV1, "program execution receipt"
    )
    artifacts = []
    for value_artifact in _tuple(
        record["output_artifacts"], "execution output artifacts"
    ):
        artifacts.append(
            _construct(
                TrustedArtifactRefV1,
                _strict_record(
                    value_artifact,
                    TrustedArtifactRefV1,
                    "trusted output artifact",
                ),
                "trusted output artifact",
            )
        )
    record["output_artifacts"] = tuple(artifacts)
    for name in ("validator_receipt_sha256s", "findings"):
        record[name] = tuple(
            str(item) for item in _tuple(record[name], f"execution {name}")
        )
    return _construct(
        ProgramExecutionReceiptV1, record, "program execution receipt"
    )


@dataclass(frozen=True)
class WorkflowNodeLaunchReservationV1:
    """Durable pre-launch fence for one exact workflow node invocation."""

    schema_version: str
    reservation_id: str
    run_id: str
    workflow_id: str
    node_id: str
    plan_sha256: str
    materialized_workflow_sha256: str
    approval_id: str
    approval_sha256: str
    invocation_sha256: str
    consumes_approval: bool
    state: str
    reserved_at: str
    reservation_sha256: str
    admission_sha256: str = ""
    data_edge_binding_sha256s: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if (
            self.schema_version
            != "chemsmart.workflow-node-launch-reservation.v1"
        ):
            raise ContractError("unsupported workflow launch reservation")
        for name, value in (
            ("reservation_id", self.reservation_id),
            ("run_id", self.run_id),
            ("workflow_id", self.workflow_id),
            ("node_id", self.node_id),
            ("approval_id", self.approval_id),
        ):
            require_identifier(value, name)
        for name, digest in (
            ("plan_sha256", self.plan_sha256),
            (
                "materialized_workflow_sha256",
                self.materialized_workflow_sha256,
            ),
            ("approval_sha256", self.approval_sha256),
            ("invocation_sha256", self.invocation_sha256),
        ):
            require_sha256(digest, name)
        if self.state != "running":
            raise ContractError("launch reservation must be running")
        if not str(self.reserved_at).strip():
            raise ContractError("launch reservation requires a timestamp")
        legacy_body = {
            "schema_version": self.schema_version,
            "reservation_id": self.reservation_id,
            "run_id": self.run_id,
            "workflow_id": self.workflow_id,
            "node_id": self.node_id,
            "plan_sha256": self.plan_sha256,
            "materialized_workflow_sha256": self.materialized_workflow_sha256,
            "approval_id": self.approval_id,
            "approval_sha256": self.approval_sha256,
            "invocation_sha256": self.invocation_sha256,
            "consumes_approval": self.consumes_approval,
            "state": self.state,
            "reserved_at": self.reserved_at,
        }
        if not self.admission_sha256:
            if self.data_edge_binding_sha256s:
                raise ContractError(
                    "legacy reservation cannot bind data edges"
                )
            if self.reservation_sha256 != canonical_sha256(legacy_body):
                raise ContractError(
                    "workflow launch reservation digest mismatch"
                )
            return
        require_sha256(self.admission_sha256, "admission_sha256")
        if self.data_edge_binding_sha256s != tuple(
            sorted(set(self.data_edge_binding_sha256s))
        ):
            raise ContractError(
                "data edge binding digests must be sorted and unique"
            )
        for digest in self.data_edge_binding_sha256s:
            require_sha256(digest, "data_edge_binding_sha256")
        body = {
            **legacy_body,
            "admission_sha256": self.admission_sha256,
            "data_edge_binding_sha256s": self.data_edge_binding_sha256s,
        }
        if self.reservation_sha256 != canonical_sha256(body):
            raise ContractError("workflow launch reservation digest mismatch")


def build_workflow_node_launch_reservation(
    *,
    run_id: str,
    plan: ScientificWorkflowPlanV2,
    materialized_workflow: MaterializedWorkflowV1,
    approval: FrozenWorkflowApprovalV1,
    invocation: ProgramExecutionInvocationV1,
    data_edge_bindings: tuple[ValidatedDataEdgeBindingV1, ...] = (),
    consumes_approval: bool,
    reserved_at: str,
) -> WorkflowNodeLaunchReservationV1:
    if materialized_workflow.workflow_id != plan.workflow_id:
        raise ContractError("materialized workflow belongs to another plan")
    if materialized_workflow.plan_sha256 != plan.plan_sha256:
        raise ContractError("materialized workflow plan digest differs")
    if approval.workflow_id != plan.workflow_id:
        raise ContractError("frozen approval belongs to another plan")
    if approval.plan_sha256 != plan.plan_sha256:
        raise ContractError("frozen approval plan digest differs")
    if approval.materialized_workflow_sha256 != (
        materialized_workflow.materialized_sha256
    ):
        raise ContractError("frozen approval materialization differs")
    if not approval.execution_admissible:
        raise ContractError(
            "legacy frozen approval is readable but cannot admit execution"
        )
    if invocation.node_id not in approval.approved_node_ids:
        raise ContractError("invocation node is absent from frozen approval")
    if invocation.resource_sha256 != approval.resource_sha256:
        raise ContractError("invocation resources differ from frozen approval")
    # Compare identity with identity. The invocation's *receipt* digest names
    # the exact observation that produced its binding, but it folds in a
    # capability receipt that changes with the active tool overlay, so it
    # never equals the digest a planning session recorded for the same
    # unchanged machine. Comparing it here rejected the very environment the
    # approval named, and no reviewer could supply an accepted digest because
    # it does not exist until execution is already authorised.
    # An approval built from identities is matched by identity; one built
    # before identities existed still holds receipt digests and is matched the
    # old way. Requiring identity outright would make every approval already
    # on disk unexecutable, which is the failure this change exists to remove.
    approved_environments = approval.environment_identity_sha256s
    if (
        invocation.environment_identity_sha256 not in approved_environments
        and invocation.environment_receipt_sha256 not in approved_environments
    ):
        raise ContractError(
            "invocation environment differs from frozen approval"
        )
    preview_binding = approval.preview_binding(invocation.node_id)
    producer_rules = approval.producer_rules_for(invocation.node_id)
    matched_data_bindings = tuple(
        sorted(
            (
                item
                for item in data_edge_bindings
                if item.target_node_id == invocation.node_id
            ),
            key=lambda item: item.consumer_input_id,
        )
    )
    if preview_binding is not None:
        if producer_rules or matched_data_bindings:
            raise ContractError(
                "frozen preview node cannot also consume a future data edge"
            )
        # Content and identity are comparable across sessions; the two raw
        # digests are not.  ``invocation_sha256`` covers argv, whose absolute
        # paths point into the planning session's own timestamped run
        # directory, and the two sides digest different record types over
        # different argv besides.  ``environment_receipt_sha256`` folds in a
        # capability receipt that moves with the tool overlay.  Comparing
        # either made an approval unexecutable by any later process.
        # ``build_program_execution_invocation`` has already matched the
        # invocation's scientific identity to the approved *node*.  The plan
        # identity below is the aggregate identity of every starting geometry
        # in the workflow, so comparing a per-node identity with it rejects
        # every valid multi-geometry workflow.
        if (
            invocation.input_sha256 != preview_binding.input_artifact_sha256
            or invocation.project_sha256
            != preview_binding.project_artifact_sha256
        ):
            raise ContractError("invocation differs from exact frozen preview")
        frozen_identity = preview_binding.invocation_identity_sha256
        if frozen_identity:
            if invocation.invocation_identity_sha256 != frozen_identity:
                raise ContractError(
                    "invocation differs from exact frozen preview"
                )
        elif invocation.invocation_sha256 != preview_binding.invocation_sha256:
            # A preview frozen before invocation identity existed can only be
            # matched the old way, which is exact but session-scoped.
            raise ContractError("invocation differs from exact frozen preview")
    else:
        if not producer_rules or len(matched_data_bindings) != len(
            producer_rules
        ):
            raise ContractError(
                "future producer-dependent invocation lacks exact data bindings"
            )
        if tuple(
            item.producer_rule_sha256 for item in matched_data_bindings
        ) != tuple(item.rule_sha256 for item in producer_rules):
            raise ContractError(
                "data bindings differ from frozen producer rules"
            )
        selected = tuple(
            item
            for item in matched_data_bindings
            if item.selected_artifact_sha256 == invocation.input_sha256
            and item.consumer_scientific_identity_sha256
            == invocation.scientific_identity_sha256
        )
        if len(selected) != 1:
            raise ContractError(
                "invocation input/state differs from selected producer artifact"
            )
    node = next(
        (item for item in plan.nodes if item.node_id == invocation.node_id),
        None,
    )
    if node is None:
        raise ContractError("invocation node is absent from scientific plan")
    if (node.program, node.engine, node.stage) != (
        invocation.program,
        invocation.engine,
        invocation.jobtype,
    ):
        raise ContractError("invocation semantics differ from scientific plan")
    reservation_id = "reservation." + invocation.idempotency_key
    body = {
        "schema_version": "chemsmart.workflow-node-launch-reservation.v1",
        "reservation_id": reservation_id,
        "run_id": require_identifier(run_id, "run_id"),
        "workflow_id": plan.workflow_id,
        "node_id": invocation.node_id,
        "plan_sha256": plan.plan_sha256,
        "materialized_workflow_sha256": (
            materialized_workflow.materialized_sha256
        ),
        "approval_id": approval.approval_id,
        "approval_sha256": approval.approval_sha256,
        "invocation_sha256": invocation.invocation_sha256,
        "consumes_approval": bool(consumes_approval),
        "state": "running",
        "reserved_at": str(reserved_at).strip(),
        "admission_sha256": approval.admission_sha256,
        "data_edge_binding_sha256s": tuple(
            sorted(item.receipt_sha256 for item in matched_data_bindings)
        ),
    }
    return WorkflowNodeLaunchReservationV1(
        **body, reservation_sha256=canonical_sha256(body)
    )


def workflow_node_launch_reservation_from_record(
    value: Mapping[str, Any],
) -> WorkflowNodeLaunchReservationV1:
    if not isinstance(value, Mapping):
        raise ContractError("workflow launch reservation must be a mapping")
    normalized = dict(value)
    legacy_fields = {
        item.name for item in fields(WorkflowNodeLaunchReservationV1)
    }.difference({"admission_sha256", "data_edge_binding_sha256s"})
    if set(normalized) == legacy_fields:
        normalized.update(admission_sha256="", data_edge_binding_sha256s=())
    record = _strict_record(
        normalized,
        WorkflowNodeLaunchReservationV1,
        "workflow launch reservation",
    )
    record["data_edge_binding_sha256s"] = tuple(
        str(item)
        for item in _tuple(
            record["data_edge_binding_sha256s"],
            "reservation data edge bindings",
        )
    )
    return _construct(
        WorkflowNodeLaunchReservationV1,
        record,
        "workflow launch reservation",
    )


@dataclass(frozen=True)
class ReconstructedWorkflowFrontierV1:
    """Typed execution frontier reconstructed solely from durable records."""

    workflow_id: str
    run_id: str
    plan: ScientificWorkflowPlanV2 | None
    materialized_workflow: MaterializedWorkflowV1 | None
    approval: FrozenWorkflowApprovalV1 | None
    run_state: WorkflowRunStateV1 | None
    reservations: tuple[WorkflowNodeLaunchReservationV1, ...]
    invocations: tuple[ProgramExecutionInvocationV1, ...]
    execution_receipts: tuple[ProgramExecutionReceiptV1, ...]
    data_edge_bindings: tuple[ValidatedDataEdgeBindingV1, ...]
    legacy_incomplete: bool

    def receipt_for_node(
        self, node_id: str
    ) -> ProgramExecutionReceiptV1 | None:
        matches = tuple(
            item for item in self.execution_receipts if item.node_id == node_id
        )
        if len(matches) > 1:
            raise ContractError(
                "multiple execution receipts exist for one node"
            )
        return matches[0] if matches else None

    def reservation_for_node(
        self, node_id: str
    ) -> WorkflowNodeLaunchReservationV1 | None:
        matches = tuple(
            item for item in self.reservations if item.node_id == node_id
        )
        if len(matches) > 1:
            raise ContractError(
                "multiple launch reservations exist for one node"
            )
        return matches[0] if matches else None


@dataclass(frozen=True)
class LaunchFenceResultV1:
    """Result of consulting and, when safe, advancing the durable fence."""

    status: str
    reservation: WorkflowNodeLaunchReservationV1 | None
    execution_receipt: ProgramExecutionReceiptV1 | None
    run_state: WorkflowRunStateV1 | None

    def __post_init__(self) -> None:
        if self.status not in {"reserved", "terminal_replay"}:
            raise ContractError("unsupported launch fence result")
        if self.status == "reserved" and (
            self.reservation is None
            or self.execution_receipt is not None
            or self.run_state is None
        ):
            raise ContractError("reserved launch fence result is incomplete")
        if self.status == "terminal_replay" and (
            self.reservation is not None or self.execution_receipt is None
        ):
            raise ContractError("terminal replay result is incomplete")


def reconstruct_workflow_frontier(
    state: Any, *, workflow_id: str, run_id: str
) -> ReconstructedWorkflowFrontierV1:
    """Rebuild one frontier and classify old digest-only streams as incomplete."""

    workflow_id = require_identifier(workflow_id, "workflow_id")
    run_id = require_identifier(run_id, "run_id")
    run_record = getattr(state, "workflow_run_records", {}).get(run_id)
    run_state = (
        workflow_run_state_from_record(run_record)
        if isinstance(run_record, Mapping)
        else None
    )
    plan_digest = run_state.plan_sha256 if run_state is not None else ""

    approval_records = getattr(state, "frozen_workflow_approval_records", {})
    approvals = tuple(
        frozen_workflow_approval_from_record(record)
        for record in approval_records.values()
        if isinstance(record, Mapping)
        and str(record.get("workflow_id") or "") == workflow_id
    )
    if run_state is not None:
        approvals = tuple(
            item
            for item in approvals
            if item.approval_sha256 == run_state.approval_sha256
        )
    if len(approvals) > 1:
        raise ContractError("multiple frozen approvals match one workflow run")
    approval = approvals[0] if approvals else None
    if approval is not None:
        plan_digest = approval.plan_sha256

    plan_record = getattr(state, "scientific_workflow_plan_records", {}).get(
        plan_digest
    )
    plan = (
        scientific_workflow_plan_from_record(plan_record)
        if isinstance(plan_record, Mapping)
        else None
    )
    materialized = None
    if approval is not None:
        materialized_record = getattr(
            state, "materialized_workflow_records", {}
        ).get(approval.materialized_workflow_sha256)
        if isinstance(materialized_record, Mapping):
            materialized = materialized_workflow_from_record(
                materialized_record
            )

    reservations = tuple(
        sorted(
            (
                workflow_node_launch_reservation_from_record(record)
                for record in getattr(
                    state, "workflow_launch_reservation_records", {}
                ).values()
                if isinstance(record, Mapping)
                and str(record.get("workflow_id") or "") == workflow_id
                and str(record.get("run_id") or "") == run_id
            ),
            key=lambda item: item.node_id,
        )
    )
    invocations = tuple(
        sorted(
            (
                program_execution_invocation_from_record(record)
                for record in getattr(
                    state, "program_execution_invocation_records", {}
                ).values()
                if isinstance(record, Mapping)
                and any(
                    reservation.invocation_sha256
                    == str(record.get("invocation_sha256") or "")
                    for reservation in reservations
                )
            ),
            key=lambda item: item.node_id,
        )
    )
    invocation_digests = {item.invocation_sha256 for item in invocations}
    receipts = tuple(
        sorted(
            (
                program_execution_receipt_from_record(record)
                for record in getattr(
                    state, "program_execution_receipt_records", {}
                ).values()
                if isinstance(record, Mapping)
                and str(record.get("invocation_sha256") or "")
                in invocation_digests
            ),
            key=lambda item: item.node_id,
        )
    )
    data_edge_bindings = tuple(
        sorted(
            (
                validated_data_edge_binding_from_record(record)
                for record in getattr(
                    state, "validated_data_edge_binding_records", {}
                ).values()
                if isinstance(record, Mapping)
                and str(record.get("workflow_id") or "") == workflow_id
                and str(record.get("run_id") or "") == run_id
            ),
            key=lambda item: (
                item.target_node_id,
                item.consumer_input_id,
            ),
        )
    )
    has_legacy_execution = bool(
        set(
            getattr(state, "legacy_incomplete_execution_node_ids", ())
        ).intersection(
            {node.node_id for node in plan.nodes}
            if plan is not None
            else set()
        )
    )
    legacy_incomplete = bool(
        has_legacy_execution
        or (run_state is not None and not reservations)
        or (
            approval is not None
            and approval.approval_id
            in set(getattr(state, "consumed_workflow_approval_ids", ()))
            and not reservations
        )
        or (approval is not None and not approval.execution_admissible)
        or any(not item.admission_sha256 for item in reservations)
    )
    return ReconstructedWorkflowFrontierV1(
        workflow_id=workflow_id,
        run_id=run_id,
        plan=plan,
        materialized_workflow=materialized,
        approval=approval,
        run_state=run_state,
        reservations=reservations,
        invocations=invocations,
        execution_receipts=receipts,
        data_edge_bindings=data_edge_bindings,
        legacy_incomplete=legacy_incomplete,
    )


def canonical_record(value: Any) -> dict[str, Any]:
    record = canonical_data(value)
    if not isinstance(record, dict):  # pragma: no cover - defensive boundary
        raise ContractError("runtime contract did not serialize to a mapping")
    return record


__all__ = [
    "LaunchFenceResultV1",
    "ReconstructedWorkflowFrontierV1",
    "WorkflowNodeLaunchReservationV1",
    "build_workflow_node_launch_reservation",
    "canonical_record",
    "frozen_workflow_approval_from_record",
    "materialized_workflow_from_record",
    "program_execution_invocation_from_record",
    "program_execution_receipt_from_record",
    "reconstruct_workflow_frontier",
    "scientific_workflow_plan_from_record",
    "workflow_node_launch_reservation_from_record",
    "workflow_run_state_from_record",
    "validated_data_edge_binding_from_record",
]
