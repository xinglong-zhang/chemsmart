"""Deterministic predecessor context for scientific workflow tasks.

The model does not decide which prior records are relevant.  This module
derives an immutable context packet from the host-owned scientific DAG and an
experiment-bound policy.  It supports three deliberately different arms:

``P0``
    The target node and workflow-scoped evidence only.  No predecessor packet.
``P1``
    Every public typed record at or before the target in topological order.
``P2``
    The exact backward dependency slice, including all fan-in branches.

Only references to public typed records are carried here.  Provider-private
reasoning, raw transcripts, paths, and native input text have no representation
in these contracts.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping, Sequence

from chemsmart.agent._contracts import (
    ContractError,
    canonical_data,
    canonical_json,
    canonical_sha256,
    require_identifier,
    require_sha256,
)
from chemsmart.agent.workflows import (
    ScientificWorkflowEdgeV2,
    ScientificWorkflowPlanV2,
)


CONTEXT_MODES = ("dependency_projected", "full_history", "none")
CONTEXT_ARM_MODES = {
    "p0": "none",
    "p1": "full_history",
    "p2": "dependency_projected",
}
CONTEXT_SELECTOR_VERSION = "chemsmart.task-dependency-context-selector.v1"


def _sorted_unique(values: Sequence[str], field_name: str) -> tuple[str, ...]:
    normalized = tuple(str(value) for value in values)
    if normalized != tuple(sorted(set(normalized))):
        raise ContractError(f"{field_name} must be sorted and unique")
    return normalized


def _unique_identifiers(
    values: Sequence[str], field_name: str
) -> tuple[str, ...]:
    normalized = tuple(str(value) for value in values)
    if len(normalized) != len(set(normalized)):
        raise ContractError(f"{field_name} must be unique")
    for value in normalized:
        require_identifier(value, field_name)
    return normalized


def _require_optional_identifier(value: str, field_name: str) -> str:
    if value:
        require_identifier(value, field_name)
    return value


@dataclass(frozen=True)
class TaskDependencyContextPolicyV1:
    """Experiment-bound rule for selecting predecessor context."""

    schema_version: str
    policy_id: str
    arm_id: str
    mode: str
    selector_version: str
    record_classes: tuple[str, ...]
    max_public_bytes: int
    policy_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.task-dependency-context-policy.v1":
            raise ContractError("unsupported task dependency context policy")
        require_identifier(self.policy_id, "policy_id")
        require_identifier(self.arm_id, "arm_id")
        if self.arm_id not in CONTEXT_ARM_MODES:
            raise ContractError("context arm must be p0, p1, or p2")
        if self.mode != CONTEXT_ARM_MODES[self.arm_id]:
            raise ContractError("context arm and mode disagree")
        if self.selector_version != CONTEXT_SELECTOR_VERSION:
            raise ContractError("unsupported dependency context selector")
        _sorted_unique(self.record_classes, "record_classes")
        if not self.record_classes:
            raise ContractError("context policy requires a record class")
        for record_class in self.record_classes:
            require_identifier(record_class, "record_class")
        if self.max_public_bytes < 1:
            raise ContractError("max_public_bytes must be positive")
        body = {
            "schema_version": self.schema_version,
            "policy_id": self.policy_id,
            "arm_id": self.arm_id,
            "mode": self.mode,
            "selector_version": self.selector_version,
            "record_classes": self.record_classes,
            "max_public_bytes": self.max_public_bytes,
        }
        if self.policy_sha256 != canonical_sha256(body):
            raise ContractError("task dependency context policy digest mismatch")


def build_task_dependency_context_policy(
    *,
    policy_id: str,
    arm_id: str,
    record_classes: Sequence[str],
    max_public_bytes: int,
) -> TaskDependencyContextPolicyV1:
    """Build a canonical P0, P1, or P2 context policy."""

    normalized_arm = require_identifier(arm_id, "arm_id")
    if normalized_arm not in CONTEXT_ARM_MODES:
        raise ContractError("context arm must be p0, p1, or p2")
    body = {
        "schema_version": "chemsmart.task-dependency-context-policy.v1",
        "policy_id": require_identifier(policy_id, "policy_id"),
        "arm_id": normalized_arm,
        "mode": CONTEXT_ARM_MODES[normalized_arm],
        "selector_version": CONTEXT_SELECTOR_VERSION,
        "record_classes": tuple(sorted(set(record_classes))),
        "max_public_bytes": int(max_public_bytes),
    }
    return TaskDependencyContextPolicyV1(
        **body, policy_sha256=canonical_sha256(body)
    )


@dataclass(frozen=True)
class PredecessorEvidenceRefV1:
    """Public typed record eligible for deterministic context selection.

    An empty ``node_id`` denotes workflow-scoped evidence.  Node-scoped
    evidence must name an actual plan node when selected.
    """

    schema_version: str
    record_id: str
    record_class: str
    record_sha256: str
    node_id: str
    size_bytes: int
    evidence_ref_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.predecessor-evidence-ref.v1":
            raise ContractError("unsupported predecessor evidence reference")
        require_identifier(self.record_id, "record_id")
        require_identifier(self.record_class, "record_class")
        require_sha256(self.record_sha256, "record_sha256")
        _require_optional_identifier(self.node_id, "node_id")
        if self.size_bytes < 0:
            raise ContractError("evidence size_bytes must be non-negative")
        body = {
            "schema_version": self.schema_version,
            "record_id": self.record_id,
            "record_class": self.record_class,
            "record_sha256": self.record_sha256,
            "node_id": self.node_id,
            "size_bytes": self.size_bytes,
        }
        if self.evidence_ref_sha256 != canonical_sha256(body):
            raise ContractError("predecessor evidence reference digest mismatch")


def build_predecessor_evidence_ref(
    *,
    record_id: str,
    record_class: str,
    record_sha256: str,
    size_bytes: int,
    node_id: str = "",
) -> PredecessorEvidenceRefV1:
    body = {
        "schema_version": "chemsmart.predecessor-evidence-ref.v1",
        "record_id": require_identifier(record_id, "record_id"),
        "record_class": require_identifier(record_class, "record_class"),
        "record_sha256": require_sha256(record_sha256, "record_sha256"),
        "node_id": _require_optional_identifier(node_id, "node_id"),
        "size_bytes": int(size_bytes),
    }
    return PredecessorEvidenceRefV1(
        **body, evidence_ref_sha256=canonical_sha256(body)
    )


def bind_selected_public_records(
    *,
    context: TaskDependencyContextV1,
    records: Mapping[str, Mapping[str, Any]],
) -> tuple[dict[str, Any], ...]:
    """Bind selected references to their exact canonical public payloads.

    ``TaskDependencyContextV1`` deliberately selects references rather than
    carrying record bodies.  The live-session composition root calls this
    function before provider work so an omitted, injected, or changed record
    fails closed.  The returned entries follow the context's canonical
    evidence-reference order and contain detached canonical payloads.
    """

    expected_ids = tuple(item.record_id for item in context.evidence_refs)
    if len(expected_ids) != len(set(expected_ids)):
        raise ContractError(
            "dependency context record IDs must be unique for payload binding"
        )

    supplied_ids = tuple(str(record_id) for record_id in records)
    if len(supplied_ids) != len(set(supplied_ids)):
        raise ContractError("dependency public record IDs must be unique")
    if set(supplied_ids) != set(expected_ids):
        missing = sorted(set(expected_ids).difference(supplied_ids))
        unexpected = sorted(set(supplied_ids).difference(expected_ids))
        raise ContractError(
            "dependency public records must exactly match selected records; "
            f"missing={missing}, unexpected={unexpected}"
        )

    bound: list[dict[str, Any]] = []
    for evidence_ref in context.evidence_refs:
        payload = records[evidence_ref.record_id]
        if not isinstance(payload, Mapping):
            raise ContractError("dependency public record must be a mapping")
        canonical_payload = canonical_data(payload)
        observed_sha256 = canonical_sha256(canonical_payload)
        if observed_sha256 != evidence_ref.record_sha256:
            raise ContractError(
                "dependency public record digest mismatch for "
                f"{evidence_ref.record_id}"
            )
        observed_size = len(canonical_json(canonical_payload).encode("utf-8"))
        if observed_size != evidence_ref.size_bytes:
            raise ContractError(
                "dependency public record size mismatch for "
                f"{evidence_ref.record_id}"
            )
        bound.append(
            {
                "record_id": evidence_ref.record_id,
                "record_class": evidence_ref.record_class,
                "node_id": evidence_ref.node_id,
                "record_sha256": evidence_ref.record_sha256,
                "public_record": canonical_payload,
            }
        )
    return tuple(bound)


@dataclass(frozen=True)
class EvidenceExclusionV1:
    """Why one public typed record was omitted from a context packet."""

    evidence_ref_sha256: str
    reason: str

    def __post_init__(self) -> None:
        require_sha256(self.evidence_ref_sha256, "evidence_ref_sha256")
        if self.reason not in {
            "context_budget_exceeded",
            "node_not_selected",
            "record_class_not_selected",
        }:
            raise ContractError("unsupported evidence exclusion reason")


@dataclass(frozen=True)
class TaskDependencyContextV1:
    """Canonical public context selected for one workflow target."""

    schema_version: str
    workflow_id: str
    plan_sha256: str
    target_node_id: str
    policy_sha256: str
    mode: str
    direct_predecessor_node_ids: tuple[str, ...]
    transitive_ancestor_node_ids: tuple[str, ...]
    included_node_ids: tuple[str, ...]
    excluded_node_ids: tuple[str, ...]
    non_dependency_node_ids: tuple[str, ...]
    selected_edges: tuple[ScientificWorkflowEdgeV2, ...]
    evidence_refs: tuple[PredecessorEvidenceRefV1, ...]
    context_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.task-dependency-context.v1":
            raise ContractError("unsupported task dependency context")
        require_identifier(self.workflow_id, "workflow_id")
        require_sha256(self.plan_sha256, "plan_sha256")
        require_identifier(self.target_node_id, "target_node_id")
        require_sha256(self.policy_sha256, "policy_sha256")
        if self.mode not in CONTEXT_MODES:
            raise ContractError("unsupported task dependency context mode")
        for name, node_ids in (
            ("direct_predecessor_node_ids", self.direct_predecessor_node_ids),
            ("transitive_ancestor_node_ids", self.transitive_ancestor_node_ids),
            ("included_node_ids", self.included_node_ids),
            ("excluded_node_ids", self.excluded_node_ids),
            ("non_dependency_node_ids", self.non_dependency_node_ids),
        ):
            _unique_identifiers(node_ids, name)
        included = set(self.included_node_ids)
        excluded = set(self.excluded_node_ids)
        if self.target_node_id not in included:
            raise ContractError("dependency context must include its target node")
        if included.intersection(excluded):
            raise ContractError("included and excluded context nodes overlap")
        ancestors = set(self.direct_predecessor_node_ids).union(
            self.transitive_ancestor_node_ids
        )
        if not ancestors.issubset(included):
            raise ContractError("selected ancestors must be included nodes")
        if not set(self.non_dependency_node_ids).issubset(included):
            raise ContractError("non-dependency nodes must be included nodes")
        if set(self.non_dependency_node_ids).intersection(
            ancestors.union({self.target_node_id})
        ):
            raise ContractError("dependency and non-dependency nodes overlap")
        if self.mode == "none" and (
            self.direct_predecessor_node_ids
            or self.transitive_ancestor_node_ids
            or self.non_dependency_node_ids
            or self.selected_edges
        ):
            raise ContractError("P0 context cannot expose predecessors")
        if self.mode == "dependency_projected" and self.non_dependency_node_ids:
            raise ContractError("P2 context cannot include unrelated nodes")
        edge_ids = tuple(edge.edge_id for edge in self.selected_edges)
        if edge_ids != tuple(sorted(set(edge_ids))):
            raise ContractError("selected context edges must be sorted and unique")
        positions = {
            node_id: index for index, node_id in enumerate(self.included_node_ids)
        }
        for edge in self.selected_edges:
            if (
                edge.source_node_id not in included
                or edge.target_node_id not in included
            ):
                raise ContractError("selected edge leaves the included context")
            if positions[edge.source_node_id] >= positions[edge.target_node_id]:
                raise ContractError("included context must be topologically ordered")
        evidence_digests = tuple(
            evidence.evidence_ref_sha256 for evidence in self.evidence_refs
        )
        if evidence_digests != tuple(sorted(set(evidence_digests))):
            raise ContractError("context evidence refs must be sorted and unique")
        record_digests = tuple(
            evidence.record_sha256 for evidence in self.evidence_refs
        )
        if len(record_digests) != len(set(record_digests)):
            raise ContractError("context cannot repeat a public record")
        for evidence in self.evidence_refs:
            if evidence.node_id and evidence.node_id not in included:
                raise ContractError("context evidence belongs to an excluded node")
        body = {
            "schema_version": self.schema_version,
            "workflow_id": self.workflow_id,
            "plan_sha256": self.plan_sha256,
            "target_node_id": self.target_node_id,
            "policy_sha256": self.policy_sha256,
            "mode": self.mode,
            "direct_predecessor_node_ids": self.direct_predecessor_node_ids,
            "transitive_ancestor_node_ids": self.transitive_ancestor_node_ids,
            "included_node_ids": self.included_node_ids,
            "excluded_node_ids": self.excluded_node_ids,
            "non_dependency_node_ids": self.non_dependency_node_ids,
            "selected_edges": self.selected_edges,
            "evidence_refs": self.evidence_refs,
        }
        if self.context_sha256 != canonical_sha256(body):
            raise ContractError("task dependency context digest mismatch")


@dataclass(frozen=True)
class ContextSelectionReceiptV1:
    """Observable result of applying a dependency-context policy."""

    schema_version: str
    selection_id: str
    workflow_id: str
    plan_sha256: str
    target_node_id: str
    policy_sha256: str
    context_sha256: str
    candidate_evidence_ref_sha256s: tuple[str, ...]
    included_evidence_ref_sha256s: tuple[str, ...]
    exclusions: tuple[EvidenceExclusionV1, ...]
    selected_public_bytes: int
    status: str
    reason: str
    receipt_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.context-selection-receipt.v1":
            raise ContractError("unsupported context selection receipt")
        for name, value in (
            ("selection_id", self.selection_id),
            ("workflow_id", self.workflow_id),
            ("target_node_id", self.target_node_id),
        ):
            require_identifier(value, name)
        for name, digest in (
            ("plan_sha256", self.plan_sha256),
            ("policy_sha256", self.policy_sha256),
        ):
            require_sha256(digest, name)
        if self.context_sha256:
            require_sha256(self.context_sha256, "context_sha256")
        for name, values in (
            (
                "candidate_evidence_ref_sha256s",
                self.candidate_evidence_ref_sha256s,
            ),
            (
                "included_evidence_ref_sha256s",
                self.included_evidence_ref_sha256s,
            ),
        ):
            _sorted_unique(values, name)
            for digest in values:
                require_sha256(digest, name)
        exclusion_digests = tuple(
            item.evidence_ref_sha256 for item in self.exclusions
        )
        if exclusion_digests != tuple(sorted(set(exclusion_digests))):
            raise ContractError("context exclusions must be sorted and unique")
        if set(self.included_evidence_ref_sha256s).intersection(
            exclusion_digests
        ):
            raise ContractError("included and excluded context evidence overlap")
        if self.selected_public_bytes < 0:
            raise ContractError("selected_public_bytes must be non-negative")
        if self.status not in {"blocked_context_budget", "selected"}:
            raise ContractError("unsupported context selection status")
        if self.status == "selected":
            require_sha256(self.context_sha256, "context_sha256")
            if self.reason:
                raise ContractError("selected context cannot carry a block reason")
        else:
            if self.context_sha256:
                raise ContractError("blocked context cannot claim a context digest")
            if not self.reason.strip():
                raise ContractError("blocked context must state a reason")
            if self.included_evidence_ref_sha256s:
                raise ContractError("blocked context cannot include records")
        body = {
            "schema_version": self.schema_version,
            "selection_id": self.selection_id,
            "workflow_id": self.workflow_id,
            "plan_sha256": self.plan_sha256,
            "target_node_id": self.target_node_id,
            "policy_sha256": self.policy_sha256,
            "context_sha256": self.context_sha256,
            "candidate_evidence_ref_sha256s": (
                self.candidate_evidence_ref_sha256s
            ),
            "included_evidence_ref_sha256s": (
                self.included_evidence_ref_sha256s
            ),
            "exclusions": self.exclusions,
            "selected_public_bytes": self.selected_public_bytes,
            "status": self.status,
            "reason": self.reason,
        }
        if self.receipt_sha256 != canonical_sha256(body):
            raise ContractError("context selection receipt digest mismatch")


def build_dependency_context_public_projection(
    *,
    context: TaskDependencyContextV1,
    selection_receipt: ContextSelectionReceiptV1,
    records: Mapping[str, Mapping[str, Any]],
) -> dict[str, Any]:
    """Build the exact dependency-context object shown to a coordinator."""

    for field_name in (
        "workflow_id",
        "plan_sha256",
        "target_node_id",
        "policy_sha256",
        "context_sha256",
    ):
        if getattr(context, field_name) != getattr(
            selection_receipt, field_name
        ):
            raise ContractError(
                f"dependency context selection receipt {field_name} mismatch"
            )
    expected_refs = tuple(
        item.evidence_ref_sha256 for item in context.evidence_refs
    )
    if selection_receipt.status != "selected" or (
        selection_receipt.included_evidence_ref_sha256s != expected_refs
    ):
        raise ContractError(
            "dependency context selection receipt does not bind selected records"
        )
    return {
        "task_dependency_context": canonical_data(context),
        "context_selection_receipt": canonical_data(selection_receipt),
        "selected_predecessor_records": bind_selected_public_records(
            context=context,
            records=records,
        ),
    }


@dataclass(frozen=True)
class ContextManifestV2:
    """Least-privilege task context with exact selected public records."""

    schema_version: str
    manifest_id: str
    workflow_id: str
    task_spec_sha256: str
    scientific_identity_sha256: str
    target_node_id: str
    dependency_context_sha256: str
    selection_receipt_sha256: str
    included_record_sha256s: tuple[str, ...]
    source_sha256s: tuple[str, ...]
    artifact_sha256s: tuple[str, ...]
    tool_schema_sha256: str
    allowed_tools: tuple[str, ...]
    token_budget: int
    tool_call_budget: int
    wall_time_seconds: int
    manifest_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.context-manifest.v2":
            raise ContractError("unsupported context manifest schema")
        for name, value in (
            ("manifest_id", self.manifest_id),
            ("workflow_id", self.workflow_id),
            ("target_node_id", self.target_node_id),
        ):
            require_identifier(value, name)
        for name, digest in (
            ("task_spec_sha256", self.task_spec_sha256),
            ("scientific_identity_sha256", self.scientific_identity_sha256),
            ("dependency_context_sha256", self.dependency_context_sha256),
            ("selection_receipt_sha256", self.selection_receipt_sha256),
            ("tool_schema_sha256", self.tool_schema_sha256),
        ):
            require_sha256(digest, name)
        for name, values in (
            ("included_record_sha256s", self.included_record_sha256s),
            ("source_sha256s", self.source_sha256s),
            ("artifact_sha256s", self.artifact_sha256s),
        ):
            _sorted_unique(values, name)
            for digest in values:
                require_sha256(digest, name)
        _sorted_unique(self.allowed_tools, "allowed_tools")
        for tool in self.allowed_tools:
            require_identifier(tool, "allowed_tool")
        if min(self.token_budget, self.tool_call_budget, self.wall_time_seconds) < 1:
            raise ContractError("context manifest budgets must be positive")
        body = {
            "schema_version": self.schema_version,
            "manifest_id": self.manifest_id,
            "workflow_id": self.workflow_id,
            "task_spec_sha256": self.task_spec_sha256,
            "scientific_identity_sha256": self.scientific_identity_sha256,
            "target_node_id": self.target_node_id,
            "dependency_context_sha256": self.dependency_context_sha256,
            "selection_receipt_sha256": self.selection_receipt_sha256,
            "included_record_sha256s": self.included_record_sha256s,
            "source_sha256s": self.source_sha256s,
            "artifact_sha256s": self.artifact_sha256s,
            "tool_schema_sha256": self.tool_schema_sha256,
            "allowed_tools": self.allowed_tools,
            "token_budget": self.token_budget,
            "tool_call_budget": self.tool_call_budget,
            "wall_time_seconds": self.wall_time_seconds,
        }
        if self.manifest_sha256 != canonical_sha256(body):
            raise ContractError("context manifest digest mismatch")


@dataclass(frozen=True)
class SpecialistTaskPacketV2:
    """Bounded specialist task tied to an exact dependency projection."""

    schema_version: str
    packet_id: str
    workflow_id: str
    role: str
    target_node_ids: tuple[str, ...]
    context_manifest_sha256: str
    context_selection_receipt_sha256: str
    input_record_sha256s: tuple[str, ...]
    expected_output_schema: str
    owner: str
    merge_key: str
    packet_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.specialist-task-packet.v2":
            raise ContractError("unsupported specialist task packet schema")
        for name, value in (
            ("packet_id", self.packet_id),
            ("workflow_id", self.workflow_id),
            ("role", self.role),
            ("expected_output_schema", self.expected_output_schema),
            ("owner", self.owner),
            ("merge_key", self.merge_key),
        ):
            require_identifier(value, name)
        _unique_identifiers(self.target_node_ids, "target_node_ids")
        if not self.target_node_ids:
            raise ContractError("specialist packet requires a target node")
        for name, digest in (
            ("context_manifest_sha256", self.context_manifest_sha256),
            (
                "context_selection_receipt_sha256",
                self.context_selection_receipt_sha256,
            ),
        ):
            require_sha256(digest, name)
        _sorted_unique(self.input_record_sha256s, "input_record_sha256s")
        for digest in self.input_record_sha256s:
            require_sha256(digest, "input_record_sha256")
        body = {
            "schema_version": self.schema_version,
            "packet_id": self.packet_id,
            "workflow_id": self.workflow_id,
            "role": self.role,
            "target_node_ids": self.target_node_ids,
            "context_manifest_sha256": self.context_manifest_sha256,
            "context_selection_receipt_sha256": (
                self.context_selection_receipt_sha256
            ),
            "input_record_sha256s": self.input_record_sha256s,
            "expected_output_schema": self.expected_output_schema,
            "owner": self.owner,
            "merge_key": self.merge_key,
        }
        if self.packet_sha256 != canonical_sha256(body):
            raise ContractError("specialist task packet digest mismatch")


def build_specialist_task_packet_v2(
    *,
    packet_id: str,
    manifest: ContextManifestV2,
    role: str,
    target_node_ids: Sequence[str],
    expected_output_schema: str,
    owner: str,
    merge_key: str,
) -> SpecialistTaskPacketV2:
    """Create a specialist packet whose inputs exactly match the manifest."""

    normalized_targets = tuple(str(item) for item in target_node_ids)
    if manifest.target_node_id not in normalized_targets:
        raise ContractError("specialist targets must include the manifest target")
    body = {
        "schema_version": "chemsmart.specialist-task-packet.v2",
        "packet_id": require_identifier(packet_id, "packet_id"),
        "workflow_id": manifest.workflow_id,
        "role": require_identifier(role, "role"),
        "target_node_ids": normalized_targets,
        "context_manifest_sha256": manifest.manifest_sha256,
        "context_selection_receipt_sha256": (
            manifest.selection_receipt_sha256
        ),
        "input_record_sha256s": manifest.included_record_sha256s,
        "expected_output_schema": require_identifier(
            expected_output_schema, "expected_output_schema"
        ),
        "owner": require_identifier(owner, "owner"),
        "merge_key": require_identifier(merge_key, "merge_key"),
    }
    return SpecialistTaskPacketV2(
        **body, packet_sha256=canonical_sha256(body)
    )


@dataclass(frozen=True)
class HarnessContextArmV1:
    """Binds P0/P1/P2 to one experiment case without ambient state."""

    schema_version: str
    experiment_id: str
    case_id: str
    arm_id: str
    mode: str
    policy_sha256: str
    prompt_sha256: str
    tool_schema_sha256: str
    arm_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.harness-context-arm.v1":
            raise ContractError("unsupported harness context arm")
        for name, value in (
            ("experiment_id", self.experiment_id),
            ("case_id", self.case_id),
            ("arm_id", self.arm_id),
        ):
            require_identifier(value, name)
        if self.arm_id not in CONTEXT_ARM_MODES:
            raise ContractError("context arm must be p0, p1, or p2")
        if self.mode != CONTEXT_ARM_MODES[self.arm_id]:
            raise ContractError("context arm and mode disagree")
        for name, digest in (
            ("policy_sha256", self.policy_sha256),
            ("prompt_sha256", self.prompt_sha256),
            ("tool_schema_sha256", self.tool_schema_sha256),
        ):
            require_sha256(digest, name)
        body = {
            "schema_version": self.schema_version,
            "experiment_id": self.experiment_id,
            "case_id": self.case_id,
            "arm_id": self.arm_id,
            "mode": self.mode,
            "policy_sha256": self.policy_sha256,
            "prompt_sha256": self.prompt_sha256,
            "tool_schema_sha256": self.tool_schema_sha256,
        }
        if self.arm_sha256 != canonical_sha256(body):
            raise ContractError("harness context arm digest mismatch")


def build_harness_context_arm(
    *,
    experiment_id: str,
    case_id: str,
    policy: TaskDependencyContextPolicyV1,
    prompt_sha256: str,
    tool_schema_sha256: str,
) -> HarnessContextArmV1:
    body = {
        "schema_version": "chemsmart.harness-context-arm.v1",
        "experiment_id": require_identifier(experiment_id, "experiment_id"),
        "case_id": require_identifier(case_id, "case_id"),
        "arm_id": policy.arm_id,
        "mode": policy.mode,
        "policy_sha256": policy.policy_sha256,
        "prompt_sha256": require_sha256(prompt_sha256, "prompt_sha256"),
        "tool_schema_sha256": require_sha256(
            tool_schema_sha256, "tool_schema_sha256"
        ),
    }
    return HarnessContextArmV1(**body, arm_sha256=canonical_sha256(body))


def _plan_graph(
    plan: ScientificWorkflowPlanV2,
) -> tuple[tuple[str, ...], dict[str, set[str]]]:
    if not isinstance(plan, ScientificWorkflowPlanV2):
        raise ContractError("dependency context requires ScientificWorkflowPlanV2")
    node_ids = tuple(node.node_id for node in plan.nodes)
    parents = {node_id: set() for node_id in node_ids}
    for edge in plan.edges:
        parents[edge.target_node_id].add(edge.source_node_id)
    return node_ids, parents


def _ancestor_sets(
    target_node_id: str, parents: dict[str, set[str]]
) -> tuple[set[str], set[str]]:
    direct = set(parents[target_node_id])
    all_ancestors: set[str] = set()
    frontier = list(direct)
    while frontier:
        node_id = frontier.pop()
        if node_id in all_ancestors:
            continue
        all_ancestors.add(node_id)
        frontier.extend(parents[node_id])
    return direct, all_ancestors.difference(direct)


def _context_payload_bytes(
    *,
    included_node_ids: tuple[str, ...],
    selected_edges: tuple[ScientificWorkflowEdgeV2, ...],
    evidence_refs: tuple[PredecessorEvidenceRefV1, ...],
) -> int:
    envelope_bytes = len(
        canonical_json(
            {
                "included_node_ids": included_node_ids,
                "selected_edges": selected_edges,
                "evidence_refs": evidence_refs,
            }
        ).encode("utf-8")
    )
    return envelope_bytes + sum(item.size_bytes for item in evidence_refs)


def select_task_dependency_context(
    *,
    selection_id: str,
    plan: ScientificWorkflowPlanV2,
    target_node_id: str,
    policy: TaskDependencyContextPolicyV1,
    evidence_refs: Sequence[PredecessorEvidenceRefV1] = (),
) -> tuple[TaskDependencyContextV1 | None, ContextSelectionReceiptV1]:
    """Select an exact P0/P1/P2 context and issue an observable receipt.

    Required records are never truncated.  If the selected packet would exceed
    the public byte budget, no context is returned and the receipt records a
    blocked selection.
    """

    normalized_selection_id = require_identifier(selection_id, "selection_id")
    normalized_target = require_identifier(target_node_id, "target_node_id")
    node_ids, parents = _plan_graph(plan)
    if normalized_target not in parents:
        raise ContractError("dependency context target is not a workflow node")
    positions = {node_id: index for index, node_id in enumerate(node_ids)}
    direct, transitive = _ancestor_sets(normalized_target, parents)
    ancestors = direct.union(transitive)

    if policy.mode == "none":
        included = {normalized_target}
        exposed_direct: set[str] = set()
        exposed_transitive: set[str] = set()
    elif policy.mode == "full_history":
        included = set(node_ids[: positions[normalized_target] + 1])
        exposed_direct = direct
        exposed_transitive = transitive
    elif policy.mode == "dependency_projected":
        included = ancestors.union({normalized_target})
        exposed_direct = direct
        exposed_transitive = transitive
    else:  # guarded by the policy contract; retained for defensive callers.
        raise ContractError("unsupported task dependency context mode")

    included_node_ids = tuple(
        node_id for node_id in node_ids if node_id in included
    )
    excluded_node_ids = tuple(
        node_id for node_id in node_ids if node_id not in included
    )
    non_dependency_node_ids = tuple(
        node_id
        for node_id in node_ids
        if node_id
        in included.difference(ancestors).difference({normalized_target})
    )
    selected_edges = (
        ()
        if policy.mode == "none"
        else tuple(
            sorted(
                (
                    edge
                    for edge in plan.edges
                    if edge.source_node_id in included
                    and edge.target_node_id in included
                ),
                key=lambda edge: edge.edge_id,
            )
        )
    )

    refs = tuple(sorted(evidence_refs, key=lambda item: item.evidence_ref_sha256))
    ref_digests = tuple(item.evidence_ref_sha256 for item in refs)
    if len(ref_digests) != len(set(ref_digests)):
        raise ContractError("predecessor evidence references must be unique")
    record_digests = tuple(item.record_sha256 for item in refs)
    if len(record_digests) != len(set(record_digests)):
        raise ContractError("predecessor evidence records must be unique")
    known_node_ids = set(node_ids)
    for item in refs:
        if item.node_id and item.node_id not in known_node_ids:
            raise ContractError("predecessor evidence names an unknown node")

    selected_refs: list[PredecessorEvidenceRefV1] = []
    excluded: list[EvidenceExclusionV1] = []
    allowed_classes = set(policy.record_classes)
    for item in refs:
        if item.record_class not in allowed_classes:
            excluded.append(
                EvidenceExclusionV1(
                    item.evidence_ref_sha256, "record_class_not_selected"
                )
            )
        elif item.node_id and item.node_id not in included:
            excluded.append(
                EvidenceExclusionV1(
                    item.evidence_ref_sha256, "node_not_selected"
                )
            )
        else:
            selected_refs.append(item)
    selected_refs_tuple = tuple(selected_refs)
    selected_public_bytes = _context_payload_bytes(
        included_node_ids=included_node_ids,
        selected_edges=selected_edges,
        evidence_refs=selected_refs_tuple,
    )

    receipt_base: dict[str, Any] = {
        "schema_version": "chemsmart.context-selection-receipt.v1",
        "selection_id": normalized_selection_id,
        "workflow_id": plan.workflow_id,
        "plan_sha256": plan.plan_sha256,
        "target_node_id": normalized_target,
        "policy_sha256": policy.policy_sha256,
        "candidate_evidence_ref_sha256s": tuple(
            item.evidence_ref_sha256 for item in selected_refs_tuple
        ),
        "selected_public_bytes": selected_public_bytes,
    }
    if selected_public_bytes > policy.max_public_bytes:
        blocked_exclusions = tuple(
            sorted(
                (
                    *excluded,
                    *(
                        EvidenceExclusionV1(
                            item.evidence_ref_sha256,
                            "context_budget_exceeded",
                        )
                        for item in selected_refs_tuple
                    ),
                ),
                key=lambda item: item.evidence_ref_sha256,
            )
        )
        receipt_body = {
            **receipt_base,
            "context_sha256": "",
            "included_evidence_ref_sha256s": (),
            "exclusions": blocked_exclusions,
            "status": "blocked_context_budget",
            "reason": (
                f"required public context is {selected_public_bytes} bytes; "
                f"policy permits {policy.max_public_bytes}"
            ),
        }
        receipt = ContextSelectionReceiptV1(
            **receipt_body, receipt_sha256=canonical_sha256(receipt_body)
        )
        return None, receipt

    context_body = {
        "schema_version": "chemsmart.task-dependency-context.v1",
        "workflow_id": plan.workflow_id,
        "plan_sha256": plan.plan_sha256,
        "target_node_id": normalized_target,
        "policy_sha256": policy.policy_sha256,
        "mode": policy.mode,
        "direct_predecessor_node_ids": tuple(
            node_id for node_id in node_ids if node_id in exposed_direct
        ),
        "transitive_ancestor_node_ids": tuple(
            node_id for node_id in node_ids if node_id in exposed_transitive
        ),
        "included_node_ids": included_node_ids,
        "excluded_node_ids": excluded_node_ids,
        "non_dependency_node_ids": non_dependency_node_ids,
        "selected_edges": selected_edges,
        "evidence_refs": selected_refs_tuple,
    }
    context = TaskDependencyContextV1(
        **context_body, context_sha256=canonical_sha256(context_body)
    )
    receipt_body = {
        **receipt_base,
        "context_sha256": context.context_sha256,
        "included_evidence_ref_sha256s": tuple(
            item.evidence_ref_sha256 for item in selected_refs_tuple
        ),
        "exclusions": tuple(
            sorted(excluded, key=lambda item: item.evidence_ref_sha256)
        ),
        "status": "selected",
        "reason": "",
    }
    receipt = ContextSelectionReceiptV1(
        **receipt_body, receipt_sha256=canonical_sha256(receipt_body)
    )
    return context, receipt


def build_context_manifest_v2(
    *,
    manifest_id: str,
    plan: ScientificWorkflowPlanV2,
    context: TaskDependencyContextV1,
    selection_receipt: ContextSelectionReceiptV1,
    source_sha256s: Sequence[str],
    artifact_sha256s: Sequence[str],
    tool_schema_sha256: str,
    allowed_tools: Sequence[str],
    token_budget: int,
    tool_call_budget: int,
    wall_time_seconds: int,
) -> ContextManifestV2:
    """Bind the selected records and context receipt to specialist budgets."""

    if (
        context.workflow_id != plan.workflow_id
        or context.plan_sha256 != plan.plan_sha256
    ):
        raise ContractError("dependency context belongs to another plan")
    if (
        selection_receipt.workflow_id != context.workflow_id
        or selection_receipt.plan_sha256 != context.plan_sha256
        or selection_receipt.target_node_id != context.target_node_id
        or selection_receipt.policy_sha256 != context.policy_sha256
        or selection_receipt.context_sha256 != context.context_sha256
    ):
        raise ContractError("selection receipt belongs to another context")
    if selection_receipt.status != "selected":
        raise ContractError("blocked context cannot produce a manifest")
    evidence_ref_sha256s = tuple(
        item.evidence_ref_sha256 for item in context.evidence_refs
    )
    if (
        selection_receipt.included_evidence_ref_sha256s
        != evidence_ref_sha256s
    ):
        raise ContractError("selection receipt records differ from context")
    record_sha256s = tuple(
        sorted(item.record_sha256 for item in context.evidence_refs)
    )
    body = {
        "schema_version": "chemsmart.context-manifest.v2",
        "manifest_id": require_identifier(manifest_id, "manifest_id"),
        "workflow_id": plan.workflow_id,
        "task_spec_sha256": plan.task_spec_sha256,
        "scientific_identity_sha256": plan.scientific_identity_sha256,
        "target_node_id": context.target_node_id,
        "dependency_context_sha256": context.context_sha256,
        "selection_receipt_sha256": selection_receipt.receipt_sha256,
        "included_record_sha256s": record_sha256s,
        "source_sha256s": tuple(sorted(set(source_sha256s))),
        "artifact_sha256s": tuple(sorted(set(artifact_sha256s))),
        "tool_schema_sha256": require_sha256(
            tool_schema_sha256, "tool_schema_sha256"
        ),
        "allowed_tools": tuple(sorted(set(allowed_tools))),
        "token_budget": int(token_budget),
        "tool_call_budget": int(tool_call_budget),
        "wall_time_seconds": int(wall_time_seconds),
    }
    return ContextManifestV2(**body, manifest_sha256=canonical_sha256(body))


__all__ = [
    "CONTEXT_ARM_MODES",
    "CONTEXT_MODES",
    "CONTEXT_SELECTOR_VERSION",
    "ContextManifestV2",
    "ContextSelectionReceiptV1",
    "EvidenceExclusionV1",
    "HarnessContextArmV1",
    "PredecessorEvidenceRefV1",
    "SpecialistTaskPacketV2",
    "TaskDependencyContextPolicyV1",
    "TaskDependencyContextV1",
    "build_context_manifest_v2",
    "build_dependency_context_public_projection",
    "build_harness_context_arm",
    "build_predecessor_evidence_ref",
    "bind_selected_public_records",
    "build_specialist_task_packet_v2",
    "build_task_dependency_context_policy",
    "select_task_dependency_context",
]
