"""Deterministic predecessor context for scientific workflow tasks.

The model does not decide which prior records are relevant.  This module
derives an immutable context packet from the host-owned scientific DAG and a
host-selected policy.  ``target_only`` exposes only the target,
``ordered_predecessors`` exposes preceding nodes in topological order, and
``dependency_projected`` exposes the exact backward dependency slice.

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

CONTEXT_MODES = (
    "dependency_projected",
    "ordered_predecessors",
    "target_only",
)
CONTEXT_SELECTOR_VERSION = "chemsmart.task-dependency-context-selector.v2"


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
class TaskDependencyContextPolicyV2:
    """Host-neutral rule for selecting predecessor context."""

    schema_version: str
    policy_id: str
    mode: str
    selector_version: str
    record_classes: tuple[str, ...]
    max_public_bytes: int
    policy_sha256: str

    def __post_init__(self) -> None:
        if (
            self.schema_version
            != "chemsmart.task-dependency-context-policy.v2"
        ):
            raise ContractError("unsupported task dependency context policy")
        require_identifier(self.policy_id, "policy_id")
        if self.mode not in CONTEXT_MODES:
            raise ContractError("unsupported task dependency context mode")
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
            "mode": self.mode,
            "selector_version": self.selector_version,
            "record_classes": self.record_classes,
            "max_public_bytes": self.max_public_bytes,
        }
        if self.policy_sha256 != canonical_sha256(body):
            raise ContractError(
                "task dependency context policy digest mismatch"
            )


def build_task_dependency_context_policy(
    *,
    policy_id: str,
    mode: str,
    record_classes: Sequence[str],
    max_public_bytes: int,
) -> TaskDependencyContextPolicyV2:
    """Build a canonical host-neutral context policy."""

    normalized_mode = str(mode).strip().lower()
    if normalized_mode not in CONTEXT_MODES:
        raise ContractError("unsupported task dependency context mode")
    body = {
        "schema_version": "chemsmart.task-dependency-context-policy.v2",
        "policy_id": require_identifier(policy_id, "policy_id"),
        "mode": normalized_mode,
        "selector_version": CONTEXT_SELECTOR_VERSION,
        "record_classes": tuple(sorted(set(record_classes))),
        "max_public_bytes": int(max_public_bytes),
    }
    return TaskDependencyContextPolicyV2(
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
            raise ContractError(
                "predecessor evidence reference digest mismatch"
            )


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
    context: TaskDependencyContextV2,
    records: Mapping[str, Mapping[str, Any]],
) -> tuple[dict[str, Any], ...]:
    """Bind selected references to their exact canonical public payloads.

    ``TaskDependencyContextV2`` deliberately selects references rather than
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
class TaskDependencyContextV2:
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
        if self.schema_version != "chemsmart.task-dependency-context.v2":
            raise ContractError("unsupported task dependency context")
        require_identifier(self.workflow_id, "workflow_id")
        require_sha256(self.plan_sha256, "plan_sha256")
        require_identifier(self.target_node_id, "target_node_id")
        require_sha256(self.policy_sha256, "policy_sha256")
        if self.mode not in CONTEXT_MODES:
            raise ContractError("unsupported task dependency context mode")
        for name, node_ids in (
            ("direct_predecessor_node_ids", self.direct_predecessor_node_ids),
            (
                "transitive_ancestor_node_ids",
                self.transitive_ancestor_node_ids,
            ),
            ("included_node_ids", self.included_node_ids),
            ("excluded_node_ids", self.excluded_node_ids),
            ("non_dependency_node_ids", self.non_dependency_node_ids),
        ):
            _unique_identifiers(node_ids, name)
        included = set(self.included_node_ids)
        excluded = set(self.excluded_node_ids)
        if self.target_node_id not in included:
            raise ContractError(
                "dependency context must include its target node"
            )
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
        if self.mode == "target_only" and (
            self.direct_predecessor_node_ids
            or self.transitive_ancestor_node_ids
            or self.non_dependency_node_ids
            or self.selected_edges
        ):
            raise ContractError(
                "target-only context cannot expose predecessors"
            )
        if (
            self.mode == "dependency_projected"
            and self.non_dependency_node_ids
        ):
            raise ContractError(
                "dependency-projected context cannot include unrelated nodes"
            )
        edge_ids = tuple(edge.edge_id for edge in self.selected_edges)
        if edge_ids != tuple(sorted(set(edge_ids))):
            raise ContractError(
                "selected context edges must be sorted and unique"
            )
        positions = {
            node_id: index
            for index, node_id in enumerate(self.included_node_ids)
        }
        for edge in self.selected_edges:
            if (
                edge.source_node_id not in included
                or edge.target_node_id not in included
            ):
                raise ContractError(
                    "selected edge leaves the included context"
                )
            if (
                positions[edge.source_node_id]
                >= positions[edge.target_node_id]
            ):
                raise ContractError(
                    "included context must be topologically ordered"
                )
        evidence_digests = tuple(
            evidence.evidence_ref_sha256 for evidence in self.evidence_refs
        )
        if evidence_digests != tuple(sorted(set(evidence_digests))):
            raise ContractError(
                "context evidence refs must be sorted and unique"
            )
        record_digests = tuple(
            evidence.record_sha256 for evidence in self.evidence_refs
        )
        if len(record_digests) != len(set(record_digests)):
            raise ContractError("context cannot repeat a public record")
        for evidence in self.evidence_refs:
            if evidence.node_id and evidence.node_id not in included:
                raise ContractError(
                    "context evidence belongs to an excluded node"
                )
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
            raise ContractError(
                "included and excluded context evidence overlap"
            )
        if self.selected_public_bytes < 0:
            raise ContractError("selected_public_bytes must be non-negative")
        if self.status not in {"blocked_context_budget", "selected"}:
            raise ContractError("unsupported context selection status")
        if self.status == "selected":
            require_sha256(self.context_sha256, "context_sha256")
            if self.reason:
                raise ContractError(
                    "selected context cannot carry a block reason"
                )
        else:
            if self.context_sha256:
                raise ContractError(
                    "blocked context cannot claim a context digest"
                )
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
    context: TaskDependencyContextV2,
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


def _plan_graph(
    plan: ScientificWorkflowPlanV2,
) -> tuple[tuple[str, ...], dict[str, set[str]]]:
    if not isinstance(plan, ScientificWorkflowPlanV2):
        raise ContractError(
            "dependency context requires ScientificWorkflowPlanV2"
        )
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
    policy: TaskDependencyContextPolicyV2,
    evidence_refs: Sequence[PredecessorEvidenceRefV1] = (),
) -> tuple[TaskDependencyContextV2 | None, ContextSelectionReceiptV1]:
    """Select an exact host-policy context and issue an observable receipt.

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

    if policy.mode == "target_only":
        included = {normalized_target}
        exposed_direct: set[str] = set()
        exposed_transitive: set[str] = set()
    elif policy.mode == "ordered_predecessors":
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
        if policy.mode == "target_only"
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

    refs = tuple(
        sorted(evidence_refs, key=lambda item: item.evidence_ref_sha256)
    )
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
        "schema_version": "chemsmart.task-dependency-context.v2",
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
    context = TaskDependencyContextV2(
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


__all__ = [
    "CONTEXT_MODES",
    "CONTEXT_SELECTOR_VERSION",
    "ContextSelectionReceiptV1",
    "EvidenceExclusionV1",
    "PredecessorEvidenceRefV1",
    "TaskDependencyContextPolicyV2",
    "TaskDependencyContextV2",
    "build_dependency_context_public_projection",
    "build_predecessor_evidence_ref",
    "bind_selected_public_records",
    "build_task_dependency_context_policy",
    "select_task_dependency_context",
]
