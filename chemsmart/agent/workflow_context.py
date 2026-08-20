"""Project the workflow state the host already knows back to the model.

The harness holds a complete, typed model of a workflow: which node feeds
which, which artifact each edge carries, and which nodes could run right now.
It used to spend that knowledge entirely on refusal -- ``derive_ready_node_ids``
existed only to reject an out-of-order execution -- while the model received
two flat lists of node IDs and had to reconstruct the branching from memory
across a long tool conversation.

That is the wrong direction for the dependency.  A model that has lost track of
its own DAG re-plans, and re-planning is where a stage carrying the answer gets
dropped.  So the host states, per node, what it is waiting for, which upstream
output feeds which input, and what it will produce.

The projection is *derived*.  It is computed from the plan plus the set of
artifacts that actually exist, and it never appears in any tool's input schema.
A model cannot assert that a node is ready; it can only be told.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable, Mapping, Sequence

from chemsmart.agent._contracts import ContractError, require_identifier

#: ``completed`` -- this exact approved node already has a validated result.
#: ``ready`` -- every input exists, so this node can be materialized now.
#: ``waiting`` -- an upstream producer has not yet produced its artifact.
#: ``blocked`` -- required, but outside the current program surface.
WORKFLOW_NODE_CONTEXT_STATES = ("blocked", "completed", "ready", "waiting")


@dataclass(frozen=True)
class UnsatisfiedInputV1:
    """One input edge whose artifact does not exist yet."""

    binding_id: str
    artifact_class: str
    producer_node_id: str
    producer_output_id: str

    def __post_init__(self) -> None:
        require_identifier(self.binding_id, "binding_id")


@dataclass(frozen=True)
class WorkflowNodeContextV1:
    """What one node is waiting for, and what depends on it."""

    node_id: str
    state: str
    waiting_on: tuple[str, ...]
    unsatisfied_inputs: tuple[UnsatisfiedInputV1, ...]
    produces: tuple[str, ...]
    dependents: tuple[str, ...]
    reason: str

    def __post_init__(self) -> None:
        require_identifier(self.node_id, "node_id")
        if self.state not in WORKFLOW_NODE_CONTEXT_STATES:
            raise ContractError(
                "state must be one of "
                f"{sorted(WORKFLOW_NODE_CONTEXT_STATES)}"
            )
        if self.state == "ready" and self.unsatisfied_inputs:
            raise ContractError("a ready node cannot have unsatisfied inputs")
        if self.state == "waiting" and not self.waiting_on:
            raise ContractError("a waiting node must name what it waits for")


@dataclass(frozen=True)
class WorkflowContextProjectionV1:
    """The frontier and the dependency shape, as the host sees them."""

    schema_version: str
    workflow_id: str
    completed_node_ids: tuple[str, ...]
    ready_node_ids: tuple[str, ...]
    waiting_node_ids: tuple[str, ...]
    blocked_node_ids: tuple[str, ...]
    nodes: tuple[WorkflowNodeContextV1, ...]

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.workflow-context.v1":
            raise ContractError("unsupported workflow context schema")
        require_identifier(self.workflow_id, "workflow_id")

    def node(self, node_id: str) -> WorkflowNodeContextV1 | None:
        for item in self.nodes:
            if item.node_id == node_id:
                return item
        return None


def _producer_edges(node: Any) -> tuple[Any, ...]:
    return tuple(getattr(node, "inputs", ()) or ())


def _input_is_satisfied(
    item: Any,
    materialized: Mapping[str, Any] | Iterable[str],
    *,
    consumer_node_id: str,
    materialized_producer_inputs: set[tuple[str, str, str, str]],
) -> bool:
    """Return whether an external artifact or observed producer handoff exists."""

    artifact_id = str(getattr(item, "artifact_id", "") or "")
    if artifact_id:
        return artifact_id in materialized
    producer_node_id = str(getattr(item, "producer_node_id", "") or "")
    producer_output_id = str(getattr(item, "producer_output_id", "") or "")
    if not producer_node_id or not producer_output_id:
        return False
    return (
        consumer_node_id,
        str(getattr(item, "binding_id", "") or ""),
        producer_node_id,
        producer_output_id,
    ) in materialized_producer_inputs


def project_workflow_context(
    *,
    workflow_id: str,
    nodes: Sequence[Any],
    materialized_artifact_ids: Mapping[str, Any] | Iterable[str] = (),
    materialized_producer_inputs: Iterable[tuple[str, str, str, str]] = (),
    completed_node_ids: Iterable[str] = (),
    blocked_reasons: Mapping[str, str] | None = None,
) -> WorkflowContextProjectionV1:
    """Derive per-node readiness and dependency context from a workflow.

    ``nodes`` may be draft nodes or plan nodes: only ``node_id``,
    ``dependencies``, ``inputs``, and ``expected_outputs`` are read, so one
    projection serves planning and execution alike.

    Args:
        workflow_id: The workflow these nodes belong to.
        nodes: The workflow's nodes, in any order.
        materialized_artifact_ids: Artifacts that actually exist right now.
        materialized_producer_inputs: Exact ``(consumer, binding, producer,
            output)`` tuples whose validated handoff artifacts exist.
        completed_node_ids: Exact approved nodes with validated execution
            results. These are reported as complete, not offered for rerun.
        blocked_reasons: Node IDs that cannot be expressed at all, and why.

    Returns:
        The projection, with nodes sorted by ID so the record is canonical.
    """

    require_identifier(workflow_id, "workflow_id")
    materialized = (
        materialized_artifact_ids
        if isinstance(materialized_artifact_ids, Mapping)
        else set(materialized_artifact_ids)
    )
    materialized_producer_inputs = set(materialized_producer_inputs)
    completed_node_ids = set(completed_node_ids)
    blocked_reasons = dict(blocked_reasons or {})

    dependents: dict[str, set[str]] = {
        str(node.node_id): set() for node in nodes
    }
    for node in nodes:
        for dependency in getattr(node, "dependencies", ()) or ():
            if dependency in dependents:
                dependents[dependency].add(str(node.node_id))
        for item in _producer_edges(node):
            producer = str(getattr(item, "producer_node_id", "") or "")
            if producer in dependents:
                dependents[producer].add(str(node.node_id))

    contexts: list[WorkflowNodeContextV1] = []
    for node in nodes:
        node_id = str(node.node_id)
        unsatisfied = tuple(
            UnsatisfiedInputV1(
                binding_id=str(item.binding_id),
                artifact_class=str(getattr(item, "artifact_class", "") or ""),
                producer_node_id=str(
                    getattr(item, "producer_node_id", "") or ""
                ),
                producer_output_id=str(
                    getattr(item, "producer_output_id", "") or ""
                ),
            )
            for item in sorted(
                _producer_edges(node),
                key=lambda entry: str(entry.binding_id),
            )
            if not _input_is_satisfied(
                item,
                materialized,
                consumer_node_id=node_id,
                materialized_producer_inputs=materialized_producer_inputs,
            )
        )
        satisfied_data_producers = {
            str(getattr(item, "producer_node_id", "") or "")
            for item in _producer_edges(node)
            if str(getattr(item, "producer_node_id", "") or "")
            and _input_is_satisfied(
                item,
                materialized,
                consumer_node_id=node_id,
                materialized_producer_inputs=materialized_producer_inputs,
            )
        }
        waiting_on = tuple(
            sorted(
                {
                    item.producer_node_id
                    for item in unsatisfied
                    if item.producer_node_id
                }.union(
                    dependency
                    for dependency in (getattr(node, "dependencies", ()) or ())
                    if dependency not in completed_node_ids
                    and dependency not in satisfied_data_producers
                )
            )
        )
        if node_id in blocked_reasons:
            state = "blocked"
            reason = blocked_reasons[node_id]
        elif node_id in completed_node_ids:
            state = "completed"
            reason = "validated execution result already exists"
        elif waiting_on:
            state = "waiting"
            producer_outputs = tuple(
                f"{item.producer_node_id}.{item.producer_output_id}"
                for item in unsatisfied
                if item.producer_node_id
            )
            control_dependencies = tuple(
                dependency
                for dependency in waiting_on
                if dependency
                not in {item.producer_node_id for item in unsatisfied}
            )
            waiting_labels = (
                *producer_outputs,
                *(
                    dependency + " completion"
                    for dependency in control_dependencies
                ),
            )
            reason = "waiting for " + ", ".join(waiting_labels)
        elif unsatisfied:
            # An input with no producer and no artifact is an external artifact
            # the plan never bound -- a different defect from waiting on a
            # producer, and one the model must repair rather than wait out.
            state = "blocked"
            reason = (
                "input "
                + ", ".join(item.binding_id for item in unsatisfied)
                + " names no artifact and no producer"
            )
        else:
            state = "ready"
            reason = "every input exists"
        contexts.append(
            WorkflowNodeContextV1(
                node_id=node_id,
                state=state,
                waiting_on=waiting_on,
                unsatisfied_inputs=unsatisfied,
                produces=tuple(
                    sorted(
                        str(item.output_id)
                        for item in getattr(node, "expected_outputs", ())
                    )
                ),
                dependents=tuple(sorted(dependents.get(node_id, ()))),
                reason=reason,
            )
        )

    contexts.sort(key=lambda item: item.node_id)
    return WorkflowContextProjectionV1(
        schema_version="chemsmart.workflow-context.v1",
        workflow_id=workflow_id,
        completed_node_ids=tuple(
            item.node_id for item in contexts if item.state == "completed"
        ),
        ready_node_ids=tuple(
            item.node_id for item in contexts if item.state == "ready"
        ),
        waiting_node_ids=tuple(
            item.node_id for item in contexts if item.state == "waiting"
        ),
        blocked_node_ids=tuple(
            item.node_id for item in contexts if item.state == "blocked"
        ),
        nodes=tuple(contexts),
    )


__all__ = [
    "UnsatisfiedInputV1",
    "WORKFLOW_NODE_CONTEXT_STATES",
    "WorkflowContextProjectionV1",
    "WorkflowNodeContextV1",
    "project_workflow_context",
]
