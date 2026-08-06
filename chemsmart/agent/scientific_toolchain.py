"""Planning contracts for connected calculation and analysis tool chains.

This module deliberately models *future* data flow.  A paper-level plan must be
able to say that an analysis consumes a quantity produced by a calculation
before that calculation has run.  Consequently these intent objects bind
producer node/output names, not artifact paths or receipt hashes.  Materialized
artifacts are attached later by the existing execution and analysis tools.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Mapping, Sequence

from chemsmart.agent._contracts import (
    ContractError,
    canonical_sha256,
    require_identifier,
)
from chemsmart.agent.workflows import CommandNodeIntentV1


ANALYSIS_INTENT_KINDS = (
    "claim_rendering",
    "quantity_expression",
    "result_extraction",
    "scientific_validation",
    "thermochemistry",
    "unsupported_external",
)
ANALYSIS_INTENT_SUPPORT_STATES = ("blocked_unsupported", "planned")


class ScientificToolchainContractError(ContractError):
    """Raised when a proposed paper-level tool chain is inconsistent."""


def _identifier(value: str, field: str) -> str:
    try:
        return require_identifier(value, field)
    except ContractError as exc:
        raise ScientificToolchainContractError(str(exc)) from exc


def _sorted_unique(values: Sequence[str], field: str) -> tuple[str, ...]:
    normalized = tuple(sorted(set(values)))
    if tuple(values) != normalized:
        raise ScientificToolchainContractError(
            f"{field} must be sorted and unique"
        )
    return normalized


@dataclass(frozen=True)
class AnalysisInputIntentV1:
    """A future program output or analysis quantity required by an analysis."""

    input_id: str
    source_kind: str
    producer_node_id: str
    producer_output_id: str

    def __post_init__(self) -> None:
        _identifier(self.input_id, "analysis input_id")
        if self.source_kind not in {"analysis_output", "program_output"}:
            raise ScientificToolchainContractError(
                "analysis source_kind must be analysis_output or program_output"
            )
        _identifier(self.producer_node_id, "analysis producer_node_id")
        _identifier(self.producer_output_id, "analysis producer_output_id")


@dataclass(frozen=True)
class AnalysisSelectorIntentV1:
    """Program-neutral quantity requested from a future result artifact."""

    quantity_id: str
    selector: str

    def __post_init__(self) -> None:
        _identifier(self.quantity_id, "analysis selector quantity_id")
        _identifier(self.selector, "analysis selector")


@dataclass(frozen=True)
class AnalysisOutputIntentV1:
    """A named quantity promised by a planned analysis node."""

    output_id: str
    quantity_kind: str
    unit: str

    def __post_init__(self) -> None:
        _identifier(self.output_id, "analysis output_id")
        _identifier(self.quantity_kind, "analysis quantity_kind")
        if not str(self.unit).strip():
            raise ScientificToolchainContractError(
                "analysis output unit must not be empty"
            )


@dataclass(frozen=True)
class AnalysisNodeIntentV1:
    """One extraction, validation, mathematics, or reporting intent."""

    node_id: str
    analysis_kind: str
    dependencies: tuple[str, ...]
    inputs: tuple[AnalysisInputIntentV1, ...]
    selectors: tuple[AnalysisSelectorIntentV1, ...]
    outputs: tuple[AnalysisOutputIntentV1, ...]
    expression_nodes: tuple[dict[str, object], ...]
    expression_output_node_ids: tuple[str, ...]
    temperature_k: float | None
    pressure_atm: float | None
    support_state: str
    blocked_reason: str

    def __post_init__(self) -> None:
        _identifier(self.node_id, "analysis node_id")
        if self.analysis_kind not in ANALYSIS_INTENT_KINDS:
            raise ScientificToolchainContractError(
                "unsupported analysis intent kind"
            )
        if self.support_state not in ANALYSIS_INTENT_SUPPORT_STATES:
            raise ScientificToolchainContractError(
                "unsupported analysis intent support_state"
            )
        object.__setattr__(self, "dependencies", tuple(self.dependencies))
        object.__setattr__(self, "inputs", tuple(self.inputs))
        object.__setattr__(self, "selectors", tuple(self.selectors))
        object.__setattr__(self, "outputs", tuple(self.outputs))
        object.__setattr__(
            self,
            "expression_nodes",
            tuple(dict(item) for item in self.expression_nodes),
        )
        object.__setattr__(
            self,
            "expression_output_node_ids",
            tuple(self.expression_output_node_ids),
        )
        _sorted_unique(self.dependencies, "analysis dependencies")
        _sorted_unique(
            tuple(item.input_id for item in self.inputs),
            "analysis input IDs",
        )
        _sorted_unique(
            tuple(item.quantity_id for item in self.selectors),
            "analysis selector IDs",
        )
        _sorted_unique(
            tuple(item.output_id for item in self.outputs),
            "analysis output IDs",
        )
        if not self.outputs:
            raise ScientificToolchainContractError(
                "analysis intent must preserve at least one output"
            )
        if self.support_state == "blocked_unsupported":
            if not str(self.blocked_reason).strip():
                raise ScientificToolchainContractError(
                    "blocked analysis intent requires a localized reason"
                )
        elif str(self.blocked_reason).strip():
            raise ScientificToolchainContractError(
                "blocked_reason applies only to blocked_unsupported"
            )
        if self.analysis_kind == "result_extraction":
            if len(self.inputs) != 1 or not self.selectors:
                raise ScientificToolchainContractError(
                    "result extraction needs one result input and selectors"
                )
        elif self.selectors:
            raise ScientificToolchainContractError(
                "selectors apply only to result extraction"
            )
        if self.analysis_kind == "thermochemistry":
            if self.temperature_k is None or self.temperature_k <= 0:
                raise ScientificToolchainContractError(
                    "thermochemistry requires a positive temperature"
                )
            if self.pressure_atm is None or self.pressure_atm <= 0:
                raise ScientificToolchainContractError(
                    "thermochemistry requires a positive pressure"
                )
        elif self.temperature_k is not None or self.pressure_atm is not None:
            raise ScientificToolchainContractError(
                "temperature and pressure apply only to thermochemistry"
            )
        if self.analysis_kind == "quantity_expression":
            if not self.expression_nodes or not self.expression_output_node_ids:
                raise ScientificToolchainContractError(
                    "quantity expression requires a typed expression DAG"
                )
        elif self.expression_nodes or self.expression_output_node_ids:
            raise ScientificToolchainContractError(
                "expression fields apply only to quantity_expression"
            )


@dataclass(frozen=True)
class ScientificToolchainPlanV1:
    """Connected calculation and postprocessing plan proposed in one tool call."""

    schema_version: str
    plan_id: str
    workflow_id: str
    command_workflow_draft_sha256: str
    calculation_node_ids: tuple[str, ...]
    calculation_observables: tuple[tuple[str, tuple[str, ...]], ...]
    analysis_nodes: tuple[AnalysisNodeIntentV1, ...]
    required_output_ids: tuple[str, ...]
    node_order: tuple[str, ...]
    status: str
    plan_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.scientific-toolchain-plan.v1":
            raise ScientificToolchainContractError(
                "unsupported scientific toolchain schema"
            )
        _identifier(self.plan_id, "toolchain plan_id")
        _identifier(self.workflow_id, "toolchain workflow_id")
        if len(self.command_workflow_draft_sha256) != 64:
            raise ScientificToolchainContractError(
                "toolchain requires a command workflow draft digest"
            )
        _sorted_unique(self.calculation_node_ids, "calculation node IDs")
        _sorted_unique(self.required_output_ids, "required output IDs")
        if self.status != "planned":
            raise ScientificToolchainContractError(
                "scientific toolchain status must be planned"
            )
        body = scientific_toolchain_plan_body(self)
        if self.plan_sha256 != canonical_sha256(body):
            raise ScientificToolchainContractError(
                "scientific toolchain plan digest mismatch"
            )


def scientific_toolchain_plan_body(
    plan: ScientificToolchainPlanV1,
) -> dict[str, object]:
    return {
        "schema_version": plan.schema_version,
        "plan_id": plan.plan_id,
        "workflow_id": plan.workflow_id,
        "command_workflow_draft_sha256": plan.command_workflow_draft_sha256,
        "calculation_node_ids": plan.calculation_node_ids,
        "calculation_observables": plan.calculation_observables,
        "analysis_nodes": plan.analysis_nodes,
        "required_output_ids": plan.required_output_ids,
        "node_order": plan.node_order,
        "status": plan.status,
    }


def build_scientific_toolchain_plan(
    *,
    plan_id: str,
    workflow_id: str,
    command_workflow_draft_sha256: str,
    calculation_nodes: Iterable[CommandNodeIntentV1],
    calculation_observables: Mapping[str, Iterable[str]],
    analysis_nodes: Iterable[AnalysisNodeIntentV1],
    required_output_ids: Iterable[str],
) -> ScientificToolchainPlanV1:
    """Validate and topologically order a connected paper-level plan."""

    calculations = tuple(calculation_nodes)
    analyses = tuple(analysis_nodes)
    calculation_by_id = {node.node_id: node for node in calculations}
    analysis_by_id = {node.node_id: node for node in analyses}
    if len(calculation_by_id) != len(calculations):
        raise ScientificToolchainContractError(
            "calculation node IDs must be unique"
        )
    if len(analysis_by_id) != len(analyses):
        raise ScientificToolchainContractError("analysis node IDs must be unique")
    overlap = set(calculation_by_id).intersection(analysis_by_id)
    if overlap:
        raise ScientificToolchainContractError(
            "calculation and analysis node IDs must be distinct"
        )
    all_ids = set(calculation_by_id) | set(analysis_by_id)
    produced: dict[str, set[str]] = {
        node.node_id: {item.output_id for item in node.expected_outputs}
        for node in calculations
    }
    produced.update(
        {
            node.node_id: {item.output_id for item in node.outputs}
            for node in analyses
        }
    )
    dependencies: dict[str, set[str]] = {
        node_id: set() for node_id in all_ids
    }
    for node in calculations:
        for dependency in node.dependencies:
            if dependency not in calculation_by_id:
                raise ScientificToolchainContractError(
                    "calculation dependency must name a calculation node"
                )
            dependencies[node.node_id].add(dependency)
    normalized_analyses: list[AnalysisNodeIntentV1] = []
    for node in analyses:
        effective_dependencies = set(node.dependencies)
        for item in node.inputs:
            producer = item.producer_node_id
            if producer not in all_ids:
                raise ScientificToolchainContractError(
                    "analysis input references an unknown producer"
                )
            expected_kind = (
                "program_output"
                if producer in calculation_by_id
                else "analysis_output"
            )
            if item.source_kind != expected_kind:
                raise ScientificToolchainContractError(
                    "analysis input source_kind differs from its producer"
                )
            if item.producer_output_id not in produced[producer]:
                raise ScientificToolchainContractError(
                    "analysis input references an unknown producer output"
                )
            effective_dependencies.add(producer)
        unknown = effective_dependencies.difference(all_ids)
        if unknown:
            raise ScientificToolchainContractError(
                "analysis dependency references an unknown node"
            )
        effective_dependencies.discard(node.node_id)
        normalized = AnalysisNodeIntentV1(
            node_id=node.node_id,
            analysis_kind=node.analysis_kind,
            dependencies=tuple(sorted(effective_dependencies)),
            inputs=node.inputs,
            selectors=node.selectors,
            outputs=node.outputs,
            expression_nodes=node.expression_nodes,
            expression_output_node_ids=node.expression_output_node_ids,
            temperature_k=node.temperature_k,
            pressure_atm=node.pressure_atm,
            support_state=node.support_state,
            blocked_reason=node.blocked_reason,
        )
        normalized_analyses.append(normalized)
        dependencies[node.node_id].update(effective_dependencies)
    children: dict[str, set[str]] = {node_id: set() for node_id in all_ids}
    for node_id, parents in dependencies.items():
        for parent in parents:
            children[parent].add(node_id)
    ready = sorted(node_id for node_id, parents in dependencies.items() if not parents)
    ordered: list[str] = []
    while ready:
        node_id = ready.pop(0)
        ordered.append(node_id)
        for child in sorted(children[node_id]):
            dependencies[child].discard(node_id)
            if not dependencies[child] and child not in ordered + ready:
                ready.append(child)
                ready.sort()
    if len(ordered) != len(all_ids):
        raise ScientificToolchainContractError(
            "scientific toolchain contains a dependency cycle"
        )
    all_output_ids = {
        output_id for output_ids in produced.values() for output_id in output_ids
    }
    required = tuple(sorted(set(required_output_ids)))
    missing = sorted(set(required).difference(all_output_ids))
    if missing:
        raise ScientificToolchainContractError(
            f"required output(s) have no producer: {missing}"
        )
    observable_rows = tuple(
        sorted(
            (
                node_id,
                tuple(sorted(set(calculation_observables.get(node_id, ())))),
            )
            for node_id in calculation_by_id
        )
    )
    analysis_by_id = {node.node_id: node for node in normalized_analyses}
    body = {
        "schema_version": "chemsmart.scientific-toolchain-plan.v1",
        "plan_id": plan_id,
        "workflow_id": workflow_id,
        "command_workflow_draft_sha256": command_workflow_draft_sha256,
        "calculation_node_ids": tuple(sorted(calculation_by_id)),
        "calculation_observables": observable_rows,
        "analysis_nodes": tuple(
            analysis_by_id[node_id]
            for node_id in ordered
            if node_id in analysis_by_id
        ),
        "required_output_ids": required,
        "node_order": tuple(ordered),
        "status": "planned",
    }
    return ScientificToolchainPlanV1(
        **body, plan_sha256=canonical_sha256(body)
    )


def project_scientific_toolchain_frontier(
    plan: ScientificToolchainPlanV1,
    *,
    actionable_calculation_node_ids: Iterable[str] = (),
    unresolved_calculation_node_ids: Iterable[str] = (),
) -> dict[str, object]:
    """Return model-useful next actions without pretending future data exist."""

    actionable = set(actionable_calculation_node_ids)
    unresolved = set(unresolved_calculation_node_ids)
    states: dict[str, str] = {}
    nodes: list[dict[str, object]] = []
    analysis_by_id = {node.node_id: node for node in plan.analysis_nodes}
    # Which nodes consume each producer.  The model submitted the edges; the
    # host can therefore state the branching instead of leaving the model to
    # recall its own DAG, which is where a stage carrying the answer is lost.
    dependents: dict[str, set[str]] = {}
    for node in plan.analysis_nodes:
        for edge in node.inputs:
            dependents.setdefault(edge.producer_node_id, set()).add(node.node_id)
    for node_id in plan.node_order:
        waiting_on: tuple[str, ...] = ()
        unsatisfied: tuple[dict[str, str], ...] = ()
        if node_id in plan.calculation_node_ids:
            state = "actionable" if node_id in actionable else "waiting"
            if node_id in unresolved:
                state = "waiting"
            reason = (
                "compile_and_preview"
                if state == "actionable"
                else "resolve_inputs_or_scientific_settings"
            )
        else:
            node = analysis_by_id[node_id]
            blocked_parents = tuple(
                sorted(
                    parent
                    for parent in node.dependencies
                    if states.get(parent)
                    in {"blocked_unsupported", "blocked_upstream"}
                )
            )
            if node.support_state == "blocked_unsupported":
                state = "blocked_unsupported"
                reason = node.blocked_reason
            elif blocked_parents:
                state = "blocked_upstream"
                waiting_on = blocked_parents
                reason = (
                    "blocked because "
                    + ", ".join(blocked_parents)
                    + (" is" if len(blocked_parents) == 1 else " are")
                    + " blocked"
                )
            elif node.dependencies:
                state = "waiting_for_artifact"
                # Name the exact producer outputs.  "await producer outputs"
                # is the same sentence for every waiting node in every
                # workflow, so it carries none of the dependency structure the
                # host already holds.
                unsatisfied = tuple(
                    {
                        "input_id": edge.input_id,
                        "producer_node_id": edge.producer_node_id,
                        "producer_output_id": edge.producer_output_id,
                        "source_kind": edge.source_kind,
                    }
                    for edge in sorted(
                        node.inputs,
                        key=lambda item: (
                            item.producer_node_id,
                            item.producer_output_id,
                            item.input_id,
                        ),
                    )
                )
                waiting_on = tuple(
                    sorted({edge.producer_node_id for edge in node.inputs})
                )
                named = ", ".join(
                    f"{edge['producer_node_id']}.{edge['producer_output_id']}"
                    for edge in unsatisfied
                ) or ", ".join(waiting_on)
                reason = f"await {named}"
            else:
                state = "actionable"
                reason = "execute the registered analysis operation"
        states[node_id] = state
        entry: dict[str, object] = {
            "node_id": node_id,
            "state": state,
            "next_action": reason,
            "feeds": tuple(sorted(dependents.get(node_id, ()))),
        }
        if waiting_on:
            entry["waiting_on"] = waiting_on
        if unsatisfied:
            entry["unsatisfied_inputs"] = unsatisfied
        nodes.append(entry)
    return {
        "schema_version": "chemsmart.scientific-toolchain-frontier.v1",
        "workflow_id": plan.workflow_id,
        "nodes": tuple(nodes),
        "actionable_node_ids": tuple(
            item["node_id"] for item in nodes if item["state"] == "actionable"
        ),
        "waiting_node_ids": tuple(
            item["node_id"]
            for item in nodes
            if item["state"] in {"waiting", "waiting_for_artifact"}
        ),
        "blocked_node_ids": tuple(
            item["node_id"]
            for item in nodes
            if item["state"] in {"blocked_unsupported", "blocked_upstream"}
        ),
    }


__all__ = [
    "ANALYSIS_INTENT_KINDS",
    "ANALYSIS_INTENT_SUPPORT_STATES",
    "AnalysisInputIntentV1",
    "AnalysisNodeIntentV1",
    "AnalysisOutputIntentV1",
    "AnalysisSelectorIntentV1",
    "ScientificToolchainContractError",
    "ScientificToolchainPlanV1",
    "build_scientific_toolchain_plan",
    "project_scientific_toolchain_frontier",
    "scientific_toolchain_plan_body",
]
