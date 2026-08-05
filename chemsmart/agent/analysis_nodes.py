"""Typed analysis nodes over trusted computational-chemistry results.

The scientific workflow records calculations.  This module adds a separate,
task-bound analysis overlay for deterministic extraction and postprocessing.
Analysis nodes therefore share workflow identity without pretending to be CLI
jobs, consuming execution approval, or accepting model-authored code and file
paths.

The parser boundary is deliberately program-neutral.  A semantic selector is
valid independently of any one engine; a registered adapter decides whether a
particular result artifact can satisfy it.  PySCF is the first adapter and
wraps the existing, provenance-checking HDF5 extraction implementation.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Iterable, Protocol, Sequence, runtime_checkable

from chemsmart.agent._contracts import (
    ContractError,
    TrustedArtifactRefV1,
    canonical_sha256,
    require_sha256,
)
from chemsmart.analysis.result_quantities import (
    SUPPORTED_PYSCF_SELECTORS,
    QuantityContractError,
    QuantityExtractionError,
    QuantityExtractionReceiptV1,
    QuantitySelectorV1,
    ResultQuantityExtractionRequestV1,
    extract_pyscf_quantities,
)


_SYMBOL_RE = re.compile(r"^[a-z][a-z0-9_.:-]{0,127}$")
_DESCRIPTOR_RE = re.compile(r"^[A-Za-z][A-Za-z0-9_.:-]{0,255}$")

ANALYSIS_KINDS = (
    "claim_rendering",
    "quantity_expression",
    "result_extraction",
    "thermochemistry",
)
ANALYSIS_SUPPORT_STATES = ("blocked_unsupported", "resolvable")
ANALYSIS_EXECUTION_STATES = (
    "blocked_unsupported",
    "derived",
    "failed",
)


class AnalysisContractError(ContractError):
    """Raised when an analysis overlay is internally inconsistent."""


class UnsupportedResultParserError(AnalysisContractError):
    """A registered deterministic parser cannot satisfy a selector request."""

    def __init__(self, message: str, *, rule_id: str) -> None:
        super().__init__(message)
        self.rule_id = rule_id


def _require_symbol(value: str, field_name: str) -> str:
    candidate = str(value or "")
    if candidate != candidate.strip().lower() or not _SYMBOL_RE.fullmatch(
        candidate
    ):
        raise AnalysisContractError(
            f"{field_name} must be a canonical public symbol"
        )
    return candidate


def _require_descriptor(value: str, field_name: str) -> str:
    candidate = str(value or "")
    if candidate != candidate.strip() or not _DESCRIPTOR_RE.fullmatch(candidate):
        raise AnalysisContractError(
            f"{field_name} must be a stable parser descriptor"
        )
    return candidate


def _require_sorted_unique(values: tuple[str, ...], field_name: str) -> None:
    if values != tuple(sorted(set(values))):
        raise AnalysisContractError(f"{field_name} must be sorted and unique")


@dataclass(frozen=True)
class ResultQuantitySelectorV1:
    """Program-neutral semantic quantity selector.

    Syntax validation here intentionally does not consult a program inventory.
    Adapter resolution later determines whether, for example, a PySCF HDF5 or
    an ORCA output can provide the requested semantic quantity.
    """

    schema_version: str
    quantity_id: str
    selector: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.result-quantity-selector.v1":
            raise AnalysisContractError("unsupported result selector schema")
        _require_symbol(self.quantity_id, "quantity_id")
        _require_symbol(self.selector, "selector")


@dataclass(frozen=True)
class AnalysisInputRefV1:
    """Exact program-artifact or upstream-analysis input binding."""

    input_id: str
    source_kind: str
    producer_node_id: str
    producer_output_id: str
    source_id: str
    source_sha256: str
    source_receipt_sha256: str
    program: str = ""

    def __post_init__(self) -> None:
        _require_symbol(self.input_id, "input_id")
        if self.source_kind not in {"analysis_quantity", "program_artifact"}:
            raise AnalysisContractError(
                "source_kind must be analysis_quantity or program_artifact"
            )
        _require_symbol(self.producer_node_id, "producer_node_id")
        _require_symbol(self.producer_output_id, "producer_output_id")
        _require_symbol(self.source_id, "source_id")
        require_sha256(self.source_sha256, "source_sha256")
        require_sha256(
            self.source_receipt_sha256, "source_receipt_sha256"
        )
        if self.source_kind == "program_artifact":
            _require_symbol(self.program, "program")
        elif self.program:
            raise AnalysisContractError(
                "program applies only to a program_artifact input"
            )


@dataclass(frozen=True)
class AnalysisOutputSpecV1:
    """Declared quantity and unit produced by one analysis node."""

    quantity_id: str
    unit: str
    dimension: tuple[int, ...]
    data_kind: str = "scalar"

    def __post_init__(self) -> None:
        _require_symbol(self.quantity_id, "quantity_id")
        if not str(self.unit).strip():
            raise AnalysisContractError("analysis output unit must not be empty")
        object.__setattr__(self, "dimension", tuple(self.dimension))
        if len(self.dimension) != 6 or not all(
            isinstance(value, int) for value in self.dimension
        ):
            raise AnalysisContractError(
                "analysis output dimension must contain six integers"
            )
        if self.data_kind not in {
            "integer",
            "matrix",
            "scalar",
            "text",
            "text_vector",
            "vector",
        }:
            raise AnalysisContractError("unsupported analysis output data kind")


@dataclass(frozen=True)
class AnalysisNodeSpecV1:
    """One deterministic analysis stage bound to a scientific workflow."""

    schema_version: str
    node_id: str
    task_spec_sha256: str
    workflow_id: str
    scientific_workflow_plan_sha256: str
    analysis_kind: str
    inputs: tuple[AnalysisInputRefV1, ...]
    selectors: tuple[ResultQuantitySelectorV1, ...]
    outputs: tuple[AnalysisOutputSpecV1, ...]
    evidence_requirements: tuple[str, ...]
    support_state: str
    blocked_reason: str
    spec_sha256: str

    def __post_init__(self) -> None:
        object.__setattr__(self, "inputs", tuple(self.inputs))
        object.__setattr__(self, "selectors", tuple(self.selectors))
        object.__setattr__(self, "outputs", tuple(self.outputs))
        object.__setattr__(
            self, "evidence_requirements", tuple(self.evidence_requirements)
        )
        if self.schema_version != "chemsmart.analysis-node-spec.v1":
            raise AnalysisContractError("unsupported analysis node schema")
        _require_symbol(self.node_id, "node_id")
        require_sha256(self.task_spec_sha256, "task_spec_sha256")
        _require_symbol(self.workflow_id, "workflow_id")
        require_sha256(
            self.scientific_workflow_plan_sha256,
            "scientific_workflow_plan_sha256",
        )
        if self.analysis_kind not in ANALYSIS_KINDS:
            raise AnalysisContractError("unsupported analysis_kind")
        if self.support_state not in ANALYSIS_SUPPORT_STATES:
            raise AnalysisContractError("unsupported analysis support_state")
        if self.support_state == "blocked_unsupported":
            if not str(self.blocked_reason).strip():
                raise AnalysisContractError(
                    "blocked analysis node must preserve its reason"
                )
        elif str(self.blocked_reason).strip():
            raise AnalysisContractError(
                "blocked_reason applies only to blocked_unsupported nodes"
            )
        if not self.outputs:
            raise AnalysisContractError("analysis node must declare outputs")
        input_ids = tuple(item.input_id for item in self.inputs)
        selector_ids = tuple(item.quantity_id for item in self.selectors)
        output_ids = tuple(item.quantity_id for item in self.outputs)
        _require_sorted_unique(input_ids, "analysis input IDs")
        _require_sorted_unique(selector_ids, "selector quantity IDs")
        _require_sorted_unique(output_ids, "analysis output quantity IDs")
        _require_sorted_unique(
            self.evidence_requirements, "evidence requirements"
        )
        for requirement in self.evidence_requirements:
            _require_symbol(requirement, "evidence requirement")
        if any(
            item.producer_node_id == self.node_id for item in self.inputs
        ):
            raise AnalysisContractError("analysis node cannot consume itself")
        if self.analysis_kind == "result_extraction":
            if len(self.inputs) != 1 or self.inputs[0].source_kind != (
                "program_artifact"
            ):
                raise AnalysisContractError(
                    "result extraction requires one program artifact input"
                )
            if not self.selectors:
                raise AnalysisContractError(
                    "result extraction requires semantic selectors"
                )
            if selector_ids != output_ids:
                raise AnalysisContractError(
                    "result extraction selectors must match declared outputs"
                )
        elif self.selectors:
            raise AnalysisContractError(
                "semantic selectors apply only to result extraction nodes"
            )
        body = analysis_node_spec_body(self)
        if self.spec_sha256 != canonical_sha256(body):
            raise AnalysisContractError("analysis node spec digest mismatch")


def analysis_node_spec_body(node: AnalysisNodeSpecV1) -> dict[str, object]:
    return {
        "schema_version": node.schema_version,
        "node_id": node.node_id,
        "task_spec_sha256": node.task_spec_sha256,
        "workflow_id": node.workflow_id,
        "scientific_workflow_plan_sha256": (
            node.scientific_workflow_plan_sha256
        ),
        "analysis_kind": node.analysis_kind,
        "inputs": node.inputs,
        "selectors": node.selectors,
        "outputs": node.outputs,
        "evidence_requirements": node.evidence_requirements,
        "support_state": node.support_state,
        "blocked_reason": node.blocked_reason,
    }


def build_analysis_node_spec(
    *,
    node_id: str,
    task_spec_sha256: str,
    workflow_id: str,
    scientific_workflow_plan_sha256: str,
    analysis_kind: str,
    inputs: Iterable[AnalysisInputRefV1],
    outputs: Iterable[AnalysisOutputSpecV1],
    selectors: Iterable[ResultQuantitySelectorV1] = (),
    evidence_requirements: Iterable[str] = (),
    support_state: str = "resolvable",
    blocked_reason: str = "",
) -> AnalysisNodeSpecV1:
    body = {
        "schema_version": "chemsmart.analysis-node-spec.v1",
        "node_id": node_id,
        "task_spec_sha256": task_spec_sha256,
        "workflow_id": workflow_id,
        "scientific_workflow_plan_sha256": (
            scientific_workflow_plan_sha256
        ),
        "analysis_kind": analysis_kind,
        "inputs": tuple(sorted(inputs, key=lambda item: item.input_id)),
        "selectors": tuple(
            sorted(selectors, key=lambda item: item.quantity_id)
        ),
        "outputs": tuple(sorted(outputs, key=lambda item: item.quantity_id)),
        "evidence_requirements": tuple(sorted(set(evidence_requirements))),
        "support_state": support_state,
        "blocked_reason": blocked_reason,
    }
    return AnalysisNodeSpecV1(
        **body, spec_sha256=canonical_sha256(body)
    )


def _canonical_analysis_node_order(
    nodes: Sequence[AnalysisNodeSpecV1],
) -> tuple[str, ...]:
    by_id = {node.node_id: node for node in nodes}
    if len(by_id) != len(nodes):
        raise AnalysisContractError("analysis node IDs must be unique")
    produced: dict[str, str] = {}
    dependencies: dict[str, set[str]] = {
        node.node_id: set() for node in nodes
    }
    children: dict[str, set[str]] = {node.node_id: set() for node in nodes}
    for node in nodes:
        for output in node.outputs:
            if output.quantity_id in produced:
                raise AnalysisContractError(
                    "analysis output quantity IDs must be globally unique"
                )
            produced[output.quantity_id] = node.node_id
        for item in node.inputs:
            if item.source_kind == "program_artifact":
                if item.producer_node_id in by_id:
                    raise AnalysisContractError(
                        "program artifact producer collides with analysis node"
                    )
                continue
            producer = by_id.get(item.producer_node_id)
            if producer is None:
                raise AnalysisContractError(
                    "analysis input references an unknown analysis producer"
                )
            producer_outputs = {
                output.quantity_id for output in producer.outputs
            }
            if item.producer_output_id not in producer_outputs:
                raise AnalysisContractError(
                    "analysis input references an unknown producer output"
                )
            if item.source_id != item.producer_output_id:
                raise AnalysisContractError(
                    "analysis quantity source must name the producer output"
                )
            dependencies[node.node_id].add(producer.node_id)
            children[producer.node_id].add(node.node_id)
    ready = sorted(
        node_id for node_id, values in dependencies.items() if not values
    )
    ordered: list[str] = []
    while ready:
        node_id = ready.pop(0)
        ordered.append(node_id)
        for child in sorted(children[node_id]):
            dependencies[child].discard(node_id)
            if not dependencies[child] and child not in ordered + ready:
                ready.append(child)
                ready.sort()
    if len(ordered) != len(nodes):
        raise AnalysisContractError("scientific analysis graph contains a cycle")
    return tuple(ordered)


@dataclass(frozen=True)
class ScientificAnalysisPlanV1:
    """Canonical analysis overlay bound to one scientific workflow plan."""

    schema_version: str
    analysis_plan_id: str
    task_spec_sha256: str
    workflow_id: str
    scientific_workflow_plan_sha256: str
    nodes: tuple[AnalysisNodeSpecV1, ...]
    required_output_quantity_ids: tuple[str, ...]
    status: str
    plan_sha256: str

    def __post_init__(self) -> None:
        object.__setattr__(self, "nodes", tuple(self.nodes))
        object.__setattr__(
            self,
            "required_output_quantity_ids",
            tuple(self.required_output_quantity_ids),
        )
        if self.schema_version != "chemsmart.scientific-analysis-plan.v1":
            raise AnalysisContractError("unsupported scientific analysis plan")
        _require_symbol(self.analysis_plan_id, "analysis_plan_id")
        require_sha256(self.task_spec_sha256, "task_spec_sha256")
        _require_symbol(self.workflow_id, "workflow_id")
        require_sha256(
            self.scientific_workflow_plan_sha256,
            "scientific_workflow_plan_sha256",
        )
        if not self.nodes:
            raise AnalysisContractError("scientific analysis plan needs nodes")
        for node in self.nodes:
            if (
                node.task_spec_sha256 != self.task_spec_sha256
                or node.workflow_id != self.workflow_id
                or node.scientific_workflow_plan_sha256
                != self.scientific_workflow_plan_sha256
            ):
                raise AnalysisContractError(
                    "analysis node differs from its bound workflow identity"
                )
        canonical_order = _canonical_analysis_node_order(self.nodes)
        if tuple(node.node_id for node in self.nodes) != canonical_order:
            raise AnalysisContractError(
                "analysis nodes must use canonical topological order"
            )
        _require_sorted_unique(
            self.required_output_quantity_ids,
            "required output quantity IDs",
        )
        for quantity_id in self.required_output_quantity_ids:
            _require_symbol(quantity_id, "required output quantity ID")
        produced = {
            output.quantity_id
            for node in self.nodes
            for output in node.outputs
        }
        missing = sorted(set(self.required_output_quantity_ids) - produced)
        if missing:
            raise AnalysisContractError(
                f"no analysis node preserves required output(s) {missing}"
            )
        if self.status != "planned":
            raise AnalysisContractError(
                "scientific analysis plan status must be planned"
            )
        if self.plan_sha256 != canonical_sha256(
            scientific_analysis_plan_body(self)
        ):
            raise AnalysisContractError("scientific analysis plan digest mismatch")


def scientific_analysis_plan_body(
    plan: ScientificAnalysisPlanV1,
) -> dict[str, object]:
    return {
        "schema_version": plan.schema_version,
        "analysis_plan_id": plan.analysis_plan_id,
        "task_spec_sha256": plan.task_spec_sha256,
        "workflow_id": plan.workflow_id,
        "scientific_workflow_plan_sha256": (
            plan.scientific_workflow_plan_sha256
        ),
        "nodes": plan.nodes,
        "required_output_quantity_ids": plan.required_output_quantity_ids,
        "status": plan.status,
    }


def build_scientific_analysis_plan(
    *,
    analysis_plan_id: str,
    task_spec_sha256: str,
    workflow_id: str,
    scientific_workflow_plan_sha256: str,
    nodes: Iterable[AnalysisNodeSpecV1],
    required_output_quantity_ids: Iterable[str] = (),
) -> ScientificAnalysisPlanV1:
    node_values = tuple(nodes)
    order = _canonical_analysis_node_order(node_values)
    by_id = {node.node_id: node for node in node_values}
    body = {
        "schema_version": "chemsmart.scientific-analysis-plan.v1",
        "analysis_plan_id": analysis_plan_id,
        "task_spec_sha256": task_spec_sha256,
        "workflow_id": workflow_id,
        "scientific_workflow_plan_sha256": (
            scientific_workflow_plan_sha256
        ),
        "nodes": tuple(by_id[node_id] for node_id in order),
        "required_output_quantity_ids": tuple(
            sorted(set(required_output_quantity_ids))
        ),
        "status": "planned",
    }
    return ScientificAnalysisPlanV1(
        **body, plan_sha256=canonical_sha256(body)
    )


@runtime_checkable
class ResultParserAdapterV1(Protocol):
    """Registered, deterministic reader for one program result family."""

    parser_id: str
    program: str
    artifact_kinds: tuple[str, ...]
    supported_selectors: tuple[str, ...]

    def extract(
        self,
        *,
        artifact: TrustedArtifactRefV1,
        selectors: tuple[ResultQuantitySelectorV1, ...],
    ) -> QuantityExtractionReceiptV1: ...


@dataclass(frozen=True)
class PySCFResultParserAdapterV1:
    """Adapter over ChemSmart's provenance-validating PySCF HDF5 reader."""

    parser_id: str = "chemsmart.io.pyscf.output.PySCFOutput"
    program: str = "pyscf"
    artifact_kinds: tuple[str, ...] = ("pyscf_hdf5",)
    supported_selectors: tuple[str, ...] = tuple(
        sorted(SUPPORTED_PYSCF_SELECTORS)
    )

    def extract(
        self,
        *,
        artifact: TrustedArtifactRefV1,
        selectors: tuple[ResultQuantitySelectorV1, ...],
    ) -> QuantityExtractionReceiptV1:
        if artifact.kind not in self.artifact_kinds:
            raise UnsupportedResultParserError(
                "PySCF parser does not support this artifact kind",
                rule_id="result-parser.artifact-kind-unsupported",
            )
        unsupported = sorted(
            {item.selector for item in selectors}
            - set(self.supported_selectors)
        )
        if unsupported:
            raise UnsupportedResultParserError(
                f"PySCF parser does not support selectors {unsupported}",
                rule_id="result-parser.selector-unsupported",
            )
        legacy_selectors = tuple(
            QuantitySelectorV1(
                quantity_id=item.quantity_id,
                selector=item.selector,
            )
            for item in selectors
        )
        request = ResultQuantityExtractionRequestV1(
            schema_version="chemsmart.quantity-extraction-request.v1",
            artifact_id=artifact.artifact_id,
            artifact_sha256=artifact.sha256,
            program=self.program,
            selectors=legacy_selectors,
        )
        receipt = extract_pyscf_quantities(
            request=request,
            artifact_path=artifact.path,
        )
        if (
            receipt.program != self.program
            or receipt.parser_id != self.parser_id
            or receipt.artifact_id != artifact.artifact_id
            or receipt.artifact_sha256 != artifact.sha256
        ):
            raise QuantityExtractionError(
                "PySCF parser receipt differs from the adapter binding"
            )
        observed = tuple(
            sorted(quantity.quantity_id for quantity in receipt.quantities)
        )
        expected = tuple(sorted(item.quantity_id for item in selectors))
        if observed != expected:
            raise QuantityExtractionError(
                "PySCF parser receipt differs from requested selectors"
            )
        return receipt


class ResultParserRegistryV1:
    """Host-owned registry of result parser adapters."""

    def __init__(
        self, adapters: Iterable[ResultParserAdapterV1] = ()
    ) -> None:
        self._adapters: dict[
            tuple[str, str], ResultParserAdapterV1
        ] = {}
        for adapter in adapters:
            self.register(adapter)

    def register(self, adapter: ResultParserAdapterV1) -> None:
        parser_id = _require_descriptor(adapter.parser_id, "parser_id")
        program = _require_symbol(adapter.program, "program")
        artifact_kinds = tuple(adapter.artifact_kinds)
        supported = tuple(adapter.supported_selectors)
        _require_sorted_unique(artifact_kinds, "adapter artifact kinds")
        _require_sorted_unique(supported, "adapter selectors")
        if not artifact_kinds or not supported:
            raise AnalysisContractError(
                "result parser adapter must declare artifacts and selectors"
            )
        for kind in artifact_kinds:
            _require_symbol(kind, "artifact kind")
        for selector in supported:
            _require_symbol(selector, "supported selector")
        for kind in artifact_kinds:
            key = (program, kind)
            if key in self._adapters:
                raise AnalysisContractError(
                    "duplicate result parser program/artifact registration"
                )
            self._adapters[key] = adapter
        if parser_id != adapter.parser_id:
            raise AnalysisContractError("parser_id must already be canonical")

    @property
    def registry_sha256(self) -> str:
        unique_adapters = {
            (adapter.program, adapter.parser_id): adapter
            for adapter in self._adapters.values()
        }
        descriptors = tuple(
            {
                "program": adapter.program,
                "artifact_kinds": tuple(adapter.artifact_kinds),
                "parser_id": adapter.parser_id,
                "supported_selectors": tuple(adapter.supported_selectors),
            }
            for adapter in (
                unique_adapters[key] for key in sorted(unique_adapters)
            )
        )
        return canonical_sha256(
            {
                "schema_version": "chemsmart.result-parser-registry.v1",
                "adapters": descriptors,
            }
        )

    def resolve(
        self,
        *,
        program: str,
        artifact_kind: str,
        selectors: Sequence[ResultQuantitySelectorV1],
    ) -> ResultParserAdapterV1:
        program_id = _require_symbol(program, "program")
        artifact_id = _require_symbol(artifact_kind, "artifact_kind")
        adapter = self._adapters.get((program_id, artifact_id))
        if adapter is None:
            raise UnsupportedResultParserError(
                "no result parser is registered for the program artifact",
                rule_id="result-parser.adapter-unavailable",
            )
        unsupported = sorted(
            {item.selector for item in selectors}
            - set(adapter.supported_selectors)
        )
        if unsupported:
            raise UnsupportedResultParserError(
                f"registered parser does not support selectors {unsupported}",
                rule_id="result-parser.selector-unsupported",
            )
        return adapter


def build_default_result_parser_registry() -> ResultParserRegistryV1:
    return ResultParserRegistryV1((PySCFResultParserAdapterV1(),))


@dataclass(frozen=True)
class AnalysisOutputQuantityRefV1:
    """Digest-bound output from a component analysis receipt."""

    quantity_id: str
    value_sha256: str
    unit: str
    dimension: tuple[int, ...]
    data_kind: str
    source_receipt_sha256: str

    def __post_init__(self) -> None:
        _require_symbol(self.quantity_id, "quantity_id")
        require_sha256(self.value_sha256, "value_sha256")
        require_sha256(
            self.source_receipt_sha256, "source_receipt_sha256"
        )
        if not str(self.unit).strip():
            raise AnalysisContractError("analysis output unit must not be empty")
        object.__setattr__(self, "dimension", tuple(self.dimension))
        if len(self.dimension) != 6 or not all(
            isinstance(value, int) for value in self.dimension
        ):
            raise AnalysisContractError(
                "analysis output dimension must contain six integers"
            )
        if self.data_kind not in {
            "integer",
            "matrix",
            "scalar",
            "text",
            "text_vector",
            "vector",
        }:
            raise AnalysisContractError("unsupported analysis output data kind")


@dataclass(frozen=True)
class AnalysisExecutionFindingV1:
    rule_id: str
    severity: str
    message: str

    def __post_init__(self) -> None:
        _require_symbol(self.rule_id, "rule_id")
        if self.severity not in {"error", "info", "warning"}:
            raise AnalysisContractError("unsupported analysis finding severity")
        if not str(self.message).strip():
            raise AnalysisContractError("analysis finding message is required")


@dataclass(frozen=True)
class AnalysisExecutionReceiptV1:
    """Node-level outcome referencing existing component receipts."""

    schema_version: str
    task_spec_sha256: str
    scientific_workflow_plan_sha256: str
    analysis_plan_sha256: str
    node_id: str
    analysis_node_spec_sha256: str
    input_source_sha256s: tuple[str, ...]
    component_receipt_sha256s: tuple[str, ...]
    outputs: tuple[AnalysisOutputQuantityRefV1, ...]
    findings: tuple[AnalysisExecutionFindingV1, ...]
    status: str
    receipt_sha256: str

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "input_source_sha256s", tuple(self.input_source_sha256s)
        )
        object.__setattr__(
            self,
            "component_receipt_sha256s",
            tuple(self.component_receipt_sha256s),
        )
        object.__setattr__(self, "outputs", tuple(self.outputs))
        object.__setattr__(self, "findings", tuple(self.findings))
        if self.schema_version != "chemsmart.analysis-execution-receipt.v1":
            raise AnalysisContractError("unsupported analysis execution receipt")
        for value, name in (
            (self.task_spec_sha256, "task_spec_sha256"),
            (
                self.scientific_workflow_plan_sha256,
                "scientific_workflow_plan_sha256",
            ),
            (self.analysis_plan_sha256, "analysis_plan_sha256"),
            (self.analysis_node_spec_sha256, "analysis_node_spec_sha256"),
        ):
            require_sha256(value, name)
        _require_symbol(self.node_id, "node_id")
        _require_sorted_unique(
            self.input_source_sha256s, "input source digests"
        )
        _require_sorted_unique(
            self.component_receipt_sha256s, "component receipt digests"
        )
        for value in (
            self.input_source_sha256s + self.component_receipt_sha256s
        ):
            require_sha256(value, "analysis receipt dependency")
        output_ids = tuple(item.quantity_id for item in self.outputs)
        _require_sorted_unique(output_ids, "analysis execution outputs")
        finding_keys = tuple(
            (item.rule_id, item.severity, item.message)
            for item in self.findings
        )
        if finding_keys != tuple(sorted(set(finding_keys))):
            raise AnalysisContractError(
                "analysis findings must be sorted and unique"
            )
        if self.status not in ANALYSIS_EXECUTION_STATES:
            raise AnalysisContractError("unsupported analysis execution status")
        error_findings = tuple(
            item for item in self.findings if item.severity == "error"
        )
        if self.status == "derived":
            if not self.outputs or not self.component_receipt_sha256s:
                raise AnalysisContractError(
                    "derived analysis requires outputs and component receipts"
                )
            if error_findings:
                raise AnalysisContractError(
                    "derived analysis cannot retain error findings"
                )
        else:
            if self.outputs:
                raise AnalysisContractError(
                    "blocked or failed analysis cannot claim outputs"
                )
            if not error_findings:
                raise AnalysisContractError(
                    "blocked or failed analysis requires an error finding"
                )
        if self.receipt_sha256 != canonical_sha256(
            analysis_execution_receipt_body(self)
        ):
            raise AnalysisContractError("analysis execution receipt mismatch")


def analysis_execution_receipt_body(
    receipt: AnalysisExecutionReceiptV1,
) -> dict[str, object]:
    return {
        "schema_version": receipt.schema_version,
        "task_spec_sha256": receipt.task_spec_sha256,
        "scientific_workflow_plan_sha256": (
            receipt.scientific_workflow_plan_sha256
        ),
        "analysis_plan_sha256": receipt.analysis_plan_sha256,
        "node_id": receipt.node_id,
        "analysis_node_spec_sha256": receipt.analysis_node_spec_sha256,
        "input_source_sha256s": receipt.input_source_sha256s,
        "component_receipt_sha256s": (
            receipt.component_receipt_sha256s
        ),
        "outputs": receipt.outputs,
        "findings": receipt.findings,
        "status": receipt.status,
    }


def build_analysis_execution_receipt(
    *,
    analysis_plan_sha256: str,
    node: AnalysisNodeSpecV1,
    status: str,
    component_receipt_sha256s: Iterable[str] = (),
    outputs: Iterable[AnalysisOutputQuantityRefV1] = (),
    findings: Iterable[AnalysisExecutionFindingV1] = (),
) -> AnalysisExecutionReceiptV1:
    output_values = tuple(sorted(outputs, key=lambda item: item.quantity_id))
    if status == "derived":
        declared = {item.quantity_id: item for item in node.outputs}
        observed = {item.quantity_id: item for item in output_values}
        if set(observed) != set(declared):
            raise AnalysisContractError(
                "derived outputs differ from the analysis node specification"
            )
        for quantity_id, value in observed.items():
            expected = declared[quantity_id]
            if (
                value.unit != expected.unit
                or value.dimension != expected.dimension
                or value.data_kind != expected.data_kind
            ):
                raise AnalysisContractError(
                    "derived output semantics differ from the node specification"
                )
    body = {
        "schema_version": "chemsmart.analysis-execution-receipt.v1",
        "task_spec_sha256": node.task_spec_sha256,
        "scientific_workflow_plan_sha256": (
            node.scientific_workflow_plan_sha256
        ),
        "analysis_plan_sha256": analysis_plan_sha256,
        "node_id": node.node_id,
        "analysis_node_spec_sha256": node.spec_sha256,
        "input_source_sha256s": tuple(
            sorted({item.source_sha256 for item in node.inputs})
        ),
        "component_receipt_sha256s": tuple(
            sorted(set(component_receipt_sha256s))
        ),
        "outputs": output_values,
        "findings": tuple(
            sorted(
                findings,
                key=lambda item: (item.rule_id, item.severity, item.message),
            )
        ),
        "status": status,
    }
    return AnalysisExecutionReceiptV1(
        **body, receipt_sha256=canonical_sha256(body)
    )


def _finding(rule_id: str, message: str) -> AnalysisExecutionFindingV1:
    return AnalysisExecutionFindingV1(
        rule_id=rule_id,
        severity="error",
        message=message,
    )


def execute_result_extraction_node(
    *,
    analysis_plan_sha256: str,
    node: AnalysisNodeSpecV1,
    artifact: TrustedArtifactRefV1,
    registry: ResultParserRegistryV1,
) -> tuple[AnalysisExecutionReceiptV1, QuantityExtractionReceiptV1 | None]:
    """Execute one extraction node through a registered host parser.

    Unsupported capability is an honest blocked result.  A registered parser
    rejecting the exact artifact is a failed observation.  Neither path can
    manufacture a numerical output.
    """

    require_sha256(analysis_plan_sha256, "analysis_plan_sha256")
    if node.analysis_kind != "result_extraction":
        raise AnalysisContractError(
            "result extraction executor requires a result_extraction node"
        )
    source = node.inputs[0]
    if node.support_state == "blocked_unsupported":
        finding = _finding(
            "analysis-node.declared-unsupported",
            node.blocked_reason,
        )
        return (
            build_analysis_execution_receipt(
                analysis_plan_sha256=analysis_plan_sha256,
                node=node,
                status="blocked_unsupported",
                findings=(finding,),
            ),
            None,
        )
    if (
        artifact.artifact_id != source.source_id
        or artifact.sha256 != source.source_sha256
    ):
        finding = _finding(
            "result-parser.artifact-binding-mismatch",
            "trusted artifact differs from the analysis input binding",
        )
        return (
            build_analysis_execution_receipt(
                analysis_plan_sha256=analysis_plan_sha256,
                node=node,
                status="failed",
                findings=(finding,),
            ),
            None,
        )
    try:
        adapter = registry.resolve(
            program=source.program,
            artifact_kind=artifact.kind,
            selectors=node.selectors,
        )
    except UnsupportedResultParserError as exc:
        finding = _finding(exc.rule_id, str(exc))
        return (
            build_analysis_execution_receipt(
                analysis_plan_sha256=analysis_plan_sha256,
                node=node,
                status="blocked_unsupported",
                findings=(finding,),
            ),
            None,
        )
    try:
        extraction = adapter.extract(
            artifact=artifact,
            selectors=node.selectors,
        )
    except UnsupportedResultParserError as exc:
        finding = _finding(exc.rule_id, str(exc))
        return (
            build_analysis_execution_receipt(
                analysis_plan_sha256=analysis_plan_sha256,
                node=node,
                status="blocked_unsupported",
                findings=(finding,),
            ),
            None,
        )
    except (QuantityContractError, QuantityExtractionError, OSError):
        finding = _finding(
            "result-parser.extraction-failed",
            "registered parser rejected the bound result artifact",
        )
        return (
            build_analysis_execution_receipt(
                analysis_plan_sha256=analysis_plan_sha256,
                node=node,
                status="failed",
                findings=(finding,),
            ),
            None,
        )
    try:
        output_refs = tuple(
            AnalysisOutputQuantityRefV1(
                quantity_id=quantity.quantity_id,
                value_sha256=quantity.value_sha256,
                unit=quantity.unit,
                dimension=quantity.dimension,
                data_kind=quantity.data_kind,
                source_receipt_sha256=extraction.receipt_sha256,
            )
            for quantity in extraction.quantities
        )
        execution = build_analysis_execution_receipt(
            analysis_plan_sha256=analysis_plan_sha256,
            node=node,
            status="derived",
            component_receipt_sha256s=(extraction.receipt_sha256,),
            outputs=output_refs,
        )
    except AnalysisContractError:
        finding = _finding(
            "result-parser.output-contract-mismatch",
            "parser output differs from the declared analysis quantities",
        )
        execution = build_analysis_execution_receipt(
            analysis_plan_sha256=analysis_plan_sha256,
            node=node,
            status="failed",
            component_receipt_sha256s=(extraction.receipt_sha256,),
            findings=(finding,),
        )
    return execution, extraction


__all__ = [
    "ANALYSIS_EXECUTION_STATES",
    "ANALYSIS_KINDS",
    "ANALYSIS_SUPPORT_STATES",
    "AnalysisContractError",
    "AnalysisExecutionFindingV1",
    "AnalysisExecutionReceiptV1",
    "AnalysisInputRefV1",
    "AnalysisNodeSpecV1",
    "AnalysisOutputQuantityRefV1",
    "AnalysisOutputSpecV1",
    "PySCFResultParserAdapterV1",
    "ResultParserAdapterV1",
    "ResultParserRegistryV1",
    "ResultQuantitySelectorV1",
    "ScientificAnalysisPlanV1",
    "UnsupportedResultParserError",
    "analysis_execution_receipt_body",
    "analysis_node_spec_body",
    "build_analysis_execution_receipt",
    "build_analysis_node_spec",
    "build_default_result_parser_registry",
    "build_scientific_analysis_plan",
    "execute_result_extraction_node",
    "scientific_analysis_plan_body",
]
