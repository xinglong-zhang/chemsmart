"""Paper-to-workflow experiments over the normal ChemSmart agent path.

This module is intentionally a thin experiment layer.  It does not search the
literature, embed benchmark answers, or introduce a second model loop.  A host
selects a source-complete paper packet, supplies an exact coordinate artifact,
and calls :func:`run_live_agent_session` separately.  The helpers here make the
public prompt, construct the P0/P1/P2 predecessor-context arms, and grade only
observable provider messages and typed tool calls.

Paper-specific source text and hidden values belong in a private run directory,
not in this module or the model tool surface.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
import json
import math
from pathlib import Path
import re
from typing import TYPE_CHECKING, Any, Mapping, Sequence

import numpy as np
import yaml

from chemsmart.agent._contracts import (
    ContractError,
    canonical_data,
    canonical_json,
    canonical_sha256,
    file_sha256,
    require_identifier,
    require_sha256,
)
from chemsmart.agent.dependency_context import (
    ContextSelectionReceiptV1,
    TaskDependencyContextV1,
    build_predecessor_evidence_ref,
    build_task_dependency_context_policy,
    select_task_dependency_context,
)
from chemsmart.agent.projects import (
    ProjectDocumentV1,
    project_effective_section_settings,
)
from chemsmart.agent.workflows import ScientificWorkflowPlanV2
from chemsmart.analysis.quantity_expressions import normalize_numeric_value

if TYPE_CHECKING:
    from chemsmart.agent.analysis_claims import AnalysisClaimRecordV1
    from chemsmart.agent.live_session import LiveAgentSessionResultV1


_ELEMENT = re.compile(r"^[A-Z][a-z]?$", re.ASCII)
_CONTEXT_ARMS = ("p0", "p1", "p2")
_GENERIC_CLAIM_TOKENS = frozenset(
    {
        "calculated",
        "claim",
        "computed",
        "display",
        "molar",
        "reported",
        "standard",
        "value",
    }
)
_SCIENTIFIC_IDENTIFIER_TOKEN_ALIASES = {
    "eh": "hartree",
    "ha": "hartree",
}


@dataclass(frozen=True)
class MultiRecordXyzObservationV1:
    """One exact coordinate record extracted from a multi-XYZ collection."""

    schema_version: str
    record_ordinal: int
    source_sha256: str
    source_comment: str
    atom_count: int
    atom_order: tuple[str, ...]
    coordinate_units: str
    output_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.multi-record-xyz-observation.v1":
            raise ContractError("unsupported multi-record XYZ observation")
        if self.record_ordinal < 1 or self.atom_count < 1:
            raise ContractError("XYZ record ordinal and atom count must be positive")
        require_sha256(self.source_sha256, "source_sha256")
        require_sha256(self.output_sha256, "output_sha256")
        if len(self.atom_order) != self.atom_count:
            raise ContractError("XYZ atom order does not match atom count")
        if self.coordinate_units != "angstrom":
            raise ContractError("multi-record XYZ coordinates must be angstrom")


def extract_multi_xyz_record(
    *, source: str | Path, record_ordinal: int, destination: str | Path
) -> MultiRecordXyzObservationV1:
    """Extract one ordinary XYZ frame without generating or changing coordinates.

    The source format is the common concatenation of ``N/comment/N atom rows``.
    Coordinate lexemes are copied verbatim; only the comment is replaced by a
    provenance-neutral record label so a source-side property cannot leak into a
    black-box task.
    """

    source_path = Path(source).expanduser().resolve()
    destination_path = Path(destination).expanduser().resolve()
    if record_ordinal < 1:
        raise ContractError("record_ordinal must be positive")
    if not source_path.is_file() or source_path.is_symlink():
        raise ContractError("multi-record XYZ source must be a regular file")
    lines = source_path.read_text(encoding="utf-8").splitlines()
    cursor = 0
    observed = 0
    selected: tuple[int, str, tuple[str, ...], tuple[str, ...]] | None = None
    while cursor < len(lines):
        if not lines[cursor].strip():
            cursor += 1
            continue
        try:
            atom_count = int(lines[cursor].strip())
        except ValueError as exc:
            raise ContractError("multi-record XYZ has a malformed atom count") from exc
        if atom_count < 1 or cursor + atom_count + 1 >= len(lines):
            raise ContractError("multi-record XYZ is truncated")
        comment = lines[cursor + 1]
        atom_rows = tuple(lines[cursor + 2 : cursor + 2 + atom_count])
        atom_order: list[str] = []
        for row in atom_rows:
            fields = row.split()
            if len(fields) < 4 or _ELEMENT.fullmatch(fields[0]) is None:
                raise ContractError("multi-record XYZ has a malformed atom row")
            try:
                coordinates = tuple(float(value) for value in fields[1:4])
            except ValueError as exc:
                raise ContractError("multi-record XYZ coordinates are not numeric") from exc
            if not all(math.isfinite(value) for value in coordinates):
                raise ContractError("multi-record XYZ coordinates must be finite")
            atom_order.append(fields[0])
        observed += 1
        if observed == record_ordinal:
            selected = (atom_count, comment, atom_rows, tuple(atom_order))
            break
        cursor += atom_count + 2
    if selected is None:
        raise ContractError("requested XYZ record is absent")

    atom_count, comment, atom_rows, atom_order = selected
    payload = "\n".join(
        (
            str(atom_count),
            f"official multi-XYZ record {record_ordinal}; coordinates in angstrom",
            *atom_rows,
            "",
        )
    ).encode("utf-8")
    destination_path.parent.mkdir(parents=True, exist_ok=True)
    destination_path.write_bytes(payload)
    return MultiRecordXyzObservationV1(
        schema_version="chemsmart.multi-record-xyz-observation.v1",
        record_ordinal=record_ordinal,
        source_sha256=file_sha256(source_path),
        source_comment=comment,
        atom_count=atom_count,
        atom_order=atom_order,
        coordinate_units="angstrom",
        output_sha256=file_sha256(destination_path),
    )


@dataclass(frozen=True)
class PaperModelInputV1:
    """Minimal black-box material visible to the benchmark model.

    Bibliographic metadata, host case identifiers, private answer material,
    graders, and earlier trajectories deliberately have no representation in
    this object.  Provenance remains host-side and is joined to the resulting
    evidence only after the model episode.
    """

    methods_excerpt: str
    coordinate_context: str
    requested_result: str
    execution_policy: str

    def __post_init__(self) -> None:
        for name in (
            "methods_excerpt",
            "coordinate_context",
            "requested_result",
            "execution_policy",
        ):
            if not str(getattr(self, name)).strip():
                raise ContractError(f"model input {name} must not be empty")

    def render(self) -> str:
        return (
            "Use ChemSmart and its normal typed tools to answer the scientific "
            "request below. The workspace contains the exact Cartesian XYZ "
            "artifact(s) described here. Derive the required project YAML, "
            "program-relative calculation workflow, result extraction, and "
            "postprocessing from the supplied Methods text and coordinates. "
            "Do not invent a missing scientific setting or numerical result; "
            "leave an unsupported or underdetermined stage explicit.\n\n"
            "METHODS\n"
            "---\n"
            f"{self.methods_excerpt.strip()}\n\n"
            "COORDINATE CONTEXT\n"
            "---\n"
            f"{self.coordinate_context.strip()}\n\n"
            "REQUESTED RESULT AND CONDITIONS\n"
            "---\n"
            f"{self.requested_result.strip()}\n\n"
            "EXECUTION BOUNDARY\n"
            "---\n"
            f"{self.execution_policy.strip()}"
        )


@dataclass(frozen=True)
class PaperBlackBoxCaseV1:
    """Public model packet plus host-only source identity.

    The case deliberately has no expected workflow or reference value field.
    Those belong to a separate private oracle so they cannot enter the prompt
    accidentally through serialization.
    """

    schema_version: str
    case_id: str
    paper_doi: str
    article_title: str
    source_locator: str
    methods_excerpt: str
    coordinate_sha256: str
    coordinate_description: str
    requested_result: str
    execution_policy: str
    system_scale: str
    case_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.paper-black-box-case.v1":
            raise ContractError("unsupported paper black-box case")
        require_identifier(self.case_id, "case_id")
        require_sha256(self.coordinate_sha256, "coordinate_sha256")
        for name in (
            "paper_doi",
            "article_title",
            "source_locator",
            "methods_excerpt",
            "coordinate_description",
            "requested_result",
            "execution_policy",
        ):
            if not str(getattr(self, name)).strip():
                raise ContractError(f"{name} must not be empty")
        if self.system_scale not in {
            "small_molecule",
            "intermediate_or_larger",
            "diagnostic_unspecified",
        }:
            raise ContractError("paper case system scale is invalid")
        if self.case_sha256 != canonical_sha256(self._body()):
            raise ContractError("paper black-box case digest mismatch")

    def _body(self) -> dict[str, Any]:
        return {
            key: value
            for key, value in self.__dict__.items()
            if key != "case_sha256"
        }

    def public_task(self) -> str:
        """Render only the minimal model-visible scientific input.

        Benchmark dispatchers must call :func:`prepare_benchmark_dispatch`
        instead.  Keeping this renderer answer-free makes it usable for
        diagnostic planning pilots without making those pilots benchmarks.
        """

        return self.model_input().render()

    def model_input(self) -> PaperModelInputV1:
        """Project host evidence onto the strict public input allowlist."""

        return PaperModelInputV1(
            methods_excerpt=self.methods_excerpt,
            coordinate_context=self.coordinate_description,
            requested_result=self.requested_result,
            execution_policy=self.execution_policy,
        )


def build_paper_black_box_case(
    *,
    case_id: str,
    paper_doi: str,
    article_title: str,
    source_locator: str,
    methods_excerpt: str,
    coordinate_sha256: str,
    coordinate_description: str,
    requested_result: str,
    execution_policy: str,
    system_scale: str = "diagnostic_unspecified",
) -> PaperBlackBoxCaseV1:
    body = {
        "schema_version": "chemsmart.paper-black-box-case.v1",
        "case_id": require_identifier(case_id, "case_id"),
        "paper_doi": str(paper_doi).strip(),
        "article_title": str(article_title).strip(),
        "source_locator": str(source_locator).strip(),
        "methods_excerpt": str(methods_excerpt).strip(),
        "coordinate_sha256": require_sha256(
            coordinate_sha256, "coordinate_sha256"
        ),
        "coordinate_description": str(coordinate_description).strip(),
        "requested_result": str(requested_result).strip(),
        "execution_policy": str(execution_policy).strip(),
        "system_scale": str(system_scale).strip(),
    }
    return PaperBlackBoxCaseV1(
        **body, case_sha256=canonical_sha256(body)
    )


@dataclass(frozen=True)
class ReferenceQuantityV1:
    """One private numerical answer and its declared acceptance tolerance."""

    schema_version: str
    quantity_id: str
    expected_value: Any
    unit: str
    absolute_tolerance: float
    relative_tolerance: float
    evidence_locator: str
    accepted_identifiers: tuple[str, ...]
    reference_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.reference-quantity.v1":
            raise ContractError("unsupported reference quantity schema")
        require_identifier(self.quantity_id, "quantity_id")
        object.__setattr__(
            self, "accepted_identifiers", tuple(self.accepted_identifiers)
        )
        if (
            not self.accepted_identifiers
            or self.quantity_id not in self.accepted_identifiers
            or self.accepted_identifiers
            != tuple(sorted(set(self.accepted_identifiers)))
        ):
            raise ContractError(
                "reference accepted identifiers must be sorted, unique, and "
                "include quantity_id"
            )
        for identifier in self.accepted_identifiers:
            require_identifier(identifier, "accepted_identifier")
        try:
            normalize_numeric_value(self.expected_value, self.unit)
            normalize_numeric_value(self.absolute_tolerance, self.unit)
        except (TypeError, ValueError) as exc:
            raise ContractError(f"invalid reference quantity: {exc}") from exc
        if (
            not math.isfinite(self.absolute_tolerance)
            or self.absolute_tolerance < 0.0
            or not math.isfinite(self.relative_tolerance)
            or self.relative_tolerance < 0.0
        ):
            raise ContractError("reference tolerances must be finite and non-negative")
        if not self.evidence_locator.strip():
            raise ContractError("reference quantity requires an evidence locator")
        if self.reference_sha256 != canonical_sha256(self._body()):
            raise ContractError("reference quantity digest mismatch")

    def _body(self) -> dict[str, Any]:
        return {
            key: value
            for key, value in self.__dict__.items()
            if key != "reference_sha256"
        }


def build_reference_quantity(
    *,
    quantity_id: str,
    expected_value: Any,
    unit: str,
    absolute_tolerance: float,
    relative_tolerance: float = 0.0,
    evidence_locator: str,
    accepted_identifiers: Sequence[str] = (),
) -> ReferenceQuantityV1:
    identifiers = tuple(
        sorted(
            {
                require_identifier(quantity_id, "quantity_id"),
                *(
                    require_identifier(item, "accepted_identifier")
                    for item in accepted_identifiers
                ),
            }
        )
    )
    body = {
        "schema_version": "chemsmart.reference-quantity.v1",
        "quantity_id": require_identifier(quantity_id, "quantity_id"),
        "expected_value": expected_value,
        "unit": str(unit).strip(),
        "absolute_tolerance": float(absolute_tolerance),
        "relative_tolerance": float(relative_tolerance),
        "evidence_locator": str(evidence_locator).strip(),
        "accepted_identifiers": identifiers,
    }
    return ReferenceQuantityV1(
        **body, reference_sha256=canonical_sha256(body)
    )


@dataclass(frozen=True)
class PaperAnswerKeyV1:
    """Private, case-bound answers that must exist before model dispatch."""

    schema_version: str
    answer_key_id: str
    case_id: str
    case_sha256: str
    answer_source: str
    quantities: tuple[ReferenceQuantityV1, ...]
    answer_key_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.paper-answer-key.v1":
            raise ContractError("unsupported paper answer-key schema")
        require_identifier(self.answer_key_id, "answer_key_id")
        require_identifier(self.case_id, "case_id")
        require_sha256(self.case_sha256, "case_sha256")
        if not self.answer_source.strip():
            raise ContractError("paper answer key requires a source description")
        object.__setattr__(self, "quantities", tuple(self.quantities))
        quantity_ids = tuple(item.quantity_id for item in self.quantities)
        if not quantity_ids or quantity_ids != tuple(sorted(set(quantity_ids))):
            raise ContractError(
                "paper answer-key quantities must be non-empty, sorted, and unique"
            )
        identifier_owners: dict[str, str] = {}
        for quantity in self.quantities:
            for identifier in quantity.accepted_identifiers:
                owner = identifier_owners.setdefault(identifier, quantity.quantity_id)
                if owner != quantity.quantity_id:
                    raise ContractError(
                        "accepted claim identifiers must identify only one "
                        "reference quantity"
                    )
        if self.answer_key_sha256 != canonical_sha256(self._body()):
            raise ContractError("paper answer-key digest mismatch")

    def _body(self) -> dict[str, Any]:
        return {
            key: value
            for key, value in self.__dict__.items()
            if key != "answer_key_sha256"
        }


def build_paper_answer_key(
    *,
    answer_key_id: str,
    case: PaperBlackBoxCaseV1,
    answer_source: str,
    quantities: Sequence[ReferenceQuantityV1],
) -> PaperAnswerKeyV1:
    ordered = tuple(sorted(quantities, key=lambda item: item.quantity_id))
    body = {
        "schema_version": "chemsmart.paper-answer-key.v1",
        "answer_key_id": require_identifier(answer_key_id, "answer_key_id"),
        "case_id": case.case_id,
        "case_sha256": case.case_sha256,
        "answer_source": str(answer_source).strip(),
        "quantities": ordered,
    }
    return PaperAnswerKeyV1(
        **body, answer_key_sha256=canonical_sha256(body)
    )


@dataclass(frozen=True)
class WorkflowProjectAnswerV1:
    """One expert-prepared project-YAML semantic document."""

    project_id: str
    document: ProjectDocumentV1

    def __post_init__(self) -> None:
        require_identifier(self.project_id, "project_id")
        if not isinstance(self.document, ProjectDocumentV1):
            raise ContractError("workflow project answer must be a typed project")


@dataclass(frozen=True)
class WorkflowDependencyAnswerV1:
    """One data/control dependency in a semantic calculation-analysis DAG."""

    source_node_id: str
    source_output: str
    target_input: str

    def __post_init__(self) -> None:
        require_identifier(self.source_node_id, "source_node_id")
        require_identifier(self.source_output, "source_output")
        require_identifier(self.target_input, "target_input")


@dataclass(frozen=True)
class WorkflowSemanticNodeV1:
    """Program-relative calculation or deterministic analysis semantics."""

    node_id: str
    node_kind: str
    program: str = ""
    jobtype: str = ""
    project_id: str = ""
    operation: str = ""
    input_geometry_sha256s: tuple[str, ...] = ()
    dependencies: tuple[WorkflowDependencyAnswerV1, ...] = ()
    output_semantics: tuple[str, ...] = ()
    semantic_parameters: tuple[tuple[str, Any], ...] = ()

    def __post_init__(self) -> None:
        require_identifier(self.node_id, "node_id")
        if self.node_kind not in {"calculation", "analysis"}:
            raise ContractError("workflow semantic node kind is invalid")
        if self.node_kind == "calculation":
            for name, value in (
                ("program", self.program),
                ("jobtype", self.jobtype),
                ("project_id", self.project_id),
            ):
                require_identifier(value, name)
            if self.operation:
                raise ContractError("calculation nodes do not carry analysis operations")
        else:
            require_identifier(self.operation, "operation")
            if self.program or self.jobtype or self.project_id:
                raise ContractError("analysis nodes do not carry program project fields")
        object.__setattr__(
            self, "input_geometry_sha256s", tuple(self.input_geometry_sha256s)
        )
        object.__setattr__(self, "dependencies", tuple(self.dependencies))
        object.__setattr__(self, "output_semantics", tuple(self.output_semantics))
        object.__setattr__(
            self, "semantic_parameters", tuple(self.semantic_parameters)
        )
        if self.input_geometry_sha256s != tuple(
            sorted(set(self.input_geometry_sha256s))
        ):
            raise ContractError("workflow geometry digests must be sorted and unique")
        for digest in self.input_geometry_sha256s:
            require_sha256(digest, "input_geometry_sha256")
        dependency_keys = tuple(
            (item.source_node_id, item.source_output, item.target_input)
            for item in self.dependencies
        )
        if dependency_keys != tuple(sorted(set(dependency_keys))):
            raise ContractError("workflow dependencies must be sorted and unique")
        if any(item.source_node_id == self.node_id for item in self.dependencies):
            raise ContractError("workflow semantic node cannot depend on itself")
        if self.output_semantics != tuple(sorted(set(self.output_semantics))):
            raise ContractError("workflow output semantics must be sorted and unique")
        for output in self.output_semantics:
            require_identifier(output, "output_semantic")
        parameter_names = tuple(item[0] for item in self.semantic_parameters)
        if parameter_names != tuple(sorted(set(parameter_names))):
            raise ContractError(
                "workflow semantic parameters must be sorted and unique"
            )
        for name, value in self.semantic_parameters:
            require_identifier(name, "semantic_parameter")
            canonical_data(value)


@dataclass(frozen=True)
class PaperWorkflowAnswerKeyV1:
    """Private expert answer for non-executed intermediate/large systems."""

    schema_version: str
    answer_key_id: str
    case_id: str
    case_sha256: str
    answer_source: str
    projects: tuple[WorkflowProjectAnswerV1, ...]
    nodes: tuple[WorkflowSemanticNodeV1, ...]
    answer_key_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.paper-workflow-answer-key.v1":
            raise ContractError("unsupported paper workflow answer-key schema")
        require_identifier(self.answer_key_id, "answer_key_id")
        require_identifier(self.case_id, "case_id")
        require_sha256(self.case_sha256, "case_sha256")
        if not self.answer_source.strip():
            raise ContractError("workflow answer key requires a source description")
        object.__setattr__(self, "projects", tuple(self.projects))
        object.__setattr__(self, "nodes", tuple(self.nodes))
        project_ids = tuple(item.project_id for item in self.projects)
        node_ids = tuple(item.node_id for item in self.nodes)
        if not project_ids or project_ids != tuple(sorted(set(project_ids))):
            raise ContractError("workflow answer projects must be sorted and unique")
        if not node_ids or node_ids != tuple(sorted(set(node_ids))):
            raise ContractError("workflow answer nodes must be sorted and unique")
        _canonical_workflow_signatures(self.projects, self.nodes)
        if self.answer_key_sha256 != canonical_sha256(self._body()):
            raise ContractError("paper workflow answer-key digest mismatch")

    def _body(self) -> dict[str, Any]:
        return {
            key: value
            for key, value in self.__dict__.items()
            if key != "answer_key_sha256"
        }


def build_paper_workflow_answer_key(
    *,
    answer_key_id: str,
    case: PaperBlackBoxCaseV1,
    answer_source: str,
    projects: Sequence[WorkflowProjectAnswerV1],
    nodes: Sequence[WorkflowSemanticNodeV1],
) -> PaperWorkflowAnswerKeyV1:
    ordered_projects = tuple(sorted(projects, key=lambda item: item.project_id))
    ordered_nodes = tuple(sorted(nodes, key=lambda item: item.node_id))
    body = {
        "schema_version": "chemsmart.paper-workflow-answer-key.v1",
        "answer_key_id": require_identifier(answer_key_id, "answer_key_id"),
        "case_id": case.case_id,
        "case_sha256": case.case_sha256,
        "answer_source": str(answer_source).strip(),
        "projects": ordered_projects,
        "nodes": ordered_nodes,
    }
    return PaperWorkflowAnswerKeyV1(
        **body, answer_key_sha256=canonical_sha256(body)
    )


def _canonical_workflow_signatures(
    projects: Sequence[WorkflowProjectAnswerV1],
    nodes: Sequence[WorkflowSemanticNodeV1],
) -> tuple[tuple[str, ...], tuple[str, ...]]:
    """Canonicalize YAML and DAG semantics independently of model node names."""

    project_by_id = {item.project_id: item.document for item in projects}
    if len(project_by_id) != len(tuple(projects)):
        raise ContractError("workflow project IDs must be unique")
    node_by_id = {item.node_id: item for item in nodes}
    if len(node_by_id) != len(tuple(nodes)):
        raise ContractError("workflow node IDs must be unique")
    def normalized_settings(
        settings_rows: Sequence[tuple[str, Any]],
    ) -> dict[str, Any]:
        settings: dict[str, Any] = {}
        for name, value in settings_rows:
            if isinstance(value, str):
                normalized: Any = value.strip()
                if name not in {
                    "custom_solvent",
                    "additional_solvent_options",
                    "input_string",
                }:
                    normalized = normalized.lower()
                if name == "additional_route_parameters":
                    normalized = " ".join(sorted(normalized.split()))
            elif name == "heavy_elements" and isinstance(value, (list, tuple)):
                normalized = tuple(
                    sorted(str(item).capitalize() for item in value)
                )
            else:
                normalized = canonical_data(value)
            settings[name] = normalized
        return settings

    project_jobtypes: dict[str, set[str]] = {
        item.project_id: set() for item in projects
    }
    for node in nodes:
        if node.node_kind == "calculation":
            project_jobtypes.setdefault(node.project_id, set()).add(
                node.jobtype
            )

    def effective_project_sha256(
        document: ProjectDocumentV1, jobtype: str
    ) -> str:
        return canonical_sha256(
            {
                "program": document.program,
                "jobtype": jobtype,
                "settings": normalized_settings(
                    project_effective_section_settings(
                        document, jobtype=jobtype
                    )
                ),
            }
        )

    project_semantics = {
        item.project_id: tuple(
            (jobtype, effective_project_sha256(item.document, jobtype))
            for jobtype in sorted(project_jobtypes.get(item.project_id) or ())
        )
        for item in projects
    }
    project_signatures = tuple(
        sorted(
            canonical_sha256(
                {
                    "program": item.document.program,
                    "effective_job_settings": project_semantics[item.project_id],
                }
            )
            for item in projects
        )
    )
    visiting: set[str] = set()
    memo: dict[str, str] = {}

    def node_signature(node_id: str) -> str:
        if node_id in memo:
            return memo[node_id]
        if node_id in visiting:
            raise ContractError("workflow answer DAG contains a cycle")
        try:
            node = node_by_id[node_id]
        except KeyError as exc:
            raise ContractError("workflow answer references an unknown producer") from exc
        visiting.add(node_id)
        if node.node_kind == "calculation":
            try:
                project = project_by_id[node.project_id]
            except KeyError as exc:
                raise ContractError("workflow calculation references an unknown project") from exc
            program_fields: Mapping[str, Any] = {
                "program": node.program,
                "jobtype": node.jobtype,
                "project_program": project.program,
                "project_effective_sha256": effective_project_sha256(
                    project, node.jobtype
                ),
            }
        else:
            program_fields = {"operation": node.operation}
        dependency_signatures = tuple(
            sorted(
                (
                    node_signature(edge.source_node_id),
                    edge.source_output,
                    edge.target_input,
                )
                for edge in node.dependencies
            )
        )

        def bind_lineage(value: Any) -> Any:
            if isinstance(value, Mapping):
                if set(value) == {"source_node_id", "source_output"}:
                    source_node_id = str(value["source_node_id"])
                    return {
                        "source_node_signature": node_signature(source_node_id),
                        "source_output": value["source_output"],
                    }
                return {
                    str(key): bind_lineage(item)
                    for key, item in sorted(value.items(), key=lambda pair: str(pair[0]))
                }
            if isinstance(value, (list, tuple)):
                return tuple(bind_lineage(item) for item in value)
            return canonical_data(value)

        signature = canonical_sha256(
            {
                "node_kind": node.node_kind,
                **program_fields,
                "input_geometry_sha256s": node.input_geometry_sha256s,
                "dependencies": dependency_signatures,
                "output_semantics": node.output_semantics,
                "semantic_parameters": tuple(
                    (name, bind_lineage(value))
                    for name, value in node.semantic_parameters
                ),
            }
        )
        visiting.remove(node_id)
        memo[node_id] = signature
        return signature

    node_signatures = tuple(sorted(node_signature(node.node_id) for node in nodes))
    return project_signatures, node_signatures


@dataclass(frozen=True)
class PaperBenchmarkEligibilityV1:
    """Proof that an eligible answer key existed before task dispatch."""

    schema_version: str
    case_id: str
    case_sha256: str
    answer_key_sha256: str
    evaluation_mode: str
    required_quantity_ids: tuple[str, ...]
    required_workflow_node_count: int
    status: str
    receipt_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.paper-benchmark-eligibility.v1":
            raise ContractError("unsupported paper benchmark eligibility schema")
        if self.status != "eligible_answer_key_registered":
            raise ContractError("paper benchmark is not eligible for dispatch")
        require_identifier(self.case_id, "case_id")
        require_sha256(self.case_sha256, "case_sha256")
        require_sha256(self.answer_key_sha256, "answer_key_sha256")
        if self.evaluation_mode not in {
            "small_molecule_numerical",
            "canonical_yaml_cli_dag",
        }:
            raise ContractError("unsupported paper benchmark evaluation mode")
        if self.required_quantity_ids != tuple(sorted(set(self.required_quantity_ids))):
            raise ContractError("benchmark quantity IDs must be sorted and unique")
        if self.evaluation_mode == "small_molecule_numerical":
            if not self.required_quantity_ids or self.required_workflow_node_count:
                raise ContractError("small-molecule benchmarks require numerical answers")
        elif self.required_quantity_ids or self.required_workflow_node_count < 1:
            raise ContractError(
                "intermediate benchmarks require a canonical YAML-CLI DAG answer"
            )
        if self.receipt_sha256 != canonical_sha256(self._body()):
            raise ContractError("paper benchmark eligibility digest mismatch")

    def _body(self) -> dict[str, Any]:
        return {
            key: value
            for key, value in self.__dict__.items()
            if key != "receipt_sha256"
        }


def prepare_benchmark_dispatch(
    *,
    case: PaperBlackBoxCaseV1,
    answer_key: "PaperAnswerKeyV1 | PaperWorkflowAnswerKeyV1",
) -> tuple[str, PaperBenchmarkEligibilityV1]:
    """Return a public task only after the appropriate private key is registered."""

    if answer_key.case_id != case.case_id or answer_key.case_sha256 != case.case_sha256:
        raise ContractError("paper answer key targets another black-box case")
    if isinstance(answer_key, PaperAnswerKeyV1):
        if case.system_scale != "small_molecule":
            raise ContractError(
                "numerical execution benchmarks are restricted to small molecules"
            )
        evaluation_mode = "small_molecule_numerical"
        required_quantity_ids = tuple(
            quantity.quantity_id for quantity in answer_key.quantities
        )
        required_workflow_node_count = 0
    elif isinstance(answer_key, PaperWorkflowAnswerKeyV1):
        if case.system_scale != "intermediate_or_larger":
            raise ContractError(
                "canonical YAML-CLI DAG grading targets intermediate or larger systems"
            )
        evaluation_mode = "canonical_yaml_cli_dag"
        required_quantity_ids = ()
        required_workflow_node_count = len(answer_key.nodes)
    else:
        raise ContractError("benchmark dispatch requires a typed private answer key")
    body = {
        "schema_version": "chemsmart.paper-benchmark-eligibility.v1",
        "case_id": case.case_id,
        "case_sha256": case.case_sha256,
        "answer_key_sha256": answer_key.answer_key_sha256,
        "evaluation_mode": evaluation_mode,
        "required_quantity_ids": required_quantity_ids,
        "required_workflow_node_count": required_workflow_node_count,
        "status": "eligible_answer_key_registered",
    }
    receipt = PaperBenchmarkEligibilityV1(
        **body, receipt_sha256=canonical_sha256(body)
    )
    return case.public_task(), receipt


@dataclass(frozen=True)
class PaperContextRecordV1:
    """One public predecessor record eligible for P0/P1/P2 selection."""

    record_id: str
    record_class: str
    node_id: str
    public_record: Mapping[str, Any]

    def __post_init__(self) -> None:
        require_identifier(self.record_id, "record_id")
        require_identifier(self.record_class, "record_class")
        if self.node_id:
            require_identifier(self.node_id, "node_id")
        if not isinstance(self.public_record, Mapping):
            raise ContractError("paper context public record must be a mapping")


@dataclass(frozen=True)
class PaperContextArmV1:
    arm_id: str
    context: TaskDependencyContextV1
    selection_receipt: ContextSelectionReceiptV1
    public_records: Mapping[str, Mapping[str, Any]]


def build_paper_context_arms(
    *,
    case_id: str,
    plan: ScientificWorkflowPlanV2,
    target_node_id: str,
    records: Sequence[PaperContextRecordV1],
    max_public_bytes: int = 1_000_000,
) -> tuple[PaperContextArmV1, ...]:
    """Create a complete P0/P1/P2 block over one fixed plan and record set."""

    normalized_case_id = require_identifier(case_id, "case_id")
    refs = tuple(
        build_predecessor_evidence_ref(
            record_id=record.record_id,
            record_class=record.record_class,
            node_id=record.node_id,
            record_sha256=canonical_sha256(record.public_record),
            size_bytes=len(canonical_json(record.public_record).encode("utf-8")),
        )
        for record in records
    )
    all_payloads = {
        record.record_id: dict(record.public_record) for record in records
    }
    arms: list[PaperContextArmV1] = []
    for arm_id in _CONTEXT_ARMS:
        policy = build_task_dependency_context_policy(
            policy_id=f"{normalized_case_id}-{arm_id}",
            arm_id=arm_id,
            record_classes=tuple(
                sorted({record.record_class for record in records})
            ),
            max_public_bytes=max_public_bytes,
        )
        context, receipt = select_task_dependency_context(
            selection_id=f"{normalized_case_id}-{arm_id}-selection",
            plan=plan,
            target_node_id=target_node_id,
            policy=policy,
            evidence_refs=refs,
        )
        if context is None:
            raise ContractError(receipt.reason)
        selected_ids = {item.record_id for item in context.evidence_refs}
        arms.append(
            PaperContextArmV1(
                arm_id=arm_id,
                context=context,
                selection_receipt=receipt,
                public_records={
                    record_id: payload
                    for record_id, payload in all_payloads.items()
                    if record_id in selected_ids
                },
            )
        )
    return tuple(arms)


@dataclass(frozen=True)
class PaperWorkflowOracleV1:
    """Paper-independent structural expectations for one private case."""

    oracle_id: str
    required_tools: tuple[str, ...]
    required_transcript_terms: tuple[str, ...]
    required_final_terms: tuple[str, ...]
    unsupported_terms: tuple[str, ...]
    require_preview: bool = True

    def __post_init__(self) -> None:
        require_identifier(self.oracle_id, "oracle_id")
        for values in (
            self.required_tools,
            self.required_transcript_terms,
            self.required_final_terms,
            self.unsupported_terms,
        ):
            if values != tuple(sorted(set(values))):
                raise ContractError("paper oracle terms must be sorted and unique")


@dataclass(frozen=True)
class PaperWorkflowGradeV1:
    """Diagnostic grade for workflow structure, never benchmark acceptance."""

    schema_version: str
    oracle_id: str
    evaluation_role: str
    terminal_state: str
    observed_tools: tuple[str, ...]
    missing_tools: tuple[str, ...]
    missing_transcript_terms: tuple[str, ...]
    missing_final_terms: tuple[str, ...]
    unsupported_terms_without_blocking: tuple[str, ...]
    preview_observed: bool
    strict_pass: bool
    grade_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.paper-workflow-grade.v1":
            raise ContractError("unsupported paper workflow grade")
        if self.evaluation_role != "diagnostic_structure_only":
            raise ContractError("paper workflow grade cannot be a primary grade")
        if self.grade_sha256 != canonical_sha256(self._body()):
            raise ContractError("paper workflow grade digest mismatch")

    def _body(self) -> dict[str, Any]:
        return {
            key: value
            for key, value in self.__dict__.items()
            if key != "grade_sha256"
        }


def grade_paper_workflow_result(
    *, result: "LiveAgentSessionResultV1", oracle: PaperWorkflowOracleV1
) -> PaperWorkflowGradeV1:
    """Grade a live run without using provider-private reasoning."""

    observed_tools: set[str] = set()
    transcript_text = canonical_json(result.public_transcript).lower()
    for message in result.public_transcript:
        for tool_call in message.get("tool_calls") or ():
            function = tool_call.get("function") or {}
            name = str(function.get("name") or "").strip()
            if name:
                observed_tools.add(name)
    final_text = result.final_text.lower()
    # Public scientific identifiers appear legitimately as ``def2SVP``,
    # ``def2-SVP``, ``excitation energies``, or ``excitation_energies`` across
    # project loaders, parser selectors, and prose.  This oracle grades the
    # semantic token, not punctuation-sensitive serialization.  Raw transcript
    # and strict tool-call grades remain available separately.
    semantic_transcript = re.sub(r"[^a-z0-9]+", "", transcript_text)
    semantic_final = re.sub(r"[^a-z0-9]+", "", final_text)

    def _present(term: str, text: str) -> bool:
        return re.sub(r"[^a-z0-9]+", "", term.lower()) in text

    missing_tools = tuple(sorted(set(oracle.required_tools) - observed_tools))
    missing_transcript = tuple(
        term
        for term in oracle.required_transcript_terms
        if not _present(term, semantic_transcript)
    )
    missing_final = tuple(
        term
        for term in oracle.required_final_terms
        if not _present(term, semantic_final)
    )
    blocking_language = any(
        term in transcript_text
        for term in ("blocked_unsupported", "unsupported", "not supported")
    )
    unsupported_without_blocking = tuple(
        term
        for term in oracle.unsupported_terms
        if _present(term, semantic_transcript) and not blocking_language
    )
    preview_observed = bool(
        {"prepare_program_node", "preview_command"}.intersection(observed_tools)
        or "safe_preview" in transcript_text
    )
    body = {
        "schema_version": "chemsmart.paper-workflow-grade.v1",
        "oracle_id": oracle.oracle_id,
        "evaluation_role": "diagnostic_structure_only",
        "terminal_state": result.terminal_state,
        "observed_tools": tuple(sorted(observed_tools)),
        "missing_tools": missing_tools,
        "missing_transcript_terms": missing_transcript,
        "missing_final_terms": missing_final,
        "unsupported_terms_without_blocking": unsupported_without_blocking,
        "preview_observed": preview_observed,
        "strict_pass": bool(
            result.terminal_state in {"complete", "planned"}
            and not missing_tools
            and not missing_transcript
            and not missing_final
            and not unsupported_without_blocking
            and (preview_observed or not oracle.require_preview)
        ),
    }
    return PaperWorkflowGradeV1(
        **body, grade_sha256=canonical_sha256(body)
    )


@dataclass(frozen=True)
class NumericalQuantityGradeV1:
    """Comparison of one host-rendered numerical claim with its private key."""

    quantity_id: str
    status: str
    expected_value: Any
    observed_value: Any
    canonical_unit: str
    maximum_absolute_error: float | None
    maximum_relative_error: float | None
    within_tolerance: bool

    def __post_init__(self) -> None:
        require_identifier(self.quantity_id, "quantity_id")
        if self.status not in {
            "within_tolerance",
            "outside_tolerance",
            "missing_claim",
            "incompatible_unit",
            "shape_mismatch",
        }:
            raise ContractError("unsupported numerical quantity grade status")
        for value in (self.maximum_absolute_error, self.maximum_relative_error):
            if value is not None and (not math.isfinite(value) or value < 0.0):
                raise ContractError("numerical errors must be finite and non-negative")
        if self.within_tolerance != (self.status == "within_tolerance"):
            raise ContractError("numerical grade status and disposition disagree")


@dataclass(frozen=True)
class PaperNumericalGradeV1:
    """Primary paper benchmark grade over answer-key-bound numerical claims."""

    schema_version: str
    case_id: str
    case_sha256: str
    answer_key_sha256: str
    claim_record_sha256: str
    quantity_grades: tuple[NumericalQuantityGradeV1, ...]
    primary_pass: bool
    grade_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.paper-numerical-grade.v1":
            raise ContractError("unsupported paper numerical grade schema")
        require_identifier(self.case_id, "case_id")
        require_sha256(self.case_sha256, "case_sha256")
        require_sha256(self.answer_key_sha256, "answer_key_sha256")
        if self.claim_record_sha256:
            require_sha256(self.claim_record_sha256, "claim_record_sha256")
        object.__setattr__(self, "quantity_grades", tuple(self.quantity_grades))
        quantity_ids = tuple(item.quantity_id for item in self.quantity_grades)
        if quantity_ids != tuple(sorted(set(quantity_ids))) or not quantity_ids:
            raise ContractError("numerical quantity grades must be sorted and unique")
        if self.primary_pass != all(
            item.within_tolerance for item in self.quantity_grades
        ):
            raise ContractError("primary pass must require every numerical target")
        if self.grade_sha256 != canonical_sha256(self._body()):
            raise ContractError("paper numerical grade digest mismatch")

    def _body(self) -> dict[str, Any]:
        return {
            key: value
            for key, value in self.__dict__.items()
            if key != "grade_sha256"
        }


def grade_paper_numerical_claims(
    *,
    case: PaperBlackBoxCaseV1,
    answer_key: PaperAnswerKeyV1,
    claim_record: "AnalysisClaimRecordV1 | None",
) -> PaperNumericalGradeV1:
    """Apply the primary answer criterion to host-rendered numerical claims.

    Structural planning, command previews, or terminal text cannot satisfy this
    grade.  Every required quantity must be present in the typed claim record,
    dimensionally compatible, and inside its preregistered absolute/relative
    tolerance envelope.
    """

    if answer_key.case_id != case.case_id or answer_key.case_sha256 != case.case_sha256:
        raise ContractError("paper answer key targets another black-box case")
    claims_by_id = {
        claim.claim_id: claim
        for claim in claim_record.claims
    } if claim_record is not None else {}
    claims_by_quantity: dict[str, list[Any]] = {}
    if claim_record is not None:
        for claim in claim_record.claims:
            claims_by_quantity.setdefault(claim.quantity_id, []).append(claim)
    grades: list[NumericalQuantityGradeV1] = []
    for reference in answer_key.quantities:
        expected, canonical_unit, expected_dimension = normalize_numeric_value(
            reference.expected_value, reference.unit
        )
        claim = claims_by_id.get(reference.quantity_id)
        if claim is None:
            quantity_matches = claims_by_quantity.get(reference.quantity_id, [])
            claim = _collapse_equivalent_claim_presentations(quantity_matches)
        if claim is None:
            accepted_token_sets = tuple(
                _semantic_quantity_tokens(identifier)
                for identifier in reference.accepted_identifiers
            )
            semantic_matches = tuple(
                candidate
                for candidate in claims_by_id.values()
                if any(
                    accepted_tokens
                    and accepted_tokens.issubset(
                        _semantic_quantity_tokens(candidate.claim_id).union(
                            _semantic_quantity_tokens(candidate.quantity_id)
                        )
                    )
                    for accepted_tokens in accepted_token_sets
                )
            )
            claim = _collapse_equivalent_claim_presentations(semantic_matches)
        if claim is None:
            grades.append(
                NumericalQuantityGradeV1(
                    quantity_id=reference.quantity_id,
                    status="missing_claim",
                    expected_value=expected,
                    observed_value=None,
                    canonical_unit=canonical_unit,
                    maximum_absolute_error=None,
                    maximum_relative_error=None,
                    within_tolerance=False,
                )
            )
            continue
        try:
            observed, observed_unit, observed_dimension = normalize_numeric_value(
                claim.display_value, claim.display_unit
            )
        except (TypeError, ValueError):
            observed = None
            observed_unit = ""
            observed_dimension = ()
        if observed_dimension != expected_dimension:
            grades.append(
                NumericalQuantityGradeV1(
                    quantity_id=reference.quantity_id,
                    status="incompatible_unit",
                    expected_value=expected,
                    observed_value=observed,
                    canonical_unit=canonical_unit,
                    maximum_absolute_error=None,
                    maximum_relative_error=None,
                    within_tolerance=False,
                )
            )
            continue
        expected_array = np.asarray(expected, dtype=float)
        observed_array = np.asarray(observed, dtype=float)
        if expected_array.shape != observed_array.shape:
            grades.append(
                NumericalQuantityGradeV1(
                    quantity_id=reference.quantity_id,
                    status="shape_mismatch",
                    expected_value=expected,
                    observed_value=observed,
                    canonical_unit=observed_unit,
                    maximum_absolute_error=None,
                    maximum_relative_error=None,
                    within_tolerance=False,
                )
            )
            continue
        absolute_error = np.abs(observed_array - expected_array)
        absolute_tolerance, _, _ = normalize_numeric_value(
            reference.absolute_tolerance, reference.unit
        )
        allowed = float(absolute_tolerance) + (
            reference.relative_tolerance * np.abs(expected_array)
        )
        with np.errstate(divide="ignore", invalid="ignore"):
            relative_error = np.where(
                expected_array != 0.0,
                absolute_error / np.abs(expected_array),
                np.where(absolute_error == 0.0, 0.0, np.inf),
            )
        within = bool(np.all(absolute_error <= allowed))
        maximum_relative_error = float(np.max(relative_error))
        if not math.isfinite(maximum_relative_error):
            # Infinite relative error at an exact-zero reference is informative
            # but canonical grade records admit finite values only.  The
            # absolute error still determines the declared tolerance result.
            maximum_relative_error = float(np.finfo(float).max)
        grades.append(
            NumericalQuantityGradeV1(
                quantity_id=reference.quantity_id,
                status="within_tolerance" if within else "outside_tolerance",
                expected_value=expected,
                observed_value=observed,
                canonical_unit=canonical_unit,
                maximum_absolute_error=float(np.max(absolute_error)),
                maximum_relative_error=maximum_relative_error,
                within_tolerance=within,
            )
        )
    ordered = tuple(sorted(grades, key=lambda item: item.quantity_id))
    body = {
        "schema_version": "chemsmart.paper-numerical-grade.v1",
        "case_id": case.case_id,
        "case_sha256": case.case_sha256,
        "answer_key_sha256": answer_key.answer_key_sha256,
        "claim_record_sha256": (
            claim_record.receipt_sha256 if claim_record is not None else ""
        ),
        "quantity_grades": ordered,
        "primary_pass": all(item.within_tolerance for item in ordered),
    }
    return PaperNumericalGradeV1(
        **body, grade_sha256=canonical_sha256(body)
    )


def _semantic_quantity_tokens(identifier: str) -> frozenset[str]:
    """Return conservative tokens for matching answer and claim semantics.

    Model-selected local IDs are presentation details.  Matching strips only
    generic reporting qualifiers and never introduces chemical aliases, so a
    Gibbs energy cannot become an electronic energy and an entropy cannot
    become a heat capacity merely because the dimensions coincide.
    """

    return frozenset(
        _SCIENTIFIC_IDENTIFIER_TOKEN_ALIASES.get(token, token)
        for token in re.split(r"[-_]+", str(identifier).strip().lower())
        if token and token not in _GENERIC_CLAIM_TOKENS
    )


def _collapse_equivalent_claim_presentations(
    candidates: Sequence[Any],
) -> Any | None:
    """Return one claim when candidates are displays of one typed quantity.

    A single expression result is often reported in hartree, kJ/mol, and
    kcal/mol.  Those are not three competing scientific answers when every
    claim points to the same host quantity value and source receipt.  Distinct
    source values remain ambiguous and therefore fail closed.
    """

    values = tuple(candidates)
    if not values:
        return None
    identities = {
        (
            claim.source_kind,
            claim.source_receipt_sha256,
            claim.quantity_id,
            claim.quantity_value_sha256,
        )
        for claim in values
    }
    if len(identities) != 1:
        return None
    return min(values, key=lambda claim: claim.claim_id)


def analysis_claim_record_from_live_result(
    result: "LiveAgentSessionResultV1",
) -> "AnalysisClaimRecordV1 | None":
    """Recover the last host-rendered claim record from a live episode.

    A paper benchmark must grade the typed value copied by the host, not a
    number repeated in model prose.  The public transcript already contains
    the sanitized ``record_analysis_claims`` tool result, so this adapter does
    not reopen private provider state or depend on Runtime V2 filesystem paths.
    When a bounded repair records a replacement claim set, the latest
    successful host record is authoritative.
    """

    from chemsmart.agent.analysis_claims import (
        analysis_claim_record_from_record,
    )

    observed: AnalysisClaimRecordV1 | None = None
    for message in result.public_transcript:
        if message.get("role") != "tool":
            continue
        content = message.get("content")
        if not isinstance(content, str):
            continue
        try:
            envelope = json.loads(content)
        except json.JSONDecodeError:
            continue
        if not isinstance(envelope, Mapping) or (
            envelope.get("schema_version") != "chemsmart.tool-result.v1"
            or envelope.get("status") != "ok"
            or envelope.get("tool") != "record_analysis_claims"
        ):
            continue
        record = envelope.get("result")
        if not isinstance(record, Mapping):
            raise ContractError(
                "successful analysis-claim tool result is not a record"
            )
        values = dict(record)
        receipt_sha256 = str(values.pop("receipt_sha256", ""))
        observed = analysis_claim_record_from_record(
            values, receipt_sha256=receipt_sha256
        )
        if observed.task_spec_sha256 != result.task_spec_sha256:
            raise ContractError(
                "analysis claim record belongs to another live task"
            )
    return observed


@dataclass(frozen=True)
class PaperWorkflowObservationV1:
    """Host-normalized YAML and DAG semantics observed from one model run."""

    schema_version: str
    projects: tuple[WorkflowProjectAnswerV1, ...]
    nodes: tuple[WorkflowSemanticNodeV1, ...]

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.paper-workflow-observation.v1":
            raise ContractError("unsupported paper workflow observation schema")
        object.__setattr__(self, "projects", tuple(self.projects))
        object.__setattr__(self, "nodes", tuple(self.nodes))
        _canonical_workflow_signatures(self.projects, self.nodes)


def _successful_live_tool_results(
    result: "LiveAgentSessionResultV1",
) -> tuple[tuple[str, Mapping[str, Any]], ...]:
    observed: list[tuple[str, Mapping[str, Any]]] = []
    for message in result.public_transcript:
        if message.get("role") != "tool":
            continue
        content = message.get("content")
        if not isinstance(content, str):
            continue
        try:
            envelope = json.loads(content)
        except json.JSONDecodeError:
            continue
        if not isinstance(envelope, Mapping) or (
            envelope.get("schema_version") != "chemsmart.tool-result.v1"
            or envelope.get("status") != "ok"
        ):
            continue
        tool = str(envelope.get("tool") or "")
        payload = envelope.get("result")
        if tool and isinstance(payload, Mapping):
            observed.append((tool, payload))
    return tuple(observed)


def _semantic_name(value: Any) -> str:
    text = re.sub(r"[^a-z0-9]+", "_", str(value).strip().lower()).strip("_")
    aliases = {
        "energy": "total_energy",
        "energies": "total_energy",
        "final_energy": "total_energy",
        "frequency": "vibrational_frequencies",
        "frequencies": "vibrational_frequencies",
        "geometry_xyz": "geometry",
        "optimized_coordinates": "optimized_geometry",
        "orca_output": "program_result",
        "orca_result": "program_result",
        "result": "program_result",
    }
    normalized = aliases.get(text, text)
    return require_identifier(normalized, "workflow semantic")


def _artifact_output_semantic(
    *, artifact_class: str, producer_jobtype: str
) -> str:
    normalized = _semantic_name(artifact_class)
    if normalized == "geometry" and producer_jobtype == "opt":
        return "optimized_geometry"
    if normalized.endswith("_result") or normalized.endswith("_output"):
        return "program_result"
    return normalized


def _effective_project_validation(
    validations: Sequence[Mapping[str, Any]],
    *,
    artifact_id: str,
    program: str,
    jobtype: str,
) -> Mapping[str, Any]:
    candidates = tuple(
        payload
        for payload in validations
        if payload.get("project_artifact_id") == artifact_id
        and payload.get("program") == program
        and payload.get("jobtype") == jobtype
        and payload.get("status") == "valid"
    )
    if not candidates:
        raise ContractError(
            f"workflow project {artifact_id!r} lacks a valid {program}/{jobtype} "
            "loader observation"
        )
    selected = candidates[-1]
    settings = dict(selected.get("settings") or ())
    method = settings.get("functional") or settings.get("ab_initio") or settings.get(
        "semiempirical"
    )
    basis = (
        settings.get("basis")
        or settings.get("light_elements_basis")
        or settings.get("gfn_version")
    )
    if program in {"gaussian", "orca", "pyscf"} and (not method or not basis):
        raise ContractError(
            f"workflow project {artifact_id!r} validated without effective "
            f"{program}/{jobtype} method and basis settings"
        )
    return selected


_COMMUTATIVE_EXPRESSION_OPERATIONS = frozenset(
    {"add", "multiply", "sum", "weighted_sum"}
)


def _canonical_expression_outputs(
    node: Mapping[str, Any],
    *,
    input_semantics: Mapping[str, Any],
) -> tuple[Mapping[str, Any], ...]:
    raw_nodes = tuple(node.get("expression_nodes") or ())
    by_id = {
        str(item.get("node_id")): item
        for item in raw_nodes
        if isinstance(item, Mapping) and item.get("node_id")
    }
    memo: dict[str, Mapping[str, Any]] = {}
    visiting: set[str] = set()

    def visit(node_id: str) -> Mapping[str, Any]:
        if node_id in memo:
            return memo[node_id]
        if node_id in visiting or node_id not in by_id:
            raise ContractError("analysis expression is cyclic or incomplete")
        visiting.add(node_id)
        raw = by_id[node_id]
        operation = _semantic_name(raw.get("operation"))
        if operation == "ref":
            reference = str(raw.get("reference") or "")
            try:
                record: dict[str, Any] = {
                    "operation": "ref",
                    "input_semantic": input_semantics[reference],
                }
            except KeyError as exc:
                raise ContractError(
                    "analysis expression references an unknown typed input"
                ) from exc
        else:
            arguments = [visit(str(item)) for item in raw.get("input_ids") or ()]
            if operation in _COMMUTATIVE_EXPRESSION_OPERATIONS:
                arguments.sort(key=canonical_json)
            record = {"operation": operation, "arguments": tuple(arguments)}
            for field in (
                "literal_value",
                "literal_unit",
                "target_unit",
                "scale_factor",
                "cardinal_numbers",
                "extrapolation_exponent",
            ):
                if raw.get(field) is not None:
                    record[field] = canonical_data(raw[field])
        visiting.remove(node_id)
        memo[node_id] = record
        return record

    outputs = tuple(
        visit(str(item)) for item in node.get("expression_output_node_ids") or ()
    )
    return tuple(sorted(outputs, key=canonical_json))


def paper_workflow_observation_from_live_result(
    result: "LiveAgentSessionResultV1",
) -> PaperWorkflowObservationV1:
    """Normalize the last typed paper plan into private-answer semantics.

    Only successful public tool records are observed.  Model prose, provider
    reasoning, expected answers, and benchmark rubrics are deliberately absent
    from this adapter.  A project counts only when the same promoted artifact
    passed the checked-out loader for the job form used by the final plan.
    """

    from chemsmart.agent.projects import project_document

    tools = _successful_live_tool_results(result)
    renders: dict[str, Mapping[str, Any]] = {}
    promotions: dict[str, str] = {}
    validations: list[Mapping[str, Any]] = []
    identities: dict[str, tuple[str, int, int]] = {}
    final_plan: Mapping[str, Any] | None = None
    for tool, payload in tools:
        if tool == "render_project_yaml":
            renders[str(payload.get("receipt_sha256"))] = payload
        elif tool == "promote_project_yaml":
            artifact = payload.get("artifact") or {}
            promotion = payload.get("promotion") or {}
            if isinstance(artifact, Mapping) and isinstance(promotion, Mapping):
                promotions[str(artifact.get("artifact_id"))] = str(
                    promotion.get("render_receipt_sha256")
                )
        elif tool == "validate_project_yaml":
            validations.append(payload)
        elif tool == "bind_scientific_identity":
            artifact_id = str(payload.get("geometry_artifact_id") or "")
            if artifact_id:
                identities[artifact_id] = (
                    require_sha256(
                        str(payload.get("geometry_artifact_sha256")),
                        "geometry_artifact_sha256",
                    ),
                    int(payload.get("charge")),
                    int(payload.get("multiplicity")),
                )
        elif tool in {
            "plan_scientific_workflow",
            "rebind_scientific_workflow_projects",
        }:
            final_plan = payload
    if final_plan is None:
        raise ContractError("live episode contains no successful scientific workflow plan")

    calculation_plan = final_plan.get("calculation_plan") or {}
    draft = calculation_plan.get("workflow_draft") or {}
    raw_calculations = tuple(draft.get("nodes") or ())
    toolchain = final_plan.get("scientific_toolchain_plan") or {}
    raw_analyses = tuple(toolchain.get("analysis_nodes") or ())
    if not raw_calculations or not raw_analyses:
        raise ContractError("scientific workflow plan is incomplete")
    calculation_by_id = {
        str(item.get("node_id")): item
        for item in raw_calculations
        if isinstance(item, Mapping)
    }
    if len(calculation_by_id) != len(raw_calculations):
        raise ContractError("scientific workflow has duplicate calculation IDs")

    artifact_sha256s = {
        str(item.get("artifact_id")): str(item.get("sha256"))
        for item in result.artifact_records
        if item.get("artifact_id") and item.get("sha256")
    }
    observable_rows = dict(toolchain.get("calculation_observables") or ())
    output_maps: dict[str, dict[str, str]] = {}
    root_geometry_memo: dict[str, tuple[str, ...]] = {}

    def calculation_roots(node_id: str, visiting: set[str] | None = None) -> tuple[str, ...]:
        if node_id in root_geometry_memo:
            return root_geometry_memo[node_id]
        active = set() if visiting is None else set(visiting)
        if node_id in active:
            raise ContractError("scientific calculation DAG contains a cycle")
        active.add(node_id)
        node = calculation_by_id[node_id]
        roots: set[str] = set()
        for item in node.get("inputs") or ():
            artifact_id = str(item.get("artifact_id") or "")
            producer = str(item.get("producer_node_id") or "")
            if artifact_id:
                try:
                    roots.add(require_sha256(artifact_sha256s[artifact_id], "geometry"))
                except KeyError as exc:
                    raise ContractError(
                        "workflow references an unobserved initial geometry"
                    ) from exc
            elif producer:
                roots.update(calculation_roots(producer, active))
        value = tuple(sorted(roots))
        root_geometry_memo[node_id] = value
        return value

    observed_projects: dict[str, WorkflowProjectAnswerV1] = {}
    calculation_nodes: list[WorkflowSemanticNodeV1] = []
    for node_id, raw in calculation_by_id.items():
        program = _semantic_name(raw.get("program"))
        jobtype = _semantic_name(raw.get("jobtype"))
        project_id = str(raw.get("project_role") or "")
        try:
            render = renders[promotions[project_id]]
        except KeyError as exc:
            raise ContractError(
                f"workflow project role {project_id!r} is not a promoted project"
            ) from exc
        if str(render.get("program")) != program:
            raise ContractError("workflow project program differs from calculation")
        parsed = yaml.safe_load(str(render.get("rendered_yaml") or ""))
        if not isinstance(parsed, Mapping) or not parsed:
            raise ContractError("workflow project render is not a YAML mapping")
        document = project_document(
            program=program,
            sections={str(key): dict(value) for key, value in parsed.items()},
        )
        _effective_project_validation(
            validations,
            artifact_id=project_id,
            program=program,
            jobtype=jobtype,
        )
        observed_projects[project_id] = WorkflowProjectAnswerV1(
            project_id=project_id,
            document=document,
        )

        output_map: dict[str, str] = {}
        output_semantics: set[str] = set()
        for output in raw.get("expected_outputs") or ():
            semantic = _artifact_output_semantic(
                artifact_class=str(output.get("artifact_class") or ""),
                producer_jobtype=jobtype,
            )
            output_map[str(output.get("output_id"))] = semantic
            output_semantics.add(semantic)
        for observable in observable_rows.get(node_id, ()):
            output_semantics.add(_semantic_name(observable))
        output_maps[node_id] = output_map

        dependencies: set[tuple[str, str, str]] = set()
        data_producers: set[str] = set()
        for item in raw.get("inputs") or ():
            producer = str(item.get("producer_node_id") or "")
            if not producer:
                continue
            data_producers.add(producer)
            source_output = output_maps.get(producer, {}).get(
                str(item.get("producer_output_id") or ""),
                _artifact_output_semantic(
                    artifact_class=str(item.get("artifact_class") or ""),
                    producer_jobtype=_semantic_name(
                        calculation_by_id[producer].get("jobtype")
                    ),
                ),
            )
            target_input = _semantic_name(item.get("artifact_class") or "artifact")
            dependencies.add((producer, source_output, target_input))
        for producer in raw.get("dependencies") or ():
            producer_id = str(producer)
            if producer_id not in data_producers:
                dependencies.add((producer_id, "completion", "control"))

        roots = calculation_roots(node_id)
        states = {
            identities[artifact_id]
            for artifact_id, digest in artifact_sha256s.items()
            if digest in roots and artifact_id in identities
        }
        semantic_parameters: list[tuple[str, Any]] = []
        if states:
            semantic_parameters.append(
                (
                    "electronic_states",
                    tuple(sorted((charge, multiplicity) for _, charge, multiplicity in states)),
                )
            )
        calculation_nodes.append(
            WorkflowSemanticNodeV1(
                node_id=node_id,
                node_kind="calculation",
                program=program,
                jobtype=jobtype,
                project_id=project_id,
                input_geometry_sha256s=roots,
                dependencies=tuple(
                    WorkflowDependencyAnswerV1(*item)
                    for item in sorted(dependencies)
                ),
                output_semantics=tuple(sorted(output_semantics)),
                semantic_parameters=tuple(sorted(semantic_parameters)),
            )
        )

    analysis_output_maps: dict[str, dict[str, str]] = {}
    for raw in raw_analyses:
        analysis_output_maps[str(raw.get("node_id"))] = {
            str(item.get("output_id")): _semantic_name(item.get("quantity_kind"))
            for item in raw.get("outputs") or ()
        }

    required_outputs = set(toolchain.get("required_output_ids") or ())
    analysis_nodes: list[WorkflowSemanticNodeV1] = []
    for raw in raw_analyses:
        node_id = str(raw.get("node_id"))
        dependencies: set[tuple[str, str, str]] = set()
        data_producers: set[str] = set()
        input_semantics: dict[str, Any] = {}
        for item in raw.get("inputs") or ():
            producer = str(item.get("producer_node_id") or "")
            producer_output = str(item.get("producer_output_id") or "")
            data_producers.add(producer)
            source_semantic = (
                output_maps.get(producer, {}).get(producer_output)
                or analysis_output_maps.get(producer, {}).get(producer_output)
                or _semantic_name(producer_output)
            )
            input_id = str(item.get("input_id") or "")
            input_semantics[input_id] = {
                "source_node_id": producer,
                "source_output": source_semantic,
            }
            dependencies.add((producer, source_semantic, source_semantic))
        for producer in raw.get("dependencies") or ():
            producer_id = str(producer)
            if producer_id not in data_producers:
                dependencies.add((producer_id, "completion", "control"))

        outputs = tuple(raw.get("outputs") or ())
        parameters: list[tuple[str, Any]] = [
            ("support_state", _semantic_name(raw.get("support_state"))),
            (
                "output_units",
                tuple(
                    sorted(
                        (
                            _semantic_name(item.get("quantity_kind")),
                            str(item.get("unit") or "").strip().lower(),
                        )
                        for item in outputs
                    )
                ),
            ),
        ]
        selectors = tuple(
            sorted(_semantic_name(item.get("selector")) for item in raw.get("selectors") or ())
        )
        if selectors:
            parameters.append(("selectors", selectors))
        validation_rules = tuple(
            sorted(
                (
                    _semantic_name(item.get("predicate")),
                    tuple(
                        sorted(
                            (
                                input_semantics.get(str(input_id), {}) or {}
                            ).get("source_output", _semantic_name(input_id))
                            for input_id in item.get("input_ids") or ()
                        )
                    ),
                    canonical_data(item.get("threshold")),
                    canonical_data(item.get("expected_count")),
                    str(item.get("unit") or "").strip().lower(),
                )
                for item in raw.get("validation_rules") or ()
            ),
            key=canonical_json,
        )
        if validation_rules:
            parameters.append(("validation_rules", validation_rules))
        expression_outputs = _canonical_expression_outputs(
            raw, input_semantics=input_semantics
        )
        if expression_outputs:
            parameters.append(("expression_outputs", expression_outputs))
        required_semantics = tuple(
            sorted(
                _semantic_name(item.get("quantity_kind"))
                for item in outputs
                if item.get("output_id") in required_outputs
            )
        )
        if required_semantics:
            parameters.append(("required_outputs", required_semantics))
        if raw.get("temperature_k") is not None:
            parameters.append(("temperature_k", float(raw["temperature_k"])))
        if raw.get("pressure_atm") is not None:
            parameters.append(("pressure_atm", float(raw["pressure_atm"])))
        analysis_nodes.append(
            WorkflowSemanticNodeV1(
                node_id=node_id,
                node_kind="analysis",
                operation=_semantic_name(raw.get("analysis_kind")),
                dependencies=tuple(
                    WorkflowDependencyAnswerV1(*item)
                    for item in sorted(dependencies)
                ),
                output_semantics=tuple(
                    sorted(
                        {
                            _semantic_name(item.get("quantity_kind"))
                            for item in outputs
                        }
                    )
                ),
                semantic_parameters=tuple(sorted(parameters)),
            )
        )

    return PaperWorkflowObservationV1(
        schema_version="chemsmart.paper-workflow-observation.v1",
        projects=tuple(sorted(observed_projects.values(), key=lambda item: item.project_id)),
        nodes=tuple(sorted((*calculation_nodes, *analysis_nodes), key=lambda item: item.node_id)),
    )


@dataclass(frozen=True)
class PaperWorkflowAnswerGradeV1:
    """Primary strict grade for a non-executed canonical YAML-CLI DAG."""

    schema_version: str
    case_id: str
    case_sha256: str
    answer_key_sha256: str
    missing_project_signatures: tuple[str, ...]
    unexpected_project_signatures: tuple[str, ...]
    missing_node_signatures: tuple[str, ...]
    unexpected_node_signatures: tuple[str, ...]
    strict_pass: bool
    grade_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.paper-workflow-answer-grade.v1":
            raise ContractError("unsupported paper workflow answer grade schema")
        require_identifier(self.case_id, "case_id")
        require_sha256(self.case_sha256, "case_sha256")
        require_sha256(self.answer_key_sha256, "answer_key_sha256")
        for values in (
            self.missing_project_signatures,
            self.unexpected_project_signatures,
            self.missing_node_signatures,
            self.unexpected_node_signatures,
        ):
            for digest in values:
                require_sha256(digest, "workflow_signature")
        expected_pass = not any(
            (
                self.missing_project_signatures,
                self.unexpected_project_signatures,
                self.missing_node_signatures,
                self.unexpected_node_signatures,
            )
        )
        if self.strict_pass != expected_pass:
            raise ContractError("workflow strict-pass disposition is inconsistent")
        if self.grade_sha256 != canonical_sha256(self._body()):
            raise ContractError("paper workflow answer grade digest mismatch")

    def _body(self) -> dict[str, Any]:
        return {
            key: value
            for key, value in self.__dict__.items()
            if key != "grade_sha256"
        }


def grade_paper_workflow_answer(
    *,
    case: PaperBlackBoxCaseV1,
    answer_key: PaperWorkflowAnswerKeyV1,
    observation: PaperWorkflowObservationV1,
) -> PaperWorkflowAnswerGradeV1:
    """Compare exact project semantics and a program-relative DAG answer.

    Model-selected node identifiers do not affect the result.  Project-YAML
    contents, target program/job form, exact input geometries, producer/consumer
    edge roles, analysis operations, and declared outputs do.
    """

    if answer_key.case_id != case.case_id or answer_key.case_sha256 != case.case_sha256:
        raise ContractError("paper workflow answer key targets another case")
    expected_projects, expected_nodes = _canonical_workflow_signatures(
        answer_key.projects, answer_key.nodes
    )
    observed_projects, observed_nodes = _canonical_workflow_signatures(
        observation.projects, observation.nodes
    )

    def difference(left: Sequence[str], right: Sequence[str]) -> tuple[str, ...]:
        return tuple(sorted((Counter(left) - Counter(right)).elements()))

    body = {
        "schema_version": "chemsmart.paper-workflow-answer-grade.v1",
        "case_id": case.case_id,
        "case_sha256": case.case_sha256,
        "answer_key_sha256": answer_key.answer_key_sha256,
        "missing_project_signatures": difference(
            expected_projects, observed_projects
        ),
        "unexpected_project_signatures": difference(
            observed_projects, expected_projects
        ),
        "missing_node_signatures": difference(expected_nodes, observed_nodes),
        "unexpected_node_signatures": difference(observed_nodes, expected_nodes),
    }
    body["strict_pass"] = not any(
        body[name]
        for name in (
            "missing_project_signatures",
            "unexpected_project_signatures",
            "missing_node_signatures",
            "unexpected_node_signatures",
        )
    )
    return PaperWorkflowAnswerGradeV1(
        **body, grade_sha256=canonical_sha256(body)
    )


def summarize_api_usage(events_file: str | Path) -> dict[str, int]:
    """Sum provider-reported usage from the public Runtime V2 event stream."""

    path = Path(events_file).expanduser().resolve()
    totals = {
        "attempts": 0,
        "input_tokens": 0,
        "output_tokens": 0,
        "reasoning_tokens": 0,
        "latency_ms": 0,
    }
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        event = json.loads(line)
        if event.get("kind") != "api_attempt_observed":
            continue
        payload = event.get("payload") or {}
        totals["attempts"] += 1
        for key in (
            "input_tokens",
            "output_tokens",
            "reasoning_tokens",
            "latency_ms",
        ):
            totals[key] += int(payload.get(key) or 0)
    return totals


__all__ = [
    "MultiRecordXyzObservationV1",
    "NumericalQuantityGradeV1",
    "PaperAnswerKeyV1",
    "PaperBenchmarkEligibilityV1",
    "PaperBlackBoxCaseV1",
    "PaperModelInputV1",
    "PaperContextArmV1",
    "PaperContextRecordV1",
    "PaperWorkflowGradeV1",
    "PaperWorkflowOracleV1",
    "PaperNumericalGradeV1",
    "PaperWorkflowAnswerGradeV1",
    "PaperWorkflowAnswerKeyV1",
    "PaperWorkflowObservationV1",
    "ReferenceQuantityV1",
    "WorkflowDependencyAnswerV1",
    "WorkflowProjectAnswerV1",
    "WorkflowSemanticNodeV1",
    "build_paper_answer_key",
    "build_paper_black_box_case",
    "build_paper_context_arms",
    "build_paper_workflow_answer_key",
    "build_reference_quantity",
    "analysis_claim_record_from_live_result",
    "paper_workflow_observation_from_live_result",
    "extract_multi_xyz_record",
    "grade_paper_workflow_result",
    "grade_paper_numerical_claims",
    "grade_paper_workflow_answer",
    "prepare_benchmark_dispatch",
    "summarize_api_usage",
]
