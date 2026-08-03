"""Host-owned contracts for approved program execution.

This module deliberately sits after command compilation and safe preview.  It
does not expose shell construction, arbitrary paths, resources, approval, or
readiness to a model.  A model-facing execution tool supplies only a stable
``node_id``; the host resolves that identifier through the contracts below.

The helpers are intentionally additive to Runtime V2.  They do not change the
meaning of existing preview receipts and they keep engine completion separate
from scientific validation, which permits older successful result artifacts
to be inspected without treating missing modern provenance as a reason to
rerun them.
"""

from __future__ import annotations

import math
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

from chemsmart.agent._contracts import (
    ContractError,
    TrustedArtifactRefV1,
    canonical_data,
    canonical_json,
    canonical_sha256,
    file_sha256,
    require_identifier,
    require_sha256,
)
from chemsmart.agent.projects import (
    ProjectRenderReceiptV1,
    ProjectValidationReceiptV1,
)


def _require_nonempty_rows(values: Iterable[str], field_name: str) -> None:
    if any(not str(value).strip() for value in values):
        raise ContractError(f"{field_name} cannot contain empty values")


def _require_sorted_unique(values: tuple[str, ...], field_name: str) -> None:
    if values != tuple(sorted(set(values))):
        raise ContractError(f"{field_name} must be sorted and unique")


def _require_current_artifact(
    artifact: TrustedArtifactRefV1, field_name: str
) -> Path:
    path = Path(artifact.path)
    if not path.is_file() or path.is_symlink():
        raise ContractError(f"{field_name} must be a current regular file")
    if path.stat().st_size != artifact.size_bytes:
        raise ContractError(f"{field_name} size differs from its binding")
    if file_sha256(path) != artifact.sha256:
        raise ContractError(f"{field_name} digest differs from its binding")
    return path.resolve()


def _absolute_workspace(value: str | Path) -> Path:
    workspace = Path(value)
    if not workspace.is_absolute():
        raise ContractError("approved workspace must be absolute")
    workspace.mkdir(parents=True, exist_ok=True)
    if workspace.is_symlink():
        raise ContractError("approved workspace cannot be a symbolic link")
    return workspace.resolve()


def _target_below(workspace: Path, *parts: str) -> Path:
    target = workspace.joinpath(*parts)
    target.parent.mkdir(parents=True, exist_ok=True)
    resolved_parent = target.parent.resolve()
    try:
        resolved_parent.relative_to(workspace)
    except ValueError as exc:
        raise ContractError(
            "artifact target escapes the approved workspace"
        ) from exc
    if target.is_symlink():
        raise ContractError("artifact target cannot be a symbolic link")
    return resolved_parent / target.name


def _write_exact_once(path: Path, payload: bytes) -> None:
    """Create ``path`` once, accepting an exact existing artifact.

    The target is always host-derived.  A conflicting existing file is never
    overwritten, while an identical file makes promotion and handoff
    naturally idempotent.
    """

    if path.exists():
        if (
            path.is_symlink()
            or not path.is_file()
            or path.read_bytes() != payload
        ):
            raise ContractError(
                "existing artifact differs from requested bytes"
            )
        return
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    descriptor = os.open(path, flags, 0o600)
    try:
        view = memoryview(payload)
        while view:
            written = os.write(descriptor, view)
            view = view[written:]
        os.fsync(descriptor)
    except Exception:
        os.close(descriptor)
        descriptor = -1
        try:
            path.unlink()
        except FileNotFoundError:
            pass
        raise
    finally:
        if descriptor >= 0:
            os.close(descriptor)


@dataclass(frozen=True)
class ScientificDecisionRecordV1:
    """Public scientific rationale, never hidden provider reasoning."""

    schema_version: str
    decision_id: str
    task_spec_sha256: str
    stage_order: tuple[str, ...]
    assumptions: tuple[str, ...]
    method_rationale: str
    alternatives: tuple[str, ...]
    uncertainties: tuple[str, ...]
    diagnostics: tuple[str, ...]
    evidence_refs: tuple[str, ...]
    record_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.scientific-decision-record.v1":
            raise ContractError("unsupported scientific decision schema")
        require_identifier(self.decision_id, "decision_id")
        require_sha256(self.task_spec_sha256, "task_spec_sha256")
        if not self.stage_order:
            raise ContractError("scientific decision requires a stage order")
        for stage in self.stage_order:
            require_identifier(stage, "stage_order")
        if len(set(self.stage_order)) != len(self.stage_order):
            raise ContractError("stage_order cannot repeat a stage")
        if not self.method_rationale.strip():
            raise ContractError("method_rationale must not be empty")
        for name, rows in (
            ("assumptions", self.assumptions),
            ("alternatives", self.alternatives),
            ("uncertainties", self.uncertainties),
            ("diagnostics", self.diagnostics),
            ("evidence_refs", self.evidence_refs),
        ):
            _require_nonempty_rows(rows, name)
        body = {
            "schema_version": self.schema_version,
            "decision_id": self.decision_id,
            "task_spec_sha256": self.task_spec_sha256,
            "stage_order": self.stage_order,
            "assumptions": self.assumptions,
            "method_rationale": self.method_rationale,
            "alternatives": self.alternatives,
            "uncertainties": self.uncertainties,
            "diagnostics": self.diagnostics,
            "evidence_refs": self.evidence_refs,
        }
        if self.record_sha256 != canonical_sha256(body):
            raise ContractError("scientific decision digest mismatch")


def build_scientific_decision_record(
    *,
    decision_id: str,
    task_spec_sha256: str,
    stage_order: Sequence[str],
    assumptions: Sequence[str],
    method_rationale: str,
    alternatives: Sequence[str] = (),
    uncertainties: Sequence[str] = (),
    diagnostics: Sequence[str] = (),
    evidence_refs: Sequence[str] = (),
) -> ScientificDecisionRecordV1:
    """Build a durable public decision summary from visible statements."""

    body = {
        "schema_version": "chemsmart.scientific-decision-record.v1",
        "decision_id": require_identifier(decision_id, "decision_id"),
        "task_spec_sha256": require_sha256(
            task_spec_sha256, "task_spec_sha256"
        ),
        "stage_order": tuple(
            str(item).strip().lower() for item in stage_order
        ),
        "assumptions": tuple(str(item).strip() for item in assumptions),
        "method_rationale": str(method_rationale).strip(),
        "alternatives": tuple(str(item).strip() for item in alternatives),
        "uncertainties": tuple(str(item).strip() for item in uncertainties),
        "diagnostics": tuple(str(item).strip() for item in diagnostics),
        "evidence_refs": tuple(str(item).strip() for item in evidence_refs),
    }
    return ScientificDecisionRecordV1(
        **body, record_sha256=canonical_sha256(body)
    )


@dataclass(frozen=True)
class ExecutionResourceSpecV1:
    """Host-owned compute allocation for every node in one workflow."""

    schema_version: str
    execution_target: str
    cores: int
    memory_gb: float
    gpu_count: int
    scratch_policy: str
    node_timeout_seconds: int
    resource_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.execution-resource-spec.v1":
            raise ContractError("unsupported execution resource schema")
        if self.execution_target not in {"run", "sub"}:
            raise ContractError("execution_target must be run or sub")
        if self.cores < 1:
            raise ContractError("cores must be positive")
        if not math.isfinite(self.memory_gb) or self.memory_gb <= 0:
            raise ContractError("memory_gb must be finite and positive")
        if self.gpu_count < 0:
            raise ContractError("gpu_count must be non-negative")
        if self.scratch_policy not in {"none", "workspace", "server"}:
            raise ContractError("unsupported scratch policy")
        if self.node_timeout_seconds < 1:
            raise ContractError("node timeout must be positive")
        body = {
            "schema_version": self.schema_version,
            "execution_target": self.execution_target,
            "cores": self.cores,
            "memory_gb": self.memory_gb,
            "gpu_count": self.gpu_count,
            "scratch_policy": self.scratch_policy,
            "node_timeout_seconds": self.node_timeout_seconds,
        }
        if self.resource_sha256 != canonical_sha256(body):
            raise ContractError("execution resource digest mismatch")


def build_execution_resource_spec(
    *,
    execution_target: str,
    cores: int,
    memory_gb: float,
    gpu_count: int,
    scratch_policy: str,
    node_timeout_seconds: int,
) -> ExecutionResourceSpecV1:
    body = {
        "schema_version": "chemsmart.execution-resource-spec.v1",
        "execution_target": str(execution_target).strip().lower(),
        "cores": int(cores),
        "memory_gb": float(memory_gb),
        "gpu_count": int(gpu_count),
        "scratch_policy": str(scratch_policy).strip().lower(),
        "node_timeout_seconds": int(node_timeout_seconds),
    }
    return ExecutionResourceSpecV1(
        **body, resource_sha256=canonical_sha256(body)
    )


@dataclass(frozen=True)
class ProjectArtifactPromotionV1:
    """Receipt for a candidate materialized in a host-owned workspace."""

    schema_version: str
    artifact_id: str
    workspace: str
    render_receipt_sha256: str
    rendered_sha256: str
    project_artifact_sha256: str
    validation_receipt_sha256: str
    validation_status: str
    status: str
    rule_ids: tuple[str, ...]
    receipt_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.project-artifact-promotion.v1":
            raise ContractError("unsupported project promotion schema")
        require_identifier(self.artifact_id, "artifact_id")
        if not Path(self.workspace).is_absolute():
            raise ContractError("promotion workspace must be absolute")
        for name, digest in (
            ("render_receipt_sha256", self.render_receipt_sha256),
            ("rendered_sha256", self.rendered_sha256),
            ("project_artifact_sha256", self.project_artifact_sha256),
        ):
            require_sha256(digest, name)
        if self.rendered_sha256 != self.project_artifact_sha256:
            raise ContractError("promoted project bytes differ from rendering")
        if self.validation_status not in {"pending", "valid", "invalid"}:
            raise ContractError("unsupported project validation state")
        if self.status not in {"materialized", "validated", "rejected"}:
            raise ContractError("unsupported project promotion state")
        expected_status = {
            "pending": "materialized",
            "valid": "validated",
            "invalid": "rejected",
        }[self.validation_status]
        if self.status != expected_status:
            raise ContractError(
                "project promotion and validation states differ"
            )
        if self.validation_status == "pending":
            if self.validation_receipt_sha256:
                raise ContractError("pending promotion cannot bind validation")
        else:
            require_sha256(
                self.validation_receipt_sha256,
                "validation_receipt_sha256",
            )
        _require_sorted_unique(self.rule_ids, "rule_ids")
        body = {
            "schema_version": self.schema_version,
            "artifact_id": self.artifact_id,
            "workspace": self.workspace,
            "render_receipt_sha256": self.render_receipt_sha256,
            "rendered_sha256": self.rendered_sha256,
            "project_artifact_sha256": self.project_artifact_sha256,
            "validation_receipt_sha256": self.validation_receipt_sha256,
            "validation_status": self.validation_status,
            "status": self.status,
            "rule_ids": self.rule_ids,
        }
        if self.receipt_sha256 != canonical_sha256(body):
            raise ContractError("project promotion receipt digest mismatch")


def promote_project_candidate(
    render_receipt: ProjectRenderReceiptV1,
    *,
    approved_workspace: str | Path,
    artifact_id: str,
) -> tuple[TrustedArtifactRefV1, ProjectArtifactPromotionV1]:
    """Materialize rendered YAML without claiming loader validation."""

    workspace = _absolute_workspace(approved_workspace)
    normalized_id = require_identifier(artifact_id, "artifact_id")
    target = _target_below(workspace, "projects", f"{normalized_id}.yaml")
    payload = render_receipt.rendered_yaml.encode("utf-8")
    _write_exact_once(target, payload)
    observed = file_sha256(target)
    if observed != render_receipt.rendered_sha256:
        raise ContractError("materialized project differs from render receipt")
    artifact = TrustedArtifactRefV1(
        artifact_id=normalized_id,
        kind="project_yaml",
        sha256=observed,
        size_bytes=target.stat().st_size,
        path=str(target),
        cli_value=str(target),
    )
    body = {
        "schema_version": "chemsmart.project-artifact-promotion.v1",
        "artifact_id": normalized_id,
        "workspace": str(workspace),
        "render_receipt_sha256": render_receipt.receipt_sha256,
        "rendered_sha256": render_receipt.rendered_sha256,
        "project_artifact_sha256": artifact.sha256,
        "validation_receipt_sha256": "",
        "validation_status": "pending",
        "status": "materialized",
        "rule_ids": ("project.promotion.host_workspace",),
    }
    return artifact, ProjectArtifactPromotionV1(
        **body, receipt_sha256=canonical_sha256(body)
    )


def bind_project_promotion_validation(
    promotion: ProjectArtifactPromotionV1,
    artifact: TrustedArtifactRefV1,
    validation: ProjectValidationReceiptV1,
) -> ProjectArtifactPromotionV1:
    """Bind the program loader's independent result to a promotion."""

    _require_current_artifact(artifact, "project artifact")
    if promotion.validation_status != "pending":
        raise ContractError("only a pending promotion can bind validation")
    if artifact.artifact_id != promotion.artifact_id:
        raise ContractError("validation artifact ID differs from promotion")
    if artifact.sha256 != promotion.project_artifact_sha256:
        raise ContractError(
            "validation artifact digest differs from promotion"
        )
    if validation.project_artifact_id != artifact.artifact_id:
        raise ContractError("loader validation refers to another artifact")
    if validation.project_sha256 != artifact.sha256:
        raise ContractError("loader validation refers to different bytes")
    validation_status = "valid" if validation.status == "valid" else "invalid"
    status = "validated" if validation_status == "valid" else "rejected"
    rule = (
        "project.promotion.loader_valid"
        if validation_status == "valid"
        else "project.promotion.loader_rejected"
    )
    body = {
        "schema_version": promotion.schema_version,
        "artifact_id": promotion.artifact_id,
        "workspace": promotion.workspace,
        "render_receipt_sha256": promotion.render_receipt_sha256,
        "rendered_sha256": promotion.rendered_sha256,
        "project_artifact_sha256": promotion.project_artifact_sha256,
        "validation_receipt_sha256": validation.receipt_sha256,
        "validation_status": validation_status,
        "status": status,
        "rule_ids": tuple(sorted(set(promotion.rule_ids + (rule,)))),
    }
    return ProjectArtifactPromotionV1(
        **body, receipt_sha256=canonical_sha256(body)
    )


@dataclass(frozen=True)
class ProgramConformanceProbeV1:
    """Observed live plumbing; never an execution-readiness assertion."""

    schema_version: str
    probe_id: str
    program: str
    engine: str
    registry_sha256: str
    live_cli_schema_sha256: str
    covered_jobtypes: tuple[str, ...]
    loader_status: str
    environment_status: str
    preview_status: str
    project_validation_receipt_sha256s: tuple[str, ...]
    environment_receipt_sha256: str
    preview_receipt_sha256s: tuple[str, ...]
    findings: tuple[str, ...]
    status: str
    receipt_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.program-conformance-probe.v1":
            raise ContractError("unsupported program conformance schema")
        require_identifier(self.probe_id, "probe_id")
        require_identifier(self.program, "program")
        require_identifier(self.engine, "engine")
        require_sha256(self.registry_sha256, "registry_sha256")
        require_sha256(self.live_cli_schema_sha256, "live_cli_schema_sha256")
        _require_sorted_unique(self.covered_jobtypes, "covered_jobtypes")
        for jobtype in self.covered_jobtypes:
            require_identifier(jobtype, "covered_jobtype")
        allowed_component_statuses = {
            "passed",
            "failed",
            "not_required",
            "not_observed",
            "available",
            "missing",
        }
        for value in (
            self.loader_status,
            self.environment_status,
            self.preview_status,
        ):
            if value not in allowed_component_statuses:
                raise ContractError("invalid conformance component status")
        for name, digests in (
            (
                "project_validation_receipt_sha256s",
                self.project_validation_receipt_sha256s,
            ),
            ("preview_receipt_sha256s", self.preview_receipt_sha256s),
        ):
            _require_sorted_unique(digests, name)
            for digest in digests:
                require_sha256(digest, name)
        if self.environment_receipt_sha256:
            require_sha256(
                self.environment_receipt_sha256,
                "environment_receipt_sha256",
            )
        if self.loader_status == "passed" and not (
            self.project_validation_receipt_sha256s
        ):
            raise ContractError(
                "passed loader probe requires validation evidence"
            )
        if self.environment_status == "available" and not (
            self.environment_receipt_sha256
        ):
            raise ContractError(
                "available environment requires an observation"
            )
        if (
            self.preview_status == "passed"
            and not self.preview_receipt_sha256s
        ):
            raise ContractError(
                "passed preview probe requires preview evidence"
            )
        if self.status not in {
            "available_for_preflight",
            "preview_only",
            "blocked",
            "incomplete",
        }:
            raise ContractError("invalid program conformance result")
        _require_nonempty_rows(self.findings, "findings")
        body = {
            "schema_version": self.schema_version,
            "probe_id": self.probe_id,
            "program": self.program,
            "engine": self.engine,
            "registry_sha256": self.registry_sha256,
            "live_cli_schema_sha256": self.live_cli_schema_sha256,
            "covered_jobtypes": self.covered_jobtypes,
            "loader_status": self.loader_status,
            "environment_status": self.environment_status,
            "preview_status": self.preview_status,
            "project_validation_receipt_sha256s": (
                self.project_validation_receipt_sha256s
            ),
            "environment_receipt_sha256": self.environment_receipt_sha256,
            "preview_receipt_sha256s": self.preview_receipt_sha256s,
            "findings": self.findings,
            "status": self.status,
        }
        if self.receipt_sha256 != canonical_sha256(body):
            raise ContractError("program conformance receipt digest mismatch")


def build_program_conformance_probe(
    *,
    probe_id: str,
    program: str,
    engine: str,
    registry_sha256: str,
    live_cli_schema_sha256: str,
    covered_jobtypes: Sequence[str],
    loader_status: str,
    environment_status: str,
    preview_status: str,
    project_validation_receipt_sha256s: Sequence[str] = (),
    environment_receipt_sha256: str = "",
    preview_receipt_sha256s: Sequence[str] = (),
    findings: Sequence[str] = (),
) -> ProgramConformanceProbeV1:
    """Summarize live observations without claiming node readiness."""

    component_values = (loader_status, environment_status, preview_status)
    if "failed" in component_values or environment_status == "missing":
        status = "blocked"
    elif (
        loader_status in {"passed", "not_required"}
        and preview_status == "passed"
        and environment_status == "available"
    ):
        status = "available_for_preflight"
    elif (
        loader_status in {"passed", "not_required"}
        and preview_status == "passed"
    ):
        status = "preview_only"
    else:
        status = "incomplete"
    body = {
        "schema_version": "chemsmart.program-conformance-probe.v1",
        "probe_id": require_identifier(probe_id, "probe_id"),
        "program": require_identifier(program, "program"),
        "engine": require_identifier(engine, "engine"),
        "registry_sha256": require_sha256(
            registry_sha256, "registry_sha256"
        ),
        "live_cli_schema_sha256": require_sha256(
            live_cli_schema_sha256, "live_cli_schema_sha256"
        ),
        "covered_jobtypes": tuple(sorted(set(covered_jobtypes))),
        "loader_status": str(loader_status),
        "environment_status": str(environment_status),
        "preview_status": str(preview_status),
        "project_validation_receipt_sha256s": tuple(
            sorted(set(project_validation_receipt_sha256s))
        ),
        "environment_receipt_sha256": str(environment_receipt_sha256),
        "preview_receipt_sha256s": tuple(
            sorted(set(preview_receipt_sha256s))
        ),
        "findings": tuple(str(item).strip() for item in findings),
        "status": status,
    }
    return ProgramConformanceProbeV1(
        **body, receipt_sha256=canonical_sha256(body)
    )


@dataclass(frozen=True)
class ProducerEdgeRuleV1:
    """Approval rule for selecting a future producer artifact."""

    producer_node_id: str
    consumer_node_id: str
    artifact_kind: str
    selection_rule: str
    edge_sha256: str

    def __post_init__(self) -> None:
        require_identifier(self.producer_node_id, "producer_node_id")
        require_identifier(self.consumer_node_id, "consumer_node_id")
        require_identifier(self.artifact_kind, "artifact_kind")
        require_identifier(self.selection_rule, "selection_rule")
        if self.producer_node_id == self.consumer_node_id:
            raise ContractError("producer and consumer nodes must differ")
        body = {
            "producer_node_id": self.producer_node_id,
            "consumer_node_id": self.consumer_node_id,
            "artifact_kind": self.artifact_kind,
            "selection_rule": self.selection_rule,
        }
        if self.edge_sha256 != canonical_sha256(body):
            raise ContractError("producer edge digest mismatch")


def build_producer_edge_rule(
    *,
    producer_node_id: str,
    consumer_node_id: str,
    artifact_kind: str,
    selection_rule: str,
) -> ProducerEdgeRuleV1:
    body = {
        "producer_node_id": require_identifier(
            producer_node_id, "producer_node_id"
        ),
        "consumer_node_id": require_identifier(
            consumer_node_id, "consumer_node_id"
        ),
        "artifact_kind": require_identifier(artifact_kind, "artifact_kind"),
        "selection_rule": require_identifier(selection_rule, "selection_rule"),
    }
    return ProducerEdgeRuleV1(**body, edge_sha256=canonical_sha256(body))


@dataclass(frozen=True)
class ApprovedNodeBindingV1:
    """Host-resolved immutable scientific settings for one approved node."""

    node_id: str
    program: str
    engine: str
    jobtype: str
    project_artifact_sha256: str
    settings_sha256: str
    charge: int
    multiplicity: int
    input_mode: str
    initial_artifact_id: str
    initial_artifact_sha256: str
    scientific_identity_sha256: str
    producer_edge_sha256: str

    def __post_init__(self) -> None:
        for name, value in (
            ("node_id", self.node_id),
            ("program", self.program),
            ("engine", self.engine),
            ("jobtype", self.jobtype),
        ):
            require_identifier(value, name)
        require_sha256(
            self.project_artifact_sha256, "project_artifact_sha256"
        )
        require_sha256(self.settings_sha256, "settings_sha256")
        if self.multiplicity < 1:
            raise ContractError("multiplicity must be positive")
        if self.input_mode not in {"initial", "producer"}:
            raise ContractError("input_mode must be initial or producer")
        if self.input_mode == "initial":
            if not self.initial_artifact_id:
                raise ContractError("initial node requires an artifact ID")
            require_sha256(
                self.initial_artifact_sha256, "initial_artifact_sha256"
            )
            require_sha256(
                self.scientific_identity_sha256,
                "scientific_identity_sha256",
            )
            if self.producer_edge_sha256:
                raise ContractError("initial node cannot bind a producer edge")
        else:
            if self.initial_artifact_id or self.initial_artifact_sha256:
                raise ContractError("producer input cannot bind initial bytes")
            if self.scientific_identity_sha256:
                raise ContractError(
                    "future geometry identity cannot be precomputed"
                )
            require_sha256(self.producer_edge_sha256, "producer_edge_sha256")


@dataclass(frozen=True)
class WorkflowExecutionApprovalV1:
    """One user approval for exact nodes and typed producer-edge rules."""

    schema_version: str
    approval_id: str
    workflow_id: str
    workflow_sha256: str
    task_spec_sha256: str
    workspace: str
    resource_sha256: str
    node_bindings: tuple[ApprovedNodeBindingV1, ...]
    producer_edges: tuple[ProducerEdgeRuleV1, ...]
    status: str
    approval_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.workflow-execution-approval.v1":
            raise ContractError("unsupported workflow approval schema")
        require_identifier(self.approval_id, "approval_id")
        require_identifier(self.workflow_id, "workflow_id")
        require_sha256(self.workflow_sha256, "workflow_sha256")
        require_sha256(self.task_spec_sha256, "task_spec_sha256")
        require_sha256(self.resource_sha256, "resource_sha256")
        if not Path(self.workspace).is_absolute():
            raise ContractError(
                "approved execution workspace must be absolute"
            )
        if not self.node_bindings:
            raise ContractError("approval requires at least one node")
        node_ids = tuple(item.node_id for item in self.node_bindings)
        if len(set(node_ids)) != len(node_ids):
            raise ContractError("approval node IDs must be unique")
        edge_digests = tuple(item.edge_sha256 for item in self.producer_edges)
        if len(set(edge_digests)) != len(edge_digests):
            raise ContractError("approval producer edges must be unique")
        edge_by_digest = {
            item.edge_sha256: item for item in self.producer_edges
        }
        for binding in self.node_bindings:
            if binding.input_mode == "producer":
                edge = edge_by_digest.get(binding.producer_edge_sha256)
                if edge is None or edge.consumer_node_id != binding.node_id:
                    raise ContractError(
                        "producer node binding has no approved edge"
                    )
                if edge.producer_node_id not in node_ids:
                    raise ContractError(
                        "producer edge names an unapproved node"
                    )
        if self.status != "approved":
            raise ContractError("workflow approval record must be approved")
        body = {
            "schema_version": self.schema_version,
            "approval_id": self.approval_id,
            "workflow_id": self.workflow_id,
            "workflow_sha256": self.workflow_sha256,
            "task_spec_sha256": self.task_spec_sha256,
            "workspace": self.workspace,
            "resource_sha256": self.resource_sha256,
            "node_bindings": self.node_bindings,
            "producer_edges": self.producer_edges,
            "status": self.status,
        }
        if self.approval_sha256 != canonical_sha256(body):
            raise ContractError("workflow approval digest mismatch")

    def node(self, node_id: str) -> ApprovedNodeBindingV1:
        normalized = require_identifier(node_id, "node_id")
        match = next(
            (
                item
                for item in self.node_bindings
                if item.node_id == normalized
            ),
            None,
        )
        if match is None:
            raise ContractError("node is not covered by workflow approval")
        return match


def build_workflow_execution_approval(
    *,
    approval_id: str,
    workflow_id: str,
    workflow_sha256: str,
    task_spec_sha256: str,
    approved_workspace: str | Path,
    resources: ExecutionResourceSpecV1,
    node_bindings: Sequence[ApprovedNodeBindingV1],
    producer_edges: Sequence[ProducerEdgeRuleV1] = (),
) -> WorkflowExecutionApprovalV1:
    """Create the host record consumed by the execution-only tool profile."""

    workspace = _absolute_workspace(approved_workspace)
    body = {
        "schema_version": "chemsmart.workflow-execution-approval.v1",
        "approval_id": require_identifier(approval_id, "approval_id"),
        "workflow_id": require_identifier(workflow_id, "workflow_id"),
        "workflow_sha256": require_sha256(workflow_sha256, "workflow_sha256"),
        "task_spec_sha256": require_sha256(
            task_spec_sha256, "task_spec_sha256"
        ),
        "workspace": str(workspace),
        "resource_sha256": resources.resource_sha256,
        "node_bindings": tuple(node_bindings),
        "producer_edges": tuple(producer_edges),
        "status": "approved",
    }
    return WorkflowExecutionApprovalV1(
        **body, approval_sha256=canonical_sha256(body)
    )


def build_locked_pyscf_sp_opt_hess_approval(
    *,
    approval_id: str,
    workflow_id: str,
    task_spec_sha256: str,
    approved_workspace: str | Path,
    resources: ExecutionResourceSpecV1,
    initial_artifact: TrustedArtifactRefV1,
    project_artifact: TrustedArtifactRefV1,
) -> WorkflowExecutionApprovalV1:
    """Build the fixed CPU water integration-workflow approval.

    The project is read through the production PySCF YAML loader.  Stage
    settings digests are derived with the same registry projection used by
    :func:`chemsmart.agent.projects.validate_project_yaml`, avoiding a
    hand-maintained digest or an approval that depends on model-selected
    artifact IDs.  The exact initial geometry is approved for SP and OPT;
    HESS is approved only through the validated OPT producer-edge rule.

    This helper does not execute, preview, or claim readiness.  It merely
    constructs the private approval consumed by the existing host gates.
    """

    workspace = _absolute_workspace(approved_workspace)
    initial_path = _require_current_artifact(
        initial_artifact, "initial geometry artifact"
    )
    project_path = _require_current_artifact(
        project_artifact, "PySCF project artifact"
    )
    for label, path in (
        ("initial geometry", initial_path),
        ("PySCF project", project_path),
    ):
        try:
            path.relative_to(workspace)
        except ValueError as exc:
            raise ContractError(
                f"{label} must be inside the approved workspace"
            ) from exc
    if initial_artifact.kind != "geometry_xyz":
        raise ContractError("locked workflow requires an exact XYZ artifact")
    if project_artifact.kind != "project_yaml":
        raise ContractError("locked workflow requires a project YAML artifact")
    if resources.execution_target != "run":
        raise ContractError("locked workflow requires local run execution")
    if (
        resources.cores != 4
        or resources.memory_gb != 4.0
        or resources.gpu_count != 0
        or resources.scratch_policy != "none"
    ):
        raise ContractError("resources differ from the locked CPU profile")

    from chemsmart.agent.commands import build_scientific_identity_binding
    from chemsmart.agent.capabilities import load_program_capabilities
    from chemsmart.settings.pyscf import YamlPySCFProjectSettings

    identity = build_scientific_identity_binding(
        task_spec_sha256=task_spec_sha256,
        geometry_artifact=initial_artifact,
        charge=0,
        multiplicity=1,
    )
    project = YamlPySCFProjectSettings.from_yaml(project_path)
    capability = load_program_capabilities().get("pyscf")
    if capability is None:
        raise ContractError("PySCF is absent from the program registry")
    allowed = set(capability.project_owned_parameters)
    allowed.update({"jobtype", "engine", "freq"})
    settings_digests: dict[str, str] = {}
    for jobtype in ("sp", "opt", "hess"):
        settings = getattr(project, f"{jobtype}_settings")()
        settings.validate()
        _require_locked_pyscf_settings(settings, jobtype=jobtype)
        effective = {
            key: canonical_data(getattr(settings, key))
            for key in sorted(allowed)
            if hasattr(settings, key)
        }
        settings_digests[jobtype] = canonical_sha256(effective)

    edge = build_producer_edge_rule(
        producer_node_id="opt-initial",
        consumer_node_id="hess-optimized",
        artifact_kind="geometry_xyz",
        selection_rule="validated_optimized_geometry",
    )
    initial_common = {
        "program": "pyscf",
        "engine": "cpu",
        "project_artifact_sha256": project_artifact.sha256,
        "charge": 0,
        "multiplicity": 1,
        "input_mode": "initial",
        "initial_artifact_id": initial_artifact.artifact_id,
        "initial_artifact_sha256": initial_artifact.sha256,
        "scientific_identity_sha256": identity.binding_sha256,
        "producer_edge_sha256": "",
    }
    nodes = (
        ApprovedNodeBindingV1(
            node_id="sp-initial",
            jobtype="sp",
            settings_sha256=settings_digests["sp"],
            **initial_common,
        ),
        ApprovedNodeBindingV1(
            node_id="opt-initial",
            jobtype="opt",
            settings_sha256=settings_digests["opt"],
            **initial_common,
        ),
        ApprovedNodeBindingV1(
            node_id="hess-optimized",
            program="pyscf",
            engine="cpu",
            jobtype="hess",
            project_artifact_sha256=project_artifact.sha256,
            settings_sha256=settings_digests["hess"],
            charge=0,
            multiplicity=1,
            input_mode="producer",
            initial_artifact_id="",
            initial_artifact_sha256="",
            scientific_identity_sha256="",
            producer_edge_sha256=edge.edge_sha256,
        ),
    )
    workflow_sha256 = canonical_sha256(
        {
            "schema_version": "chemsmart.pyscf-sp-opt-hess-workflow.v1",
            "task_spec_sha256": task_spec_sha256,
            "nodes": nodes,
            "producer_edges": (edge,),
        }
    )
    return build_workflow_execution_approval(
        approval_id=approval_id,
        workflow_id=workflow_id,
        workflow_sha256=workflow_sha256,
        task_spec_sha256=task_spec_sha256,
        approved_workspace=workspace,
        resources=resources,
        node_bindings=nodes,
        producer_edges=(edge,),
    )


def workflow_execution_approval_json(
    approval: WorkflowExecutionApprovalV1,
) -> str:
    """Return the exact private JSON envelope accepted by ``agent run``."""

    return canonical_json({"workflow_approval": approval}) + "\n"


def _require_locked_pyscf_settings(settings: object, *, jobtype: str) -> None:
    expected = {
        "functional": "b3lyp",
        "basis": "def2-svp",
        "defgrid": "defgrid2",
        "density_fit": False,
        "dispersion": None,
        "solvent_model": None,
        "solvent_id": None,
        "scf_tol": 1e-9,
        "scf_maxiter": 100,
        "engine": "cpu",
        "jobtype": jobtype,
    }
    if jobtype == "opt":
        expected.update({"opt_solver": "geometric", "opt_maxsteps": 100})
    if jobtype == "hess":
        expected["freq"] = True
    for name, required in expected.items():
        observed = getattr(settings, name, None)
        if isinstance(required, str) and isinstance(observed, str):
            matches = observed.strip().lower() == required
        else:
            matches = observed == required
        if not matches:
            raise ContractError(
                f"locked PySCF {jobtype} setting {name} differs: "
                f"expected {required!r}, observed {observed!r}"
            )
    if getattr(settings, "charge", None) is not None or getattr(
        settings, "multiplicity", None
    ) is not None:
        raise ContractError(
            "locked project must inherit charge and multiplicity from XYZ"
        )


@dataclass(frozen=True)
class ProgramExecutionInvocationV1:
    """Exact host-compiled argv for one approved node; never shell text."""

    schema_version: str
    node_id: str
    approval_sha256: str
    program: str
    engine: str
    jobtype: str
    project_sha256: str
    input_artifact_id: str
    input_sha256: str
    scientific_identity_sha256: str
    environment_receipt_sha256: str
    resource_sha256: str
    workspace: str
    argv: tuple[str, ...]
    idempotency_key: str
    status: str
    invocation_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.program-execution-invocation.v1":
            raise ContractError("unsupported program execution invocation")
        for name, value in (
            ("node_id", self.node_id),
            ("program", self.program),
            ("engine", self.engine),
            ("jobtype", self.jobtype),
        ):
            require_identifier(value, name)
        for name, digest in (
            ("approval_sha256", self.approval_sha256),
            ("project_sha256", self.project_sha256),
            ("input_sha256", self.input_sha256),
            ("scientific_identity_sha256", self.scientific_identity_sha256),
            ("environment_receipt_sha256", self.environment_receipt_sha256),
            ("resource_sha256", self.resource_sha256),
            ("idempotency_key", self.idempotency_key),
        ):
            require_sha256(digest, name)
        if not self.input_artifact_id:
            raise ContractError(
                "execution input artifact ID must not be empty"
            )
        if not Path(self.workspace).is_absolute():
            raise ContractError("execution workspace must be absolute")
        if not self.argv or any(
            not str(item) or "\x00" in item for item in self.argv
        ):
            raise ContractError(
                "execution argv must contain non-empty safe tokens"
            )
        if self.status != "ready":
            raise ContractError("execution invocation must be ready")
        body = {
            "schema_version": self.schema_version,
            "node_id": self.node_id,
            "approval_sha256": self.approval_sha256,
            "program": self.program,
            "engine": self.engine,
            "jobtype": self.jobtype,
            "project_sha256": self.project_sha256,
            "input_artifact_id": self.input_artifact_id,
            "input_sha256": self.input_sha256,
            "scientific_identity_sha256": self.scientific_identity_sha256,
            "environment_receipt_sha256": self.environment_receipt_sha256,
            "resource_sha256": self.resource_sha256,
            "workspace": self.workspace,
            "argv": self.argv,
            "idempotency_key": self.idempotency_key,
            "status": self.status,
        }
        if self.invocation_sha256 != canonical_sha256(body):
            raise ContractError("execution invocation digest mismatch")


def build_program_execution_invocation(
    *,
    node_id: str,
    approval: WorkflowExecutionApprovalV1,
    project_artifact: TrustedArtifactRefV1,
    input_artifact: TrustedArtifactRefV1,
    scientific_identity_sha256: str,
    environment_receipt_sha256: str,
    resources: ExecutionResourceSpecV1,
    argv: Sequence[str],
    handoff: OptimizedGeometryHandoffV1 | None = None,
) -> ProgramExecutionInvocationV1:
    """Resolve a ``node_id`` to an exact invocation using host bindings."""

    node = approval.node(node_id)
    _require_current_artifact(project_artifact, "project artifact")
    _require_current_artifact(input_artifact, "input artifact")
    if resources.resource_sha256 != approval.resource_sha256:
        raise ContractError("resources differ from workflow approval")
    if project_artifact.sha256 != node.project_artifact_sha256:
        raise ContractError("project bytes differ from workflow approval")
    identity_sha256 = require_sha256(
        scientific_identity_sha256, "scientific_identity_sha256"
    )
    if node.input_mode == "initial":
        if input_artifact.artifact_id != node.initial_artifact_id:
            raise ContractError(
                "initial input ID differs from workflow approval"
            )
        if input_artifact.sha256 != node.initial_artifact_sha256:
            raise ContractError(
                "initial input bytes differ from workflow approval"
            )
        if identity_sha256 != node.scientific_identity_sha256:
            raise ContractError(
                "initial scientific identity differs from approval"
            )
        if handoff is not None:
            raise ContractError("initial input node cannot consume a handoff")
    else:
        if handoff is None:
            raise ContractError(
                "producer input node requires a validated handoff"
            )
        if handoff.consumer_node_id != node.node_id:
            raise ContractError("handoff targets another consumer node")
        if handoff.producer_edge_sha256 != node.producer_edge_sha256:
            raise ContractError("handoff edge differs from workflow approval")
        if handoff.geometry_artifact_id != input_artifact.artifact_id:
            raise ContractError(
                "handoff artifact ID differs from execution input"
            )
        if handoff.geometry_artifact_sha256 != input_artifact.sha256:
            raise ContractError("handoff bytes differ from execution input")
        if (handoff.charge, handoff.multiplicity) != (
            node.charge,
            node.multiplicity,
        ):
            raise ContractError(
                "handoff electronic state differs from approval"
            )
    environment_sha256 = require_sha256(
        environment_receipt_sha256, "environment_receipt_sha256"
    )
    argv_tuple = tuple(str(item) for item in argv)
    idempotency_key = canonical_sha256(
        {
            "approval_sha256": approval.approval_sha256,
            "node_id": node.node_id,
            "project_sha256": project_artifact.sha256,
            "input_sha256": input_artifact.sha256,
            "scientific_identity_sha256": identity_sha256,
            "environment_receipt_sha256": environment_sha256,
            "resource_sha256": resources.resource_sha256,
            "argv": argv_tuple,
        }
    )
    body = {
        "schema_version": "chemsmart.program-execution-invocation.v1",
        "node_id": node.node_id,
        "approval_sha256": approval.approval_sha256,
        "program": node.program,
        "engine": node.engine,
        "jobtype": node.jobtype,
        "project_sha256": project_artifact.sha256,
        "input_artifact_id": input_artifact.artifact_id,
        "input_sha256": input_artifact.sha256,
        "scientific_identity_sha256": identity_sha256,
        "environment_receipt_sha256": environment_sha256,
        "resource_sha256": resources.resource_sha256,
        "workspace": approval.workspace,
        "argv": argv_tuple,
        "idempotency_key": idempotency_key,
        "status": "ready",
    }
    return ProgramExecutionInvocationV1(
        **body, invocation_sha256=canonical_sha256(body)
    )


@dataclass(frozen=True)
class ProgramExecutionReceiptV1:
    """Observable engine outcome and independent validation state."""

    schema_version: str
    invocation_sha256: str
    node_id: str
    idempotency_key: str
    execution_state: str
    exit_status: int | None
    engine_complete: bool
    validated: bool
    environment_receipt_sha256: str
    output_artifacts: tuple[TrustedArtifactRefV1, ...]
    validator_receipt_sha256s: tuple[str, ...]
    findings: tuple[str, ...]
    started_at: str
    finished_at: str
    receipt_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.program-execution-receipt.v1":
            raise ContractError("unsupported program execution receipt")
        require_sha256(self.invocation_sha256, "invocation_sha256")
        require_identifier(self.node_id, "node_id")
        require_sha256(self.idempotency_key, "idempotency_key")
        require_sha256(
            self.environment_receipt_sha256,
            "environment_receipt_sha256",
        )
        if self.execution_state not in {
            "not_started",
            "running",
            "engine_complete",
            "failed",
            "ambiguous",
        }:
            raise ContractError("invalid execution state")
        if self.engine_complete and self.exit_status != 0:
            raise ContractError("engine completion requires exit status zero")
        if (
            self.execution_state == "engine_complete"
            and not self.engine_complete
        ):
            raise ContractError(
                "engine_complete state requires engine completion"
            )
        if self.validated and not self.engine_complete:
            raise ContractError("validation requires engine completion")
        if self.validated and not self.validator_receipt_sha256s:
            raise ContractError(
                "validated execution requires validator evidence"
            )
        if self.validated and self.findings:
            raise ContractError(
                "validated execution cannot retain unresolved findings"
            )
        if self.execution_state in {"not_started", "running", "ambiguous"}:
            if self.engine_complete or self.validated:
                raise ContractError("unfinished execution cannot be complete")
        artifact_ids = tuple(
            item.artifact_id for item in self.output_artifacts
        )
        if len(set(artifact_ids)) != len(artifact_ids):
            raise ContractError("execution output artifact IDs must be unique")
        _require_sorted_unique(
            self.validator_receipt_sha256s,
            "validator_receipt_sha256s",
        )
        for digest in self.validator_receipt_sha256s:
            require_sha256(digest, "validator_receipt_sha256")
        _require_nonempty_rows(self.findings, "findings")
        if not self.started_at.strip():
            raise ContractError("execution receipt requires a start time")
        if self.execution_state not in {"running", "not_started"} and not (
            self.finished_at.strip()
        ):
            raise ContractError(
                "terminal execution receipt requires a finish time"
            )
        body = {
            "schema_version": self.schema_version,
            "invocation_sha256": self.invocation_sha256,
            "node_id": self.node_id,
            "idempotency_key": self.idempotency_key,
            "execution_state": self.execution_state,
            "exit_status": self.exit_status,
            "engine_complete": self.engine_complete,
            "validated": self.validated,
            "environment_receipt_sha256": self.environment_receipt_sha256,
            "output_artifacts": self.output_artifacts,
            "validator_receipt_sha256s": self.validator_receipt_sha256s,
            "findings": self.findings,
            "started_at": self.started_at,
            "finished_at": self.finished_at,
        }
        if self.receipt_sha256 != canonical_sha256(body):
            raise ContractError("program execution receipt digest mismatch")


def build_program_execution_receipt(
    invocation: ProgramExecutionInvocationV1,
    *,
    execution_state: str,
    exit_status: int | None,
    engine_complete: bool,
    validated: bool,
    output_artifacts: Sequence[TrustedArtifactRefV1] = (),
    validator_receipt_sha256s: Sequence[str] = (),
    findings: Sequence[str] = (),
    started_at: str,
    finished_at: str = "",
) -> ProgramExecutionReceiptV1:
    """Record an observable outcome without conflating validation."""

    for artifact in output_artifacts:
        _require_current_artifact(artifact, "execution output")
    body = {
        "schema_version": "chemsmart.program-execution-receipt.v1",
        "invocation_sha256": invocation.invocation_sha256,
        "node_id": invocation.node_id,
        "idempotency_key": invocation.idempotency_key,
        "execution_state": str(execution_state),
        "exit_status": exit_status,
        "engine_complete": bool(engine_complete),
        "validated": bool(validated),
        "environment_receipt_sha256": invocation.environment_receipt_sha256,
        "output_artifacts": tuple(output_artifacts),
        "validator_receipt_sha256s": tuple(
            sorted(set(validator_receipt_sha256s))
        ),
        "findings": tuple(str(item).strip() for item in findings),
        "started_at": str(started_at).strip(),
        "finished_at": str(finished_at).strip(),
    }
    return ProgramExecutionReceiptV1(
        **body, receipt_sha256=canonical_sha256(body)
    )


@dataclass(frozen=True)
class OptimizedGeometryHandoffV1:
    """Validated PySCF OPT result converted to one immutable XYZ frame."""

    schema_version: str
    producer_node_id: str
    consumer_node_id: str
    producer_edge_sha256: str
    producer_execution_receipt_sha256: str
    result_artifact_id: str
    result_artifact_sha256: str
    geometry_artifact_id: str
    geometry_artifact_sha256: str
    atom_count: int
    symbols: tuple[str, ...]
    positions_sha256: str
    charge: int
    multiplicity: int
    status: str
    receipt_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.optimized-geometry-handoff.v1":
            raise ContractError("unsupported optimized geometry handoff")
        require_identifier(self.producer_node_id, "producer_node_id")
        require_identifier(self.consumer_node_id, "consumer_node_id")
        for name, digest in (
            ("producer_edge_sha256", self.producer_edge_sha256),
            (
                "producer_execution_receipt_sha256",
                self.producer_execution_receipt_sha256,
            ),
            ("result_artifact_sha256", self.result_artifact_sha256),
            ("geometry_artifact_sha256", self.geometry_artifact_sha256),
            ("positions_sha256", self.positions_sha256),
        ):
            require_sha256(digest, name)
        if not self.result_artifact_id or not self.geometry_artifact_id:
            raise ContractError("handoff artifact IDs must not be empty")
        if self.atom_count < 1 or self.atom_count != len(self.symbols):
            raise ContractError("handoff atom count differs from symbols")
        _require_nonempty_rows(self.symbols, "symbols")
        if self.multiplicity < 1:
            raise ContractError("multiplicity must be positive")
        if self.status != "validated_handoff":
            raise ContractError("optimized geometry handoff must be validated")
        body = {
            "schema_version": self.schema_version,
            "producer_node_id": self.producer_node_id,
            "consumer_node_id": self.consumer_node_id,
            "producer_edge_sha256": self.producer_edge_sha256,
            "producer_execution_receipt_sha256": (
                self.producer_execution_receipt_sha256
            ),
            "result_artifact_id": self.result_artifact_id,
            "result_artifact_sha256": self.result_artifact_sha256,
            "geometry_artifact_id": self.geometry_artifact_id,
            "geometry_artifact_sha256": self.geometry_artifact_sha256,
            "atom_count": self.atom_count,
            "symbols": self.symbols,
            "positions_sha256": self.positions_sha256,
            "charge": self.charge,
            "multiplicity": self.multiplicity,
            "status": self.status,
        }
        if self.receipt_sha256 != canonical_sha256(body):
            raise ContractError("optimized geometry handoff digest mismatch")


def handoff_optimized_pyscf_geometry(
    *,
    producer_receipt: ProgramExecutionReceiptV1,
    result_artifact: TrustedArtifactRefV1,
    producer_edge: ProducerEdgeRuleV1,
    approved_workspace: str | Path,
    geometry_artifact_id: str,
    expected_charge: int,
    expected_multiplicity: int,
) -> tuple[TrustedArtifactRefV1, OptimizedGeometryHandoffV1]:
    """Extract the validated final OPT frame and materialize exact XYZ text."""

    if not producer_receipt.validated:
        raise ContractError("optimized geometry requires a validated producer")
    if producer_receipt.node_id != producer_edge.producer_node_id:
        raise ContractError("execution receipt does not match producer edge")
    if producer_edge.selection_rule != "validated_optimized_geometry":
        raise ContractError("unsupported optimized geometry selection rule")
    source = _require_current_artifact(
        result_artifact, "PySCF result artifact"
    )
    if source.suffix.lower() != ".h5" or result_artifact.kind != "pyscf_hdf5":
        raise ContractError("optimized handoff requires a PySCF HDF5 artifact")
    if not any(
        item.artifact_id == result_artifact.artifact_id
        and item.sha256 == result_artifact.sha256
        for item in producer_receipt.output_artifacts
    ):
        raise ContractError("result artifact is not bound to producer receipt")

    before = file_sha256(source)
    try:
        from chemsmart.io.pyscf.output import read_pyscf_h5

        spec, _provenance, status, results = read_pyscf_h5(source)
    except (ImportError, KeyError, OSError, TypeError, ValueError) as exc:
        raise ContractError("PySCF HDF5 result is not readable") from exc
    after = file_sha256(source)
    if before != result_artifact.sha256 or after != before:
        raise ContractError("PySCF result changed while extracting geometry")
    stages = tuple(str(item) for item in (spec.get("stages") or ()))
    opt_status = (status.get("stages") or {}).get("opt") or {}
    if "opt" not in stages or not bool(opt_status.get("converged", False)):
        raise ContractError("PySCF artifact has no converged optimization")
    if not bool(status.get("normal_termination", False)):
        raise ContractError("PySCF optimization did not terminate normally")

    symbols = tuple(str(item) for item in (spec.get("symbols") or ()))
    positions_value = results.get("positions")
    if positions_value is None:
        raise ContractError("PySCF optimization has no final positions")
    try:
        positions = tuple(
            tuple(float(component) for component in row)
            for row in positions_value
        )
    except (TypeError, ValueError) as exc:
        raise ContractError("PySCF final positions are not numeric") from exc
    if not symbols or len(symbols) != len(positions):
        raise ContractError("PySCF symbols and positions differ in length")
    if any(
        len(row) != 3 or any(not math.isfinite(value) for value in row)
        for row in positions
    ):
        raise ContractError("PySCF final positions must be finite Nx3 values")
    charge = int(spec.get("charge"))
    multiplicity = int(spec.get("multiplicity"))
    if (charge, multiplicity) != (
        int(expected_charge),
        int(expected_multiplicity),
    ):
        raise ContractError("optimized electronic state differs from approval")

    normalized_id = require_identifier(geometry_artifact_id, "artifact_id")
    workspace = _absolute_workspace(approved_workspace)
    target = _target_below(
        workspace,
        "artifacts",
        (
            f"{producer_edge.producer_node_id}--"
            f"{producer_edge.consumer_node_id}.xyz"
        ),
    )
    comment = (
        "ChemSmart validated PySCF OPT handoff; "
        f"charge={charge}; multiplicity={multiplicity}; "
        f"source_sha256={result_artifact.sha256}"
    )
    lines = [str(len(symbols)), comment]
    for symbol, (x, y, z) in zip(symbols, positions):
        lines.append(f"{symbol:<3} {x:.17g} {y:.17g} {z:.17g}")
    payload = ("\n".join(lines) + "\n").encode("utf-8")
    _write_exact_once(target, payload)
    geometry_artifact = TrustedArtifactRefV1(
        artifact_id=normalized_id,
        kind="geometry_xyz",
        sha256=file_sha256(target),
        size_bytes=target.stat().st_size,
        path=str(target),
        cli_value=str(target),
    )
    positions_sha256 = canonical_sha256(
        {"symbols": symbols, "positions_angstrom": positions}
    )
    body = {
        "schema_version": "chemsmart.optimized-geometry-handoff.v1",
        "producer_node_id": producer_edge.producer_node_id,
        "consumer_node_id": producer_edge.consumer_node_id,
        "producer_edge_sha256": producer_edge.edge_sha256,
        "producer_execution_receipt_sha256": producer_receipt.receipt_sha256,
        "result_artifact_id": result_artifact.artifact_id,
        "result_artifact_sha256": result_artifact.sha256,
        "geometry_artifact_id": geometry_artifact.artifact_id,
        "geometry_artifact_sha256": geometry_artifact.sha256,
        "atom_count": len(symbols),
        "symbols": symbols,
        "positions_sha256": positions_sha256,
        "charge": charge,
        "multiplicity": multiplicity,
        "status": "validated_handoff",
    }
    return geometry_artifact, OptimizedGeometryHandoffV1(
        **body, receipt_sha256=canonical_sha256(body)
    )


__all__ = [
    "ApprovedNodeBindingV1",
    "ExecutionResourceSpecV1",
    "OptimizedGeometryHandoffV1",
    "ProducerEdgeRuleV1",
    "ProgramConformanceProbeV1",
    "ProgramExecutionInvocationV1",
    "ProgramExecutionReceiptV1",
    "ProjectArtifactPromotionV1",
    "ScientificDecisionRecordV1",
    "WorkflowExecutionApprovalV1",
    "bind_project_promotion_validation",
    "build_execution_resource_spec",
    "build_locked_pyscf_sp_opt_hess_approval",
    "build_producer_edge_rule",
    "build_program_conformance_probe",
    "build_program_execution_invocation",
    "build_program_execution_receipt",
    "build_scientific_decision_record",
    "build_workflow_execution_approval",
    "handoff_optimized_pyscf_geometry",
    "promote_project_candidate",
    "workflow_execution_approval_json",
]
