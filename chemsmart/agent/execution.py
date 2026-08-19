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

import json
import math
import os
import sys
from dataclasses import dataclass
from pathlib import Path, PureWindowsPath
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

from chemsmart.agent._contracts import (
    AuxiliaryArtifactBindingV1,
    ContractError,
    TrustedArtifactRefV1,
    canonical_data,
    canonical_json,
    canonical_sha256,
    file_sha256,
    require_identifier,
    require_auxiliary_artifact_bindings,
    require_sha256,
)
from chemsmart.agent.scientific_toolchain import ScientificToolchainPlanV1
from chemsmart.agent.projects import (
    ProjectRenderReceiptV1,
    ProjectValidationReceiptV1,
)
from chemsmart.agent.workflows import (
    PREVIEW_RESOURCE_SHA256,
    MaterializedWorkflowV1,
    ScientificWorkflowEdgeV2,
    ScientificWorkflowPlanV2,
    StationaryPointValidationPolicyV1,
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


@dataclass(frozen=True)
class MolecularCompositionReceiptV1:
    """Host-owned lineage of one two-fragment molecular arrangement.

    The host owns the placement mathematics; the model owns the choice of
    fragments, contact atoms, and distance. Composition never infers an
    electronic state: charge and multiplicity of the arrangement are bound
    separately and explicitly, and the stage that consumes the composed
    geometry is a new workflow for human review.
    """

    schema_version: str
    composed_artifact_id: str
    composed_artifact_sha256: str
    fragment_a_artifact_id: str
    fragment_a_sha256: str
    fragment_a_identity_sha256: str
    fragment_b_artifact_id: str
    fragment_b_sha256: str
    fragment_b_identity_sha256: str
    placement: dict[str, Any]
    achieved_contact_distance_angstrom: float
    min_interfragment_distance_angstrom: float
    atom_count: int
    formula: str
    atom_order_note: str
    status: str
    receipt_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.molecular-composition.v1":
            raise ContractError("unsupported molecular composition receipt")
        for name, digest in (
            ("composed_artifact_sha256", self.composed_artifact_sha256),
            ("fragment_a_sha256", self.fragment_a_sha256),
            ("fragment_a_identity_sha256", self.fragment_a_identity_sha256),
            ("fragment_b_sha256", self.fragment_b_sha256),
            ("fragment_b_identity_sha256", self.fragment_b_identity_sha256),
        ):
            require_sha256(digest, name)
        if self.atom_count < 2:
            raise ContractError("a composition carries at least two atoms")
        if self.status != "composed":
            raise ContractError("molecular composition must be composed")
        object.__setattr__(
            self, "placement", canonical_data(dict(self.placement))
        )
        body = {
            key: value
            for key, value in self.__dict__.items()
            if key != "receipt_sha256"
        }
        if self.receipt_sha256 != canonical_sha256(body):
            raise ContractError("molecular composition digest mismatch")


def compose_trusted_molecular_arrangement(
    *,
    approved_workspace: str | Path,
    composed_artifact_id: str,
    fragment_a: TrustedArtifactRefV1,
    fragment_a_identity_sha256: str,
    fragment_b: TrustedArtifactRefV1,
    fragment_b_identity_sha256: str,
    fragment_a_atom: int,
    fragment_b_atom: int,
    distance_angstrom: float,
) -> tuple[TrustedArtifactRefV1, MolecularCompositionReceiptV1]:
    """Place fragment B against fragment A at one explicit atomic contact.

    Deterministic host-owned placement: the named contact pair is held at
    the requested distance while every other interfragment pair stays
    outside covalent-radii-plus-buffer, and the remaining freedom maximises
    separation (the iterate module's SLSQP machinery, reused with the
    model's contact distance instead of a covalent bond length). Fragment A
    keeps its coordinates; the composed file lists fragment A's atoms first.
    """

    import numpy as np

    from chemsmart.io.molecules.structure import Molecule
    from chemsmart.jobs.iterate.iterate import (
        DEFAULT_BUFFER,
        IterateAnalyzer,
    )
    from chemsmart.utils.periodictable import covalent_radii

    if fragment_a.kind != "geometry_xyz" or fragment_b.kind != "geometry_xyz":
        raise ContractError(
            "composition consumes two geometry_xyz artifacts"
        )
    distance = float(distance_angstrom)
    if not 0.5 <= distance <= 10.0:
        raise ContractError(
            "contact distance must lie in [0.5, 10.0] angstrom; got "
            f"{distance}"
        )
    path_a = _require_current_artifact(fragment_a, "composition fragment A")
    path_b = _require_current_artifact(fragment_b, "composition fragment B")
    molecule_a = Molecule.from_filepath(str(path_a))
    molecule_b = Molecule.from_filepath(str(path_b))
    count_a = len(molecule_a.chemical_symbols)
    count_b = len(molecule_b.chemical_symbols)
    if not 1 <= int(fragment_a_atom) <= count_a:
        raise ContractError(
            f"fragment_a_atom must be 1..{count_a}; got {fragment_a_atom}"
        )
    if not 1 <= int(fragment_b_atom) <= count_b:
        raise ContractError(
            f"fragment_b_atom must be 1..{count_b}; got {fragment_b_atom}"
        )
    link_a = int(fragment_a_atom) - 1
    link_b = int(fragment_b_atom) - 1

    array_a = IterateAnalyzer._molecule_to_array(molecule_a)
    array_b = IterateAnalyzer._molecule_to_array(molecule_b)
    relative = IterateAnalyzer._calc_relative_coords(array_b, link_b)
    elements_a = array_a[:, 0].astype(int)
    elements_b = array_b[:, 0].astype(int)
    radii_a = np.array([covalent_radii[z] for z in elements_a])
    radii_b = np.array([covalent_radii[z] for z in elements_b])
    min_dist_matrix = (
        radii_b[:, np.newaxis] + radii_a[np.newaxis, :] + DEFAULT_BUFFER
    )
    ineq_mask = np.ones((count_b, count_a), dtype=bool)
    ineq_mask[link_b, link_a] = False
    placed = IterateAnalyzer._optimize_lagrange(
        array_b,
        array_a[:, 1:4],
        96,
        6,
        link_a,
        relative[:, 1:4],
        distance,
        min_dist_matrix,
        link_b,
        ineq_mask,
    )
    if placed is None:
        raise ContractError(
            "no clash-free arrangement satisfies the requested contact: "
            "raise distance_angstrom or choose different contact atoms"
        )
    positions_a = array_a[:, 1:4]
    positions_b = placed[:, 1:4]
    achieved = float(
        np.linalg.norm(positions_b[link_b] - positions_a[link_a])
    )
    pair_distances = np.linalg.norm(
        positions_b[:, np.newaxis, :] - positions_a[np.newaxis, :, :],
        axis=2,
    )
    minimum_separation = float(pair_distances.min())
    symbols = tuple(molecule_a.chemical_symbols) + tuple(
        molecule_b.chemical_symbols
    )
    composed = Molecule(
        symbols=list(symbols),
        positions=np.vstack((positions_a, positions_b)),
    )
    lines = [str(len(symbols)), (
        "ChemSmart composed arrangement; fragment A atoms first; "
        f"contact {fragment_a_atom}(A)-{fragment_b_atom}(B) at "
        f"{achieved:.4f} angstrom; electronic state deliberately unbound"
    )]
    for symbol, position in zip(symbols, composed.positions):
        lines.append(
            f"{symbol:<3} {position[0]:.10f} {position[1]:.10f} "
            f"{position[2]:.10f}"
        )
    payload = ("\n".join(lines) + "\n").encode("utf-8")
    target = _target_below(
        _absolute_workspace(approved_workspace),
        "artifacts",
        f"{require_identifier(composed_artifact_id, 'artifact_id')}.xyz",
    )
    _write_exact_once(target, payload)
    artifact = TrustedArtifactRefV1(
        artifact_id=composed_artifact_id,
        kind="geometry_xyz",
        sha256=file_sha256(target),
        size_bytes=target.stat().st_size,
        path=str(target),
        cli_value=str(target),
    )
    placement = {
        "schema_version": "chemsmart.placement-spec.v1",
        "mode": "contact",
        "fragment_a_atom": int(fragment_a_atom),
        "fragment_b_atom": int(fragment_b_atom),
        "distance_angstrom": distance,
        "buffer_angstrom": DEFAULT_BUFFER,
        "sphere_direction_samples": 96,
        "axial_rotation_samples": 6,
    }
    body = {
        "schema_version": "chemsmart.molecular-composition.v1",
        "composed_artifact_id": artifact.artifact_id,
        "composed_artifact_sha256": artifact.sha256,
        "fragment_a_artifact_id": fragment_a.artifact_id,
        "fragment_a_sha256": fragment_a.sha256,
        "fragment_a_identity_sha256": fragment_a_identity_sha256,
        "fragment_b_artifact_id": fragment_b.artifact_id,
        "fragment_b_sha256": fragment_b.sha256,
        "fragment_b_identity_sha256": fragment_b_identity_sha256,
        "placement": placement,
        "achieved_contact_distance_angstrom": round(achieved, 6),
        "min_interfragment_distance_angstrom": round(minimum_separation, 6),
        "atom_count": len(symbols),
        "formula": composed.get_chemical_formula(),
        "atom_order_note": "fragment A atoms first, then fragment B",
        "status": "composed",
    }
    return artifact, MolecularCompositionReceiptV1(
        **body, receipt_sha256=canonical_sha256(body)
    )


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
        "registry_sha256": require_sha256(registry_sha256, "registry_sha256"),
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
        "preview_receipt_sha256s": tuple(sorted(set(preview_receipt_sha256s))),
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
    auxiliary_input_bindings: tuple[AuxiliaryArtifactBindingV1, ...] = ()
    #: Which internal coordinates this node drives or holds. The executor
    #: rebuilds each invocation from this binding alone, so a coordinate that
    #: lives only in the planning draft is not there to rebuild from: two scans
    #: of the same molecule over different ranges then synthesise to the same
    #: coordinate-free argv, which correctly fails the comparison against what
    #: a human approved. It is also the honest place for it -- the range that
    #: was displayed is part of what was approved.
    #:
    #: Optional and omitted from the canonical body when absent, so every
    #: already-recorded approval keeps its digest and stays replayable.
    internal_coordinates: Mapping[str, Any] | None = None

    def __post_init__(self) -> None:
        for name, value in (
            ("node_id", self.node_id),
            ("program", self.program),
            ("engine", self.engine),
            ("jobtype", self.jobtype),
        ):
            require_identifier(value, name)
        require_sha256(self.project_artifact_sha256, "project_artifact_sha256")
        require_sha256(self.settings_sha256, "settings_sha256")
        require_auxiliary_artifact_bindings(self.auxiliary_input_bindings)
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


def _approved_node_binding_body(
    binding: ApprovedNodeBindingV1,
) -> dict[str, Any]:
    """Canonical approval projection with additive auxiliary inputs."""

    body = {
        "node_id": binding.node_id,
        "program": binding.program,
        "engine": binding.engine,
        "jobtype": binding.jobtype,
        "project_artifact_sha256": binding.project_artifact_sha256,
        "settings_sha256": binding.settings_sha256,
        "charge": binding.charge,
        "multiplicity": binding.multiplicity,
        "input_mode": binding.input_mode,
        "initial_artifact_id": binding.initial_artifact_id,
        "initial_artifact_sha256": binding.initial_artifact_sha256,
        "scientific_identity_sha256": binding.scientific_identity_sha256,
        "producer_edge_sha256": binding.producer_edge_sha256,
    }
    if binding.auxiliary_input_bindings:
        body["auxiliary_input_bindings"] = binding.auxiliary_input_bindings
    if binding.internal_coordinates:
        body["internal_coordinates"] = canonical_data(
            binding.internal_coordinates
        )
    return body


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
        for edge in self.producer_edges:
            if (
                edge.producer_node_id not in node_ids
                or edge.consumer_node_id not in node_ids
            ):
                raise ContractError(
                    "approval producer edge endpoints must both be approved"
                )
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
            "node_bindings": tuple(
                _approved_node_binding_body(item)
                for item in self.node_bindings
            ),
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
    normalized_node_bindings = tuple(node_bindings)
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
        "node_bindings": normalized_node_bindings,
        "producer_edges": tuple(producer_edges),
        "status": "approved",
    }
    digest_body = {
        **body,
        "node_bindings": tuple(
            _approved_node_binding_body(item)
            for item in normalized_node_bindings
        ),
    }
    return WorkflowExecutionApprovalV1(
        **body, approval_sha256=canonical_sha256(digest_body)
    )


def workflow_execution_approval_json(
    approval: WorkflowExecutionApprovalV1,
) -> str:
    """Return the exact private JSON envelope accepted by ``agent run``."""

    return canonical_json({"workflow_approval": approval}) + "\n"


@dataclass(frozen=True)
class WorkflowApprovalRequestV1:
    """Exact, inert workflow request presented for human review."""

    schema_version: str
    request_id: str
    workflow_id: str
    workflow_sha256: str
    task_spec_sha256: str
    workspace: str
    resource_sha256: str
    node_bindings: tuple[ApprovedNodeBindingV1, ...]
    producer_edges: tuple[ProducerEdgeRuleV1, ...]
    status: str
    request_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.workflow-approval-request.v1":
            raise ContractError("unsupported workflow approval request schema")
        if self.status != "unapproved":
            raise ContractError(
                "a workflow approval request is never approved; approve it "
                "explicitly to obtain a workflow execution approval"
            )
        require_identifier(self.request_id, "request_id")
        require_identifier(self.workflow_id, "workflow_id")
        require_sha256(self.workflow_sha256, "workflow_sha256")
        require_sha256(self.task_spec_sha256, "task_spec_sha256")
        require_sha256(self.resource_sha256, "resource_sha256")
        if not Path(self.workspace).is_absolute():
            raise ContractError("approval request workspace must be absolute")
        object.__setattr__(self, "node_bindings", tuple(self.node_bindings))
        object.__setattr__(self, "producer_edges", tuple(self.producer_edges))
        if not self.node_bindings:
            raise ContractError("approval request requires at least one node")
        node_ids = tuple(item.node_id for item in self.node_bindings)
        if len(set(node_ids)) != len(node_ids):
            raise ContractError("approval request node IDs must be unique")
        if self.request_sha256 != canonical_sha256(self._body()):
            raise ContractError("approval request digest mismatch")

    def _body(self) -> dict[str, Any]:
        return {
            key: value
            for key, value in self.__dict__.items()
            if key != "request_sha256"
        }


def build_workflow_approval_request(
    *,
    request_id: str,
    workflow_id: str,
    workflow_sha256: str,
    task_spec_sha256: str,
    workspace: str | Path,
    resources: ExecutionResourceSpecV1,
    node_bindings: Sequence[ApprovedNodeBindingV1],
    producer_edges: Sequence[ProducerEdgeRuleV1] = (),
) -> WorkflowApprovalRequestV1:
    """Assemble the reviewable body of a plan without approving anything."""

    body = {
        "schema_version": "chemsmart.workflow-approval-request.v1",
        "request_id": require_identifier(request_id, "request_id"),
        "workflow_id": require_identifier(workflow_id, "workflow_id"),
        "workflow_sha256": require_sha256(workflow_sha256, "workflow_sha256"),
        "task_spec_sha256": require_sha256(
            task_spec_sha256, "task_spec_sha256"
        ),
        "workspace": str(_absolute_workspace(workspace)),
        "resource_sha256": resources.resource_sha256,
        "node_bindings": tuple(node_bindings),
        "producer_edges": tuple(producer_edges),
        "status": "unapproved",
    }
    return WorkflowApprovalRequestV1(
        **body, request_sha256=canonical_sha256(body)
    )


def approve_workflow_request(
    request: WorkflowApprovalRequestV1,
    *,
    approval_id: str,
    approved_request_sha256: str,
    resources: ExecutionResourceSpecV1,
) -> WorkflowExecutionApprovalV1:
    """Convert a reviewed request into an approval, one explicit act.

    ``approved_request_sha256`` is the digest the reviewer saw.  Requiring it
    means approval cannot be granted to a body that changed after review, and
    that the caller had to look at something to supply it.  ``resources`` are
    the host's current locked allocation; if they have moved since the plan,
    the request is stale and must be re-reviewed rather than re-signed.
    """

    if approved_request_sha256 != request.request_sha256:
        raise ContractError(
            "the approved request digest does not match the request; "
            f"reviewed {approved_request_sha256[:16]}..., holding "
            f"{request.request_sha256[:16]}..."
        )
    if resources.resource_sha256 != request.resource_sha256:
        raise ContractError(
            "the host's execution resources have changed since this plan was "
            "made; re-plan and review again rather than approving a request "
            "whose compute allocation no longer exists"
        )
    return build_workflow_execution_approval(
        approval_id=approval_id,
        workflow_id=request.workflow_id,
        workflow_sha256=request.workflow_sha256,
        task_spec_sha256=request.task_spec_sha256,
        approved_workspace=request.workspace,
        resources=resources,
        node_bindings=request.node_bindings,
        producer_edges=request.producer_edges,
    )


def workflow_approval_request_json(
    request: WorkflowApprovalRequestV1,
) -> str:
    """Return the reviewable envelope, deliberately not the run envelope."""

    return canonical_json({"workflow_approval_request": request}) + "\n"


def argv_shape(argv: Sequence[str]) -> tuple[str, ...]:
    """Replace absolute-path tokens with a placeholder.

    Paths in a compiled command are where files happened to sit during one
    session: a plan session writes its project YAML and server profile under a
    private run directory whose name carries a timestamp and a nonce, so those
    tokens can never recur.  What the command *does* -- the subcommand, the
    flags, their order and their literal values -- is identical wherever it
    runs, and the file contents are pinned separately by digest.
    """

    return tuple(
        "<path>" if str(token).startswith("/") else str(token)
        for token in argv
    )


def invocation_identity_sha256(
    *,
    program: str,
    engine: str,
    jobtype: str,
    project_sha256: str,
    input_sha256: str,
    scientific_identity_sha256: str,
    argv: Sequence[str],
    auxiliary_input_bindings: tuple[AuxiliaryArtifactBindingV1, ...] = (),
) -> str:
    """Identify a compiled command by what it computes, not where it ran.

    An approval that pinned ``invocation_sha256`` could not be executed by any
    later process.  That digest covers argv, and argv carries absolute paths
    into the planning session's own timestamped run directory; it also covers
    a different record type on each side of the check, since a preview pins a
    ``CanonicalCommandInvocationV1`` over the preview argv while execution
    builds a ``ProgramExecutionInvocationV1`` over a host-rewritten real argv.
    Two schemas, two argvs, one field name.

    This identity is stable across directories and machines, and still fails
    on a changed program, engine, job type, project bytes, input bytes,
    molecular identity, or command shape -- which is everything the check was
    defending.
    """

    require_auxiliary_artifact_bindings(auxiliary_input_bindings)
    body = {
        "schema_version": "chemsmart.invocation-identity.v1",
        "program": str(program),
        "engine": str(engine),
        "jobtype": str(jobtype),
        "project_sha256": str(project_sha256),
        "input_sha256": str(input_sha256),
        "scientific_identity_sha256": str(scientific_identity_sha256),
        "argv_shape": argv_shape(argv),
    }
    if auxiliary_input_bindings:
        body["auxiliary_input_bindings"] = auxiliary_input_bindings
    return canonical_sha256(body)


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
    #: Identity of the environment this invocation observed.  Carried
    #: separately from ``environment_receipt_sha256`` because a receipt digest
    #: folds in a capability receipt and therefore differs between the session
    #: that planned and the session that executes, on one unchanged machine.
    #: Admission compares this; the receipt digest still names the exact
    #: observation.  Defaulted so records built before the field existed keep
    #: their digests.
    environment_identity_sha256: str = ""
    #: Path-independent identity of the compiled command this execution
    #: realises.  See ``invocation_identity_sha256``.
    invocation_identity_sha256: str = ""
    auxiliary_input_bindings: tuple[AuxiliaryArtifactBindingV1, ...] = ()

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
        require_auxiliary_artifact_bindings(self.auxiliary_input_bindings)
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
        if self.environment_identity_sha256:
            require_sha256(
                self.environment_identity_sha256,
                "environment_identity_sha256",
            )
            body["environment_identity_sha256"] = (
                self.environment_identity_sha256
            )
        if self.invocation_identity_sha256:
            require_sha256(
                self.invocation_identity_sha256,
                "invocation_identity_sha256",
            )
            body["invocation_identity_sha256"] = (
                self.invocation_identity_sha256
            )
        if self.auxiliary_input_bindings:
            body["auxiliary_input_bindings"] = self.auxiliary_input_bindings
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
    environment_identity_sha256: str = "",
    invocation_identity_sha256: str = "",
    auxiliary_input_artifacts: Mapping[str, TrustedArtifactRefV1]
    | None = None,
    auxiliary_handoffs: Mapping[str, ORCAHessianHandoffV1] | None = None,
) -> ProgramExecutionInvocationV1:
    """Resolve a ``node_id`` to an exact invocation using host bindings.

    ``environment_identity_sha256`` is what admission compares against the
    frozen approval.  The receipt digest still records which exact observation
    produced the binding, but it cannot be compared across sessions: it folds
    in a capability receipt that changes with the active tool overlay.
    """

    node = approval.node(node_id)
    _require_current_artifact(project_artifact, "project artifact")
    _require_current_artifact(input_artifact, "input artifact")
    auxiliary_input_artifacts = dict(auxiliary_input_artifacts or {})
    auxiliary_handoffs = dict(auxiliary_handoffs or {})
    auxiliary_input_bindings = tuple(
        AuxiliaryArtifactBindingV1(
            parameter_name=parameter_name,
            artifact_id=artifact.artifact_id,
            artifact_sha256=artifact.sha256,
        )
        for parameter_name, artifact in sorted(
            auxiliary_input_artifacts.items()
        )
    )
    require_auxiliary_artifact_bindings(auxiliary_input_bindings)
    for parameter_name, artifact in sorted(auxiliary_input_artifacts.items()):
        _require_current_artifact(
            artifact, f"auxiliary input {parameter_name}"
        )
    if auxiliary_input_bindings != node.auxiliary_input_bindings:
        # A future ORCA Hessian has no bytes at review time.  It is authorized
        # by the typed producer edge, then checked here against the semantic
        # handoff selected from the validated TS receipt.
        hessian_roles = {
            "hess_filename": "validated_final_orca_ts_hessian",
            "inhess_filename": "validated_producer_orca_hessian",
        }
        aux_keys = tuple(auxiliary_input_artifacts)
        future_auxiliary_ok = bool(
            node.input_mode == "producer"
            and len(aux_keys) == 1
            and aux_keys[0] in hessian_roles
            and tuple(auxiliary_handoffs) == aux_keys
        )
        if future_auxiliary_ok:
            role = aux_keys[0]
            artifact = auxiliary_input_artifacts[role]
            auxiliary_handoff = auxiliary_handoffs[role]
            matching_edges = tuple(
                edge
                for edge in approval.producer_edges
                if edge.consumer_node_id == node.node_id
                and edge.artifact_kind == "orca_hessian"
                and edge.selection_rule == hessian_roles[role]
            )
            future_auxiliary_ok = bool(
                len(matching_edges) == 1
                and auxiliary_handoff.producer_edge_sha256
                == matching_edges[0].edge_sha256
                and auxiliary_handoff.consumer_node_id == node.node_id
                and auxiliary_handoff.selected_artifact_id
                == artifact.artifact_id
                and auxiliary_handoff.selected_artifact_sha256
                == artifact.sha256
                and auxiliary_handoff.consumer_state
                == (node.charge, node.multiplicity)
            )
        if not future_auxiliary_ok:
            raise ContractError("auxiliary inputs differ from workflow approval")
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
        if handoff.consumer_state != (
            node.charge,
            node.multiplicity,
        ):
            raise ContractError(
                "handoff electronic state differs from approval"
            )
        # A geometry edge transfers coordinates, not electronic state.  Bind
        # the downstream target state to these exact bytes before an invocation
        # can claim the corresponding scientific identity.
        from chemsmart.agent.commands import build_scientific_identity_binding

        expected_identity = build_scientific_identity_binding(
            task_spec_sha256=approval.task_spec_sha256,
            geometry_artifact=input_artifact,
            charge=node.charge,
            multiplicity=node.multiplicity,
        )
        if identity_sha256 != expected_identity.binding_sha256:
            raise ContractError(
                "handoff scientific identity differs from target state and geometry"
            )
    environment_sha256 = require_sha256(
        environment_receipt_sha256, "environment_receipt_sha256"
    )
    argv_tuple = tuple(str(item) for item in argv)
    idempotency_body = {
        "approval_sha256": approval.approval_sha256,
        "node_id": node.node_id,
        "project_sha256": project_artifact.sha256,
        "input_sha256": input_artifact.sha256,
        "scientific_identity_sha256": identity_sha256,
        "environment_receipt_sha256": environment_sha256,
        "resource_sha256": resources.resource_sha256,
        "argv": argv_tuple,
    }
    if auxiliary_input_bindings:
        idempotency_body["auxiliary_input_bindings"] = auxiliary_input_bindings
    idempotency_key = canonical_sha256(idempotency_body)
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
    if environment_identity_sha256:
        body["environment_identity_sha256"] = require_sha256(
            environment_identity_sha256, "environment_identity_sha256"
        )
    if invocation_identity_sha256:
        body["invocation_identity_sha256"] = require_sha256(
            invocation_identity_sha256, "invocation_identity_sha256"
        )
    if auxiliary_input_bindings:
        body["auxiliary_input_bindings"] = auxiliary_input_bindings
    return ProgramExecutionInvocationV1(
        **body, invocation_sha256=canonical_sha256(body)
    )


@dataclass(frozen=True)
class ProgramResultValidationReceiptV1:
    """Canonical, resolvable evidence for one post-execution validation.

    This receipt is deliberately separate from the execution envelope.  It
    records the complete deterministic observation used to decide scientific
    validity, including the two different environment evidence layers: the
    capability receipt approved before launch and the per-run receipt emitted
    by the program runner.  Those source digests are never compared as if they
    represented the same object.
    """

    schema_version: str
    validator_id: str
    validator_schema_version: str
    validator_version: str
    invocation_sha256: str
    node_id: str
    program: str
    engine: str
    jobtype: str
    input_artifact_sha256: str
    project_artifact_sha256: str
    capability_environment_receipt_sha256: str
    run_environment_receipt_sha256: str
    environment_validation_sha256: str
    stationary_point_policy_sha256: str
    output_artifacts: tuple[TrustedArtifactRefV1, ...]
    observations: Mapping[str, Any]
    findings: tuple[str, ...]
    state: str
    receipt_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != (
            "chemsmart.program-result-validation-receipt.v1"
        ):
            raise ContractError(
                "unsupported program result validation receipt"
            )
        for name, value in (
            ("validator_id", self.validator_id),
            ("node_id", self.node_id),
            ("program", self.program),
            ("engine", self.engine),
            ("jobtype", self.jobtype),
        ):
            require_identifier(value, name)
        if not self.validator_schema_version.strip():
            raise ContractError("validator schema version must not be empty")
        if not self.validator_version.strip():
            raise ContractError("validator version must not be empty")
        for name, digest in (
            ("invocation_sha256", self.invocation_sha256),
            ("input_artifact_sha256", self.input_artifact_sha256),
            ("project_artifact_sha256", self.project_artifact_sha256),
            (
                "capability_environment_receipt_sha256",
                self.capability_environment_receipt_sha256,
            ),
        ):
            require_sha256(digest, name)
        for name, digest in (
            (
                "run_environment_receipt_sha256",
                self.run_environment_receipt_sha256,
            ),
            (
                "environment_validation_sha256",
                self.environment_validation_sha256,
            ),
            (
                "stationary_point_policy_sha256",
                self.stationary_point_policy_sha256,
            ),
        ):
            if digest:
                require_sha256(digest, name)
        artifact_ids = tuple(
            item.artifact_id for item in self.output_artifacts
        )
        if len(artifact_ids) != len(set(artifact_ids)):
            raise ContractError(
                "validation output artifact IDs must be unique"
            )
        if self.findings != tuple(sorted(set(self.findings))):
            raise ContractError(
                "validation findings must be sorted and unique"
            )
        _require_nonempty_rows(self.findings, "findings")
        if self.state not in {"valid", "invalid"}:
            raise ContractError(
                "result validation state must be valid or invalid"
            )
        derived_findings = _typed_result_validation_findings(
            canonical_data(dict(self.observations))
        )
        expected_state = (
            "invalid" if self.findings or derived_findings else "valid"
        )
        if self.state != expected_state:
            raise ContractError(
                "result validation state differs from deterministic findings"
            )
        body = {
            "schema_version": self.schema_version,
            "validator_id": self.validator_id,
            "validator_schema_version": self.validator_schema_version,
            "validator_version": self.validator_version,
            "invocation_sha256": self.invocation_sha256,
            "node_id": self.node_id,
            "program": self.program,
            "engine": self.engine,
            "jobtype": self.jobtype,
            "input_artifact_sha256": self.input_artifact_sha256,
            "project_artifact_sha256": self.project_artifact_sha256,
            "capability_environment_receipt_sha256": (
                self.capability_environment_receipt_sha256
            ),
            "run_environment_receipt_sha256": (
                self.run_environment_receipt_sha256
            ),
            "environment_validation_sha256": (
                self.environment_validation_sha256
            ),
            "stationary_point_policy_sha256": (
                self.stationary_point_policy_sha256
            ),
            "output_artifacts": self.output_artifacts,
            "observations": canonical_data(dict(self.observations)),
            "findings": self.findings,
            "state": self.state,
        }
        if self.receipt_sha256 != canonical_sha256(body):
            raise ContractError("program result validation digest mismatch")


def _typed_result_validation_findings(
    observations: Mapping[str, Any],
) -> tuple[str, ...]:
    """Derive red findings from an embedded typed validation observation.

    Program-specific validators commonly expose their typed result under
    ``result_validation``.  A few compatibility callers still supply that
    object at the observation root.  In both forms an explicit non-green
    state or an embedded finding is authoritative; a caller cannot erase it
    by passing an empty top-level ``findings`` sequence.
    """

    if "result_validation" in observations:
        result_validation = observations["result_validation"]
        if not isinstance(result_validation, Mapping):
            raise ContractError(
                "result_validation observation must be a typed mapping"
            )
    elif "state" in observations:
        result_validation = observations
    else:
        # Some program validators (currently xTB) report a typed observation
        # and their normalized findings separately.  Do not invent a state
        # for those records; the caller findings remain authoritative until
        # that validator gains a nested result-validation contract.
        return ()

    state = str(result_validation.get("state") or "").strip().lower()
    downstream_hessian_classification = bool(
        observations.get("runner_validation_delegation")
        == "downstream_scientific_analysis"
        and state == "unclassified"
    )
    derived: list[str] = []
    if state not in {"valid", "validated"} and not (
        downstream_hessian_classification
    ):
        derived.append("result_validation.state_not_validated")

    embedded = result_validation.get("findings", ())
    if embedded is None:
        embedded = ()
    if isinstance(embedded, (str, bytes, bytearray)) or not isinstance(
        embedded, Sequence
    ):
        raise ContractError(
            "result_validation findings must be a canonical sequence"
        )
    for item in embedded:
        if isinstance(item, Mapping):
            rule_id = str(item.get("rule_id") or "").strip()
            derived.append(
                rule_id or "result_validation.finding_missing_rule_id"
            )
        else:
            rule_id = str(getattr(item, "rule_id", item)).strip()
            derived.append(
                rule_id or "result_validation.finding_missing_rule_id"
            )
    return tuple(sorted(set(derived)))


def build_program_result_validation_receipt(
    *,
    invocation: ProgramExecutionInvocationV1,
    validator_id: str,
    validator_schema_version: str,
    validator_version: str,
    input_artifact: TrustedArtifactRefV1,
    project_artifact: TrustedArtifactRefV1,
    capability_environment_receipt_sha256: str,
    output_artifacts: Sequence[TrustedArtifactRefV1],
    observations: Mapping[str, Any],
    findings: Sequence[str],
    run_environment_receipt_sha256: str = "",
    environment_validation_sha256: str = "",
    stationary_point_policy_sha256: str = "",
) -> ProgramResultValidationReceiptV1:
    """Build exact post-run validation evidence after rechecking artifacts."""

    _require_current_artifact(input_artifact, "validation input artifact")
    _require_current_artifact(project_artifact, "validation project artifact")
    if input_artifact.sha256 != invocation.input_sha256:
        raise ContractError(
            "validation input differs from execution invocation"
        )
    if project_artifact.sha256 != invocation.project_sha256:
        raise ContractError(
            "validation project differs from execution invocation"
        )
    for artifact in output_artifacts:
        _require_current_artifact(artifact, "validation output artifact")
    environment_sha256 = require_sha256(
        capability_environment_receipt_sha256,
        "capability_environment_receipt_sha256",
    )
    if environment_sha256 != invocation.environment_receipt_sha256:
        raise ContractError(
            "validation environment differs from execution invocation"
        )
    canonical_observations = canonical_data(dict(observations))
    derived_findings = _typed_result_validation_findings(
        canonical_observations
    )
    normalized_findings = tuple(
        sorted(
            {
                *(str(item).strip() for item in findings),
                *derived_findings,
            }
        )
    )
    body = {
        "schema_version": "chemsmart.program-result-validation-receipt.v1",
        "validator_id": require_identifier(validator_id, "validator_id"),
        "validator_schema_version": str(validator_schema_version).strip(),
        "validator_version": str(validator_version).strip(),
        "invocation_sha256": invocation.invocation_sha256,
        "node_id": invocation.node_id,
        "program": invocation.program,
        "engine": invocation.engine,
        "jobtype": invocation.jobtype,
        "input_artifact_sha256": input_artifact.sha256,
        "project_artifact_sha256": project_artifact.sha256,
        "capability_environment_receipt_sha256": environment_sha256,
        "run_environment_receipt_sha256": str(
            run_environment_receipt_sha256
        ).strip(),
        "environment_validation_sha256": str(
            environment_validation_sha256
        ).strip(),
        "stationary_point_policy_sha256": str(
            stationary_point_policy_sha256
        ).strip(),
        "output_artifacts": tuple(output_artifacts),
        "observations": canonical_observations,
        "findings": normalized_findings,
        "state": "valid" if not normalized_findings else "invalid",
    }
    return ProgramResultValidationReceiptV1(
        **body, receipt_sha256=canonical_sha256(body)
    )


@dataclass(frozen=True)
class ProgramExecutionReceiptV1:
    """Wrapper, child-engine, and independent validation observations.

    ``wrapper_exit_status`` belongs to the outer ``chemsmart`` CLI process.
    It can be non-zero after the PySCF child has completed because the wrapper
    deliberately raises when deterministic scientific validation is red.
    ``child_exit_status`` and ``engine_complete`` therefore come from the
    program-specific, digest-valid run receipt rather than the wrapper code.
    """

    schema_version: str
    invocation_sha256: str
    node_id: str
    idempotency_key: str
    execution_state: str
    wrapper_exit_status: int | None
    child_exit_status: int | None
    engine_complete: bool
    validated: bool
    engine_receipt_sha256: str
    environment_receipt_sha256: str
    result_validation_receipt_sha256: str
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
        if self.result_validation_receipt_sha256:
            require_sha256(
                self.result_validation_receipt_sha256,
                "result_validation_receipt_sha256",
            )
        if self.execution_state not in {
            "not_started",
            "running",
            "engine_complete",
            "validated",
            "failed",
            "ambiguous",
        }:
            raise ContractError("invalid execution state")
        for name, value in (
            ("wrapper_exit_status", self.wrapper_exit_status),
            ("child_exit_status", self.child_exit_status),
        ):
            if value is not None and (
                isinstance(value, bool) or not isinstance(value, int)
            ):
                raise ContractError(f"{name} must be an integer or null")
        if self.engine_complete and self.child_exit_status != 0:
            raise ContractError(
                "engine completion requires child exit status zero"
            )
        if self.engine_receipt_sha256:
            require_sha256(self.engine_receipt_sha256, "engine_receipt_sha256")
        if (
            self.execution_state == "engine_complete"
            and not self.engine_complete
        ):
            raise ContractError(
                "engine_complete state requires engine completion"
            )
        if self.validated and self.execution_state != "validated":
            raise ContractError(
                "validated execution requires validated execution state"
            )
        if self.execution_state == "validated" and not self.validated:
            raise ContractError(
                "validated execution state requires validation"
            )
        if self.validated and not self.engine_complete:
            raise ContractError("validation requires engine completion")
        if self.validated and not self.validator_receipt_sha256s:
            raise ContractError(
                "validated execution requires validator evidence"
            )
        if self.validated and not self.result_validation_receipt_sha256:
            raise ContractError(
                "validated execution requires a resolvable result validation receipt"
            )
        if self.result_validation_receipt_sha256 and (
            self.result_validation_receipt_sha256
            not in self.validator_receipt_sha256s
        ):
            raise ContractError(
                "result validation receipt must be included in validator evidence"
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
            "wrapper_exit_status": self.wrapper_exit_status,
            "child_exit_status": self.child_exit_status,
            "engine_complete": self.engine_complete,
            "validated": self.validated,
            "engine_receipt_sha256": self.engine_receipt_sha256,
            "environment_receipt_sha256": self.environment_receipt_sha256,
            "result_validation_receipt_sha256": (
                self.result_validation_receipt_sha256
            ),
            "output_artifacts": self.output_artifacts,
            "validator_receipt_sha256s": self.validator_receipt_sha256s,
            "findings": self.findings,
            "started_at": self.started_at,
            "finished_at": self.finished_at,
        }
        if self.receipt_sha256 != canonical_sha256(body):
            raise ContractError("program execution receipt digest mismatch")

    @property
    def exit_status(self) -> int | None:
        """Compatibility alias for the outer ChemSmart wrapper status."""

        return self.wrapper_exit_status


def build_program_execution_receipt(
    invocation: ProgramExecutionInvocationV1,
    *,
    execution_state: str,
    exit_status: int | None,
    engine_complete: bool,
    validated: bool,
    child_exit_status: int | None = None,
    engine_receipt_sha256: str = "",
    result_validation_receipt_sha256: str = "",
    output_artifacts: Sequence[TrustedArtifactRefV1] = (),
    validator_receipt_sha256s: Sequence[str] = (),
    findings: Sequence[str] = (),
    started_at: str,
    finished_at: str = "",
) -> ProgramExecutionReceiptV1:
    """Record wrapper transport, child completion, and validation separately.

    ``exit_status`` remains the input spelling for source compatibility and
    is materialized as ``wrapper_exit_status`` in the versioned receipt.
    Existing non-nested engines default their child status to the wrapper
    status when completion was independently asserted.
    """

    for artifact in output_artifacts:
        _require_current_artifact(artifact, "execution output")
    resolved_child_exit_status = child_exit_status
    if resolved_child_exit_status is None and engine_complete:
        resolved_child_exit_status = exit_status
    body = {
        "schema_version": "chemsmart.program-execution-receipt.v1",
        "invocation_sha256": invocation.invocation_sha256,
        "node_id": invocation.node_id,
        "idempotency_key": invocation.idempotency_key,
        "execution_state": str(execution_state),
        "wrapper_exit_status": exit_status,
        "child_exit_status": resolved_child_exit_status,
        "engine_complete": bool(engine_complete),
        "validated": bool(validated),
        "engine_receipt_sha256": str(engine_receipt_sha256).strip(),
        "environment_receipt_sha256": invocation.environment_receipt_sha256,
        "result_validation_receipt_sha256": str(
            result_validation_receipt_sha256
        ).strip(),
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
    """Validated program OPT result converted to one immutable XYZ frame."""

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
    #: The source ``charge``/``multiplicity`` above remain the state on which
    #: the geometry was optimized.  These optional fields name a deliberately
    #: different state for the consumer while preserving the exact coordinates
    #: and atom order.  They are omitted for legacy same-state handoffs.
    consumer_charge: int | None = None
    consumer_multiplicity: int | None = None

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
        if (self.consumer_charge is None) != (
            self.consumer_multiplicity is None
        ):
            raise ContractError(
                "handoff consumer charge and multiplicity must be declared together"
            )
        if self.consumer_charge is not None:
            if isinstance(self.consumer_charge, bool) or not isinstance(
                self.consumer_charge, int
            ):
                raise ContractError(
                    "handoff consumer charge must be an integer"
                )
            if isinstance(self.consumer_multiplicity, bool) or not isinstance(
                self.consumer_multiplicity, int
            ):
                raise ContractError(
                    "handoff consumer multiplicity must be an integer"
                )
            if self.consumer_multiplicity < 1:
                raise ContractError(
                    "handoff consumer multiplicity must be positive"
                )
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
        if self.consumer_charge is not None:
            body["consumer_charge"] = self.consumer_charge
            body["consumer_multiplicity"] = self.consumer_multiplicity
        if self.receipt_sha256 != canonical_sha256(body):
            raise ContractError("optimized geometry handoff digest mismatch")

    @property
    def consumer_state(self) -> tuple[int, int]:
        """Electronic state applied by the downstream calculation."""

        if self.consumer_charge is None:
            return self.charge, self.multiplicity
        return self.consumer_charge, self.consumer_multiplicity

    @property
    def selected_artifact_id(self) -> str:
        return self.geometry_artifact_id

    @property
    def selected_artifact_sha256(self) -> str:
        return self.geometry_artifact_sha256


@dataclass(frozen=True)
class ORCAHessianHandoffV1:
    """Unique final TS Hessian selected by scientific content.

    ORCA may leave more than one ``.hess`` file during an OptTS calculation.
    The suffix is not a semantic role: this record binds the Hessian whose
    atoms and frequencies agree with the validated final TS result.
    """

    schema_version: str
    producer_node_id: str
    consumer_node_id: str
    producer_edge_sha256: str
    producer_execution_receipt_sha256: str
    result_artifact_id: str
    result_artifact_sha256: str
    hessian_artifact_id: str
    hessian_artifact_sha256: str
    atom_count: int
    symbols: tuple[str, ...]
    positions_sha256: str
    frequencies_sha256: str
    charge: int
    multiplicity: int
    consequential_imaginary_mode_count: int
    status: str
    receipt_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.orca-hessian-handoff.v1":
            raise ContractError("unsupported ORCA Hessian handoff")
        require_identifier(self.producer_node_id, "producer_node_id")
        require_identifier(self.consumer_node_id, "consumer_node_id")
        for name, digest in (
            ("producer_edge_sha256", self.producer_edge_sha256),
            (
                "producer_execution_receipt_sha256",
                self.producer_execution_receipt_sha256,
            ),
            ("result_artifact_sha256", self.result_artifact_sha256),
            ("hessian_artifact_sha256", self.hessian_artifact_sha256),
            ("positions_sha256", self.positions_sha256),
            ("frequencies_sha256", self.frequencies_sha256),
        ):
            require_sha256(digest, name)
        if not self.result_artifact_id or not self.hessian_artifact_id:
            raise ContractError("ORCA Hessian handoff artifact IDs are empty")
        if self.atom_count < 1 or self.atom_count != len(self.symbols):
            raise ContractError("ORCA Hessian atom count differs from symbols")
        _require_nonempty_rows(self.symbols, "symbols")
        if self.multiplicity < 1:
            raise ContractError("ORCA Hessian multiplicity must be positive")
        if self.consequential_imaginary_mode_count != 1:
            raise ContractError("final TS Hessian requires one imaginary mode")
        if self.status != "validated_handoff":
            raise ContractError("ORCA Hessian handoff must be validated")
        body = {
            key: value
            for key, value in self.__dict__.items()
            if key != "receipt_sha256"
        }
        if self.receipt_sha256 != canonical_sha256(body):
            raise ContractError("ORCA Hessian handoff digest mismatch")

    @property
    def consumer_state(self) -> tuple[int, int]:
        return self.charge, self.multiplicity

    @property
    def selected_artifact_id(self) -> str:
        return self.hessian_artifact_id

    @property
    def selected_artifact_sha256(self) -> str:
        return self.hessian_artifact_sha256


def _handoff_consumer_fields(
    *,
    producer_charge: int,
    producer_multiplicity: int,
    consumer_charge: int | None,
    consumer_multiplicity: int | None,
) -> dict[str, int]:
    """Validate and canonically encode an optional state rebind."""

    if (consumer_charge is None) != (consumer_multiplicity is None):
        raise ContractError(
            "handoff consumer charge and multiplicity must be declared together"
        )
    if consumer_charge is None:
        return {}
    if isinstance(consumer_charge, bool) or not isinstance(
        consumer_charge, int
    ):
        raise ContractError("handoff consumer charge must be an integer")
    if isinstance(consumer_multiplicity, bool) or not isinstance(
        consumer_multiplicity, int
    ):
        raise ContractError("handoff consumer multiplicity must be an integer")
    if consumer_multiplicity < 1:
        raise ContractError("handoff consumer multiplicity must be positive")
    target = (consumer_charge, consumer_multiplicity)
    if target == (producer_charge, producer_multiplicity):
        # An explicit no-op has the same canonical record as legacy inheritance.
        return {}
    return {
        "consumer_charge": consumer_charge,
        "consumer_multiplicity": consumer_multiplicity,
    }


def handoff_optimized_pyscf_geometry(
    *,
    producer_receipt: ProgramExecutionReceiptV1,
    result_artifact: TrustedArtifactRefV1,
    producer_edge: ProducerEdgeRuleV1,
    approved_workspace: str | Path,
    geometry_artifact_id: str,
    expected_charge: int,
    expected_multiplicity: int,
    consumer_charge: int | None = None,
    consumer_multiplicity: int | None = None,
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
    consumer_fields = _handoff_consumer_fields(
        producer_charge=charge,
        producer_multiplicity=multiplicity,
        consumer_charge=consumer_charge,
        consumer_multiplicity=consumer_multiplicity,
    )

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
        **consumer_fields,
    }
    return geometry_artifact, OptimizedGeometryHandoffV1(
        **body, receipt_sha256=canonical_sha256(body)
    )


def _read_exact_xyz_geometry(
    artifact: TrustedArtifactRefV1,
    *,
    label: str,
) -> tuple[tuple[str, ...], tuple[tuple[float, float, float], ...]]:
    """Read one immutable XYZ frame without accepting trailing frames."""

    source = _require_current_artifact(artifact, label)
    if source.suffix.lower() != ".xyz" or artifact.kind != "geometry_xyz":
        raise ContractError(f"{label} must be an XYZ geometry artifact")
    before = file_sha256(source)
    try:
        lines = source.read_text(encoding="utf-8").splitlines()
        atom_count = int(lines[0].strip())
    except (IndexError, OSError, UnicodeDecodeError, ValueError) as exc:
        raise ContractError(f"{label} is not readable XYZ") from exc
    if atom_count < 1 or len(lines) != atom_count + 2:
        raise ContractError(f"{label} must contain exactly one XYZ frame")
    symbols = []
    positions = []
    try:
        for line in lines[2:]:
            fields = line.split()
            if len(fields) != 4:
                raise ValueError
            symbol = str(fields[0]).strip()
            position = tuple(float(value) for value in fields[1:])
            if not symbol or any(
                not math.isfinite(value) for value in position
            ):
                raise ValueError
            symbols.append(symbol)
            positions.append(position)
    except ValueError as exc:
        raise ContractError(
            f"{label} has invalid symbols or coordinates"
        ) from exc
    after = file_sha256(source)
    if before != artifact.sha256 or after != before:
        raise ContractError(f"{label} changed while it was read")
    return tuple(symbols), tuple(positions)


def handoff_optimized_xtb_geometry(
    *,
    producer_receipt: ProgramExecutionReceiptV1,
    result_artifact: TrustedArtifactRefV1,
    input_artifact: TrustedArtifactRefV1,
    producer_edge: ProducerEdgeRuleV1,
    approved_workspace: str | Path,
    geometry_artifact_id: str,
    expected_charge: int,
    expected_multiplicity: int,
    consumer_charge: int | None = None,
    consumer_multiplicity: int | None = None,
) -> tuple[TrustedArtifactRefV1, OptimizedGeometryHandoffV1]:
    """Bind validated xTB ``xtbopt.xyz`` to the normal geometry edge."""

    if not producer_receipt.validated:
        raise ContractError("optimized geometry requires a validated producer")
    if producer_receipt.node_id != producer_edge.producer_node_id:
        raise ContractError("execution receipt does not match producer edge")
    if producer_edge.selection_rule != "validated_optimized_geometry":
        raise ContractError("unsupported optimized geometry selection rule")
    if Path(result_artifact.path).name != "xtbopt.xyz":
        raise ContractError("xTB optimized handoff requires exact xtbopt.xyz")
    if not any(
        item.artifact_id == result_artifact.artifact_id
        and item.sha256 == result_artifact.sha256
        for item in producer_receipt.output_artifacts
    ):
        raise ContractError("xTB geometry is not bound to producer receipt")
    expected_symbols, _initial_positions = _read_exact_xyz_geometry(
        input_artifact, label="xTB optimization input"
    )
    symbols, positions = _read_exact_xyz_geometry(
        result_artifact, label="xTB optimized geometry"
    )
    if symbols != expected_symbols:
        raise ContractError(
            "xTB optimized geometry changed atom identity or atom order"
        )

    charge = int(expected_charge)
    multiplicity = int(expected_multiplicity)
    if multiplicity < 1:
        raise ContractError(
            "optimized electronic-state multiplicity is invalid"
        )
    consumer_fields = _handoff_consumer_fields(
        producer_charge=charge,
        producer_multiplicity=multiplicity,
        consumer_charge=consumer_charge,
        consumer_multiplicity=consumer_multiplicity,
    )
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
        "ChemSmart validated xTB OPT handoff; "
        f"charge={charge}; multiplicity={multiplicity}; "
        f"source_sha256={result_artifact.sha256}"
    )
    lines = [str(len(symbols)), comment]
    for symbol, (x, y, z) in zip(symbols, positions):
        lines.append(f"{symbol:<3} {x:.17g} {y:.17g} {z:.17g}")
    _write_exact_once(target, ("\n".join(lines) + "\n").encode("utf-8"))
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
        **consumer_fields,
    }
    return geometry_artifact, OptimizedGeometryHandoffV1(
        **body, receipt_sha256=canonical_sha256(body)
    )


def handoff_optimized_native_geometry(
    *,
    program: str,
    producer_receipt: ProgramExecutionReceiptV1,
    result_artifact: TrustedArtifactRefV1,
    input_artifact: TrustedArtifactRefV1,
    producer_edge: ProducerEdgeRuleV1,
    approved_workspace: str | Path,
    geometry_artifact_id: str,
    expected_charge: int,
    expected_multiplicity: int,
    consumer_charge: int | None = None,
    consumer_multiplicity: int | None = None,
) -> tuple[TrustedArtifactRefV1, OptimizedGeometryHandoffV1]:
    """Extract a validated ORCA/Gaussian OPT or TS final geometry."""

    normalized_program = require_identifier(program, "program")
    if normalized_program not in {"gaussian", "orca"}:
        raise ContractError(
            "native optimized handoff supports Gaussian or ORCA"
        )
    if not producer_receipt.validated:
        raise ContractError("optimized geometry requires a validated producer")
    if producer_receipt.node_id != producer_edge.producer_node_id:
        raise ContractError("execution receipt does not match producer edge")
    if producer_edge.selection_rule != "validated_optimized_geometry":
        raise ContractError("unsupported optimized geometry selection rule")
    expected_kind = f"{normalized_program}_output"
    source = _require_current_artifact(
        result_artifact, f"{normalized_program} result artifact"
    )
    if result_artifact.kind != expected_kind:
        raise ContractError(
            f"optimized handoff requires a {expected_kind} artifact"
        )
    if not any(
        item.artifact_id == result_artifact.artifact_id
        and item.sha256 == result_artifact.sha256
        for item in producer_receipt.output_artifacts
    ):
        raise ContractError("result artifact is not bound to producer receipt")
    expected_symbols, _initial_positions = _read_exact_xyz_geometry(
        input_artifact, label=f"{normalized_program} optimization input"
    )
    before = file_sha256(source)
    try:
        if normalized_program == "gaussian":
            from chemsmart.io.gaussian.output import Gaussian16Output

            output = Gaussian16Output(str(source))
            jobtype = str(output.jobtype or "").lower()
            normal_termination = bool(output.normal_termination)
            converged = (
                any(
                    "Optimization completed." in line
                    for line in output.contents
                )
                and not output.convergence_criterion_not_met
            )
            molecule = output.optimized_structure
            charge = output.charge
            multiplicity = output.multiplicity
        else:
            from chemsmart.io.orca.output import ORCAOutput

            output = ORCAOutput(str(source))
            jobtype = str(output.jobtype or "").lower()
            normal_termination = bool(output.normal_termination)
            converged = output.converged is True
            molecule = output.molecule
            charge = output.charge
            multiplicity = output.multiplicity
    except (AttributeError, IndexError, OSError, TypeError, ValueError) as exc:
        raise ContractError(
            f"{normalized_program} optimized result is not readable"
        ) from exc
    after = file_sha256(source)
    if before != result_artifact.sha256 or after != before:
        raise ContractError(
            f"{normalized_program} result changed while extracting geometry"
        )
    if not normal_termination or not converged or jobtype not in {"opt", "ts"}:
        raise ContractError(
            f"{normalized_program} result is not a converged OPT or TS"
        )
    if (charge, multiplicity) != (
        int(expected_charge),
        int(expected_multiplicity),
    ):
        raise ContractError("optimized electronic state differs from approval")
    consumer_fields = _handoff_consumer_fields(
        producer_charge=int(charge),
        producer_multiplicity=int(multiplicity),
        consumer_charge=consumer_charge,
        consumer_multiplicity=consumer_multiplicity,
    )
    symbols = tuple(str(item) for item in molecule.chemical_symbols)
    try:
        positions = tuple(
            tuple(float(component) for component in row)
            for row in molecule.positions
        )
    except (TypeError, ValueError) as exc:
        raise ContractError(
            "optimized final positions are not numeric"
        ) from exc
    if symbols != expected_symbols:
        raise ContractError(
            f"{normalized_program} optimized geometry changed atom identity "
            "or atom order"
        )
    if any(
        len(row) != 3 or any(not math.isfinite(value) for value in row)
        for row in positions
    ):
        raise ContractError(
            "optimized final positions must be finite Nx3 values"
        )

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
        f"ChemSmart validated {normalized_program} {jobtype.upper()} handoff; "
        f"charge={charge}; multiplicity={multiplicity}; "
        f"source_sha256={result_artifact.sha256}"
    )
    lines = [str(len(symbols)), comment]
    for symbol, (x, y, z) in zip(symbols, positions):
        lines.append(f"{symbol:<3} {x:.17g} {y:.17g} {z:.17g}")
    _write_exact_once(target, ("\n".join(lines) + "\n").encode("utf-8"))
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
        "charge": int(charge),
        "multiplicity": int(multiplicity),
        "status": "validated_handoff",
        **consumer_fields,
    }
    return geometry_artifact, OptimizedGeometryHandoffV1(
        **body, receipt_sha256=canonical_sha256(body)
    )


@dataclass(frozen=True)
class _ParsedORCAHessian:
    dimension: int
    matrix: np.ndarray
    symbols: tuple[str, ...]
    positions_bohr: np.ndarray
    frequencies_cm1: tuple[float, ...]


def _orca_hessian_sections(path: Path) -> dict[str, tuple[str, ...]]:
    try:
        lines = path.read_text(encoding="utf-8", errors="strict").splitlines()
    except (OSError, UnicodeError) as exc:
        raise ContractError("ORCA Hessian is not readable text") from exc
    sections: dict[str, list[str]] = {}
    active = ""
    for raw_line in lines:
        line = raw_line.strip()
        if line.startswith("$"):
            active = line.casefold()
            if active in sections:
                raise ContractError("ORCA Hessian repeats a required section")
            sections[active] = []
            continue
        if active:
            sections[active].append(raw_line)
    return {name: tuple(values) for name, values in sections.items()}


def _nonempty_section_lines(
    sections: Mapping[str, tuple[str, ...]], name: str
) -> list[str]:
    values = sections.get(name)
    if values is None:
        raise ContractError(f"ORCA Hessian lacks {name}")
    return [line.strip() for line in values if line.strip()]


def _parse_orca_hessian(path: Path, *, symmetry_atol: float = 1e-10) -> _ParsedORCAHessian:
    """Parse the native Hessian, atoms and frequency blocks fail-closed."""

    sections = _orca_hessian_sections(path)
    hessian_lines = _nonempty_section_lines(sections, "$hessian")
    try:
        dimension = int(hessian_lines[0])
    except (IndexError, ValueError) as exc:
        raise ContractError("ORCA Hessian dimension is invalid") from exc
    if dimension < 3:
        raise ContractError("ORCA Hessian dimension is too small")
    matrix = np.full((dimension, dimension), np.nan, dtype=float)
    cursor = 1
    try:
        while cursor < len(hessian_lines):
            columns = tuple(int(value) for value in hessian_lines[cursor].split())
            cursor += 1
            if not columns or any(
                column < 0 or column >= dimension for column in columns
            ):
                raise ValueError
            for _ in range(dimension):
                row = hessian_lines[cursor].split()
                cursor += 1
                row_index = int(row[0])
                values = tuple(float(value) for value in row[1:])
                if (
                    row_index < 0
                    or row_index >= dimension
                    or len(values) != len(columns)
                ):
                    raise ValueError
                matrix[row_index, columns] = values
    except (IndexError, ValueError) as exc:
        raise ContractError("ORCA Hessian matrix is malformed") from exc
    if not np.isfinite(matrix).all():
        raise ContractError("ORCA Hessian matrix is incomplete or non-finite")
    if not np.allclose(matrix, matrix.T, rtol=1e-9, atol=symmetry_atol):
        raise ContractError("ORCA Hessian matrix is not symmetric")

    atom_lines = _nonempty_section_lines(sections, "$atoms")
    frequency_lines = _nonempty_section_lines(
        sections, "$vibrational_frequencies"
    )
    try:
        atom_count = int(atom_lines[0])
        frequency_count = int(frequency_lines[0])
        atom_rows = tuple(line.split() for line in atom_lines[1:])
        frequency_rows = tuple(line.split() for line in frequency_lines[1:])
        if len(atom_rows) != atom_count or len(frequency_rows) != frequency_count:
            raise ValueError
        symbols = tuple(row[0] for row in atom_rows)
        positions_bohr = np.asarray(
            [[float(value) for value in row[2:5]] for row in atom_rows],
            dtype=float,
        )
        frequencies = tuple(float(row[1]) for row in frequency_rows)
        if any(int(row[0]) != index for index, row in enumerate(frequency_rows)):
            raise ValueError
    except (IndexError, TypeError, ValueError) as exc:
        raise ContractError("ORCA Hessian atoms or frequencies are malformed") from exc
    if dimension != 3 * atom_count or frequency_count != dimension:
        raise ContractError("ORCA Hessian blocks do not have finite 3N shape")
    if positions_bohr.shape != (atom_count, 3) or not np.isfinite(
        positions_bohr
    ).all():
        raise ContractError("ORCA Hessian atom coordinates are non-finite")
    if not all(math.isfinite(value) for value in frequencies):
        raise ContractError("ORCA Hessian frequencies are non-finite")
    return _ParsedORCAHessian(
        dimension=dimension,
        matrix=matrix,
        symbols=symbols,
        positions_bohr=positions_bohr,
        frequencies_cm1=frequencies,
    )


def _same_orca_geometry_up_to_rotation(
    positions_bohr: np.ndarray, positions_angstrom: np.ndarray
) -> bool:
    """Compare one atom-ordered frame while permitting ORCA rigid rotation."""

    from ase import units

    candidate = positions_bohr * float(units.Bohr)
    reference = np.asarray(positions_angstrom, dtype=float)
    if candidate.shape != reference.shape or candidate.ndim != 2:
        return False
    candidate_centered = candidate - candidate.mean(axis=0)
    reference_centered = reference - reference.mean(axis=0)
    try:
        left, _singular, right = np.linalg.svd(
            candidate_centered.T @ reference_centered
        )
    except np.linalg.LinAlgError:
        return False
    rotation = left @ right
    if np.linalg.det(rotation) < 0:
        left[:, -1] *= -1
        rotation = left @ right
    aligned = candidate_centered @ rotation
    # ORCA's output geometry is printed to 1e-6 Angstrom, whereas $atoms
    # retains substantially more digits in Bohr.  This tolerance covers only
    # that presentation precision, not a changed molecular geometry.
    return bool(np.max(np.abs(aligned - reference_centered)) <= 2.0e-6)


def handoff_final_orca_ts_hessian(
    *,
    producer_receipt: ProgramExecutionReceiptV1,
    result_artifact: TrustedArtifactRefV1,
    hessian_candidates: Sequence[TrustedArtifactRefV1],
    producer_edge: ProducerEdgeRuleV1,
    approved_workspace: str | Path,
    hessian_artifact_id: str,
    expected_charge: int,
    expected_multiplicity: int,
    imaginary_mode_cutoff_cm1: float = 20.0,
) -> tuple[TrustedArtifactRefV1, ORCAHessianHandoffV1]:
    """Select the unique Hessian belonging to a validated final ORCA TS.

    Selection uses the native scientific content.  A numbered filename is
    never assumed to mean initial or final.
    """

    if not producer_receipt.validated:
        raise ContractError("ORCA TS Hessian requires a validated producer")
    if producer_receipt.node_id != producer_edge.producer_node_id:
        raise ContractError("ORCA Hessian receipt differs from producer edge")
    if producer_edge.selection_rule != "validated_final_orca_ts_hessian":
        raise ContractError("unsupported ORCA Hessian selection rule")
    if producer_edge.artifact_kind != "orca_hessian":
        raise ContractError("ORCA Hessian edge has the wrong artifact class")
    source = _require_current_artifact(result_artifact, "ORCA TS output")
    if result_artifact.kind != "orca_output" or not any(
        item.artifact_id == result_artifact.artifact_id
        and item.sha256 == result_artifact.sha256
        for item in producer_receipt.output_artifacts
    ):
        raise ContractError("ORCA TS output is not bound to the producer")
    try:
        from chemsmart.io.orca.output import ORCAOutput

        output = ORCAOutput(str(source))
        if (
            not output.normal_termination
            or output.converged is not True
            or str(output.jobtype or "").casefold() != "ts"
        ):
            raise ContractError("ORCA output is not a converged final TS")
        if (output.charge, output.multiplicity) != (
            int(expected_charge),
            int(expected_multiplicity),
        ):
            raise ContractError("ORCA TS state differs from approval")
        final_symbols = tuple(str(value) for value in output.molecule.chemical_symbols)
        final_positions = np.asarray(output.molecule.positions, dtype=float)
        final_frequencies = tuple(
            float(value) for value in (output.all_vibrational_frequencies or ())
        )
    except ContractError:
        raise
    except (AttributeError, IndexError, OSError, TypeError, ValueError) as exc:
        raise ContractError("ORCA final TS output is unreadable") from exc
    if not final_frequencies or not np.isfinite(final_positions).all():
        raise ContractError("ORCA final TS lacks geometry or frequencies")

    consequential = tuple(
        value
        for value in final_frequencies
        if value < -abs(float(imaginary_mode_cutoff_cm1))
    )
    if len(consequential) != 1:
        raise ContractError("ORCA final TS does not have one imaginary mode")

    matches: list[tuple[TrustedArtifactRefV1, _ParsedORCAHessian]] = []
    for candidate in hessian_candidates:
        if candidate.kind != "orca_hessian" or not any(
            item.artifact_id == candidate.artifact_id
            and item.sha256 == candidate.sha256
            for item in producer_receipt.output_artifacts
        ):
            continue
        path = _require_current_artifact(candidate, "ORCA Hessian candidate")
        try:
            parsed = _parse_orca_hessian(path)
        except ContractError:
            continue
        if parsed.symbols != final_symbols or not _same_orca_geometry_up_to_rotation(
            parsed.positions_bohr, final_positions
        ):
            continue
        if len(parsed.frequencies_cm1) != len(final_frequencies) or any(
            abs(observed - expected) > 0.02
            for observed, expected in zip(
                parsed.frequencies_cm1, final_frequencies
            )
        ):
            continue
        if (
            sum(
                value < -abs(float(imaginary_mode_cutoff_cm1))
                for value in parsed.frequencies_cm1
            )
            != 1
        ):
            continue
        matches.append((candidate, parsed))
    if len(matches) != 1:
        raise ContractError(
            "ORCA final TS Hessian selection is not unique: "
            f"{len(matches)} semantic matches among {len(hessian_candidates)} candidates"
        )
    candidate, parsed = matches[0]
    candidate_path = _require_current_artifact(
        candidate, "selected ORCA Hessian"
    )
    payload = candidate_path.read_bytes()
    target = _target_below(
        _absolute_workspace(approved_workspace),
        "artifacts",
        f"{producer_edge.producer_node_id}--{producer_edge.consumer_node_id}.hess",
    )
    _write_exact_once(target, payload)
    artifact = TrustedArtifactRefV1(
        artifact_id=require_identifier(hessian_artifact_id, "artifact_id"),
        kind="orca_hessian",
        sha256=file_sha256(target),
        size_bytes=target.stat().st_size,
        path=str(target),
        cli_value=str(target),
    )
    positions_sha256 = canonical_sha256(
        {
            "symbols": parsed.symbols,
            "positions_bohr": parsed.positions_bohr.tolist(),
        }
    )
    frequencies_sha256 = canonical_sha256(
        {"frequencies_cm1": parsed.frequencies_cm1}
    )
    body = {
        "schema_version": "chemsmart.orca-hessian-handoff.v1",
        "producer_node_id": producer_edge.producer_node_id,
        "consumer_node_id": producer_edge.consumer_node_id,
        "producer_edge_sha256": producer_edge.edge_sha256,
        "producer_execution_receipt_sha256": producer_receipt.receipt_sha256,
        "result_artifact_id": result_artifact.artifact_id,
        "result_artifact_sha256": result_artifact.sha256,
        "hessian_artifact_id": artifact.artifact_id,
        "hessian_artifact_sha256": artifact.sha256,
        "atom_count": len(parsed.symbols),
        "symbols": parsed.symbols,
        "positions_sha256": positions_sha256,
        "frequencies_sha256": frequencies_sha256,
        "charge": int(output.charge),
        "multiplicity": int(output.multiplicity),
        "consequential_imaginary_mode_count": 1,
        "status": "validated_handoff",
    }
    return artifact, ORCAHessianHandoffV1(
        **body, receipt_sha256=canonical_sha256(body)
    )


@dataclass(frozen=True)
class ORCAProducerHessianHandoffV1:
    """A producer's validated Hessian carried into an ORCA TS search.

    Unlike the final-TS handoff, this is starting curvature, not a
    classification: a starting Hessian for a transition-state search may
    show any imaginary-mode count, and the observed count is recorded as a
    fact rather than enforced as a requirement.
    """

    schema_version: str
    producer_node_id: str
    consumer_node_id: str
    producer_edge_sha256: str
    producer_execution_receipt_sha256: str
    result_artifact_id: str
    result_artifact_sha256: str
    hessian_artifact_id: str
    hessian_artifact_sha256: str
    atom_count: int
    symbols: tuple[str, ...]
    positions_sha256: str
    frequencies_sha256: str
    charge: int
    multiplicity: int
    observed_imaginary_mode_count: int
    status: str
    receipt_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.orca-producer-hessian-handoff.v1":
            raise ContractError("unsupported ORCA producer Hessian handoff")
        require_identifier(self.producer_node_id, "producer_node_id")
        require_identifier(self.consumer_node_id, "consumer_node_id")
        for name, digest in (
            ("producer_edge_sha256", self.producer_edge_sha256),
            (
                "producer_execution_receipt_sha256",
                self.producer_execution_receipt_sha256,
            ),
            ("result_artifact_sha256", self.result_artifact_sha256),
            ("hessian_artifact_sha256", self.hessian_artifact_sha256),
            ("positions_sha256", self.positions_sha256),
            ("frequencies_sha256", self.frequencies_sha256),
        ):
            require_sha256(digest, name)
        if not self.result_artifact_id or not self.hessian_artifact_id:
            raise ContractError(
                "ORCA producer Hessian handoff artifact IDs are empty"
            )
        if self.atom_count < 1 or self.atom_count != len(self.symbols):
            raise ContractError(
                "ORCA producer Hessian atom count differs from symbols"
            )
        _require_nonempty_rows(self.symbols, "symbols")
        if self.multiplicity < 1:
            raise ContractError(
                "ORCA producer Hessian multiplicity must be positive"
            )
        if self.observed_imaginary_mode_count < 0:
            raise ContractError(
                "observed imaginary mode count cannot be negative"
            )
        if self.status != "validated_handoff":
            raise ContractError(
                "ORCA producer Hessian handoff must be validated"
            )
        body = {
            key: value
            for key, value in self.__dict__.items()
            if key != "receipt_sha256"
        }
        if self.receipt_sha256 != canonical_sha256(body):
            raise ContractError(
                "ORCA producer Hessian handoff digest mismatch"
            )

    @property
    def consumer_state(self) -> tuple[int, int]:
        return self.charge, self.multiplicity

    @property
    def selected_artifact_id(self) -> str:
        return self.hessian_artifact_id

    @property
    def selected_artifact_sha256(self) -> str:
        return self.hessian_artifact_sha256


def handoff_validated_orca_producer_hessian(
    *,
    producer_receipt: ProgramExecutionReceiptV1,
    result_artifact: TrustedArtifactRefV1,
    hessian_candidates: Sequence[TrustedArtifactRefV1],
    producer_edge: ProducerEdgeRuleV1,
    approved_workspace: str | Path,
    hessian_artifact_id: str,
    expected_charge: int,
    expected_multiplicity: int,
    imaginary_mode_cutoff_cm1: float = 20.0,
) -> tuple[TrustedArtifactRefV1, ORCAProducerHessianHandoffV1]:
    """Select the Hessian of a validated ORCA frequency-bearing producer.

    Selection matches the native scientific content of the producer's own
    validated result -- symbols, rotation-invariant geometry, frequency
    agreement -- exactly like the final-TS handoff. What it deliberately
    does NOT require: a ``ts`` jobtype or a one-imaginary-mode count. A
    starting Hessian for a TS search may have any imaginary-mode count;
    the observed count is recorded in the receipt.
    """

    if not producer_receipt.validated:
        raise ContractError(
            "ORCA producer Hessian requires a validated producer"
        )
    if producer_receipt.node_id != producer_edge.producer_node_id:
        raise ContractError(
            "ORCA producer Hessian receipt differs from producer edge"
        )
    if producer_edge.selection_rule != "validated_producer_orca_hessian":
        raise ContractError(
            "unsupported ORCA producer Hessian selection rule"
        )
    if producer_edge.artifact_kind != "orca_hessian":
        raise ContractError(
            "ORCA producer Hessian edge has the wrong artifact class"
        )
    source = _require_current_artifact(result_artifact, "ORCA producer output")
    if result_artifact.kind != "orca_output" or not any(
        item.artifact_id == result_artifact.artifact_id
        and item.sha256 == result_artifact.sha256
        for item in producer_receipt.output_artifacts
    ):
        raise ContractError(
            "ORCA producer output is not bound to the producer"
        )
    try:
        from chemsmart.io.orca.output import ORCAOutput

        output = ORCAOutput(str(source))
        if not output.normal_termination:
            raise ContractError(
                "ORCA producer output did not terminate normally"
            )
        if (output.charge, output.multiplicity) != (
            int(expected_charge),
            int(expected_multiplicity),
        ):
            raise ContractError(
                "ORCA producer state differs from approval"
            )
        final_symbols = tuple(
            str(value) for value in output.molecule.chemical_symbols
        )
        final_positions = np.asarray(output.molecule.positions, dtype=float)
        final_frequencies = tuple(
            float(value)
            for value in (output.all_vibrational_frequencies or ())
        )
    except ContractError:
        raise
    except (AttributeError, IndexError, OSError, TypeError, ValueError) as exc:
        raise ContractError("ORCA producer output is unreadable") from exc
    if not final_frequencies or not np.isfinite(final_positions).all():
        raise ContractError(
            "a starting Hessian needs a frequency-bearing producer: the "
            "validated result carries no vibrational frequencies"
        )
    observed_imaginary = sum(
        value < -abs(float(imaginary_mode_cutoff_cm1))
        for value in final_frequencies
    )
    matches: list[tuple[TrustedArtifactRefV1, _ParsedORCAHessian]] = []
    for candidate in hessian_candidates:
        if candidate.kind != "orca_hessian" or not any(
            item.artifact_id == candidate.artifact_id
            and item.sha256 == candidate.sha256
            for item in producer_receipt.output_artifacts
        ):
            continue
        path = _require_current_artifact(candidate, "ORCA Hessian candidate")
        try:
            # A numerically assembled opt/freq Hessian carries ~1e-3
            # asymmetry that an analytic TS Hessian does not. As a STARTING
            # Hessian it is curvature to seed the search with, not a
            # certified classification, so the strict symmetry gate of the
            # final-TS handoff would refuse exactly the artifacts this edge
            # exists to carry. Observed on a real opt+freq fixture.
            parsed = _parse_orca_hessian(path, symmetry_atol=5.0e-3)
        except ContractError:
            continue
        if (
            parsed.symbols != final_symbols
            or not _same_orca_geometry_up_to_rotation(
                parsed.positions_bohr, final_positions
            )
        ):
            continue
        if len(parsed.frequencies_cm1) != len(final_frequencies) or any(
            abs(observed - expected) > 0.02
            for observed, expected in zip(
                parsed.frequencies_cm1, final_frequencies
            )
        ):
            continue
        matches.append((candidate, parsed))
    if len(matches) != 1:
        raise ContractError(
            "ORCA producer Hessian selection is not unique: "
            f"{len(matches)} semantic matches among "
            f"{len(hessian_candidates)} candidates"
        )
    candidate, parsed = matches[0]
    candidate_path = _require_current_artifact(
        candidate, "selected ORCA Hessian"
    )
    payload = candidate_path.read_bytes()
    target = _target_below(
        _absolute_workspace(approved_workspace),
        "artifacts",
        f"{producer_edge.producer_node_id}--"
        f"{producer_edge.consumer_node_id}.hess",
    )
    _write_exact_once(target, payload)
    artifact = TrustedArtifactRefV1(
        artifact_id=require_identifier(hessian_artifact_id, "artifact_id"),
        kind="orca_hessian",
        sha256=file_sha256(target),
        size_bytes=target.stat().st_size,
        path=str(target),
        cli_value=str(target),
    )
    positions_sha256 = canonical_sha256(
        {
            "symbols": parsed.symbols,
            "positions_bohr": parsed.positions_bohr.tolist(),
        }
    )
    frequencies_sha256 = canonical_sha256(
        {"frequencies_cm1": parsed.frequencies_cm1}
    )
    body = {
        "schema_version": "chemsmart.orca-producer-hessian-handoff.v1",
        "producer_node_id": producer_edge.producer_node_id,
        "consumer_node_id": producer_edge.consumer_node_id,
        "producer_edge_sha256": producer_edge.edge_sha256,
        "producer_execution_receipt_sha256": producer_receipt.receipt_sha256,
        "result_artifact_id": result_artifact.artifact_id,
        "result_artifact_sha256": result_artifact.sha256,
        "hessian_artifact_id": artifact.artifact_id,
        "hessian_artifact_sha256": artifact.sha256,
        "atom_count": len(parsed.symbols),
        "symbols": parsed.symbols,
        "positions_sha256": positions_sha256,
        "frequencies_sha256": frequencies_sha256,
        "charge": int(output.charge),
        "multiplicity": int(output.multiplicity),
        "observed_imaginary_mode_count": int(observed_imaginary),
        "status": "validated_handoff",
    }
    return artifact, ORCAProducerHessianHandoffV1(
        **body, receipt_sha256=canonical_sha256(body)
    )


@dataclass(frozen=True)
class FrozenMaterializedNodePreviewV1:
    """Exact green preview admitted for one currently executable node."""

    schema_version: str
    node_id: str
    input_artifact_sha256: str
    project_artifact_sha256: str
    project_validation_receipt_sha256: str
    environment_receipt_sha256: str
    invocation_sha256: str
    preflight_receipt_sha256: str
    binding_sha256: str
    #: What admission actually compares.  ``invocation_sha256`` is retained as
    #: the exact record of what this session previewed, but it is path-scoped
    #: and cannot be matched by a later process.
    invocation_identity_sha256: str = ""
    auxiliary_input_bindings: tuple[AuxiliaryArtifactBindingV1, ...] = ()

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.frozen-node-preview.v1":
            raise ContractError("unsupported frozen node preview schema")
        require_identifier(self.node_id, "node_id")
        for name, digest in (
            ("input_artifact_sha256", self.input_artifact_sha256),
            ("project_artifact_sha256", self.project_artifact_sha256),
            (
                "project_validation_receipt_sha256",
                self.project_validation_receipt_sha256,
            ),
            ("environment_receipt_sha256", self.environment_receipt_sha256),
            ("invocation_sha256", self.invocation_sha256),
            ("preflight_receipt_sha256", self.preflight_receipt_sha256),
        ):
            require_sha256(digest, name)
        require_auxiliary_artifact_bindings(self.auxiliary_input_bindings)
        body = _frozen_preview_binding_body(self)
        if self.binding_sha256 != canonical_sha256(body):
            raise ContractError("frozen node preview digest mismatch")


def _frozen_preview_binding_body(
    binding: FrozenMaterializedNodePreviewV1,
) -> dict[str, Any]:
    body = {
        "schema_version": binding.schema_version,
        "node_id": binding.node_id,
        "input_artifact_sha256": binding.input_artifact_sha256,
        "project_artifact_sha256": binding.project_artifact_sha256,
        "project_validation_receipt_sha256": (
            binding.project_validation_receipt_sha256
        ),
        "environment_receipt_sha256": binding.environment_receipt_sha256,
        "invocation_sha256": binding.invocation_sha256,
        "preflight_receipt_sha256": binding.preflight_receipt_sha256,
    }
    if binding.invocation_identity_sha256:
        require_sha256(
            binding.invocation_identity_sha256,
            "invocation_identity_sha256",
        )
        body["invocation_identity_sha256"] = binding.invocation_identity_sha256
    if binding.auxiliary_input_bindings:
        body["auxiliary_input_bindings"] = binding.auxiliary_input_bindings
    return body


@dataclass(frozen=True)
class FrozenProducerEdgeRuleV1:
    """Exact future-artifact selection rule admitted before execution."""

    schema_version: str
    scientific_edge_sha256: str
    source_node_id: str
    target_node_id: str
    artifact_class: str
    producer_output_id: str
    consumer_input_id: str
    selection_rule: str
    environment_receipt_sha256: str
    preserve_atom_order: bool
    preserve_electronic_state: bool
    rule_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.frozen-producer-edge-rule.v1":
            raise ContractError("unsupported frozen producer edge rule")
        require_sha256(self.scientific_edge_sha256, "scientific_edge_sha256")
        for name, value in (
            ("source_node_id", self.source_node_id),
            ("target_node_id", self.target_node_id),
            ("artifact_class", self.artifact_class),
            ("producer_output_id", self.producer_output_id),
            ("consumer_input_id", self.consumer_input_id),
            ("selection_rule", self.selection_rule),
        ):
            require_identifier(value, name)
        require_sha256(
            self.environment_receipt_sha256,
            "environment_receipt_sha256",
        )
        if not self.preserve_atom_order or not self.preserve_electronic_state:
            raise ContractError(
                "molecular producer edges must preserve atom order and state"
            )
        body = {
            "schema_version": self.schema_version,
            "scientific_edge_sha256": self.scientific_edge_sha256,
            "source_node_id": self.source_node_id,
            "target_node_id": self.target_node_id,
            "artifact_class": self.artifact_class,
            "producer_output_id": self.producer_output_id,
            "consumer_input_id": self.consumer_input_id,
            "selection_rule": self.selection_rule,
            "environment_receipt_sha256": self.environment_receipt_sha256,
            "preserve_atom_order": self.preserve_atom_order,
            "preserve_electronic_state": self.preserve_electronic_state,
        }
        if self.rule_sha256 != canonical_sha256(body):
            raise ContractError("frozen producer edge rule digest mismatch")


def _frozen_preview_binding(node: Any) -> FrozenMaterializedNodePreviewV1:
    if node.state != "previewed":
        raise ContractError(
            "every currently materialized approved node requires exact preview"
        )
    body = {
        "schema_version": "chemsmart.frozen-node-preview.v1",
        "node_id": node.node_id,
        "input_artifact_sha256": node.input_artifact_sha256,
        "project_artifact_sha256": node.project_artifact_sha256,
        "project_validation_receipt_sha256": (
            node.project_validation_receipt_sha256
        ),
        "environment_receipt_sha256": node.environment_receipt_sha256,
        "invocation_sha256": node.invocation_sha256,
        "preflight_receipt_sha256": node.preflight_receipt_sha256,
    }
    identity = getattr(node, "invocation_identity_sha256", "")
    if identity:
        body["invocation_identity_sha256"] = identity
    auxiliary = getattr(node, "auxiliary_input_bindings", ())
    if auxiliary:
        body["auxiliary_input_bindings"] = auxiliary
    return FrozenMaterializedNodePreviewV1(
        **body, binding_sha256=canonical_sha256(body)
    )


def _frozen_producer_edge_rule(
    plan: ScientificWorkflowPlanV2,
    edge: ScientificWorkflowEdgeV2,
    *,
    environment_receipt_sha256: str,
) -> FrozenProducerEdgeRuleV1:
    if is_validated_optimized_geometry_edge(plan, edge):
        selection_rule = "validated_optimized_geometry"
    elif is_validated_orca_ts_hessian_edge(plan, edge):
        selection_rule = "validated_final_orca_ts_hessian"
    elif is_validated_producer_orca_hessian_edge(plan, edge):
        selection_rule = "validated_producer_orca_hessian"
    else:
        raise ContractError(
            "future data edge has no registered exact artifact selection rule: "
            "supported edges are an opt/ts geometry_xyz input, an ORCA TS "
            "orca_hessian input to an IRC hess_filename role, or an ORCA "
            "producer orca_hessian input to a TS inhess_filename role"
        )
    body = {
        "schema_version": "chemsmart.frozen-producer-edge-rule.v1",
        "scientific_edge_sha256": canonical_sha256(edge),
        "source_node_id": edge.source_node_id,
        "target_node_id": edge.target_node_id,
        "artifact_class": edge.artifact_class,
        "producer_output_id": edge.producer_output_id,
        "consumer_input_id": edge.consumer_input_id,
        "selection_rule": selection_rule,
        "environment_receipt_sha256": require_sha256(
            environment_receipt_sha256,
            "environment_receipt_sha256",
        ),
        "preserve_atom_order": True,
        "preserve_electronic_state": True,
    }
    return FrozenProducerEdgeRuleV1(**body, rule_sha256=canonical_sha256(body))


#: Producer stages whose geometry a consumer may wait on inside one approval.
#: Each of these ends at a single stationary structure that ChemSmart's parser
#: selects and validates without anyone choosing between candidates.
#:
#: A relaxed scan is deliberately absent. It ends at a surface, not a structure,
#: and picking which point to carry forward is a scientific judgement that the
#: computed surface has to inform -- so it cannot be settled in advance, and the
#: host must not settle it on the scientist's behalf.
DEFERRABLE_GEOMETRY_PRODUCER_STAGES = frozenset({"opt", "ts"})


def is_validated_optimized_geometry_edge(
    plan: ScientificWorkflowPlanV2,
    edge: ScientificWorkflowEdgeV2,
) -> bool:
    """Whether one future edge has ChemSmart's exact geometry selector.

    Both a minimum optimization and a transition-state optimization produce a
    stationary geometry that the normal parser can select and validate.  The
    older approval predicate admitted only ``opt``.  A real TS -> IRC/SP DAG
    was consequently rejected even though the producer output was named
    exactly and the same native handoff path supports both stages.

    The output identifier is a workflow-local presentation label.  It does
    not select a native file: ChemSmart's registered rule below selects the
    validated final geometry from an ``opt`` or ``ts`` result.  Requiring the
    model to repeat the literal ``optimized_geometry`` in that local label or
    observable list created a second authority over the typed stage and edge.
    The scientific meaning is therefore established by the producer stage and
    ``geometry_xyz`` artifact class alone.
    """

    if edge.edge_kind != "data" or edge.artifact_class != "geometry_xyz":
        return False
    source = next(
        (
            node
            for node in plan.nodes
            if node.node_id == edge.source_node_id
        ),
        None,
    )
    return bool(
        source is not None
        and source.stage in DEFERRABLE_GEOMETRY_PRODUCER_STAGES
    )


def is_validated_orca_ts_hessian_edge(
    plan: ScientificWorkflowPlanV2,
    edge: ScientificWorkflowEdgeV2,
) -> bool:
    """Whether an edge requests the native Hessian of a final ORCA TS."""

    if (
        edge.edge_kind != "data"
        or edge.artifact_class != "orca_hessian"
        or edge.consumer_input_id != "hess_filename"
    ):
        return False
    source = next(
        (node for node in plan.nodes if node.node_id == edge.source_node_id),
        None,
    )
    target = next(
        (node for node in plan.nodes if node.node_id == edge.target_node_id),
        None,
    )
    return bool(
        source is not None
        and target is not None
        and source.program == "orca"
        and source.stage == "ts"
        and target.program == "orca"
        and target.stage == "irc"
    )


def is_validated_producer_orca_hessian_edge(
    plan: ScientificWorkflowPlanV2,
    edge: ScientificWorkflowEdgeV2,
) -> bool:
    """Whether an edge feeds a producer's Hessian into an ORCA TS search.

    A starting Hessian for a transition-state SEARCH is a different fact
    from the final Hessian of a converged TS: it may come from any
    frequency-bearing ORCA stage and may show any imaginary-mode count --
    curvature information to start from, not a classification to certify.
    The consumer role is the live ``--inhess-filename`` option on
    ``run/orca/ts``; the ``inhess: true`` switch is method configuration
    and lives in the project YAML ts section.
    """

    if (
        edge.edge_kind != "data"
        or edge.artifact_class != "orca_hessian"
        or edge.consumer_input_id != "inhess_filename"
    ):
        return False
    source = next(
        (node for node in plan.nodes if node.node_id == edge.source_node_id),
        None,
    )
    target = next(
        (node for node in plan.nodes if node.node_id == edge.target_node_id),
        None,
    )
    return bool(
        source is not None
        and target is not None
        and source.program == "orca"
        and target.program == "orca"
        and target.stage == "ts"
    )


@dataclass(frozen=True)
class FrozenWorkflowApprovalV1:
    """Exact host-owned approval boundary for a scientific workflow.

    Future inputs are authorized only through the producer data edges already
    frozen in ``producer_edge_sha256s``.  Consumption is an event-store action;
    the immutable approval record never rewrites its own state.
    """

    schema_version: str
    approval_id: str
    workflow_id: str
    plan_sha256: str
    materialized_workflow_sha256: str
    task_spec_sha256: str
    scientific_identity_sha256: str
    resource_sha256: str
    environment_identity_sha256s: tuple[str, ...]
    approved_node_ids: tuple[str, ...]
    producer_edge_sha256s: tuple[str, ...]
    stationary_point_policy_sha256: str
    status: str
    approval_sha256: str
    materialized_preview_bindings: tuple[
        FrozenMaterializedNodePreviewV1, ...
    ] = ()
    producer_edge_rules: tuple[FrozenProducerEdgeRuleV1, ...] = ()
    admission_sha256: str = ""
    #: Digest of the analysis chain the same approval covers; empty when the
    #: review carried none, preserving every pre-existing approval digest.
    scientific_toolchain_plan_sha256: str = ""

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.frozen-workflow-approval.v1":
            raise ContractError("unsupported frozen workflow approval schema")
        require_identifier(self.approval_id, "approval_id")
        require_identifier(self.workflow_id, "workflow_id")
        for name, digest in (
            ("plan_sha256", self.plan_sha256),
            (
                "materialized_workflow_sha256",
                self.materialized_workflow_sha256,
            ),
            ("task_spec_sha256", self.task_spec_sha256),
            ("scientific_identity_sha256", self.scientific_identity_sha256),
            ("resource_sha256", self.resource_sha256),
        ):
            require_sha256(digest, name)
        for name, values in (
            (
                "environment_identity_sha256s",
                self.environment_identity_sha256s,
            ),
            ("producer_edge_sha256s", self.producer_edge_sha256s),
        ):
            if values != tuple(sorted(set(values))):
                raise ContractError(f"{name} must be sorted and unique")
            for digest in values:
                require_sha256(digest, name)
        if not self.environment_identity_sha256s:
            raise ContractError("approval requires environment evidence")
        if self.approved_node_ids != tuple(
            sorted(set(self.approved_node_ids))
        ):
            raise ContractError("approved_node_ids must be sorted and unique")
        if not self.approved_node_ids:
            raise ContractError("approval requires at least one node")
        for node_id in self.approved_node_ids:
            require_identifier(node_id, "approved_node_id")
        if self.stationary_point_policy_sha256:
            require_sha256(
                self.stationary_point_policy_sha256,
                "stationary_point_policy_sha256",
            )
        if self.scientific_toolchain_plan_sha256:
            require_sha256(
                self.scientific_toolchain_plan_sha256,
                "scientific_toolchain_plan_sha256",
            )
        if self.status != "approved":
            raise ContractError("frozen workflow approval must be approved")
        legacy_body = {
            "schema_version": self.schema_version,
            "approval_id": self.approval_id,
            "workflow_id": self.workflow_id,
            "plan_sha256": self.plan_sha256,
            "materialized_workflow_sha256": self.materialized_workflow_sha256,
            "task_spec_sha256": self.task_spec_sha256,
            "scientific_identity_sha256": self.scientific_identity_sha256,
            "resource_sha256": self.resource_sha256,
            "environment_identity_sha256s": self.environment_identity_sha256s,
            "approved_node_ids": self.approved_node_ids,
            "producer_edge_sha256s": self.producer_edge_sha256s,
            "stationary_point_policy_sha256": (
                self.stationary_point_policy_sha256
            ),
            "status": self.status,
        }
        if self.scientific_toolchain_plan_sha256:
            legacy_body["scientific_toolchain_plan_sha256"] = (
                self.scientific_toolchain_plan_sha256
            )
        has_admission = bool(
            self.materialized_preview_bindings
            or self.producer_edge_rules
            or self.admission_sha256
        )
        if not has_admission:
            if self.approval_sha256 != canonical_sha256(legacy_body):
                raise ContractError("frozen workflow approval digest mismatch")
            return
        preview_ids = tuple(
            item.node_id for item in self.materialized_preview_bindings
        )
        if preview_ids != tuple(sorted(set(preview_ids))):
            raise ContractError(
                "materialized preview bindings must be sorted and unique"
            )
        rule_keys = tuple(
            (item.target_node_id, item.consumer_input_id)
            for item in self.producer_edge_rules
        )
        if rule_keys != tuple(sorted(set(rule_keys))):
            raise ContractError(
                "producer edge rules must be sorted and unique"
            )
        edge_digests = tuple(
            sorted(
                item.scientific_edge_sha256
                for item in self.producer_edge_rules
            )
        )
        if edge_digests != self.producer_edge_sha256s:
            raise ContractError(
                "producer edge rules differ from approved edges"
            )
        covered = set(preview_ids).union(
            item.target_node_id for item in self.producer_edge_rules
        )
        if covered != set(self.approved_node_ids):
            raise ContractError("workflow admission does not cover every node")
        preview_binding_bodies = tuple(
            _frozen_preview_binding_body(item)
            for item in self.materialized_preview_bindings
        )
        admission_body = {
            "schema_version": "chemsmart.frozen-workflow-admission.v1",
            "plan_sha256": self.plan_sha256,
            "materialized_workflow_sha256": self.materialized_workflow_sha256,
            "materialized_preview_bindings": preview_binding_bodies,
            "producer_edge_rules": self.producer_edge_rules,
        }
        if self.admission_sha256 != canonical_sha256(admission_body):
            raise ContractError("frozen workflow admission digest mismatch")
        body = {
            **legacy_body,
            "materialized_preview_bindings": preview_binding_bodies,
            "producer_edge_rules": self.producer_edge_rules,
            "admission_sha256": self.admission_sha256,
        }
        if self.approval_sha256 != canonical_sha256(body):
            raise ContractError("frozen workflow approval digest mismatch")

    @property
    def execution_admissible(self) -> bool:
        """Whether this approval carries the complete V2 admission plane."""

        return bool(self.admission_sha256)

    def preview_binding(
        self, node_id: str
    ) -> FrozenMaterializedNodePreviewV1 | None:
        matches = tuple(
            item
            for item in self.materialized_preview_bindings
            if item.node_id == node_id
        )
        if len(matches) > 1:  # pragma: no cover - constructor prevents this
            raise ContractError("multiple frozen previews exist for one node")
        return matches[0] if matches else None

    def producer_rules_for(
        self, node_id: str
    ) -> tuple[FrozenProducerEdgeRuleV1, ...]:
        return tuple(
            item
            for item in self.producer_edge_rules
            if item.target_node_id == node_id
        )


def build_frozen_workflow_approval(
    *,
    approval_id: str,
    plan: ScientificWorkflowPlanV2,
    materialized_workflow: MaterializedWorkflowV1,
    resources: ExecutionResourceSpecV1,
    environment_identity_sha256s: Sequence[str],
    future_node_environment_identity_sha256s: Mapping[str, str] | None = None,
    environment_identity_by_receipt: Mapping[str, str] | None = None,
    stationary_point_policy: StationaryPointValidationPolicyV1 | None = None,
    non_executable_node_ids: Sequence[str] = (),
    scientific_toolchain_plan: ScientificToolchainPlanV1 | None = None,
) -> FrozenWorkflowApprovalV1:
    """Freeze one plan without pretending future data artifacts exist.

    ``non_executable_node_ids`` names plan stages this release cannot execute.
    They receive no execution admission or frozen future input, even when a
    preview-only root stage was materialized before its release support was
    classified.
    """

    if materialized_workflow.workflow_id != plan.workflow_id:
        raise ContractError("materialized workflow belongs to another plan")
    if materialized_workflow.plan_sha256 != plan.plan_sha256:
        raise ContractError("materialized workflow plan digest differs")
    if materialized_workflow.task_spec_sha256 != plan.task_spec_sha256:
        raise ContractError("materialized workflow task identity differs")
    if materialized_workflow.scientific_identity_sha256 != (
        plan.scientific_identity_sha256
    ):
        raise ContractError(
            "materialized workflow scientific identity differs"
        )
    if materialized_workflow.resource_sha256 not in (
        resources.resource_sha256,
        PREVIEW_RESOURCE_SHA256,
    ):
        raise ContractError("materialized workflow resource binding differs")
    # Preview materialization may carry the preview sentinel.  Human review
    # binds the concrete resources, while a real-profile materialization must
    # already match them exactly.
    planned_by_id = {node.node_id: node for node in plan.nodes}
    non_executable = set(str(item) for item in non_executable_node_ids)
    if not non_executable.issubset(planned_by_id):
        raise ContractError(
            "a non-executable node must belong to the approved plan"
        )
    for node_id in non_executable:
        if planned_by_id[node_id].support_state != "blocked_unsupported":
            raise ContractError(
                "a non-executable approval node must be blocked unsupported"
            )
    executed = set(planned_by_id) - non_executable
    if not executed:
        raise ContractError("an approval requires an executable node")
    if any(
        edge.source_node_id in non_executable
        and edge.target_node_id in executed
        for edge in plan.edges
    ):
        raise ContractError(
            "an executed node cannot consume a non-executable node's output"
        )
    # A non-executable stage is not a future producer target the approval has
    # to freeze an input rule for.
    unresolved = (
        set(materialized_workflow.unresolved_node_ids) - non_executable
    )
    data_edges = tuple(
        edge
        for edge in plan.edges
        if edge.edge_kind == "data"
        and edge.source_node_id not in non_executable
        and edge.target_node_id not in non_executable
    )
    data_targets = {edge.target_node_id for edge in data_edges}
    if len(data_targets) != len(data_edges):
        raise ContractError(
            "one CLI invocation cannot admit multiple future molecular inputs"
        )
    if unresolved != data_targets:
        raise ContractError(
            "future producer-dependent nodes must be exactly the unresolved set"
        )
    preview_bindings = tuple(
        sorted(
            (
                _frozen_preview_binding(node)
                for node in materialized_workflow.nodes
                if node.node_id not in non_executable
            ),
            key=lambda item: item.node_id,
        )
    )
    normalized_environment_receipts = tuple(
        sorted(set(environment_identity_sha256s))
    )
    # The materialization records the environment *receipt* each preview
    # observed, while the approval pins environment *identities*.  The two are
    # different digests for the same machine, so the completeness checks below
    # need a translation the caller supplies from the receipt bodies it holds.
    # Defaulting to the digest itself keeps callers that pin receipts working
    # unchanged, since there identity and receipt are the same value.
    identity_by_receipt = dict(environment_identity_by_receipt or {})

    def _identity_of(receipt_sha256: str) -> str:
        return identity_by_receipt.get(receipt_sha256, receipt_sha256)

    future_environment_bindings = dict(
        future_node_environment_identity_sha256s or {}
    )
    if data_targets and not future_environment_bindings:
        if len(normalized_environment_receipts) != 1:
            raise ContractError(
                "future nodes require exact node-specific environment evidence"
            )
        future_environment_bindings = {
            node_id: normalized_environment_receipts[0]
            for node_id in data_targets
        }
    if set(future_environment_bindings) != data_targets:
        raise ContractError(
            "future environment bindings must cover exactly unresolved nodes"
        )
    if {
        _identity_of(value) for value in future_environment_bindings.values()
    }.difference(normalized_environment_receipts):
        raise ContractError(
            "future node environment is absent from approval evidence"
        )
    producer_rules = tuple(
        sorted(
            (
                _frozen_producer_edge_rule(
                    plan,
                    edge,
                    environment_receipt_sha256=(
                        future_environment_bindings[edge.target_node_id]
                    ),
                )
                for edge in data_edges
            ),
            key=lambda item: (item.target_node_id, item.consumer_input_id),
        )
    )
    if {
        _identity_of(item.environment_receipt_sha256)
        for item in preview_bindings
    }.difference(normalized_environment_receipts):
        raise ContractError(
            "frozen preview environment is absent from approval evidence"
        )
    if stationary_point_policy is not None:
        if stationary_point_policy.task_spec_sha256 != plan.task_spec_sha256:
            raise ContractError(
                "stationary point policy belongs to another task"
            )
        if stationary_point_policy.hessian_node_id not in planned_by_id:
            raise ContractError(
                "stationary point policy names an unknown node"
            )
        if stationary_point_policy.hessian_node_id not in executed:
            raise ContractError(
                "execution stationary point policy requires an executable "
                "Hessian node"
            )
    producer_digests = tuple(
        sorted(item.scientific_edge_sha256 for item in producer_rules)
    )
    preview_binding_bodies = tuple(
        _frozen_preview_binding_body(item) for item in preview_bindings
    )
    admission_body = {
        "schema_version": "chemsmart.frozen-workflow-admission.v1",
        "plan_sha256": plan.plan_sha256,
        "materialized_workflow_sha256": (
            materialized_workflow.materialized_sha256
        ),
        "materialized_preview_bindings": preview_binding_bodies,
        "producer_edge_rules": producer_rules,
    }
    body = {
        "schema_version": "chemsmart.frozen-workflow-approval.v1",
        "approval_id": require_identifier(approval_id, "approval_id"),
        "workflow_id": plan.workflow_id,
        "plan_sha256": plan.plan_sha256,
        "materialized_workflow_sha256": (
            materialized_workflow.materialized_sha256
        ),
        "task_spec_sha256": plan.task_spec_sha256,
        "scientific_identity_sha256": plan.scientific_identity_sha256,
        "resource_sha256": resources.resource_sha256,
        "environment_identity_sha256s": normalized_environment_receipts,
        "approved_node_ids": tuple(
            sorted(
                node.node_id for node in plan.nodes if node.node_id in executed
            )
        ),
        "producer_edge_sha256s": producer_digests,
        "stationary_point_policy_sha256": (
            stationary_point_policy.policy_sha256
            if stationary_point_policy is not None
            else ""
        ),
        "status": "approved",
        "materialized_preview_bindings": preview_bindings,
        "producer_edge_rules": producer_rules,
        "admission_sha256": canonical_sha256(admission_body),
    }
    if scientific_toolchain_plan is not None:
        body["scientific_toolchain_plan_sha256"] = (
            scientific_toolchain_plan.plan_sha256
        )
    digest_body = {
        **body,
        "materialized_preview_bindings": preview_binding_bodies,
    }
    return FrozenWorkflowApprovalV1(
        **body, approval_sha256=canonical_sha256(digest_body)
    )


@dataclass(frozen=True)
class WorkflowEnvironmentBindingV1:
    """The exact observed execution environment assigned to one DAG node."""

    node_id: str
    program: str
    engine: str
    environment_receipt_sha256: str
    environment_identity_sha256: str

    def __post_init__(self) -> None:
        require_identifier(self.node_id, "node_id")
        require_identifier(self.program, "program")
        require_identifier(self.engine, "engine")
        require_sha256(
            self.environment_receipt_sha256,
            "environment_receipt_sha256",
        )
        require_sha256(
            self.environment_identity_sha256,
            "environment_identity_sha256",
        )


_REVIEW_SECRET_FIELD_PARTS = frozenset(
    {
        "api_key",
        "credential",
        "password",
        "private_key",
        "secret",
        "token",
    }
)


def execution_path_placeholder(role: str, sha256: str) -> str:
    """Return a human-readable, digest-bound stand-in for one approved path."""

    normalized_role = require_identifier(role, "execution path role")
    digest = require_sha256(sha256, f"{normalized_role} sha256")
    return f"<{normalized_role}:sha256={digest}>"


def execution_server_profile_sha256(
    *, resources: ExecutionResourceSpecV1, scratch_root: str
) -> str:
    """Bind the generated server profile to reviewed resources and scratch."""

    scratch = str(scratch_root).strip()
    if not scratch:
        raise ContractError("execution server profile requires a scratch root")
    return canonical_sha256(
        {
            "schema_version": "chemsmart.execution-server-profile.v1",
            "resource_sha256": resources.resource_sha256,
            "scratch_root_sha256": canonical_sha256(
                {"scratch_root": scratch}
            ),
        }
    )


def _looks_absolute_path(value: str) -> bool:
    return Path(value).is_absolute() or PureWindowsPath(value).is_absolute()


def _path_free_review_value(value: Any, *, field_name: str = "") -> Any:
    """Preserve reviewed science while withholding secrets and host paths."""

    normalized_field = str(field_name).strip().lower().replace("-", "_")
    if any(part in normalized_field for part in _REVIEW_SECRET_FIELD_PARTS):
        return execution_path_placeholder(
            "redacted-value",
            canonical_sha256({"field": normalized_field, "redacted": True}),
        )
    if isinstance(value, Mapping):
        return {
            str(key): _path_free_review_value(item, field_name=str(key))
            for key, item in sorted(value.items(), key=lambda pair: str(pair[0]))
        }
    if isinstance(value, (list, tuple)):
        return tuple(
            _path_free_review_value(item, field_name=field_name)
            for item in value
        )
    if isinstance(value, str) and _looks_absolute_path(value):
        return execution_path_placeholder(
            "private-path",
            canonical_sha256({"field": normalized_field, "value": value}),
        )
    return canonical_data(value)


def path_free_review_record(value: Mapping[str, Any]) -> dict[str, Any]:
    """Create the canonical path-free public projection of a reviewed record."""

    projected = _path_free_review_value(value)
    if not isinstance(projected, dict):  # pragma: no cover - Mapping narrows this
        raise ContractError("review record must project to an object")
    return projected


def normalized_project_settings_review(
    settings: Mapping[str, Any] | Sequence[tuple[str, Any]],
) -> tuple[str, str]:
    """Render effective loader settings without secrets or private paths."""

    projected = path_free_review_record(dict(settings))
    text = canonical_json(projected)
    return text, canonical_sha256(
        {"normalized_project_settings_text": text}
    )


def environment_review_summary(receipt: Any) -> dict[str, Any]:
    """Project execution-relevant environment facts without host paths."""

    status = getattr(receipt, "status", "")
    summary = {
        "program": getattr(receipt, "program", ""),
        "engine": getattr(receipt, "engine", ""),
        "status": getattr(status, "value", str(status)),
        "target_kind": getattr(receipt, "target_kind", ""),
        "locator": getattr(receipt, "locator", ""),
        "observed_version": getattr(receipt, "observed_version", ""),
        "observed_location_sha256": getattr(
            receipt, "observed_location_sha256", ""
        ),
        "compute_interpreter_sha256": getattr(
            receipt, "compute_interpreter_sha256", ""
        ),
        "compute_evidence_sha256": getattr(
            receipt, "compute_evidence_sha256", ""
        ),
        "dependency_versions": tuple(
            tuple(item)
            for item in getattr(receipt, "dependency_versions", ())
        ),
        "solver_evidence": tuple(
            tuple(item) for item in getattr(receipt, "solver_evidence", ())
        ),
        "gpu_evidence": tuple(
            tuple(item) for item in getattr(receipt, "gpu_evidence", ())
        ),
        "observation_method": getattr(receipt, "observation_method", ""),
    }
    # Reviews cross a JSON boundary before the provider-free executor reloads
    # them.  Normalize tuples here as JSON arrays so the review builder and
    # the execution-time comparison see the same value shape as well as the
    # same scientific environment facts.
    return canonical_data(path_free_review_record(summary))


def build_real_execution_argv(
    *,
    compiled_argv: Sequence[str],
    command_path: Sequence[str],
    resources: ExecutionResourceSpecV1,
    server: str = "",
) -> tuple[str, ...]:
    """Materialize the exact OS process argv from reviewed input."""

    path = tuple(str(item) for item in command_path)
    argv = tuple(str(item) for item in compiled_argv)
    if len(path) < 2 or path[0] != "run":
        raise ContractError("local execution requires a run command")
    program = path[1]
    try:
        program_index = argv.index(program, 2)
    except ValueError as exc:
        raise ContractError("compiled argv does not contain its program") from exc
    root = [
        sys.executable,
        "-m",
        "chemsmart",
        "run",
        "--no-fake",
        "--delete-scratch",
    ]
    root.append(
        "--no-scratch"
        if resources.scratch_policy == "none"
        else "--scratch"
    )
    root.extend(("--num-cores", str(resources.cores)))
    root.extend(("--num-gpus", str(resources.gpu_count)))
    memory = (
        str(int(resources.memory_gb))
        if resources.memory_gb.is_integer()
        else str(resources.memory_gb)
    )
    root.extend(("--mem-gb", memory))
    if server:
        root.extend(("--server", str(server)))
    return tuple(root + list(argv[program_index:]))


def project_real_execution_argv(
    argv: Sequence[str],
    *,
    path_bindings: Mapping[str, tuple[str, str]],
) -> tuple[str, ...]:
    """Replace every approved path token and reject any unreviewed path."""

    replacements: dict[str, str] = {}
    for raw_path, (role, digest) in path_bindings.items():
        token = str(raw_path)
        replacement = execution_path_placeholder(role, digest)
        previous = replacements.setdefault(token, replacement)
        if previous != replacement:
            raise ContractError("one execution path has conflicting approvals")
    projected = tuple(replacements.get(str(item), str(item)) for item in argv)
    unbound = tuple(item for item in projected if _looks_absolute_path(item))
    if unbound:
        raise ContractError(
            "real execution argv contains an unreviewed absolute path"
        )
    return projected


def real_execution_argv_sha256(argv: Sequence[str]) -> str:
    """Digest the exact path-neutral execution argv shown to the reviewer."""

    return canonical_sha256({"real_execution_argv": tuple(argv)})


def _require_path_free_review_value(value: Any, field_name: str) -> None:
    if isinstance(value, Mapping):
        for key, item in value.items():
            _require_path_free_review_value(item, f"{field_name}.{key}")
    elif isinstance(value, (list, tuple)):
        for item in value:
            _require_path_free_review_value(item, field_name)
    elif isinstance(value, str) and _looks_absolute_path(value):
        raise ContractError(f"{field_name} cannot expose an absolute host path")


@dataclass(frozen=True)
class WorkflowExecutionNodeReviewV1:
    """Self-contained human projection of one proposed real engine launch."""

    schema_version: str
    node_id: str
    stage: str
    program: str
    engine: str
    molecular_identity: dict[str, Any]
    molecular_identity_sha256: str
    project_artifact_sha256: str
    project_settings_sha256: str
    project_settings_text: str
    project_settings_text_sha256: str
    real_execution_argv: tuple[str, ...]
    real_execution_argv_sha256: str
    execution_resources: dict[str, Any]
    resource_sha256: str
    environment_summary: dict[str, Any]
    server_profile_sha256: str
    environment_receipt_sha256: str
    environment_identity_sha256: str
    projection_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.workflow-execution-node-review.v1":
            raise ContractError("unsupported workflow node review schema")
        for name, value in (
            ("node_id", self.node_id),
            ("stage", self.stage),
            ("program", self.program),
            ("engine", self.engine),
        ):
            require_identifier(value, name)
        object.__setattr__(
            self, "molecular_identity", canonical_data(self.molecular_identity)
        )
        object.__setattr__(
            self, "execution_resources", canonical_data(self.execution_resources)
        )
        object.__setattr__(
            self, "environment_summary", canonical_data(self.environment_summary)
        )
        object.__setattr__(
            self, "real_execution_argv", tuple(self.real_execution_argv)
        )
        for name, digest in (
            ("molecular_identity_sha256", self.molecular_identity_sha256),
            ("project_artifact_sha256", self.project_artifact_sha256),
            ("project_settings_sha256", self.project_settings_sha256),
            ("project_settings_text_sha256", self.project_settings_text_sha256),
            ("real_execution_argv_sha256", self.real_execution_argv_sha256),
            ("resource_sha256", self.resource_sha256),
            ("server_profile_sha256", self.server_profile_sha256),
            ("environment_receipt_sha256", self.environment_receipt_sha256),
            ("environment_identity_sha256", self.environment_identity_sha256),
        ):
            require_sha256(digest, name)
        _require_path_free_review_value(
            self.molecular_identity, "molecular_identity"
        )
        _require_path_free_review_value(
            self.execution_resources, "execution_resources"
        )
        _require_path_free_review_value(
            self.environment_summary, "environment_summary"
        )
        _require_path_free_review_value(
            self.real_execution_argv, "real_execution_argv"
        )
        if self.molecular_identity_sha256 != canonical_sha256(
            self.molecular_identity
        ):
            raise ContractError("molecular identity review digest mismatch")
        expected_settings_text_sha256 = canonical_sha256(
            {"normalized_project_settings_text": self.project_settings_text}
        )
        try:
            parsed_settings = json.loads(self.project_settings_text)
        except json.JSONDecodeError as exc:
            raise ContractError(
                "project settings review must be canonical JSON"
            ) from exc
        if not isinstance(parsed_settings, Mapping) or canonical_json(
            parsed_settings
        ) != self.project_settings_text:
            raise ContractError("project settings review must be a canonical object")
        _require_path_free_review_value(
            parsed_settings, "project_settings_text"
        )
        if self.project_settings_text_sha256 != expected_settings_text_sha256:
            raise ContractError("project settings review digest mismatch")
        if self.real_execution_argv_sha256 != real_execution_argv_sha256(
            self.real_execution_argv
        ):
            raise ContractError("real execution argv review digest mismatch")
        if self.projection_sha256 != canonical_sha256(self._body()):
            raise ContractError("workflow node review digest mismatch")

    def _body(self) -> dict[str, Any]:
        return {
            key: value
            for key, value in self.__dict__.items()
            if key != "projection_sha256"
        }


def build_workflow_execution_node_review(
    *,
    node_id: str,
    stage: str,
    program: str,
    engine: str,
    molecular_identity: Mapping[str, Any],
    project_artifact_sha256: str,
    project_settings_sha256: str,
    project_settings: Mapping[str, Any] | Sequence[tuple[str, Any]],
    real_execution_argv: Sequence[str],
    execution_resources: ExecutionResourceSpecV1,
    environment_summary: Mapping[str, Any],
    server_profile_sha256: str,
    environment_receipt_sha256: str,
    environment_identity_sha256: str,
) -> WorkflowExecutionNodeReviewV1:
    """Build and digest one path-free node projection for human review."""

    identity = path_free_review_record(dict(molecular_identity))
    settings_text, settings_text_sha256 = normalized_project_settings_review(
        project_settings
    )
    argv = tuple(str(item) for item in real_execution_argv)
    body = {
        "schema_version": "chemsmart.workflow-execution-node-review.v1",
        "node_id": node_id,
        "stage": stage,
        "program": program,
        "engine": engine,
        "molecular_identity": identity,
        "molecular_identity_sha256": canonical_sha256(identity),
        "project_artifact_sha256": project_artifact_sha256,
        "project_settings_sha256": project_settings_sha256,
        "project_settings_text": settings_text,
        "project_settings_text_sha256": settings_text_sha256,
        "real_execution_argv": argv,
        "real_execution_argv_sha256": real_execution_argv_sha256(argv),
        "execution_resources": canonical_data(execution_resources),
        "resource_sha256": execution_resources.resource_sha256,
        "environment_summary": path_free_review_record(environment_summary),
        "server_profile_sha256": server_profile_sha256,
        "environment_receipt_sha256": environment_receipt_sha256,
        "environment_identity_sha256": environment_identity_sha256,
    }
    return WorkflowExecutionNodeReviewV1(
        **body, projection_sha256=canonical_sha256(body)
    )


@dataclass(frozen=True)
class WorkflowExecutionReviewV1:
    """Inert, exact review packet emitted after planning and safe preview.

    This object intentionally contains no approval.  It binds the scientific
    plan, its command materialization, resources, operating bounds, and the
    environment identities observed during preview so that a human can review
    one digest and later approve exactly those bytes.
    """

    schema_version: str
    request: WorkflowApprovalRequestV1
    scientific_plan: ScientificWorkflowPlanV2
    materialized_workflow: MaterializedWorkflowV1
    execution_resources: ExecutionResourceSpecV1
    execution_envelope: dict[str, Any]
    environment_bindings: tuple[WorkflowEnvironmentBindingV1, ...]
    node_reviews: tuple[WorkflowExecutionNodeReviewV1, ...]
    stationary_point_policy: StationaryPointValidationPolicyV1 | None
    status: str
    review_sha256: str
    #: Plan nodes this release cannot execute.  They stay in the scientific
    #: plan, are displayed to the human, and are never approved or launched.
    non_executable_node_ids: tuple[str, ...] = ()
    #: The typed analysis chain the same /approve covers.  Until this field
    #: existed the analysis nodes were RAM-only in the planning session: the
    #: human approved a packet that never displayed them and the executor had
    #: nothing to walk, so every validation rule a session wrote remained
    #: intent.  Present, it is digest-bound with the rest of the review.
    scientific_toolchain_plan: ScientificToolchainPlanV1 | None = None

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.workflow-execution-review.v1":
            raise ContractError("unsupported workflow execution review schema")
        if self.scientific_toolchain_plan is not None:
            toolchain = self.scientific_toolchain_plan
            if toolchain.workflow_id != self.scientific_plan.workflow_id:
                raise ContractError(
                    "review toolchain plan belongs to another workflow"
                )
            if set(toolchain.calculation_node_ids) != {
                node.node_id for node in self.scientific_plan.nodes
            }:
                raise ContractError(
                    "review toolchain plan covers different calculation nodes"
                )
        object.__setattr__(
            self,
            "non_executable_node_ids",
            tuple(self.non_executable_node_ids),
        )
        if self.status != "unapproved":
            raise ContractError("an execution review is inert and unapproved")
        object.__setattr__(
            self, "execution_envelope", canonical_data(self.execution_envelope)
        )
        object.__setattr__(
            self, "environment_bindings", tuple(self.environment_bindings)
        )
        object.__setattr__(self, "node_reviews", tuple(self.node_reviews))
        if self.request.workflow_id != self.scientific_plan.workflow_id:
            raise ContractError("review request belongs to another workflow")
        if self.request.workflow_sha256 != self.scientific_plan.plan_sha256:
            raise ContractError("review request differs from scientific plan")
        if self.request.task_spec_sha256 != self.scientific_plan.task_spec_sha256:
            raise ContractError("review request differs from task specification")
        if self.materialized_workflow.workflow_id != self.scientific_plan.workflow_id:
            raise ContractError("review materialization belongs to another workflow")
        if self.materialized_workflow.plan_sha256 != self.scientific_plan.plan_sha256:
            raise ContractError("review materialization differs from scientific plan")
        if self.materialized_workflow.task_spec_sha256 != (
            self.scientific_plan.task_spec_sha256
        ):
            raise ContractError("review materialization differs from task specification")
        if self.request.resource_sha256 != self.execution_resources.resource_sha256:
            raise ContractError("review request differs from execution resources")
        if self.materialized_workflow.resource_sha256 not in {
            PREVIEW_RESOURCE_SHA256,
            self.execution_resources.resource_sha256,
        }:
            raise ContractError("review materialization differs from resources")
        planned_by_id = {
            node.node_id: node for node in self.scientific_plan.nodes
        }
        non_executable = self.non_executable_node_ids
        if non_executable != tuple(sorted(set(non_executable))):
            raise ContractError(
                "non-executable node ids must be sorted and unique"
            )
        if not set(non_executable).issubset(planned_by_id):
            raise ContractError("a non-executable node must belong to the plan")
        for node_id in non_executable:
            if planned_by_id[node_id].support_state != "blocked_unsupported":
                raise ContractError(
                    "a non-executable review node must declare why it is unsupported"
                )
        executed_ids = set(planned_by_id) - set(non_executable)
        if not executed_ids:
            raise ContractError("an execution review requires an executable node")
        for edge in self.scientific_plan.edges:
            if (
                edge.source_node_id in set(non_executable)
                and edge.target_node_id in executed_ids
            ):
                raise ContractError(
                    "an executed node cannot consume a non-executable node's output"
                )
        node_ids = tuple(item.node_id for item in self.environment_bindings)
        if node_ids != tuple(sorted(set(node_ids))):
            raise ContractError(
                "review environment bindings must be sorted and node-unique"
            )
        if set(node_ids) != executed_ids:
            raise ContractError("review environments must cover every executed node")
        review_node_ids = tuple(item.node_id for item in self.node_reviews)
        if review_node_ids != tuple(sorted(set(review_node_ids))):
            raise ContractError("node reviews must be sorted and node-unique")
        if set(review_node_ids) != executed_ids:
            raise ContractError("node reviews must cover every executed node")
        binding_by_id = {
            item.node_id: item for item in self.request.node_bindings
        }
        environment_by_id = {
            item.node_id: item for item in self.environment_bindings
        }
        if set(binding_by_id) != executed_ids:
            raise ContractError("review node bindings are incomplete")
        for item in self.node_reviews:
            planned = planned_by_id[item.node_id]
            binding = binding_by_id[item.node_id]
            environment = environment_by_id[item.node_id]
            if (item.stage, item.program, item.engine) != (
                planned.stage,
                planned.program,
                planned.engine,
            ):
                raise ContractError("node review differs from scientific plan")
            if (
                item.project_artifact_sha256
                != binding.project_artifact_sha256
                or item.project_settings_sha256 != binding.settings_sha256
                or item.resource_sha256
                != self.execution_resources.resource_sha256
            ):
                raise ContractError("node review differs from approved inputs")
            if (
                item.environment_receipt_sha256
                != environment.environment_receipt_sha256
                or item.environment_identity_sha256
                != environment.environment_identity_sha256
            ):
                raise ContractError("node review differs from environment binding")
            if (
                item.environment_summary.get("program") != item.program
                or item.environment_summary.get("engine") != item.engine
            ):
                raise ContractError("node review environment summary differs")
            molecular = item.molecular_identity
            if (
                molecular.get("charge") != binding.charge
                or molecular.get("multiplicity") != binding.multiplicity
            ):
                raise ContractError("node review molecular state differs")
            coordinate = molecular.get("coordinate_identity")
            if not isinstance(coordinate, Mapping):
                raise ContractError("node review lacks coordinate identity")
            if binding.input_mode == "initial" and (
                coordinate.get("kind") != "exact-input-artifact"
                or coordinate.get("geometry_artifact_sha256")
                != binding.initial_artifact_sha256
            ):
                raise ContractError("node review initial geometry differs")
            if binding.input_mode == "producer" and (
                coordinate.get("kind") != "validated-producer-output"
                or coordinate.get("producer_edge_sha256")
                != binding.producer_edge_sha256
            ):
                raise ContractError("node review producer geometry differs")
        scratch_root = str(self.execution_envelope.get("scratch_root", ""))
        expected_server_profile = execution_server_profile_sha256(
            resources=self.execution_resources,
            scratch_root=scratch_root,
        )
        if any(
            item.server_profile_sha256 != expected_server_profile
            for item in self.node_reviews
        ):
            raise ContractError("node review server profile differs")
        if self.stationary_point_policy is not None:
            if self.stationary_point_policy.task_spec_sha256 != (
                self.scientific_plan.task_spec_sha256
            ):
                raise ContractError("stationary-point policy belongs to another task")
            if self.stationary_point_policy.hessian_node_id not in executed_ids:
                raise ContractError(
                    "execution stationary-point policy requires an executable "
                    "Hessian node"
                )
        if self.review_sha256 != canonical_sha256(self._body()):
            raise ContractError("workflow execution review digest mismatch")

    def _body(self) -> dict[str, Any]:
        body = {
            key: value
            for key, value in self.__dict__.items()
            if key != "review_sha256"
        }
        # These additive v1 fields were introduced after review packets
        # existed.  Preserve their canonical body when the semantic is absent.
        if not self.non_executable_node_ids:
            body.pop("non_executable_node_ids", None)
        if self.scientific_toolchain_plan is None:
            body.pop("scientific_toolchain_plan", None)
        return body


def build_workflow_execution_review(
    *,
    request: WorkflowApprovalRequestV1,
    scientific_plan: ScientificWorkflowPlanV2,
    materialized_workflow: MaterializedWorkflowV1,
    execution_resources: ExecutionResourceSpecV1,
    execution_envelope: Mapping[str, Any],
    environment_bindings: Sequence[WorkflowEnvironmentBindingV1],
    node_reviews: Sequence[WorkflowExecutionNodeReviewV1],
    stationary_point_policy: StationaryPointValidationPolicyV1 | None = None,
    non_executable_node_ids: Sequence[str] = (),
    scientific_toolchain_plan: ScientificToolchainPlanV1 | None = None,
) -> WorkflowExecutionReviewV1:
    """Assemble one self-verifying review packet without granting authority.

    ``non_executable_node_ids`` names plan stages this release cannot execute.
    They are displayed with the workflow and never approved.  A nonempty set
    is digest-bound; its absence retains the original v1 canonical body.
    """

    non_executable = tuple(
        sorted(set(str(item) for item in non_executable_node_ids))
    )
    body: dict[str, Any] = {
        "schema_version": "chemsmart.workflow-execution-review.v1",
        "request": request,
        "scientific_plan": scientific_plan,
        "materialized_workflow": materialized_workflow,
        "execution_resources": execution_resources,
        "execution_envelope": canonical_data(dict(execution_envelope)),
        "environment_bindings": tuple(environment_bindings),
        "node_reviews": tuple(node_reviews),
        "stationary_point_policy": stationary_point_policy,
        "status": "unapproved",
    }
    if non_executable:
        body["non_executable_node_ids"] = non_executable
    if scientific_toolchain_plan is not None:
        body["scientific_toolchain_plan"] = scientific_toolchain_plan
    return WorkflowExecutionReviewV1(
        **body, review_sha256=canonical_sha256(body)
    )


@dataclass(frozen=True)
class WorkflowReviewResolutionV1:
    """Durable human resolution of one exact execution review digest."""

    schema_version: str
    resolution_id: str
    review_sha256: str
    decision: str
    actor: str
    approval_id: str
    one_shot: bool
    resolution_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.workflow-review-resolution.v1":
            raise ContractError("unsupported workflow review resolution schema")
        require_identifier(self.resolution_id, "resolution_id")
        require_sha256(self.review_sha256, "review_sha256")
        if self.decision not in {"approve", "deny", "revise", "quit"}:
            raise ContractError("unknown workflow review decision")
        if not str(self.actor).strip():
            raise ContractError("workflow review resolution requires an actor")
        if self.decision == "approve":
            require_identifier(self.approval_id, "approval_id")
            if not self.one_shot:
                raise ContractError("workflow execution approval must be one-shot")
        elif self.approval_id or self.one_shot:
            raise ContractError("a non-approval decision cannot grant authority")
        if self.resolution_sha256 != canonical_sha256(self._body()):
            raise ContractError("workflow review resolution digest mismatch")

    def _body(self) -> dict[str, Any]:
        return {
            key: value
            for key, value in self.__dict__.items()
            if key != "resolution_sha256"
        }


def build_workflow_review_resolution(
    *,
    resolution_id: str,
    review_sha256: str,
    decision: str,
    actor: str,
    approval_id: str = "",
) -> WorkflowReviewResolutionV1:
    normalized_decision = str(decision).strip().lower()
    body = {
        "schema_version": "chemsmart.workflow-review-resolution.v1",
        "resolution_id": require_identifier(resolution_id, "resolution_id"),
        "review_sha256": require_sha256(review_sha256, "review_sha256"),
        "decision": normalized_decision,
        "actor": str(actor).strip(),
        "approval_id": str(approval_id).strip(),
        "one_shot": normalized_decision == "approve",
    }
    return WorkflowReviewResolutionV1(
        **body, resolution_sha256=canonical_sha256(body)
    )


@dataclass(frozen=True)
class WorkflowExecutionApprovalBundleV1:
    """Compound, one-shot input consumed only by ``agent execute``."""

    schema_version: str
    review_sha256: str
    resolution: WorkflowReviewResolutionV1
    workflow_approval: WorkflowExecutionApprovalV1
    frozen_workflow_approval: FrozenWorkflowApprovalV1
    approved_scientific_plan: ScientificWorkflowPlanV2
    approved_materialized_workflow: MaterializedWorkflowV1
    execution_resources: ExecutionResourceSpecV1
    execution_envelope: dict[str, Any]
    approved_environment_identities: tuple[str, ...]
    node_reviews: tuple[WorkflowExecutionNodeReviewV1, ...]
    stationary_point_policy: StationaryPointValidationPolicyV1 | None
    one_shot: bool
    status: str
    bundle_sha256: str
    #: Approved-plan stages this release cannot execute.  The executor never
    #: launches them; they remain in the plan as declared scientific intent.
    non_executable_node_ids: tuple[str, ...] = ()
    #: The typed analysis chain the same approval covers, verbatim from the
    #: reviewed packet; None for every bundle approved before it existed.
    scientific_toolchain_plan: ScientificToolchainPlanV1 | None = None

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.workflow-execution-approval-bundle.v1":
            raise ContractError("unsupported execution approval bundle schema")
        object.__setattr__(
            self,
            "non_executable_node_ids",
            tuple(self.non_executable_node_ids),
        )
        require_sha256(self.review_sha256, "review_sha256")
        object.__setattr__(
            self, "execution_envelope", canonical_data(self.execution_envelope)
        )
        object.__setattr__(
            self,
            "approved_environment_identities",
            tuple(self.approved_environment_identities),
        )
        object.__setattr__(self, "node_reviews", tuple(self.node_reviews))
        if self.resolution.decision != "approve":
            raise ContractError("execution bundle requires an approve resolution")
        if self.resolution.review_sha256 != self.review_sha256:
            raise ContractError("execution bundle resolution targets another review")
        if self.resolution.approval_id != self.workflow_approval.approval_id:
            raise ContractError("execution bundle approval IDs differ")
        frozen = self.frozen_workflow_approval
        approval = self.workflow_approval
        if frozen.approval_id != approval.approval_id:
            raise ContractError("frozen and workflow approval IDs differ")
        if frozen.workflow_id != approval.workflow_id:
            raise ContractError("frozen and workflow approvals differ")
        if frozen.plan_sha256 != self.approved_scientific_plan.plan_sha256:
            raise ContractError("execution bundle plan differs from frozen approval")
        if frozen.materialized_workflow_sha256 != (
            self.approved_materialized_workflow.materialized_sha256
        ):
            raise ContractError(
                "execution bundle materialization differs from frozen approval"
            )
        if approval.resource_sha256 != self.execution_resources.resource_sha256:
            raise ContractError("execution bundle resources differ from approval")
        identities = tuple(sorted(set(self.approved_environment_identities)))
        if identities != self.approved_environment_identities or not identities:
            raise ContractError(
                "approved environment identities must be sorted, unique, and nonempty"
            )
        if identities != frozen.environment_identity_sha256s:
            raise ContractError("execution bundle environments differ from approval")
        planned_by_id = {
            item.node_id: item for item in self.approved_scientific_plan.nodes
        }
        non_executable = self.non_executable_node_ids
        if non_executable != tuple(sorted(set(non_executable))):
            raise ContractError(
                "non-executable node ids must be sorted and unique"
            )
        if not set(non_executable).issubset(planned_by_id):
            raise ContractError(
                "a non-executable node must belong to the approved plan"
            )
        executed_ids = set(planned_by_id) - set(non_executable)
        if not executed_ids:
            raise ContractError("an execution bundle requires an executable node")
        for node_id in non_executable:
            if planned_by_id[node_id].support_state != "blocked_unsupported":
                raise ContractError(
                    "a non-executable bundle node must be blocked unsupported"
                )
        if any(
            edge.source_node_id in set(non_executable)
            and edge.target_node_id in executed_ids
            for edge in self.approved_scientific_plan.edges
        ):
            raise ContractError(
                "an executed node cannot consume a non-executable node's output"
            )
        if set(frozen.approved_node_ids) != executed_ids:
            raise ContractError(
                "frozen approval differs from executable workflow nodes"
            )
        approval_edge_projection = tuple(
            sorted(
                (
                    edge.producer_node_id,
                    edge.consumer_node_id,
                    edge.artifact_kind,
                    edge.selection_rule,
                )
                for edge in approval.producer_edges
            )
        )
        frozen_edge_projection = tuple(
            sorted(
                (
                    edge.source_node_id,
                    edge.target_node_id,
                    edge.artifact_class,
                    edge.selection_rule,
                )
                for edge in frozen.producer_edge_rules
            )
        )
        if approval_edge_projection != frozen_edge_projection:
            raise ContractError(
                "workflow and frozen approval producer edges differ"
            )
        policy_sha256 = (
            self.stationary_point_policy.policy_sha256
            if self.stationary_point_policy is not None
            else ""
        )
        if frozen.stationary_point_policy_sha256 != policy_sha256:
            raise ContractError(
                "execution bundle stationary-point policy differs from approval"
            )
        toolchain_sha256 = (
            self.scientific_toolchain_plan.plan_sha256
            if self.scientific_toolchain_plan is not None
            else ""
        )
        if frozen.scientific_toolchain_plan_sha256 != toolchain_sha256:
            raise ContractError(
                "execution bundle analysis chain differs from approval"
            )
        if self.scientific_toolchain_plan is not None:
            toolchain = self.scientific_toolchain_plan
            if toolchain.workflow_id != self.approved_scientific_plan.workflow_id:
                raise ContractError(
                    "bundle toolchain plan belongs to another workflow"
                )
            if set(toolchain.calculation_node_ids) != set(planned_by_id):
                raise ContractError(
                    "bundle toolchain plan covers different calculation nodes"
                )
        if (
            self.stationary_point_policy is not None
            and self.stationary_point_policy.hessian_node_id not in executed_ids
        ):
            raise ContractError(
                "execution stationary-point policy requires an executable "
                "Hessian node"
            )
        node_ids = tuple(item.node_id for item in self.node_reviews)
        if node_ids != tuple(sorted(set(node_ids))):
            raise ContractError("execution bundle node reviews must be unique")
        if set(node_ids) != executed_ids:
            raise ContractError("execution bundle node reviews are incomplete")
        binding_by_id = {
            item.node_id: item for item in approval.node_bindings
        }
        if set(binding_by_id) != executed_ids:
            raise ContractError("execution bundle node bindings are incomplete")
        for item in self.node_reviews:
            binding = binding_by_id[item.node_id]
            planned = planned_by_id[item.node_id]
            if (
                item.project_artifact_sha256
                != binding.project_artifact_sha256
                or item.project_settings_sha256 != binding.settings_sha256
                or item.resource_sha256
                != self.execution_resources.resource_sha256
                or item.environment_identity_sha256 not in identities
            ):
                raise ContractError("execution bundle node review differs")
            if (item.stage, item.program, item.engine) != (
                planned.stage,
                planned.program,
                planned.engine,
            ):
                raise ContractError("execution bundle command stage differs")
            molecular = item.molecular_identity
            if (
                molecular.get("charge") != binding.charge
                or molecular.get("multiplicity") != binding.multiplicity
            ):
                raise ContractError("execution bundle molecular state differs")
            if (
                item.environment_summary.get("program") != item.program
                or item.environment_summary.get("engine") != item.engine
            ):
                raise ContractError("execution bundle environment summary differs")
            coordinate = molecular.get("coordinate_identity")
            if not isinstance(coordinate, Mapping):
                raise ContractError("execution bundle lacks coordinate identity")
            if binding.input_mode == "initial" and (
                coordinate.get("kind") != "exact-input-artifact"
                or coordinate.get("geometry_artifact_sha256")
                != binding.initial_artifact_sha256
            ):
                raise ContractError("execution bundle initial geometry differs")
            if binding.input_mode == "producer" and (
                coordinate.get("kind") != "validated-producer-output"
                or coordinate.get("producer_edge_sha256")
                != binding.producer_edge_sha256
            ):
                raise ContractError("execution bundle producer geometry differs")
        expected_server_profile = execution_server_profile_sha256(
            resources=self.execution_resources,
            scratch_root=str(self.execution_envelope.get("scratch_root", "")),
        )
        if any(
            item.server_profile_sha256 != expected_server_profile
            for item in self.node_reviews
        ):
            raise ContractError("execution bundle server profile differs")
        if not self.one_shot or self.status != "approved":
            raise ContractError("execution bundle must be approved and one-shot")
        if self.bundle_sha256 != canonical_sha256(self._body()):
            raise ContractError("workflow execution approval bundle digest mismatch")

    def _body(self) -> dict[str, Any]:
        body = {
            key: value
            for key, value in self.__dict__.items()
            if key != "bundle_sha256"
        }
        # These additive v1 fields were introduced after approval bundles
        # existed.  Preserve their canonical body when the semantic is absent.
        if not self.non_executable_node_ids:
            body.pop("non_executable_node_ids", None)
        if self.scientific_toolchain_plan is None:
            body.pop("scientific_toolchain_plan", None)
        return body

    def node_review(self, node_id: str) -> WorkflowExecutionNodeReviewV1:
        """Return the exact human-reviewed projection for one workflow node."""

        matches = tuple(item for item in self.node_reviews if item.node_id == node_id)
        if len(matches) != 1:
            raise ContractError("execution bundle has no unique node review")
        return matches[0]


def approve_workflow_execution_review(
    review: WorkflowExecutionReviewV1,
    *,
    approval_id: str,
    approved_review_sha256: str,
    actor: str,
    resolution_id: str,
) -> WorkflowExecutionApprovalBundleV1:
    """Convert the one reviewed digest into one exact executable bundle."""

    if approved_review_sha256 != review.review_sha256:
        raise ContractError(
            "the reviewed digest does not match this execution review"
        )
    approval = approve_workflow_request(
        review.request,
        approval_id=approval_id,
        approved_request_sha256=review.request.request_sha256,
        resources=review.execution_resources,
    )
    identity_by_receipt = {
        item.environment_receipt_sha256: item.environment_identity_sha256
        for item in review.environment_bindings
    }
    non_executable = set(review.non_executable_node_ids)
    data_targets = {
        edge.target_node_id
        for edge in review.scientific_plan.edges
        if edge.edge_kind == "data"
        and edge.source_node_id not in non_executable
        and edge.target_node_id not in non_executable
    }
    future_environments = {
        item.node_id: item.environment_receipt_sha256
        for item in review.environment_bindings
        if item.node_id in data_targets
    }
    identities = tuple(
        sorted(
            {
                item.environment_identity_sha256
                for item in review.environment_bindings
            }
        )
    )
    frozen = build_frozen_workflow_approval(
        approval_id=approval_id,
        plan=review.scientific_plan,
        materialized_workflow=review.materialized_workflow,
        resources=review.execution_resources,
        environment_identity_sha256s=identities,
        future_node_environment_identity_sha256s=future_environments,
        environment_identity_by_receipt=identity_by_receipt,
        stationary_point_policy=review.stationary_point_policy,
        non_executable_node_ids=review.non_executable_node_ids,
        scientific_toolchain_plan=review.scientific_toolchain_plan,
    )
    resolution = build_workflow_review_resolution(
        resolution_id=resolution_id,
        review_sha256=review.review_sha256,
        decision="approve",
        actor=actor,
        approval_id=approval_id,
    )
    body: dict[str, Any] = {
        "schema_version": "chemsmart.workflow-execution-approval-bundle.v1",
        "review_sha256": review.review_sha256,
        "resolution": resolution,
        "workflow_approval": approval,
        "frozen_workflow_approval": frozen,
        "approved_scientific_plan": review.scientific_plan,
        "approved_materialized_workflow": review.materialized_workflow,
        "execution_resources": review.execution_resources,
        "execution_envelope": review.execution_envelope,
        "approved_environment_identities": identities,
        "node_reviews": review.node_reviews,
        "stationary_point_policy": review.stationary_point_policy,
        "one_shot": True,
        "status": "approved",
    }
    if review.non_executable_node_ids:
        body["non_executable_node_ids"] = review.non_executable_node_ids
    if review.scientific_toolchain_plan is not None:
        body["scientific_toolchain_plan"] = review.scientific_toolchain_plan
    return WorkflowExecutionApprovalBundleV1(
        **body, bundle_sha256=canonical_sha256(body)
    )


def workflow_execution_review_json(review: WorkflowExecutionReviewV1) -> str:
    return canonical_json({"workflow_execution_review": review}) + "\n"


def workflow_execution_approval_bundle_json(
    bundle: WorkflowExecutionApprovalBundleV1,
) -> str:
    return canonical_json({"workflow_execution_approval_bundle": bundle}) + "\n"


@dataclass(frozen=True)
class ValidatedDataEdgeBindingV1:
    """Exact validated producer output selected for one scientific data edge."""

    schema_version: str
    run_id: str
    workflow_id: str
    plan_sha256: str
    approval_sha256: str
    scientific_edge_sha256: str
    producer_rule_sha256: str
    source_node_id: str
    target_node_id: str
    artifact_class: str
    producer_output_id: str
    consumer_input_id: str
    selection_rule: str
    producer_execution_receipt_sha256: str
    producer_validator_receipt_sha256s: tuple[str, ...]
    source_artifact_id: str
    source_artifact_sha256: str
    selected_artifact_id: str
    selected_artifact_sha256: str
    producer_scientific_identity_sha256: str
    consumer_scientific_identity_sha256: str
    atom_order_sha256: str
    positions_sha256: str
    charge: int
    multiplicity: int
    handoff_receipt_sha256: str
    status: str
    receipt_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.validated-data-edge-binding.v1":
            raise ContractError("unsupported validated data edge binding")
        for name, value in (
            ("run_id", self.run_id),
            ("workflow_id", self.workflow_id),
            ("source_node_id", self.source_node_id),
            ("target_node_id", self.target_node_id),
            ("artifact_class", self.artifact_class),
            ("producer_output_id", self.producer_output_id),
            ("consumer_input_id", self.consumer_input_id),
            ("selection_rule", self.selection_rule),
        ):
            require_identifier(value, name)
        for name, digest in (
            ("plan_sha256", self.plan_sha256),
            ("approval_sha256", self.approval_sha256),
            ("scientific_edge_sha256", self.scientific_edge_sha256),
            ("producer_rule_sha256", self.producer_rule_sha256),
            (
                "producer_execution_receipt_sha256",
                self.producer_execution_receipt_sha256,
            ),
            ("source_artifact_sha256", self.source_artifact_sha256),
            ("selected_artifact_sha256", self.selected_artifact_sha256),
            (
                "producer_scientific_identity_sha256",
                self.producer_scientific_identity_sha256,
            ),
            (
                "consumer_scientific_identity_sha256",
                self.consumer_scientific_identity_sha256,
            ),
            ("atom_order_sha256", self.atom_order_sha256),
            ("positions_sha256", self.positions_sha256),
            ("handoff_receipt_sha256", self.handoff_receipt_sha256),
        ):
            require_sha256(digest, name)
        if not self.source_artifact_id or not self.selected_artifact_id:
            raise ContractError("data edge artifacts require stable IDs")
        if (
            self.producer_validator_receipt_sha256s
            != tuple(sorted(set(self.producer_validator_receipt_sha256s)))
            or not self.producer_validator_receipt_sha256s
        ):
            raise ContractError(
                "data edge requires sorted producer validator receipts"
            )
        for digest in self.producer_validator_receipt_sha256s:
            require_sha256(digest, "producer_validator_receipt_sha256")
        if self.multiplicity < 1:
            raise ContractError("data edge multiplicity must be positive")
        if self.status != "validated":
            raise ContractError("data edge binding must be validated")
        body = {
            "schema_version": self.schema_version,
            "run_id": self.run_id,
            "workflow_id": self.workflow_id,
            "plan_sha256": self.plan_sha256,
            "approval_sha256": self.approval_sha256,
            "scientific_edge_sha256": self.scientific_edge_sha256,
            "producer_rule_sha256": self.producer_rule_sha256,
            "source_node_id": self.source_node_id,
            "target_node_id": self.target_node_id,
            "artifact_class": self.artifact_class,
            "producer_output_id": self.producer_output_id,
            "consumer_input_id": self.consumer_input_id,
            "selection_rule": self.selection_rule,
            "producer_execution_receipt_sha256": (
                self.producer_execution_receipt_sha256
            ),
            "producer_validator_receipt_sha256s": (
                self.producer_validator_receipt_sha256s
            ),
            "source_artifact_id": self.source_artifact_id,
            "source_artifact_sha256": self.source_artifact_sha256,
            "selected_artifact_id": self.selected_artifact_id,
            "selected_artifact_sha256": self.selected_artifact_sha256,
            "producer_scientific_identity_sha256": (
                self.producer_scientific_identity_sha256
            ),
            "consumer_scientific_identity_sha256": (
                self.consumer_scientific_identity_sha256
            ),
            "atom_order_sha256": self.atom_order_sha256,
            "positions_sha256": self.positions_sha256,
            "charge": self.charge,
            "multiplicity": self.multiplicity,
            "handoff_receipt_sha256": self.handoff_receipt_sha256,
            "status": self.status,
        }
        if self.receipt_sha256 != canonical_sha256(body):
            raise ContractError("validated data edge binding digest mismatch")


def build_validated_data_edge_binding(
    *,
    run_id: str,
    plan: ScientificWorkflowPlanV2,
    approval: FrozenWorkflowApprovalV1,
    scientific_edge: ScientificWorkflowEdgeV2,
    producer_edge: ProducerEdgeRuleV1,
    producer_invocation: ProgramExecutionInvocationV1,
    producer_receipt: ProgramExecutionReceiptV1,
    handoff: OptimizedGeometryHandoffV1 | ORCAHessianHandoffV1,
    producer_scientific_identity_sha256: str,
    consumer_scientific_identity_sha256: str,
) -> ValidatedDataEdgeBindingV1:
    """Bind one semantically selected producer artifact to its DAG edge."""

    if not approval.execution_admissible:
        raise ContractError("legacy frozen approval cannot admit data edges")
    if approval.plan_sha256 != plan.plan_sha256:
        raise ContractError("data edge plan differs from frozen approval")
    if scientific_edge.edge_kind != "data":
        raise ContractError("validated data binding requires a data edge")
    edge_sha256 = canonical_sha256(scientific_edge)
    rules = tuple(
        item
        for item in approval.producer_edge_rules
        if item.scientific_edge_sha256 == edge_sha256
    )
    if len(rules) != 1:
        raise ContractError("data edge lacks one exact frozen selection rule")
    rule = rules[0]
    observed_roles = (
        producer_edge.producer_node_id,
        producer_edge.consumer_node_id,
        producer_edge.artifact_kind,
        producer_edge.selection_rule,
    )
    expected_roles = (
        rule.source_node_id,
        rule.target_node_id,
        rule.artifact_class,
        rule.selection_rule,
    )
    if observed_roles != expected_roles:
        raise ContractError("producer handoff differs from frozen data edge")
    if (
        handoff.producer_edge_sha256 != producer_edge.edge_sha256
        or handoff.producer_execution_receipt_sha256
        != producer_receipt.receipt_sha256
        or handoff.producer_node_id != rule.source_node_id
        or handoff.consumer_node_id != rule.target_node_id
    ):
        raise ContractError("producer handoff identity differs from producer")
    if not producer_receipt.validated or (
        producer_receipt.node_id != rule.source_node_id
    ):
        raise ContractError("data edge producer is not validated")
    if (
        producer_invocation.node_id != rule.source_node_id
        or producer_receipt.invocation_sha256
        != producer_invocation.invocation_sha256
    ):
        raise ContractError(
            "data edge producer invocation differs from validated execution"
        )
    if not producer_receipt.validator_receipt_sha256s:
        raise ContractError("data edge producer lacks validator evidence")
    if not any(
        artifact.artifact_id == handoff.result_artifact_id
        and artifact.sha256 == handoff.result_artifact_sha256
        for artifact in producer_receipt.output_artifacts
    ):
        raise ContractError("data edge source artifact differs from producer")
    # ``plan.scientific_identity_sha256`` is an aggregate over every external
    # starting geometry in a multi-structure workflow.  It is deliberately not
    # the identity of any one producer node.  Bind the handoff to the exact
    # producer invocation instead; the validated receipt above already binds
    # that invocation to this source node and its output artifact.
    if (
        producer_scientific_identity_sha256
        != producer_invocation.scientific_identity_sha256
    ):
        raise ContractError(
            "producer scientific identity differs from execution invocation"
        )
    atom_order_sha256 = canonical_sha256({"symbols": handoff.symbols})
    consumer_charge, consumer_multiplicity = handoff.consumer_state
    body = {
        "schema_version": "chemsmart.validated-data-edge-binding.v1",
        "run_id": require_identifier(run_id, "run_id"),
        "workflow_id": plan.workflow_id,
        "plan_sha256": plan.plan_sha256,
        "approval_sha256": approval.approval_sha256,
        "scientific_edge_sha256": edge_sha256,
        "producer_rule_sha256": rule.rule_sha256,
        "source_node_id": rule.source_node_id,
        "target_node_id": rule.target_node_id,
        "artifact_class": rule.artifact_class,
        "producer_output_id": rule.producer_output_id,
        "consumer_input_id": rule.consumer_input_id,
        "selection_rule": rule.selection_rule,
        "producer_execution_receipt_sha256": producer_receipt.receipt_sha256,
        "producer_validator_receipt_sha256s": (
            producer_receipt.validator_receipt_sha256s
        ),
        "source_artifact_id": handoff.result_artifact_id,
        "source_artifact_sha256": handoff.result_artifact_sha256,
        "selected_artifact_id": handoff.selected_artifact_id,
        "selected_artifact_sha256": handoff.selected_artifact_sha256,
        "producer_scientific_identity_sha256": (
            producer_scientific_identity_sha256
        ),
        "consumer_scientific_identity_sha256": require_sha256(
            consumer_scientific_identity_sha256,
            "consumer_scientific_identity_sha256",
        ),
        "atom_order_sha256": atom_order_sha256,
        "positions_sha256": handoff.positions_sha256,
        "charge": consumer_charge,
        "multiplicity": consumer_multiplicity,
        "handoff_receipt_sha256": handoff.receipt_sha256,
        "status": "validated",
    }
    return ValidatedDataEdgeBindingV1(
        **body, receipt_sha256=canonical_sha256(body)
    )


@dataclass(frozen=True)
class WorkflowNodeRunStateV1:
    node_id: str
    state: str
    invocation_sha256: str
    execution_receipt_sha256: str
    validator_receipt_sha256s: tuple[str, ...]
    output_artifact_sha256s: tuple[str, ...]
    failure_rule_ids: tuple[str, ...]

    def __post_init__(self) -> None:
        require_identifier(self.node_id, "node_id")
        if self.state not in {
            "pending",
            "deferred",
            "running",
            "engine_complete",
            "validated",
            "failed",
            "blocked",
            "ambiguous",
        }:
            raise ContractError("unsupported workflow node run state")
        for name, digest in (
            ("invocation_sha256", self.invocation_sha256),
            ("execution_receipt_sha256", self.execution_receipt_sha256),
        ):
            if digest:
                require_sha256(digest, name)
        for name, values in (
            ("validator_receipt_sha256s", self.validator_receipt_sha256s),
            ("output_artifact_sha256s", self.output_artifact_sha256s),
        ):
            if values != tuple(sorted(set(values))):
                raise ContractError(f"{name} must be sorted and unique")
            for digest in values:
                require_sha256(digest, name)
        if self.failure_rule_ids != tuple(sorted(set(self.failure_rule_ids))):
            raise ContractError("failure_rule_ids must be sorted and unique")
        for rule_id in self.failure_rule_ids:
            require_identifier(rule_id, "failure_rule_id")
        if self.state in {"pending", "deferred"} and any(
            (
                self.invocation_sha256,
                self.execution_receipt_sha256,
                self.validator_receipt_sha256s,
                self.output_artifact_sha256s,
                self.failure_rule_ids,
            )
        ):
            raise ContractError(
                f"{self.state} node cannot carry execution evidence"
            )
        if self.state == "running" and not self.invocation_sha256:
            raise ContractError("running node requires an invocation")
        if self.state in {"engine_complete", "validated"}:
            require_sha256(
                self.execution_receipt_sha256,
                "execution_receipt_sha256",
            )
        if self.state == "validated" and not self.validator_receipt_sha256s:
            raise ContractError("validated node requires validator receipts")
        if self.state in {"failed", "blocked", "ambiguous"} and not (
            self.failure_rule_ids
        ):
            raise ContractError("unsuccessful node requires a failure rule")


@dataclass(frozen=True)
class WorkflowRunStateV1:
    """Replayable execution snapshot derived only from host events."""

    schema_version: str
    run_id: str
    workflow_id: str
    plan_sha256: str
    approval_id: str
    approval_sha256: str
    approval_consumed: bool
    state: str
    nodes: tuple[WorkflowNodeRunStateV1, ...]
    started_at: str
    finished_at: str
    run_state_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.workflow-run-state.v1":
            raise ContractError("unsupported workflow run state schema")
        for name, value in (
            ("run_id", self.run_id),
            ("workflow_id", self.workflow_id),
            ("approval_id", self.approval_id),
        ):
            require_identifier(value, name)
        require_sha256(self.plan_sha256, "plan_sha256")
        require_sha256(self.approval_sha256, "approval_sha256")
        if self.state not in {
            "waiting_for_approval",
            "approved",
            "running",
            "validated",
            "failed",
            "blocked",
            "ambiguous",
        }:
            raise ContractError("unsupported workflow run state")
        if self.state != "waiting_for_approval" and not self.approval_consumed:
            raise ContractError("active workflow requires consumed approval")
        if self.approval_consumed and self.state == "waiting_for_approval":
            raise ContractError("consumed workflow cannot wait for approval")
        node_ids = tuple(node.node_id for node in self.nodes)
        if node_ids != tuple(sorted(set(node_ids))):
            raise ContractError("workflow run nodes must be sorted and unique")
        if self.state in {
            "running",
            "validated",
            "failed",
            "ambiguous",
        } and not (self.started_at.strip()):
            raise ContractError("started workflow requires a timestamp")
        if (
            self.state in {"validated", "failed"}
            and not self.finished_at.strip()
        ):
            raise ContractError(
                "terminal workflow requires a finish timestamp"
            )
        body = {
            "schema_version": self.schema_version,
            "run_id": self.run_id,
            "workflow_id": self.workflow_id,
            "plan_sha256": self.plan_sha256,
            "approval_id": self.approval_id,
            "approval_sha256": self.approval_sha256,
            "approval_consumed": self.approval_consumed,
            "state": self.state,
            "nodes": self.nodes,
            "started_at": self.started_at,
            "finished_at": self.finished_at,
        }
        if self.run_state_sha256 != canonical_sha256(body):
            raise ContractError("workflow run state digest mismatch")


def build_workflow_run_state(
    *,
    run_id: str,
    plan: ScientificWorkflowPlanV2,
    approval: FrozenWorkflowApprovalV1,
    approval_consumed: bool = False,
) -> WorkflowRunStateV1:
    if approval.workflow_id != plan.workflow_id:
        raise ContractError("approval belongs to another workflow")
    if approval.plan_sha256 != plan.plan_sha256:
        raise ContractError("approval plan digest differs")
    planned_by_id = {node.node_id: node for node in plan.nodes}
    node_ids = tuple(sorted(planned_by_id))
    approved = set(approval.approved_node_ids)
    if not approved or not approved.issubset(planned_by_id):
        raise ContractError("approval does not match the reviewed workflow")
    deferred = set(planned_by_id) - approved
    for node_id in deferred:
        if planned_by_id[node_id].support_state != "blocked_unsupported":
            raise ContractError(
                "an unapproved workflow node must be blocked unsupported"
            )
    if any(
        edge.source_node_id in deferred and edge.target_node_id in approved
        for edge in plan.edges
    ):
        raise ContractError(
            "an approved node cannot consume a deferred node's output"
        )
    body = {
        "schema_version": "chemsmart.workflow-run-state.v1",
        "run_id": require_identifier(run_id, "run_id"),
        "workflow_id": plan.workflow_id,
        "plan_sha256": plan.plan_sha256,
        "approval_id": approval.approval_id,
        "approval_sha256": approval.approval_sha256,
        "approval_consumed": bool(approval_consumed),
        "state": "approved" if approval_consumed else "waiting_for_approval",
        "nodes": tuple(
            WorkflowNodeRunStateV1(
                node_id=node_id,
                state=("pending" if node_id in approved else "deferred"),
                invocation_sha256="",
                execution_receipt_sha256="",
                validator_receipt_sha256s=(),
                output_artifact_sha256s=(),
                failure_rule_ids=(),
            )
            for node_id in node_ids
        ),
        "started_at": "",
        "finished_at": "",
    }
    return WorkflowRunStateV1(**body, run_state_sha256=canonical_sha256(body))


def derive_ready_node_ids(
    plan: ScientificWorkflowPlanV2,
    run_state: WorkflowRunStateV1,
    data_edge_bindings: Sequence[ValidatedDataEdgeBindingV1] = (),
) -> tuple[str, ...]:
    """Derive runnable nodes from validated predecessors and exact data edges."""

    if run_state.workflow_id != plan.workflow_id:
        raise ContractError("run state belongs to another workflow")
    if run_state.plan_sha256 != plan.plan_sha256:
        raise ContractError("run state plan digest differs")
    if not run_state.approval_consumed or run_state.state in {
        "waiting_for_approval",
        "validated",
        "blocked",
        "ambiguous",
    }:
        return ()
    if run_state.state == "failed" and not any(
        node.state == "pending" for node in run_state.nodes
    ):
        return ()
    state_by_node = {node.node_id: node for node in run_state.nodes}
    binding_by_edge: dict[str, ValidatedDataEdgeBindingV1] = {}
    for binding in data_edge_bindings:
        if (
            binding.run_id != run_state.run_id
            or binding.workflow_id != plan.workflow_id
            or binding.plan_sha256 != plan.plan_sha256
            or binding.approval_sha256 != run_state.approval_sha256
        ):
            raise ContractError(
                "data edge binding belongs to another workflow run"
            )
        if binding.scientific_edge_sha256 in binding_by_edge:
            raise ContractError(
                "multiple bindings exist for one scientific edge"
            )
        binding_by_edge[binding.scientific_edge_sha256] = binding
    ready: list[str] = []
    for node in plan.nodes:
        observed = state_by_node.get(node.node_id)
        if observed is None:
            raise ContractError("run state omits a workflow node")
        if observed.state != "pending":
            continue
        incoming = tuple(
            edge for edge in plan.edges if edge.target_node_id == node.node_id
        )
        can_run = True
        for edge in incoming:
            source = state_by_node[edge.source_node_id]
            if source.state != "validated":
                can_run = False
                break
            if edge.edge_kind == "data":
                binding = binding_by_edge.get(canonical_sha256(edge))
                if binding is None:
                    can_run = False
                    break
                exact_roles = (
                    edge.source_node_id,
                    edge.target_node_id,
                    edge.artifact_class,
                    edge.producer_output_id,
                    edge.consumer_input_id,
                )
                observed_roles = (
                    binding.source_node_id,
                    binding.target_node_id,
                    binding.artifact_class,
                    binding.producer_output_id,
                    binding.consumer_input_id,
                )
                if exact_roles != observed_roles:
                    raise ContractError(
                        "data edge binding roles differ from plan"
                    )
                if (
                    binding.producer_execution_receipt_sha256
                    != source.execution_receipt_sha256
                    or binding.producer_validator_receipt_sha256s
                    != source.validator_receipt_sha256s
                    or binding.source_artifact_sha256
                    not in source.output_artifact_sha256s
                    or binding.selected_artifact_sha256
                    not in source.output_artifact_sha256s
                ):
                    can_run = False
                    break
        if can_run:
            ready.append(node.node_id)
    return tuple(ready)


def transition_workflow_node(
    run_state: WorkflowRunStateV1,
    *,
    node_id: str,
    new_state: str,
    plan: ScientificWorkflowPlanV2 | None = None,
    invocation_sha256: str = "",
    execution_receipt_sha256: str = "",
    validator_receipt_sha256s: Sequence[str] = (),
    result_validation_receipt: ProgramResultValidationReceiptV1 | None = None,
    output_artifact_sha256s: Sequence[str] = (),
    failure_rule_ids: Sequence[str] = (),
    timestamp: str,
) -> WorkflowRunStateV1:
    """Apply one deterministic host transition to an immutable run snapshot."""

    normalized_id = require_identifier(node_id, "node_id")
    by_id = {node.node_id: node for node in run_state.nodes}
    current = by_id.get(normalized_id)
    if current is None:
        raise ContractError("workflow run state has no such node")
    allowed = {
        "pending": {"running", "blocked"},
        "deferred": set(),
        "running": {"engine_complete", "failed", "ambiguous"},
        "engine_complete": {"validated", "failed"},
        "validated": set(),
        "failed": set(),
        "blocked": set(),
        "ambiguous": set(),
    }
    if new_state not in allowed[current.state]:
        raise ContractError("invalid workflow node state transition")
    if plan is not None:
        if (
            plan.workflow_id != run_state.workflow_id
            or plan.plan_sha256 != run_state.plan_sha256
        ):
            raise ContractError("workflow transition plan differs from run state")
    if new_state == "validated":
        if result_validation_receipt is None:
            raise ContractError(
                "validated transition requires a typed result validation receipt"
            )
        if result_validation_receipt.state != "valid":
            raise ContractError(
                "validated transition requires a valid result validation receipt"
            )
        if result_validation_receipt.node_id != normalized_id:
            raise ContractError(
                "result validation receipt belongs to another workflow node"
            )
        if (
            not current.invocation_sha256
            or result_validation_receipt.invocation_sha256
            != current.invocation_sha256
        ):
            raise ContractError(
                "result validation receipt differs from the node invocation"
            )
        if (
            result_validation_receipt.receipt_sha256
            not in validator_receipt_sha256s
        ):
            raise ContractError(
                "typed result validation receipt must be included in validator evidence"
            )
    elif result_validation_receipt is not None:
        raise ContractError(
            "result validation receipt is only accepted for validated transitions"
        )
    replacement = WorkflowNodeRunStateV1(
        node_id=normalized_id,
        state=new_state,
        invocation_sha256=(invocation_sha256 or current.invocation_sha256),
        execution_receipt_sha256=(
            execution_receipt_sha256 or current.execution_receipt_sha256
        ),
        validator_receipt_sha256s=tuple(
            sorted(
                set(validator_receipt_sha256s)
                or set(current.validator_receipt_sha256s)
            )
        ),
        output_artifact_sha256s=tuple(
            sorted(
                set(output_artifact_sha256s)
                or set(current.output_artifact_sha256s)
            )
        ),
        failure_rule_ids=tuple(
            sorted(set(failure_rule_ids) or set(current.failure_rule_ids))
        ),
    )
    nodes = tuple(
        replacement if node.node_id == normalized_id else node
        for node in run_state.nodes
    )
    if plan is not None:
        nodes = _block_failed_descendants(plan, nodes)
    node_states = {node.state for node in nodes}
    executable_states = node_states - {"deferred"}
    if "ambiguous" in executable_states:
        workflow_state = "ambiguous"
    elif executable_states.intersection(
        {"pending", "running", "engine_complete"}
    ):
        # A failed branch is not a failed workflow while an independent branch
        # remains runnable.  The exact DAG below has already converted only
        # causal descendants into ``blocked``; pending siblings therefore keep
        # the workflow active until the frontier is exhausted.
        workflow_state = "running"
    elif executable_states == {"validated"}:
        workflow_state = "validated"
    elif "failed" in executable_states:
        workflow_state = "failed"
    elif "blocked" in executable_states:
        workflow_state = "blocked"
    else:
        workflow_state = "running"
    started_at = run_state.started_at
    if not started_at and (
        new_state == "running" or workflow_state == "running"
    ):
        started_at = str(timestamp).strip()
    finished_at = run_state.finished_at
    if workflow_state == "running":
        # Old streams could mark the run failed as soon as the first branch
        # failed.  If replay now exposes an independent pending sibling, clear
        # that premature timestamp when the sibling advances.
        finished_at = ""
    elif workflow_state in {"validated", "failed", "blocked"}:
        finished_at = str(timestamp).strip()
    body = {
        "schema_version": run_state.schema_version,
        "run_id": run_state.run_id,
        "workflow_id": run_state.workflow_id,
        "plan_sha256": run_state.plan_sha256,
        "approval_id": run_state.approval_id,
        "approval_sha256": run_state.approval_sha256,
        "approval_consumed": run_state.approval_consumed,
        "state": workflow_state,
        "nodes": nodes,
        "started_at": started_at,
        "finished_at": finished_at,
    }
    return WorkflowRunStateV1(**body, run_state_sha256=canonical_sha256(body))


def _block_failed_descendants(
    plan: ScientificWorkflowPlanV2,
    nodes: tuple[WorkflowNodeRunStateV1, ...],
) -> tuple[WorkflowNodeRunStateV1, ...]:
    """Block only pending nodes made impossible by a failed predecessor.

    Every scientific edge is a required dependency.  Once its source is
    failed, blocked, or ambiguous, its pending target cannot run; that fact then
    propagates transitively.  Nodes outside that descendant closure are left
    untouched, which is the key difference from a workflow-global failure.
    """

    by_id = {node.node_id: node for node in nodes}
    changed = True
    while changed:
        changed = False
        for planned in plan.nodes:
            observed = by_id[planned.node_id]
            if observed.state != "pending":
                continue
            failed_sources = tuple(
                sorted(
                    {
                        edge.source_node_id
                        for edge in plan.edges
                        if edge.target_node_id == planned.node_id
                        and by_id[edge.source_node_id].state
                        in {"failed", "blocked", "ambiguous"}
                    }
                )
            )
            if not failed_sources:
                continue
            rules = tuple(
                sorted(
                    "workflow.dependency."
                    + by_id[source].state
                    + "."
                    + source
                    for source in failed_sources
                )
            )
            by_id[planned.node_id] = WorkflowNodeRunStateV1(
                node_id=planned.node_id,
                state="blocked",
                invocation_sha256="",
                execution_receipt_sha256="",
                validator_receipt_sha256s=(),
                output_artifact_sha256s=(),
                failure_rule_ids=rules,
            )
            changed = True
    return tuple(by_id[node.node_id] for node in nodes)


__all__ = [
    "ApprovedNodeBindingV1",
    "ExecutionResourceSpecV1",
    "FrozenMaterializedNodePreviewV1",
    "FrozenProducerEdgeRuleV1",
    "FrozenWorkflowApprovalV1",
    "ORCAHessianHandoffV1",
    "OptimizedGeometryHandoffV1",
    "ProducerEdgeRuleV1",
    "ProgramConformanceProbeV1",
    "ProgramExecutionInvocationV1",
    "ProgramExecutionReceiptV1",
    "ProgramResultValidationReceiptV1",
    "ProjectArtifactPromotionV1",
    "ScientificDecisionRecordV1",
    "WorkflowExecutionApprovalV1",
    "WorkflowExecutionApprovalBundleV1",
    "WorkflowExecutionNodeReviewV1",
    "WorkflowExecutionReviewV1",
    "WorkflowEnvironmentBindingV1",
    "WorkflowReviewResolutionV1",
    "WorkflowNodeRunStateV1",
    "WorkflowRunStateV1",
    "ValidatedDataEdgeBindingV1",
    "bind_project_promotion_validation",
    "build_execution_resource_spec",
    "build_real_execution_argv",
    "build_frozen_workflow_approval",
    "WorkflowApprovalRequestV1",
    "approve_workflow_request",
    "approve_workflow_execution_review",
    "build_producer_edge_rule",
    "build_workflow_approval_request",
    "build_workflow_execution_review",
    "build_workflow_execution_node_review",
    "build_workflow_review_resolution",
    "workflow_approval_request_json",
    "build_program_conformance_probe",
    "build_program_execution_invocation",
    "build_program_execution_receipt",
    "build_program_result_validation_receipt",
    "build_scientific_decision_record",
    "build_validated_data_edge_binding",
    "build_workflow_execution_approval",
    "build_workflow_run_state",
    "derive_ready_node_ids",
    "environment_review_summary",
    "handoff_optimized_native_geometry",
    "handoff_optimized_pyscf_geometry",
    "handoff_optimized_xtb_geometry",
    "handoff_final_orca_ts_hessian",
    "is_validated_orca_ts_hessian_edge",
    "is_validated_optimized_geometry_edge",
    "execution_path_placeholder",
    "execution_server_profile_sha256",
    "normalized_project_settings_review",
    "path_free_review_record",
    "project_real_execution_argv",
    "real_execution_argv_sha256",
    "promote_project_candidate",
    "transition_workflow_node",
    "workflow_execution_approval_json",
    "workflow_execution_approval_bundle_json",
    "workflow_execution_review_json",
]
