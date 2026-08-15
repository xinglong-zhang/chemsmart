"""Host-owned dispatcher for the command-compiled model tool surface."""

from __future__ import annotations

from dataclasses import dataclass, replace
from datetime import datetime, timezone
import json
import logging
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import time
from typing import Any, Mapping, Sequence

from chemsmart.agent._contracts import (
    ContractError,
    TrustedArtifactRefV1,
    canonical_data,
    canonical_sha256,
    file_sha256,
    require_identifier,
    require_sha256,
)
from chemsmart.agent.analysis_completion import (
    AnalysisIncompleteError,
    AnalysisCompletionReceiptV1,
    AnalysisCompletionPolicyV1,
    build_analysis_completion_receipt,
)
from chemsmart.agent.analysis_claims import (
    AnalysisReportedQuantityV1,
    analysis_claim_record_from_record,
    build_analysis_claim_record,
)
from chemsmart.agent.skills import resolve_skill
from chemsmart.agent.capabilities import (
    CapabilityQueryReceiptV1,
    CapabilityQueryV1,
    EnvironmentCapabilityReceiptV1,
    EnvironmentTargetV1,
    ProgramCapabilityRegistryV1,
    ProgramComponentConformanceReceiptV1,
    ProgramSupportOverlayV1,
    ResolvedEngineBindingV1,
    ResolvedProgramBindingV1,
    SupportLevel,
    TrustedComputeEnvironmentReceiptV1,
    build_approved_execution_overlay,
    build_command_compiled_preview_overlay,
    consume_pyscf_compute_environment_receipt,
    environment_identity_sha256,
    load_program_capabilities,
    query_capability,
    query_environment,
    resolve_engine_binding,
    resolve_program_binding,
)
from chemsmart.agent.cli_schema import (
    LiveClickSchemaV1,
    build_live_click_schema,
)
from chemsmart.agent.commands import (
    CanonicalCommandInvocationV1,
    CommandCounterexampleV1,
    CommandInspectionReceiptV1,
    CommandProposalV1,
    ScientificIdentityBindingV1,
    build_scientific_identity_binding,
    compile_command,
    inspect_command,
)
from chemsmart.agent.inspection import (
    GeneratedArtifactInspectionReceiptV1,
    inspect_generated_artifact,
)
from chemsmart.agent.identity import ApprovedMolecularIdentityV1
from chemsmart.agent.execution import (
    ApprovedNodeBindingV1,
    ExecutionResourceSpecV1,
    FrozenWorkflowApprovalV1,
    ORCAHessianHandoffV1,
    OptimizedGeometryHandoffV1,
    ProgramExecutionReceiptV1,
    ProgramResultValidationReceiptV1,
    ProjectArtifactPromotionV1,
    ScientificDecisionRecordV1,
    WorkflowExecutionApprovalV1,
    WorkflowExecutionReviewV1,
    WorkflowExecutionNodeReviewV1,
    WorkflowEnvironmentBindingV1,
    WorkflowNodeRunStateV1,
    WorkflowRunStateV1,
    bind_project_promotion_validation,
    build_frozen_workflow_approval,
    build_producer_edge_rule,
    build_program_execution_invocation,
    build_program_execution_receipt,
    build_program_result_validation_receipt,
    build_scientific_decision_record,
    build_validated_data_edge_binding,
    build_workflow_execution_approval,
    build_workflow_approval_request,
    build_workflow_execution_review,
    build_workflow_execution_node_review,
    build_real_execution_argv,
    derive_ready_node_ids,
    environment_review_summary,
    execution_server_profile_sha256,
    handoff_optimized_native_geometry,
    handoff_optimized_pyscf_geometry,
    handoff_optimized_xtb_geometry,
    handoff_final_orca_ts_hessian,
    invocation_identity_sha256,
    is_validated_orca_ts_hessian_edge,
    is_validated_optimized_geometry_edge,
    promote_project_candidate,
    execution_path_placeholder,
    project_real_execution_argv,
)
from chemsmart.agent.execution_envelope import BoundedExecutionEnvelopeV1
from chemsmart.agent.knowledge import (
    FunctionalEquivalenceReceiptV1,
    ProgramSubstitutionReceiptV1,
    ProgramSubstitutionApprovalV1,
    ScientificClaimEvidenceV1,
    assess_typed_program_substitution,
    build_program_substitution_request,
)
from chemsmart.agent.preflight import (
    ProgramNodePreflightReceiptV1,
    ProgramValidatorReceiptV1,
    build_program_node_preflight_request,
    evaluate_program_node_preflight,
    validator_receipt_from_safe_preview,
)
from chemsmart.agent.postprocessing import (
    derive_trusted_thermochemistry,
    evaluate_typed_quantity_expression,
    extract_trusted_result_quantities,
)
from chemsmart.agent.preview import SafePreviewReceiptV1, execute_safe_preview
from chemsmart.agent.program_verifiers import build_preview_expectation
from chemsmart.agent.projects import (
    ProjectDocumentV1,
    ProjectRenderReceiptV1,
    ProjectValidationReceiptV1,
    PySCFFunctionalResolutionReceiptV1,
    project_document,
    project_section_application_observation,
    project_scientific_materializations,
    read_project_yaml,
    render_project_yaml,
    validate_project_yaml,
)
from chemsmart.agent.runtime.event_store import RuntimeEventStore
from chemsmart.agent.runtime.events import EventKind
from chemsmart.agent.scientific_toolchain import (
    AnalysisInputIntentV1,
    AnalysisNodeIntentV1,
    AnalysisOutputIntentV1,
    AnalysisSelectorIntentV1,
    AnalysisValidationRuleIntentV1,
    RegisteredResultInputIntentV1,
    ScientificToolchainPlanV1,
    build_scientific_toolchain_plan,
    project_scientific_toolchain_frontier,
)
from chemsmart.agent.scientific_validation import (
    ScientificValidationReceiptV1,
    evaluate_planned_scientific_validation,
    scientific_validation_receipt_from_record,
)
from chemsmart.agent.tool_specs import (
    REGISTRY_PRODUCERS,
    AgentToolSurfaceV1,
    build_approved_execution_tool_surface,
    build_command_compiled_tool_surface,
)
from chemsmart.analysis.quantity_expressions import (
    QuantityExpressionNodeV1,
    QuantityExpressionRequestV1,
    convert_normalized_value,
    normalize_numeric_value,
    quantity_expression_receipt_from_record,
)
from chemsmart.analysis.result_quantities import (
    QuantityExtractionReceiptV1,
    QuantitySelectorV1,
    QuantityValueV1,
    ThermochemistryReceiptV1,
    canonical_thermochemistry_quantity,
    make_quantity_value,
    quantity_extraction_receipt_from_record,
    thermochemistry_receipt_from_record,
)
from chemsmart.agent.workflow_context import (
    project_workflow_context,
)
from chemsmart.agent.dependency_context import (
    ContextSelectionReceiptV1,
    TaskDependencyContextV2,
)
from chemsmart.agent.workflows import (
    AGGREGATE_NODE_PROGRAM,
    ArtifactInputIntentV1,
    ArtifactOutputIntentV1,
    CommandNodeIntentV1,
    CommandWorkflowDraftV1,
    MaterializedNodeV1,
    MaterializedWorkflowV1,
    ScientificWorkflowEdgeV2,
    ScientificWorkflowNodeV2,
    ScientificWorkflowPlanV2,
    StationaryPointValidationPolicyV1,
    build_command_workflow_draft,
    build_materialized_workflow,
    build_scientific_workflow_plan,
)
from chemsmart.utils.process_observation import (
    ProcessObservationV1,
    ProcessSignalGuard,
    launch_failure_observation,
    observe_process,
)



logger = logging.getLogger(__name__)

@dataclass(frozen=True)
class _CommandContext:
    proposal: CommandProposalV1
    capability: CapabilityQueryReceiptV1
    program_binding: ResolvedProgramBindingV1
    engine_binding: ResolvedEngineBindingV1
    project_artifact: TrustedArtifactRefV1 | None
    project_validation: ProjectValidationReceiptV1 | None
    input_artifact: TrustedArtifactRefV1
    scientific_identity: ScientificIdentityBindingV1
    job_artifact_options: tuple[tuple[str, TrustedArtifactRefV1], ...] = ()


def _formula_from_symbols(symbols: Sequence[str]) -> str:
    """Return a compact Hill formula for a human execution review."""

    counts: dict[str, int] = {}
    for symbol in symbols:
        canonical = str(symbol)
        counts[canonical] = counts.get(canonical, 0) + 1
    if "C" in counts:
        order = ["C"]
        if "H" in counts:
            order.append("H")
        order.extend(sorted(item for item in counts if item not in {"C", "H"}))
    else:
        order = sorted(counts)
    return "".join(
        symbol + (str(counts[symbol]) if counts[symbol] != 1 else "")
        for symbol in order
    )


def _review_molecule_identity(
    artifact: TrustedArtifactRefV1,
) -> dict[str, Any]:
    """Read path-free molecular facts from ChemSmart's trusted geometry."""

    from chemsmart.io.molecules.structure import Molecule

    molecule = Molecule.from_filepath(artifact.path)
    if molecule is None:
        raise ContractError("execution review cannot read the molecular input")
    symbols = tuple(str(symbol) for symbol in molecule.symbols)
    if not symbols:
        raise ContractError("execution review molecular input has no atoms")
    return {
        "atom_count": len(symbols),
        "atom_order": symbols,
        "formula": _formula_from_symbols(symbols),
    }


@dataclass(frozen=True)
class _ResolvedProgramWorkflow:
    """Exact host-owned view used by workflow program tools.

    A calculation-only command plan and a calculation-plus-analysis toolchain
    use different public schemas.  Program tools need their shared command
    draft and task-bound V2 plan, so retain that binding without requiring an
    unrelated analysis plan to exist.
    """

    draft: CommandWorkflowDraftV1
    scientific_plan: ScientificWorkflowPlanV2 | None
    command_result: Mapping[str, Any]
    scientific_toolchain_plan: ScientificToolchainPlanV1 | None = None


@dataclass(frozen=True)
class _PySCFEngineObservation:
    """Digest-bound child-engine evidence independent of wrapper status."""

    child_exit_status: int | None
    engine_complete: bool
    run_receipt_sha256: str
    run_receipt: Mapping[str, Any] | None
    result_artifact: TrustedArtifactRefV1 | None
    findings: tuple[str, ...]


@dataclass(frozen=True)
class _ExecutionValidationEvaluation:
    """Complete deterministic observation before receipt materialization."""

    validator_id: str
    validator_schema_version: str
    validator_version: str
    observations: Mapping[str, Any]
    findings: tuple[str, ...]
    run_environment_receipt_sha256: str = ""
    environment_validation_sha256: str = ""

    @property
    def validated(self) -> bool:
        return not self.findings


def _scientific_decision_binding_requirement(
    materializations: tuple[PySCFFunctionalResolutionReceiptV1, ...],
) -> dict[str, Any]:
    """Expose the exact durable action needed before rendering XC semantics.

    Project validation is the first point at which a PySCF functional literal
    has a host-owned applied-XC interpretation.  A model may have recorded a
    useful task-level decision earlier, but that record cannot ground a later
    implementation-specific alias or correlation-convention claim.  This
    small typed action keeps the evidence requirement visible in both full and
    causal feedback without making the host author the scientific rationale.
    """

    evidence_refs = tuple(
        sorted(
            item.evidence_ref for item in materializations if item.evidence_ref
        )
    )
    body = {
        "schema_version": "chemsmart.scientific-decision-binding.v1",
        "status": "required_if_rendering_implementation_semantics",
        "rule_id": "scientific.functional_resolution.decision_binding",
        "next_tool": "record_scientific_decision",
        "evidence_refs": evidence_refs,
        "message": (
            "After project validation, record a final scientific decision "
            "citing these exact evidence_refs before rendering an applied XC "
            "alias or correlation convention. An earlier task-level decision "
            "does not satisfy this evidence binding."
            if evidence_refs
            else "No functional materialization requires a decision binding."
        ),
    }
    return {**body, "receipt_sha256": canonical_sha256(body)}


def _postprocessing_evidence_reference(
    value: str,
) -> tuple[str, str] | None:
    """Parse canonical or hyphenated typed postprocessing references.

    Provider prose frequently varies ``-`` and ``_`` in public labels.  The
    reference grammar therefore normalizes only the type prefix, while the
    exact receipt digest remains unchanged and must still resolve in host
    state.  This is a presentation normalization, not evidence relaxation.
    """

    reference = str(value).strip()
    if ":" not in reference:
        return None
    raw_kind, digest = reference.split(":", 1)
    normalized_kind = raw_kind.strip().lower().replace("-", "_")
    aliases = {
        "quantity_extraction_receipt": "quantity_extraction",
        "thermochemistry_receipt": "thermochemistry",
        "quantity_expression_receipt": "quantity_expression",
        "scientific_validation_receipt": "scientific_validation",
        "analysis_claim_receipt": "analysis_claim",
        "receipt": "generic",
    }
    kind = aliases.get(normalized_kind)
    if kind is None:
        return None
    require_sha256(digest, "postprocessing evidence receipt")
    return kind, digest


def _current_artifact_path(
    artifact: TrustedArtifactRefV1,
    *,
    field_name: str,
) -> Path:
    """Rehash an artifact immediately before a validator opens it."""

    path = Path(artifact.path)
    if not path.is_file() or path.is_symlink():
        raise ContractError(f"{field_name} is not a current regular file")
    before = path.stat()
    if before.st_size != artifact.size_bytes:
        raise ContractError(f"{field_name} size differs from its binding")
    observed_sha256 = file_sha256(path)
    after = path.stat()
    if (
        after.st_size != before.st_size
        or after.st_mtime_ns != before.st_mtime_ns
        or observed_sha256 != artifact.sha256
    ):
        raise ContractError(f"{field_name} digest differs from its binding")
    return path.resolve()


def _environment_semantic_facts(
    receipt: (
        EnvironmentCapabilityReceiptV1 | TrustedComputeEnvironmentReceiptV1
    ),
) -> Mapping[str, Any]:
    """Return only stable facts shared by capability and per-run probes."""

    return {
        "program": receipt.program,
        "engine": receipt.engine,
        "compute_interpreter_sha256": receipt.compute_interpreter_sha256,
        "dependency_versions": receipt.dependency_versions,
        "solver_evidence": receipt.solver_evidence,
        "gpu_evidence": receipt.gpu_evidence,
    }


def _project_environment_semantic_facts(
    receipt: TrustedComputeEnvironmentReceiptV1,
    *,
    approved_facts: Mapping[str, Any],
) -> Mapping[str, Any]:
    """Project a detailed run probe onto the approved capability facts.

    A per-run PySCF probe is intentionally richer than the capability receipt
    used during planning.  Additional packages, solver probes, or GPU facts do
    not change the approved environment; only a missing or changed approved
    fact does.  Comparing the two complete records therefore rejects safe
    supersets produced by the same interpreter.
    """

    observed = dict(_environment_semantic_facts(receipt))
    for field in ("dependency_versions", "solver_evidence", "gpu_evidence"):
        available = dict(observed[field])
        observed[field] = tuple(
            (name, available.get(name))
            for name, _value in approved_facts[field]
        )
    return observed


def _pyscf_environment_evidence(
    *,
    output_artifacts: tuple[TrustedArtifactRefV1, ...],
    run_receipt: Mapping[str, Any] | None,
    capability_environment: EnvironmentCapabilityReceiptV1 | None,
) -> tuple[Mapping[str, Any], tuple[str, ...]]:
    """Compare different environment receipt types by stable semantics."""

    findings: list[str] = []
    candidates: list[
        tuple[TrustedArtifactRefV1, Mapping[str, Any] | None]
    ] = []
    for artifact in output_artifacts:
        if artifact.kind != "json":
            continue
        try:
            path = _current_artifact_path(
                artifact, field_name="PySCF environment receipt"
            )
            raw = json.loads(path.read_text(encoding="utf-8"))
        except (
            ContractError,
            OSError,
            TypeError,
            ValueError,
            json.JSONDecodeError,
        ):
            continue
        if isinstance(raw, dict) and raw.get("schema_version") == (
            "chemsmart.pyscf-environment.v1"
        ):
            candidates.append((artifact, _digest_valid_json_receipt(path)))

    observation: dict[str, Any] = {
        "capability_environment_receipt_sha256": (
            capability_environment.receipt_sha256
            if capability_environment is not None
            else ""
        ),
        "run_environment_receipt_sha256": "",
        "approved_semantic_fingerprint_sha256": "",
        "observed_semantic_fingerprint_sha256": "",
    }
    if capability_environment is None:
        findings.append("pyscf.environment.capability_receipt_unavailable")
    elif capability_environment.status.value != "available":
        findings.append("pyscf.environment.capability_not_available")
    if len(candidates) != 1:
        findings.append("pyscf.environment.run_receipt_count")
        observation["state"] = "invalid"
        observation["validation_sha256"] = canonical_sha256(observation)
        return observation, tuple(sorted(set(findings)))
    _artifact, raw_receipt = candidates[0]
    if raw_receipt is None:
        findings.append("pyscf.environment.run_receipt_digest_invalid")
        observation["state"] = "invalid"
        observation["validation_sha256"] = canonical_sha256(observation)
        return observation, tuple(sorted(set(findings)))

    run_environment_sha256 = str(raw_receipt["receipt_sha256"])
    observation["run_environment_receipt_sha256"] = run_environment_sha256
    if (
        run_receipt is None
        or run_receipt.get("environment_receipt_sha256")
        != run_environment_sha256
    ):
        findings.append("pyscf.environment.run_receipt_link_mismatch")
    if capability_environment is not None:
        try:
            adapted = consume_pyscf_compute_environment_receipt(
                raw_receipt, engine=capability_environment.engine
            )
        except ContractError:
            findings.append("pyscf.environment.semantic_adaptation_failed")
        else:
            approved_facts = _environment_semantic_facts(
                capability_environment
            )
            observed_facts = _project_environment_semantic_facts(
                adapted, approved_facts=approved_facts
            )
            approved_fingerprint = canonical_sha256(approved_facts)
            observed_fingerprint = canonical_sha256(observed_facts)
            observation.update(
                {
                    "approved_semantic_fingerprint_sha256": (
                        approved_fingerprint
                    ),
                    "observed_semantic_fingerprint_sha256": (
                        observed_fingerprint
                    ),
                    "program": adapted.program,
                    "engine": adapted.engine,
                }
            )
            if approved_fingerprint != observed_fingerprint:
                findings.append("pyscf.environment.semantic_mismatch")
    observation["state"] = "valid" if not findings else "invalid"
    observation["validation_sha256"] = canonical_sha256(observation)
    return observation, tuple(sorted(set(findings)))


def _digest_valid_json_receipt(path: str | Path) -> Mapping[str, Any] | None:
    """Load a JSON receipt only when its embedded digest is exact."""

    try:
        payload = json.loads(Path(path).read_text(encoding="utf-8"))
        if not isinstance(payload, dict):
            return None
        observed = payload.get("receipt_sha256")
        body = dict(payload)
        body.pop("receipt_sha256", None)
        if observed and observed == canonical_sha256(body):
            return payload
    except (OSError, TypeError, ValueError, json.JSONDecodeError):
        pass
    return None


def _inspect_pyscf_engine_observation(
    output_artifacts: tuple[TrustedArtifactRefV1, ...],
    *,
    launch_ambiguous: bool,
) -> _PySCFEngineObservation:
    """Resolve PySCF child completion from its exact run receipt.

    The outer ChemSmart process may return non-zero after writing a complete
    child receipt because post-run scientific validation intentionally raises.
    A timeout remains ambiguous even if partial artifacts happen to exist.
    """

    if launch_ambiguous:
        return _PySCFEngineObservation(
            child_exit_status=None,
            engine_complete=False,
            run_receipt_sha256="",
            run_receipt=None,
            result_artifact=None,
            findings=("execution.process.ambiguous",),
        )

    schema_candidates: list[
        tuple[TrustedArtifactRefV1, Mapping[str, Any] | None]
    ] = []
    for artifact in output_artifacts:
        if artifact.kind != "json":
            continue
        try:
            path = _current_artifact_path(
                artifact, field_name="PySCF run receipt"
            )
            raw = json.loads(path.read_text(encoding="utf-8"))
        except (
            ContractError,
            OSError,
            TypeError,
            ValueError,
            json.JSONDecodeError,
        ):
            continue
        if isinstance(raw, dict) and raw.get("schema_version") == (
            "chemsmart.pyscf-run.v1"
        ):
            schema_candidates.append(
                (artifact, _digest_valid_json_receipt(path))
            )

    findings: list[str] = []
    if len(schema_candidates) != 1:
        findings.append("pyscf.run_receipt.count")
        return _PySCFEngineObservation(
            child_exit_status=None,
            engine_complete=False,
            run_receipt_sha256="",
            run_receipt=None,
            result_artifact=None,
            findings=tuple(findings),
        )
    _artifact, run_receipt = schema_candidates[0]
    if run_receipt is None:
        findings.append("pyscf.run_receipt.digest_invalid")
        return _PySCFEngineObservation(
            child_exit_status=None,
            engine_complete=False,
            run_receipt_sha256="",
            run_receipt=None,
            result_artifact=None,
            findings=tuple(findings),
        )

    child_exit_status = run_receipt.get("child_returncode")
    if isinstance(child_exit_status, bool) or not isinstance(
        child_exit_status, int
    ):
        findings.append("pyscf.run_receipt.child_status_invalid")
        child_exit_status = None
    if run_receipt.get("fake") is not False:
        findings.append("pyscf.run_receipt.not_real_execution")

    results = tuple(
        artifact
        for artifact in output_artifacts
        if artifact.kind == "pyscf_hdf5"
    )
    result_artifact = results[0] if len(results) == 1 else None
    if result_artifact is None:
        findings.append("pyscf.result.hdf5_count")
    else:
        try:
            _current_artifact_path(
                result_artifact, field_name="PySCF result artifact"
            )
        except ContractError:
            findings.append("pyscf.result.artifact_binding_mismatch")
            result_artifact = None
        if (
            result_artifact is not None
            and run_receipt.get("result_sha256") != result_artifact.sha256
        ):
            findings.append("pyscf.run_receipt.result_digest_mismatch")

    receipt_declares_complete = run_receipt.get("engine_complete") is True
    if not receipt_declares_complete:
        findings.append("pyscf.run_receipt.engine_incomplete")
    engine_complete = bool(
        receipt_declares_complete
        and child_exit_status == 0
        and result_artifact is not None
        and run_receipt.get("result_sha256") == result_artifact.sha256
        and run_receipt.get("fake") is False
    )
    return _PySCFEngineObservation(
        child_exit_status=child_exit_status,
        engine_complete=engine_complete,
        run_receipt_sha256=str(run_receipt["receipt_sha256"]),
        run_receipt=run_receipt,
        result_artifact=result_artifact,
        findings=tuple(findings),
    )


def _pyscf_input_geometry(
    artifact: TrustedArtifactRefV1 | None,
) -> tuple[tuple[str, ...], tuple[tuple[float, float, float], ...]]:
    """Read exact atom order and positions from the approved input artifact."""

    if artifact is None:
        return (), ()
    try:
        path = _current_artifact_path(
            artifact, field_name="PySCF approved input artifact"
        )
        if artifact.kind == "pyscf_hdf5":
            from chemsmart.io.pyscf.output import read_pyscf_h5

            spec, _provenance, _status, results = read_pyscf_h5(path)
            symbols = tuple(str(value) for value in spec.get("symbols") or ())
            raw_positions = results.get("positions")
            positions = tuple(
                tuple(float(component) for component in row)
                for row in raw_positions
            )
            if len(symbols) != len(positions) or any(
                len(row) != 3 for row in positions
            ):
                return (), ()
            return symbols, positions
        lines = path.read_text(encoding="utf-8").splitlines()
        atom_count = int(lines[0].strip())
        atom_lines = lines[2 : atom_count + 2]
        if atom_count < 1 or len(atom_lines) != atom_count:
            return (), ()
        columns = tuple(line.split() for line in atom_lines)
        if any(len(row) < 4 for row in columns):
            return (), ()
        symbols = tuple(row[0] for row in columns)
        positions = tuple(
            (float(row[1]), float(row[2]), float(row[3])) for row in columns
        )
        return symbols, positions
    except (
        ContractError,
        OSError,
        TypeError,
        UnicodeDecodeError,
        ValueError,
        IndexError,
        KeyError,
    ):
        return (), ()


def _pyscf_input_symbols(
    artifact: TrustedArtifactRefV1 | None,
) -> tuple[str, ...]:
    """Compatibility projection of the exact approved input geometry."""

    symbols, _positions = _pyscf_input_geometry(artifact)
    return symbols


def _pyscf_input_geometry_sha256(
    artifact: TrustedArtifactRefV1 | None,
    *,
    charge: int,
    multiplicity: int,
) -> str:
    """Return the PySCF writer's canonical identity for an approved geometry."""

    symbols, positions = _pyscf_input_geometry(artifact)
    if not symbols or not positions:
        return ""
    return canonical_sha256(
        {
            "symbols": symbols,
            "positions": positions,
            "unit": "Angstrom",
            "charge": charge,
            "multiplicity": multiplicity,
        }
    )


def _pyscf_result_receipt_expectation(
    run_receipt: Mapping[str, Any] | None,
) -> Mapping[str, Any] | None:
    """Translate the digest-valid outer run receipt to HDF5 bindings."""

    if run_receipt is None:
        return None
    expected = {
        field: run_receipt.get(field)
        for field in (
            "run_id",
            "run_nonce",
            "script_sha256",
            "input_receipt_sha256",
            "environment_receipt_sha256",
            "input_geometry_sha256",
            "input_artifact_kind",
            "input_artifact_sha256",
            "requested_settings_sha256",
            "applied_settings_sha256",
        )
    }
    expected.update(
        {
            "project_yaml_digest": run_receipt.get("project_yaml_sha256"),
            "require_applied_settings_sha256": True,
            "require_engine_complete": True,
        }
    )
    return expected


def _validate_stationary_point_policy_binding(
    frozen_approval: FrozenWorkflowApprovalV1 | None,
    policy: StationaryPointValidationPolicyV1 | None,
    *,
    plan: ScientificWorkflowPlanV2 | None = None,
    hessian_node_id: str = "",
    require_for_hessian: bool = False,
) -> None:
    """Require a task-, plan-, and node-bound Hessian policy."""

    if frozen_approval is not None and not isinstance(
        frozen_approval, FrozenWorkflowApprovalV1
    ):
        raise ContractError(
            "stationary-point binding requires a real approval"
        )
    if policy is not None and not isinstance(
        policy, StationaryPointValidationPolicyV1
    ):
        raise ContractError("stationary-point binding requires a real policy")
    if frozen_approval is None:
        if policy is not None:
            raise ContractError(
                "stationary-point policy requires a frozen approval"
            )
        if require_for_hessian:
            raise ContractError(
                "Hessian execution requires a frozen workflow approval"
            )
        return
    approved_sha256 = frozen_approval.stationary_point_policy_sha256
    if not approved_sha256:
        if policy is not None:
            raise ContractError(
                "unapproved stationary-point policy was supplied"
            )
        if require_for_hessian:
            raise ContractError(
                "Hessian execution requires an approved stationary-point policy"
            )
        return
    if policy is None or policy.policy_sha256 != approved_sha256:
        raise ContractError("frozen stationary-point policy is unavailable")
    if plan is None or not isinstance(plan, ScientificWorkflowPlanV2):
        raise ContractError(
            "stationary-point policy requires its exact scientific plan"
        )
    if plan.plan_sha256 != frozen_approval.plan_sha256:
        raise ContractError(
            "stationary-point policy plan differs from approval"
        )
    if policy.task_spec_sha256 != frozen_approval.task_spec_sha256 or (
        policy.task_spec_sha256 != plan.task_spec_sha256
    ):
        raise ContractError(
            "stationary-point policy task differs from approval"
        )
    if policy.hessian_node_id not in frozen_approval.approved_node_ids:
        raise ContractError("stationary-point Hessian node is not approved")
    matching_nodes = tuple(
        node for node in plan.nodes if node.node_id == policy.hessian_node_id
    )
    if len(matching_nodes) != 1 or matching_nodes[0].stage != "hess":
        raise ContractError("stationary-point policy must bind a Hessian node")
    if hessian_node_id and policy.hessian_node_id != hessian_node_id:
        raise ContractError(
            "stationary-point policy targets another Hessian node"
        )
    if policy.require_finite_modes is not True:
        raise ContractError(
            "stationary-point policy must require finite modes"
        )
    if policy.require_symmetric_hessian is not True:
        raise ContractError(
            "stationary-point policy must require a symmetric Hessian"
        )


def _runner_defers_hessian_classification(
    *,
    run_receipt: Mapping[str, Any],
    jobtype: str,
    hessian_node_id: str,
    engine_complete: bool,
    stationary_point_policy: StationaryPointValidationPolicyV1 | None,
    approved_stationary_point_policy_sha256: str,
) -> bool:
    """Admit a complete HESS while leaving stationary-point meaning explicit.

    The CLI runner owns engine/artifact invariants but does not own an
    imaginary-mode expectation.  A Hessian is still a valid computed artifact
    when its minimum/transition-state classification has not been prescribed
    in advance.  In that case the scientific DAG must classify the modes
    downstream.  If a policy *is* supplied, it must remain exactly bound.
    """

    if jobtype != "hess" or not hessian_node_id:
        return False
    if stationary_point_policy is None:
        if approved_stationary_point_policy_sha256:
            return False
    else:
        if not isinstance(
            stationary_point_policy, StationaryPointValidationPolicyV1
        ):
            return False
        if stationary_point_policy.hessian_node_id != hessian_node_id:
            return False
        if (
            not approved_stationary_point_policy_sha256
            or stationary_point_policy.policy_sha256
            != approved_stationary_point_policy_sha256
        ):
            return False
    result_validation = run_receipt.get("result_validation")
    if not isinstance(result_validation, Mapping):
        return False
    frequency_validation = result_validation.get("frequency_validation")
    if not isinstance(frequency_validation, Mapping):
        return False
    return bool(
        engine_complete
        and run_receipt.get("engine_complete") is True
        and run_receipt.get("child_returncode") == 0
        and run_receipt.get("state") == "engine_complete"
        and run_receipt.get("scientifically_validated") is False
        and run_receipt.get("scientific_validation_state") == "unclassified"
        and not (run_receipt.get("findings") or ())
        and result_validation.get("state") == "unclassified"
        and not (result_validation.get("findings") or ())
        and frequency_validation.get("stationary_point_classification")
        == "unclassified"
        and (
            stationary_point_policy is None
            or (
                stationary_point_policy.require_finite_modes is True
                and stationary_point_policy.require_symmetric_hessian is True
            )
        )
    )


class CommandCompiledToolHostV1:
    """Resolve every model ID against immutable host-held objects."""

    def __init__(
        self,
        *,
        event_store: RuntimeEventStore,
        artifacts: Mapping[str, TrustedArtifactRefV1] = {},
        scientific_identities: Mapping[str, ScientificIdentityBindingV1] = {},
        approved_molecular_identities: Mapping[
            str, ApprovedMolecularIdentityV1
        ] = {},
        environment_targets: tuple[EnvironmentTargetV1, ...] = (),
        compute_environment_receipts: tuple[
            TrustedComputeEnvironmentReceiptV1, ...
        ] = (),
        component_conformance_receipts: tuple[
            ProgramComponentConformanceReceiptV1, ...
        ] = (),
        support_overlay: ProgramSupportOverlayV1 | None = None,
        tool_surface: AgentToolSurfaceV1 | None = None,
        settings_objects: Mapping[str, Any] = {},
        run_receipts: Mapping[str, Mapping[str, Any]] = {},
        scientific_claim_evidence: Mapping[
            str, ScientificClaimEvidenceV1
        ] = {},
        functional_equivalence_receipts: Mapping[
            str, FunctionalEquivalenceReceiptV1
        ] = {},
        substitution_approvals: Mapping[
            str, ProgramSubstitutionApprovalV1
        ] = {},
        capability_receipts: Mapping[str, CapabilityQueryReceiptV1] = {},
        environment_receipts: Mapping[
            str, EnvironmentCapabilityReceiptV1
        ] = {},
        program_binding_receipts: Mapping[str, ResolvedProgramBindingV1] = {},
        engine_binding_receipts: Mapping[str, ResolvedEngineBindingV1] = {},
        project_validation_receipts: Mapping[
            str, ProjectValidationReceiptV1
        ] = {},
        result_functional_evidence: Mapping[str, Mapping[str, Any]] = {},
        analysis_completion_policy: AnalysisCompletionPolicyV1 | None = None,
        registry: ProgramCapabilityRegistryV1 | None = None,
        live_schema: LiveClickSchemaV1 | None = None,
        task_spec_sha256s: tuple[str, ...] = (),
        approved_workspace: str | Path | None = None,
        execution_resources: ExecutionResourceSpecV1 | None = None,
        workflow_execution_approval: WorkflowExecutionApprovalV1 | None = None,
        frozen_workflow_approval: FrozenWorkflowApprovalV1 | None = None,
        bounded_execution_envelope: BoundedExecutionEnvelopeV1 | None = None,
        approved_environment_identities: tuple[str, ...] = (),
        materialized_workflow: MaterializedWorkflowV1 | None = None,
        stationary_point_policy: (
            StationaryPointValidationPolicyV1 | None
        ) = None,
        scientific_workflow_plan: ScientificWorkflowPlanV2 | None = None,
        preview_server: str = "",
        execution_server: str = "",
        execution_server_file_sha256: str = "",
        execution_environment: Mapping[str, str] = {},
        execution_environment_remove: tuple[str, ...] = (),
        dependency_context: TaskDependencyContextV2 | None = None,
        dependency_context_selection_receipt: (
            ContextSelectionReceiptV1 | None
        ) = None,
    ) -> None:
        self.event_store = event_store
        self.registry = registry or load_program_capabilities()
        self.live_schema = live_schema or build_live_click_schema()
        preview_overlay = build_command_compiled_preview_overlay(
            self.registry,
            conformance_receipts=component_conformance_receipts,
            live_schema=self.live_schema,
        )
        command_surface = build_command_compiled_tool_surface(self.registry)
        execution_surface = build_approved_execution_tool_surface(
            self.registry
        )
        if (
            tool_surface is not None
            and tool_surface.tool_schema_sha256
            not in {
                command_surface.tool_schema_sha256,
                execution_surface.tool_schema_sha256,
            }
        ):
            raise ContractError(
                "injected tool surface is not a canonical profile"
            )
        self.surface = tool_surface or command_surface
        self.artifacts = dict(artifacts)
        self.task_spec_sha256s = frozenset(task_spec_sha256s)
        self.approved_environment_identities = tuple(
            approved_environment_identities
        )
        self.approved_workspace = (
            Path(approved_workspace).resolve()
            if approved_workspace is not None
            else None
        )
        self.execution_resources = execution_resources
        self.workflow_execution_approval = workflow_execution_approval
        self.frozen_workflow_approval = frozen_workflow_approval
        self.bounded_execution_envelope = bounded_execution_envelope
        self._bounded_execution_started_at = time.monotonic()
        self.stationary_point_policy = stationary_point_policy
        _validate_stationary_point_policy_binding(
            self.frozen_workflow_approval,
            self.stationary_point_policy,
            plan=scientific_workflow_plan,
        )
        self.preview_server = str(preview_server)
        self.execution_server = str(execution_server)
        self.execution_server_file_sha256 = str(
            execution_server_file_sha256
        )
        if self.execution_server_file_sha256:
            require_sha256(
                self.execution_server_file_sha256,
                "execution_server_file_sha256",
            )
        self.execution_environment = {
            str(key): str(value)
            for key, value in execution_environment.items()
        }
        self.execution_environment_remove = tuple(execution_environment_remove)
        if self.execution_environment_remove != tuple(
            sorted(set(self.execution_environment_remove))
        ) or any(
            not str(key).strip() for key in self.execution_environment_remove
        ):
            raise ContractError(
                "execution environment removals must be sorted, unique labels"
            )
        if (dependency_context is None) != (
            dependency_context_selection_receipt is None
        ):
            raise ContractError(
                "dependency context and selection receipt must be supplied together"
            )
        if dependency_context is not None:
            if (
                dependency_context.context_sha256
                != dependency_context_selection_receipt.context_sha256
                or dependency_context.policy_sha256
                != dependency_context_selection_receipt.policy_sha256
                or dependency_context.plan_sha256
                != dependency_context_selection_receipt.plan_sha256
            ):
                raise ContractError(
                    "dependency context differs from its selection receipt"
                )
        self.dependency_context = dependency_context
        self.dependency_context_selection_receipt = (
            dependency_context_selection_receipt
        )
        if self.surface.profile == "command_compiled_approved_execution":
            if self.approved_workspace is None:
                raise ContractError(
                    "execution profile requires an approved workspace"
                )
            if self.execution_resources is None:
                raise ContractError(
                    "execution profile requires host-owned resources"
                )
            if (
                self.workflow_execution_approval is None
                and self.bounded_execution_envelope is None
            ):
                raise ContractError(
                    "execution profile requires workflow approval or bounded envelope"
                )
            if (
                self.workflow_execution_approval is not None
                and self.bounded_execution_envelope is not None
            ):
                raise ContractError(
                    "workflow approval and bounded envelope are mutually exclusive"
                )
            if self.bounded_execution_envelope is not None and (
                self.execution_resources.resource_sha256
                != self.bounded_execution_envelope.resources.resource_sha256
            ):
                raise ContractError(
                    "execution resources differ from bounded envelope"
                )
            if self.frozen_workflow_approval is not None:
                if self.workflow_execution_approval is None:
                    raise ContractError(
                        "frozen workflow approval requires its V1 approval"
                    )
                if (
                    self.frozen_workflow_approval.workflow_id
                    != self.workflow_execution_approval.workflow_id
                    or self.frozen_workflow_approval.task_spec_sha256
                    != self.workflow_execution_approval.task_spec_sha256
                    or self.frozen_workflow_approval.resource_sha256
                    != self.workflow_execution_approval.resource_sha256
                ):
                    raise ContractError(
                        "frozen workflow approval differs from V1 approval"
                    )
            if self.workflow_execution_approval is not None and (
                Path(self.workflow_execution_approval.workspace).resolve()
                != self.approved_workspace
            ):
                raise ContractError(
                    "workflow approval targets another workspace"
                )
            execution_evidence_sha256 = canonical_sha256(
                {
                    "approval_sha256": (
                        self.workflow_execution_approval.approval_sha256
                        if self.workflow_execution_approval is not None
                        else ""
                    ),
                    "bounded_execution": (
                        self.bounded_execution_envelope.public_record()
                        if self.bounded_execution_envelope is not None
                        else {}
                    ),
                    "resource_sha256": self.execution_resources.resource_sha256,
                    "compute_environment_receipts": tuple(
                        sorted(
                            item.evidence_sha256
                            for item in compute_environment_receipts
                        )
                    ),
                    "execution_server": self.execution_server,
                    "execution_environment": self.execution_environment,
                }
            )
            approved_nodes = (
                (
                    (item.program, item.jobtype, item.engine)
                    for item in self.workflow_execution_approval.node_bindings
                )
                if self.workflow_execution_approval is not None
                else self._bounded_overlay_nodes(
                    preview_overlay=preview_overlay
                )
            )
            evidence_overlay = build_approved_execution_overlay(
                registry=self.registry,
                preview_overlay=preview_overlay,
                approved_nodes=approved_nodes,
                execution_evidence_sha256=execution_evidence_sha256,
            )
        else:
            evidence_overlay = preview_overlay
        if support_overlay is not None:
            if (
                support_overlay.base_registry_sha256
                != self.registry.registry_sha256
            ):
                raise ContractError(
                    "injected support overlay uses another registry"
                )
            if (
                support_overlay.overlay_sha256
                != evidence_overlay.overlay_sha256
            ):
                raise ContractError(
                    "injected support overlay lacks matching host evidence"
                )
        self.overlay = support_overlay or evidence_overlay
        self.scientific_identities = dict(scientific_identities)
        #: Advisory skill documents the model read this session, keyed by
        #: document digest so replay can reconstruct the exact text.
        self.consulted_skills = {}
        self.approved_molecular_identities = dict(
            approved_molecular_identities
        )
        if any(
            key != value.identity_sha256
            for key, value in self.approved_molecular_identities.items()
        ):
            raise ContractError(
                "approved molecular identity registry key mismatch"
            )
        self.environment_targets = tuple(environment_targets)
        self.compute_environment_receipts = tuple(compute_environment_receipts)
        self.settings_objects = dict(settings_objects)
        self.run_receipts = {
            key: dict(value) for key, value in run_receipts.items()
        }
        self.scientific_claim_evidence = dict(scientific_claim_evidence)
        self.functional_equivalence_receipts = dict(
            functional_equivalence_receipts
        )
        self.substitution_approvals = dict(substitution_approvals)
        if any(
            key != value.claim_sha256
            for key, value in self.scientific_claim_evidence.items()
        ):
            raise ContractError("claim evidence registry key mismatch")
        if any(
            key != value.receipt_sha256
            for key, value in self.functional_equivalence_receipts.items()
        ):
            raise ContractError("equivalence receipt registry key mismatch")
        if any(
            key != value.substitution_request_sha256
            for key, value in self.substitution_approvals.items()
        ):
            raise ContractError("substitution approval registry key mismatch")
        self.capabilities: dict[str, CapabilityQueryReceiptV1] = dict(
            capability_receipts
        )
        self.environments: dict[str, EnvironmentCapabilityReceiptV1] = dict(
            environment_receipts
        )
        self.program_bindings: dict[str, ResolvedProgramBindingV1] = dict(
            program_binding_receipts
        )
        self.engine_bindings: dict[str, ResolvedEngineBindingV1] = dict(
            engine_binding_receipts
        )
        self.substitutions: dict[str, ProgramSubstitutionReceiptV1] = {}
        self.project_documents: dict[str, ProjectDocumentV1] = {}
        self.project_renders: dict[str, ProjectRenderReceiptV1] = {}
        self.project_validations: dict[str, ProjectValidationReceiptV1] = dict(
            project_validation_receipts
        )
        self.functional_resolutions: dict[
            str, PySCFFunctionalResolutionReceiptV1
        ] = {}
        for validation in self.project_validations.values():
            for resolution in project_scientific_materializations(validation):
                self.functional_resolutions[resolution.receipt_sha256] = (
                    resolution
                )
        self.result_functional_evidence = {
            str(key): dict(value)
            for key, value in result_functional_evidence.items()
        }
        for key in self.result_functional_evidence:
            require_sha256(key, "result functional evidence receipt")
        self.analysis_completion_policy = analysis_completion_policy
        if (
            self.analysis_completion_policy is not None
            and self.analysis_completion_policy.task_spec_sha256
            not in self.task_spec_sha256s
        ):
            raise ContractError(
                "analysis completion policy targets another task spec"
            )
        self.invocations: dict[str, CanonicalCommandInvocationV1] = {}
        self.command_inspections: dict[str, CommandInspectionReceiptV1] = {}
        self.safe_previews: dict[str, SafePreviewReceiptV1] = {}
        self.validators: dict[str, ProgramValidatorReceiptV1] = {}
        self.preflights: dict[str, ProgramNodePreflightReceiptV1] = {}
        self.result_inspections: dict[
            str, GeneratedArtifactInspectionReceiptV1
        ] = {}
        self.quantity_extractions: dict[str, QuantityExtractionReceiptV1] = {}
        self.quantity_extraction_selectors: dict[str, tuple[str, ...]] = {}
        self.quantity_extraction_bindings: dict[str, dict[str, str]] = {}
        self.thermochemistry_receipts: dict[str, ThermochemistryReceiptV1] = {}
        self.quantity_expression_receipts: dict[str, Any] = {}
        self.quantity_expression_requests: dict[
            str, QuantityExpressionRequestV1
        ] = {}
        self.scientific_validation_receipts: dict[
            str, ScientificValidationReceiptV1
        ] = {}
        self.analysis_claim_records: dict[str, Any] = {}
        self.analysis_completion_receipts: dict[str, Any] = {}
        self.counterexamples: dict[str, CommandCounterexampleV1] = {}
        self.workflow_drafts: dict[str, CommandWorkflowDraftV1] = {}
        self.scientific_toolchain_plans: dict[
            str, ScientificToolchainPlanV1
        ] = {}
        self._scientific_toolchain_command_results: dict[
            str, dict[str, Any]
        ] = {}
        self._latest_program_workflows: dict[
            str, _ResolvedProgramWorkflow
        ] = {}
        #: Last accepted scientific plan per workflow, so a later repair can be
        #: checked against the question and not only against its own findings.
        self.scientific_plans: dict[str, Any] = {}
        self.scientific_workflow_plans: dict[
            str, ScientificWorkflowPlanV2
        ] = {}
        self.materialized_workflows: dict[str, MaterializedWorkflowV1] = {}
        if materialized_workflow is not None:
            self.materialized_workflows[
                materialized_workflow.materialized_sha256
            ] = materialized_workflow
        if scientific_workflow_plan is not None:
            self.scientific_workflow_plans[
                scientific_workflow_plan.plan_sha256
            ] = scientific_workflow_plan
        self.project_promotions: dict[str, ProjectArtifactPromotionV1] = {}
        self.scientific_decisions: dict[str, ScientificDecisionRecordV1] = {}
        self.execution_receipts: dict[str, ProgramExecutionReceiptV1] = {}
        self.result_validation_receipts: dict[
            str, ProgramResultValidationReceiptV1
        ] = {}
        self.handoffs: dict[str, OptimizedGeometryHandoffV1] = {}
        self.hessian_handoffs: dict[str, ORCAHessianHandoffV1] = {}
        self._command_contexts: dict[str, _CommandContext] = {}
        # A node ID is unique only inside one workflow.  Bind every command
        # prepared through the workflow surface to that exact plan so a later
        # diagnostic workflow with the same local node name cannot borrow an
        # older invocation or silently change approval ownership.
        self._invocation_workflow_plan_sha256s: dict[str, str] = {}
        self._completion_sets: dict[str, tuple[str, ...]] = {}
        self._latest_environment_by_capability: dict[
            str, EnvironmentCapabilityReceiptV1
        ] = {}
        self._preflight_by_node: dict[str, ProgramNodePreflightReceiptV1] = {}
        _require_registry_keys(
            self.capabilities, "receipt_sha256", "capability receipt"
        )
        _require_registry_keys(
            self.environments, "receipt_sha256", "environment receipt"
        )
        _require_registry_keys(
            self.program_bindings, "binding_sha256", "program binding"
        )
        _require_registry_keys(
            self.engine_bindings, "binding_sha256", "engine binding"
        )
        _require_registry_keys(
            self.project_validations,
            "receipt_sha256",
            "project validation",
        )
        for receipt in self.capabilities.values():
            if (
                receipt.registry_sha256 != self.registry.registry_sha256
                or receipt.live_cli_schema_sha256
                != self.live_schema.schema_sha256
                or receipt.overlay_sha256 != self.overlay.overlay_sha256
            ):
                raise ContractError("seeded capability receipt is stale")
        for environment in self.environments.values():
            self._latest_environment_by_capability[
                environment.capability_receipt_sha256
            ] = environment
        self._rehydrate_analysis_event_records()

    def _bounded_overlay_nodes(
        self,
        *,
        preview_overlay: ProgramSupportOverlayV1 | None = None,
    ):
        """Return preview-conformant executable pairs allowed by the envelope.

        The envelope is an operating ceiling, not evidence that every allowed
        program compiled and previewed successfully in this session.  A broad
        target-host envelope may therefore name programs whose bootstrap
        conformance is currently red.  Keep those programs reference-only and
        expose execution only for the allowed pairs that already have green
        preview evidence.
        """

        envelope = getattr(self, "bounded_execution_envelope", None)
        if envelope is None:  # pragma: no cover - caller narrows this
            return ()
        nodes = []
        for program, engines in envelope.allowed_program_engines:
            capability = self.registry.get(program)
            if capability is None:
                raise ContractError(
                    f"bounded execution allows unknown program {program!r}"
                )
            allowed_engines = set(engines)
            preview_pairs = None
            if preview_overlay is not None:
                preview_rule = preview_overlay.get(program)
                if (
                    preview_rule is None
                    or preview_rule.support_level
                    is not SupportLevel.PREVIEW_ONLY
                ):
                    continue
                preview_pairs = set(
                    preview_rule.allowed_engine_job_pairs
                    or (
                        (engine, jobtype)
                        for engine in preview_rule.allowed_engines
                        for jobtype in preview_rule.allowed_jobtypes
                    )
                )
            nodes.extend(
                (program, jobtype, engine)
                for engine, jobtype in capability.execution_engine_job_pairs
                if engine in allowed_engines
                and (
                    preview_pairs is None or (engine, jobtype) in preview_pairs
                )
            )
        if not nodes:
            raise ContractError(
                "bounded execution allowlist contains no preview-conformant "
                "executable engine/job pair"
            )
        return tuple(sorted(set(nodes)))

    def _rehydrate_analysis_event_records(self) -> None:
        """Restore typed postprocessing state from canonical Runtime V2 events.

        Historical lightweight events remain replayable but cannot be resumed
        as typed analysis state.  New events carry canonical records and are
        reconstructed through the same validating dataclasses used at write
        time; malformed or substituted records fail closed during host startup.
        """

        for event in self.event_store.read_events():
            record = event.payload.get("record")
            if not isinstance(record, Mapping):
                continue
            receipt_sha256 = str(event.payload.get("receipt_sha256") or "")
            if event.kind == EventKind.RESULT_QUANTITIES_EXTRACTED.value:
                receipt = quantity_extraction_receipt_from_record(
                    record, receipt_sha256=receipt_sha256
                )
                self.quantity_extractions[receipt_sha256] = receipt
                bindings = event.payload.get("selector_bindings") or {}
                if not isinstance(bindings, Mapping):
                    raise ContractError(
                        "persisted extraction selector bindings are invalid"
                    )
                normalized_bindings = {
                    str(quantity_id): str(selector)
                    for quantity_id, selector in bindings.items()
                }
                self.quantity_extraction_bindings[receipt_sha256] = (
                    normalized_bindings
                )
                self.quantity_extraction_selectors[receipt_sha256] = tuple(
                    sorted(set(normalized_bindings.values()))
                )
            elif event.kind == EventKind.THERMOCHEMISTRY_DERIVED.value:
                receipt = thermochemistry_receipt_from_record(
                    record, receipt_sha256=receipt_sha256
                )
                self.thermochemistry_receipts[receipt_sha256] = receipt
            elif event.kind == EventKind.QUANTITY_EXPRESSION_EVALUATED.value:
                receipt = quantity_expression_receipt_from_record(
                    record, receipt_sha256=receipt_sha256
                )
                self.quantity_expression_receipts[receipt_sha256] = receipt
            elif event.kind == EventKind.SCIENTIFIC_VALIDATION_EVALUATED.value:
                receipt = scientific_validation_receipt_from_record(
                    record, receipt_sha256=receipt_sha256
                )
                self.scientific_validation_receipts[receipt_sha256] = receipt
            elif event.kind == EventKind.ANALYSIS_CLAIMS_RECORDED.value:
                receipt = analysis_claim_record_from_record(
                    dict(record), receipt_sha256=receipt_sha256
                )
                self.analysis_claim_records[receipt_sha256] = receipt
            elif event.kind == EventKind.ANALYSIS_COMPLETION_EVALUATED.value:
                values = dict(record)
                values["source_receipt_sha256s"] = tuple(
                    values.get("source_receipt_sha256s") or ()
                )
                values["findings"] = tuple(values.get("findings") or ())
                receipt = AnalysisCompletionReceiptV1(
                    **values, receipt_sha256=receipt_sha256
                )
                self.analysis_completion_receipts[receipt_sha256] = receipt
            elif event.kind == EventKind.SCIENTIFIC_DECISION_RECORDED.value:
                values = dict(record)
                for field in (
                    "stage_order",
                    "assumptions",
                    "alternatives",
                    "uncertainties",
                    "diagnostics",
                    "evidence_refs",
                ):
                    values[field] = tuple(values.get(field) or ())
                decision = ScientificDecisionRecordV1(
                    **values, record_sha256=receipt_sha256
                )
                self.scientific_decisions[receipt_sha256] = decision

    def record_seeded_evidence(self, turn_id: str) -> None:
        """Persist host-prebound evidence before any model action."""

        for receipt in self.capabilities.values():
            self._emit(
                turn_id,
                EventKind.CAPABILITY_QUERIED,
                receipt.receipt_sha256,
                status=receipt.status.value,
                program=receipt.query.program,
                jobtype=receipt.query.jobtype,
                engine=receipt.query.engine,
            )
        for receipt in self.environments.values():
            self._emit(
                turn_id,
                EventKind.ENVIRONMENT_QUERIED,
                receipt.receipt_sha256,
                status=receipt.status.value,
                program=receipt.program,
                engine=receipt.engine,
            )
        for binding in self.program_bindings.values():
            self._emit_binding(turn_id, EventKind.PROGRAM_BOUND, binding)
        for binding in self.engine_bindings.values():
            self._emit_binding(turn_id, EventKind.ENGINE_BOUND, binding)
        for receipt in self.project_validations.values():
            self._emit(
                turn_id,
                EventKind.PROJECT_VALIDATED,
                receipt.receipt_sha256,
                status=receipt.status,
                program=receipt.program,
                jobtype=receipt.jobtype,
            )
        if self.dependency_context_selection_receipt is not None:
            receipt = self.dependency_context_selection_receipt
            self._emit(
                turn_id,
                EventKind.TASK_DEPENDENCY_CONTEXT_SELECTED,
                receipt.receipt_sha256,
                status=receipt.status,
                workflow_id=receipt.workflow_id,
                target_node_id=receipt.target_node_id,
                policy_sha256=receipt.policy_sha256,
                context_sha256=receipt.context_sha256,
            )

    def register_counterexample(self, value: CommandCounterexampleV1) -> None:
        """Register a deterministic counterexample; never model-authored here."""

        self.counterexamples[value.counterexample_id] = value

    def dispatch(
        self, *, turn_id: str, tool_name: str, arguments: Mapping[str, Any]
    ) -> dict[str, Any]:
        """Validate a call and invoke exactly one approved host operation."""

        values = dict(arguments)
        _validate_tool_arguments(self.surface, tool_name, values)
        handlers = {
            "inspect_program_capability": self._inspect_program_capability,
            "inspect_program_environment": self._inspect_program_environment,
            "assess_program_candidate": self._assess_program_candidate,
            "render_project_yaml": self._render_project_yaml,
            "promote_project_yaml": self._promote_project_yaml,
            "bind_scientific_identity": self._bind_scientific_identity,
            "read_project_yaml": self._read_project_yaml,
            "validate_project_yaml": self._validate_project_yaml,
            "plan_command_workflow": self._plan_command_workflow,
            "plan_scientific_workflow": self._plan_scientific_workflow,
            "rebind_scientific_workflow_projects": (
                self._rebind_scientific_workflow_projects
            ),
            "inspect_workflow_frontier": self._inspect_workflow_frontier,
            "prepare_program_node": self._prepare_program_node,
            "synthesize_command": self._synthesize_command,
            "repair_command": self._repair_command,
            "preview_command": self._preview_command,
            "preflight_program_node": self._preflight_program_node,
            "inspect_calculation_artifact": self._inspect_calculation_artifact,
            "extract_result_quantities": self._extract_result_quantities,
            "derive_thermochemistry": self._derive_thermochemistry,
            "evaluate_quantity_expression": self._evaluate_quantity_expression,
            "evaluate_scientific_validation": (
                self._evaluate_scientific_validation
            ),
            "record_analysis_claims": self._record_analysis_claims,
            "record_scientific_decision": self._record_scientific_decision,
            "execute_approved_program_node": self._execute_approved_program_node,
            "consult_domain_skill": self._consult_domain_skill,
        }
        handler = handlers.get(tool_name)
        if handler is None:
            raise ContractError("tool is absent from command-compiled profile")
        result = handler(turn_id, values)
        return {
            "schema_version": "chemsmart.tool-result.v1",
            "tool": tool_name,
            "status": "ok",
            "result": _model_visible_data(canonical_data(result)),
        }

    def execution_wait_timeout_seconds(self) -> float:
        """Return the bounded wait advertised before an engine launch."""

        if self.execution_resources is None:
            raise ContractError("execution wait requires host-owned resources")
        return self._require_bounded_launch_budget()

    def _resolve_task_spec_reference(
        self, values: Mapping[str, Any], field_name: str
    ) -> str:
        """Bind an omitted reference only when the host task is unambiguous."""

        if field_name in values:
            supplied = str(values[field_name]).strip()
            if not supplied:
                raise ContractError(f"{field_name} must not be empty")
            return supplied
        if len(self.task_spec_sha256s) == 1:
            return next(iter(self.task_spec_sha256s))
        if not self.task_spec_sha256s:
            raise ContractError("the host has no active task spec")
        raise ContractError(
            f"{field_name} is required when multiple task specs are active"
        )

    def _consult_domain_skill(self, turn_id: str, values: dict) -> Any:
        """Return one advisory skill body with its digests.

        The body is knowledge, not evidence: it carries no readiness, no
        approval, no terminal state, and no accuracy claim.  The digests make
        the exact text the model read reconstructible from replay.
        """

        del turn_id
        skill_id = require_identifier(values["skill_id"], "skill_id")
        document = resolve_skill(skill_id)
        if document is None:
            raise ContractError(f"unknown domain skill: {skill_id}")
        self.consulted_skills[document.document_sha256] = document
        return {
            "schema_version": "chemsmart.domain-skill-consultation.v1",
            "skill_id": document.skill_id,
            "skill_version": document.skill_version,
            "origin": document.origin,
            "description": document.description,
            "body": document.body,
            "body_sha256": document.body_sha256,
            "document_sha256": document.document_sha256,
            "advisory_only": True,
            "readiness_authority": False,
            "accuracy_authority": False,
        }

    def _bind_scientific_identity(self, turn_id: str, values: dict) -> Any:
        task_spec_sha256 = self._resolve_task_spec_reference(
            values, "task_spec_sha256"
        )
        if task_spec_sha256 not in self.task_spec_sha256s:
            raise ContractError(
                "scientific identity targets an unknown task spec"
            )
        artifact = self._artifact(values["input_artifact_id"])
        binding = build_scientific_identity_binding(
            task_spec_sha256=task_spec_sha256,
            geometry_artifact=artifact,
            charge=values["charge"],
            multiplicity=values["multiplicity"],
        )
        self.scientific_identities[binding.binding_sha256] = binding
        return binding

    def _inspect_program_capability(self, turn_id: str, values: dict) -> Any:
        receipt = query_capability(
            CapabilityQueryV1(
                program=values["program"],
                jobtype=values["jobtype"],
                engine=values["engine"],
            ),
            registry=self.registry,
            overlay=self.overlay,
            live_schema=self.live_schema,
        )
        self.capabilities[receipt.receipt_sha256] = receipt
        self._emit(
            turn_id,
            EventKind.CAPABILITY_QUERIED,
            receipt.receipt_sha256,
            status=receipt.status.value,
            program=receipt.query.program,
            jobtype=receipt.query.jobtype,
            engine=receipt.query.engine,
        )
        return receipt

    def _inspect_program_environment(self, turn_id: str, values: dict) -> Any:
        capability = self._get(
            self.capabilities,
            values["capability_receipt_sha256"],
            "capability receipt",
        )
        environment = query_environment(
            capability,
            targets=self.environment_targets,
            compute_receipts=self.compute_environment_receipts,
        )
        program_binding = resolve_program_binding(capability)
        engine_binding = resolve_engine_binding(program_binding, environment)
        self.environments[environment.receipt_sha256] = environment
        self.program_bindings[program_binding.binding_sha256] = program_binding
        self.engine_bindings[engine_binding.binding_sha256] = engine_binding
        self._latest_environment_by_capability[capability.receipt_sha256] = (
            environment
        )
        self._emit(
            turn_id,
            EventKind.ENVIRONMENT_QUERIED,
            environment.receipt_sha256,
            status=environment.status.value,
            program=environment.program,
            engine=environment.engine,
        )
        self._emit_binding(turn_id, EventKind.PROGRAM_BOUND, program_binding)
        self._emit_binding(turn_id, EventKind.ENGINE_BOUND, engine_binding)
        return {
            "environment": environment,
            "program_binding": program_binding,
            "engine_binding": engine_binding,
        }

    def _assess_program_candidate(self, turn_id: str, values: dict) -> Any:
        capability = self._get(
            self.capabilities,
            values["capability_receipt_sha256"],
            "capability receipt",
        )
        claim_evidence = tuple(
            self._get(
                self.scientific_claim_evidence,
                digest,
                "scientific claim evidence",
            )
            for digest in values["source_claim_sha256s"]
        )
        verified_claims = tuple(
            item for item in claim_evidence if item.status == "verified"
        )
        identity = values["requested_program"] == values["selected_program"]
        equivalence_digest = values.get(
            "functional_equivalence_receipt_sha256", ""
        )
        equivalence = (
            self._get(
                self.functional_equivalence_receipts,
                equivalence_digest,
                "functional equivalence receipt",
            )
            if equivalence_digest
            else None
        )
        if not identity and equivalence is not None:
            if (
                equivalence.requested_program != values["requested_program"]
                or equivalence.selected_program != values["selected_program"]
                or equivalence.method_name != values["method_name"]
            ):
                raise ContractError(
                    "functional equivalence targets another request"
                )
            observed_claim_receipts = tuple(
                sorted(item.receipt_sha256 for item in verified_claims)
            )
            if (
                equivalence.status == "verified"
                and equivalence.claim_evidence_receipt_sha256s
                != observed_claim_receipts
            ):
                raise ContractError(
                    "functional equivalence uses other claim evidence"
                )
        request = build_program_substitution_request(
            request_id=values["request_id"],
            requested_program=values["requested_program"],
            selected_program=values["selected_program"],
            requested_engine=values["requested_engine"],
            selected_engine=values["selected_engine"],
            job_families=values["job_families"],
            method_family=values["method_family"],
            method_name=values["method_name"],
            basis_mode=values["basis_mode"],
            constraint_kinds=values["constraint_kinds"],
            requires_post_hf=values["requires_post_hf"],
            requires_double_hybrid=values["requires_double_hybrid"],
            functional_semantics_confirmed=bool(
                identity or (equivalence and equivalence.status == "verified")
            ),
            source_claim_sha256s=tuple(
                item.claim_sha256 for item in verified_claims
            ),
        )
        approval = self.substitution_approvals.get(request.request_sha256)
        approval_ref = ""
        if approval is not None:
            if (
                approval.substitution_request_sha256 != request.request_sha256
                or approval.decision != "approved"
            ):
                raise ContractError("substitution approval is stale or red")
            approval_ref = approval.receipt_sha256
        receipt = assess_typed_program_substitution(
            request, capability, approval_ref=approval_ref
        )
        self.substitutions[receipt.receipt_sha256] = receipt
        self._emit(
            turn_id,
            EventKind.SUBSTITUTION_ASSESSED,
            receipt.receipt_sha256,
            decision=receipt.decision,
            requested_program=receipt.requested_program,
            selected_program=receipt.selected_program,
        )
        if receipt.decision in {"exact", "approved"}:
            program_binding = resolve_program_binding(
                capability,
                requested_program=request.requested_program,
                substitution_receipt_sha256=receipt.receipt_sha256,
            )
            self.program_bindings[program_binding.binding_sha256] = (
                program_binding
            )
            self._emit_binding(
                turn_id, EventKind.PROGRAM_BOUND, program_binding
            )
            environment = self._latest_environment_by_capability.get(
                capability.receipt_sha256
            )
            if environment is not None:
                engine_binding = resolve_engine_binding(
                    program_binding, environment
                )
                self.engine_bindings[engine_binding.binding_sha256] = (
                    engine_binding
                )
                self._emit_binding(
                    turn_id, EventKind.ENGINE_BOUND, engine_binding
                )
        return receipt

    def _render_project_yaml(self, turn_id: str, values: dict) -> Any:
        document = project_document(
            program=values["program"], sections=values["sections"]
        )
        receipt = render_project_yaml(document, registry=self.registry)
        self.project_documents[document.document_sha256] = document
        self.project_renders[receipt.receipt_sha256] = receipt
        return receipt

    def _promote_project_yaml(self, turn_id: str, values: dict) -> Any:
        if self.approved_workspace is None:
            raise ContractError("project promotion requires a task workspace")
        if values["artifact_id"] in self.artifacts:
            # Naming the taken IDs is what makes this actionable. Without them
            # a caller can only guess, and a live run collided five times in a
            # row on the same message.
            taken = sorted(self.artifacts)
            requested = values["artifact_id"]
            raise ContractError(
                f"artifact ID {requested!r} is already registered; "
                f"choose one not in {taken}"
            )
        render = self._get(
            self.project_renders,
            values["render_receipt_sha256"],
            "project render receipt",
        )
        artifact, promotion = promote_project_candidate(
            render,
            approved_workspace=self.approved_workspace,
            artifact_id=values["artifact_id"],
        )
        self.artifacts[artifact.artifact_id] = artifact
        self.project_promotions[artifact.artifact_id] = promotion
        self._emit(
            turn_id,
            EventKind.PROJECT_PROMOTED,
            promotion.receipt_sha256,
            status=promotion.status,
            artifact_id=artifact.artifact_id,
        )
        return {"artifact": artifact, "promotion": promotion}

    def _read_project_yaml(self, turn_id: str, values: dict) -> Any:
        artifact = self._artifact(values["project_artifact_id"])
        document = read_project_yaml(artifact, program=values["program"])
        self.project_documents[document.document_sha256] = document
        return document

    def _validate_project_yaml(self, turn_id: str, values: dict) -> Any:
        artifact = self._artifact(values["project_artifact_id"])
        capability = self._get(
            self.capabilities,
            values["capability_receipt_sha256"],
            "capability receipt",
        )
        receipt = validate_project_yaml(artifact, capability=capability)
        document = read_project_yaml(
            artifact, program=capability.query.program
        )
        section_application = project_section_application_observation(
            document,
            jobtype=capability.query.jobtype,
            applied_settings=dict(receipt.settings),
        )
        self.project_validations[receipt.receipt_sha256] = receipt
        materializations = project_scientific_materializations(receipt)
        for materialization in materializations:
            self.functional_resolutions[materialization.receipt_sha256] = (
                materialization
            )
        promotion = self.project_promotions.get(artifact.artifact_id)
        if promotion is not None and promotion.validation_status == "pending":
            self.project_promotions[artifact.artifact_id] = (
                bind_project_promotion_validation(promotion, artifact, receipt)
            )
            bound = self.project_promotions[artifact.artifact_id]
            self._emit(
                turn_id,
                EventKind.PROJECT_PROMOTED,
                bound.receipt_sha256,
                status=bound.status,
                artifact_id=artifact.artifact_id,
            )
        self._emit(
            turn_id,
            EventKind.PROJECT_VALIDATED,
            receipt.receipt_sha256,
            status=receipt.status,
            program=receipt.program,
            jobtype=receipt.jobtype,
        )
        effective_settings = dict(receipt.settings)
        frequency_modes = tuple(
            name
            for name in ("freq", "numfreq", "vpt2")
            if bool(effective_settings.get(name))
        )
        if frequency_modes:
            frequency_semantics = {
                "requested": True,
                "modes": frequency_modes,
                "produces_observable": "vibrational_frequencies",
                "planning_guidance": (
                    f"This {receipt.program}/{receipt.jobtype} project already "
                    "requests a frequency calculation. Declare "
                    "vibrational_frequencies on this same scientific node; "
                    "do not add another frequency node at the same method and "
                    "geometry unless an independent Hessian is scientifically "
                    "intended."
                ),
            }
        else:
            frequency_semantics = {
                "requested": False,
                "modes": (),
                "produces_observable": "",
                "planning_guidance": (
                    f"This {receipt.program}/{receipt.jobtype} project does "
                    "not request a frequency calculation. A distinct frequency "
                    "node is needed when vibrational evidence is required."
                ),
            }
        return {
            **canonical_data(receipt),
            "scientific_materializations": tuple(
                item.public_record() for item in materializations
            ),
            "decision_binding": _scientific_decision_binding_requirement(
                materializations
            ),
            "section_application": section_application,
            "effective_frequency_semantics": frequency_semantics,
            "workflow_binding": self._project_workflow_binding_observation(
                artifact.artifact_id
            ),
        }

    def _project_workflow_binding_observation(
        self, project_artifact_id: str
    ) -> dict[str, Any]:
        """Show whether a validated project participates in the latest DAG.

        Project promotion is intentionally append-only, so a repaired YAML gets
        a new artifact ID.  Without this small relational observation a model
        can validate the replacement yet finish with an older project still
        named by the workflow.  The observation contains only public host state
        and contains no predetermined scientific answer.
        """

        latest_by_workflow: dict[str, ScientificToolchainPlanV1] = {}
        for plan in self.scientific_toolchain_plans.values():
            latest_by_workflow[plan.workflow_id] = plan
        bindings: list[dict[str, Any]] = []
        active_roles: dict[str, tuple[str, ...]] = {}
        for workflow_id, plan in sorted(latest_by_workflow.items()):
            command_result = self._scientific_toolchain_command_results[
                plan.plan_sha256
            ]
            draft = command_result["workflow_draft"]
            roles = tuple(sorted({node.project_role for node in draft.nodes}))
            active_roles[workflow_id] = roles
            for node in draft.nodes:
                if node.project_role == project_artifact_id:
                    bindings.append(
                        {
                            "workflow_id": workflow_id,
                            "node_id": node.node_id,
                        }
                    )
        if bindings:
            return {
                "status": "bound",
                "bindings": tuple(bindings),
            }
        if active_roles:
            return {
                "status": "unbound",
                "active_project_roles": active_roles,
                "next_action": (
                    "if this project repairs an active node, call "
                    "rebind_scientific_workflow_projects"
                ),
            }
        return {"status": "no_scientific_workflow_planned"}

    def _record_scientific_decision(self, turn_id: str, values: dict) -> Any:
        task_spec_sha256 = self._resolve_task_spec_reference(
            values, "task_spec_sha256"
        )
        if task_spec_sha256 not in self.task_spec_sha256s:
            raise ContractError(
                "scientific decision targets an unknown task spec"
            )
        postprocessing_registries = (
            self.quantity_extractions,
            self.thermochemistry_receipts,
            self.quantity_expression_receipts,
            self.scientific_validation_receipts,
            self.analysis_claim_records,
        )
        postprocessing_receipt_sha256s = tuple(
            str(item)
            for item in values.get("postprocessing_receipt_sha256s", ())
        )
        for receipt_sha256 in postprocessing_receipt_sha256s:
            require_sha256(receipt_sha256, "postprocessing_receipt_sha256")
            if not any(
                receipt_sha256 in receipts
                for receipts in postprocessing_registries
            ):
                raise ContractError(
                    "scientific decision cites an unknown postprocessing receipt"
                )
        evidence_refs = tuple(values["evidence_refs"]) + tuple(
            f"receipt:{receipt_sha256}"
            for receipt_sha256 in postprocessing_receipt_sha256s
        )
        evidence_refs = tuple(dict.fromkeys(evidence_refs))
        functional_resolution_refs = set()
        for evidence_ref in evidence_refs:
            prefix = "molecular_identity:"
            reference = str(evidence_ref)
            if reference.startswith(prefix):
                identity_sha256 = reference[len(prefix) :]
                if identity_sha256 not in self.approved_molecular_identities:
                    raise ContractError(
                        "scientific decision cites an unapproved molecular identity"
                    )
                continue
            prefix = "functional_resolution:"
            if reference.startswith(prefix):
                receipt_sha256 = reference[len(prefix) :]
                if receipt_sha256 not in self.functional_resolutions:
                    raise ContractError(
                        "scientific decision cites an unknown functional resolution"
                    )
                functional_resolution_refs.add(receipt_sha256)
                continue
            parsed_postprocessing_ref = _postprocessing_evidence_reference(
                reference
            )
            if parsed_postprocessing_ref is not None:
                receipt_kind, receipt_sha256 = parsed_postprocessing_ref
                registries = {
                    "quantity_extraction": self.quantity_extractions,
                    "thermochemistry": self.thermochemistry_receipts,
                    "quantity_expression": self.quantity_expression_receipts,
                    "scientific_validation": (
                        self.scientific_validation_receipts
                    ),
                    "analysis_claim": self.analysis_claim_records,
                }
                if receipt_kind == "generic":
                    known = any(
                        receipt_sha256 in receipts
                        for receipts in registries.values()
                    )
                else:
                    known = receipt_sha256 in registries[receipt_kind]
                if not known:
                    raise ContractError(
                        "scientific decision cites an unknown "
                        "postprocessing receipt"
                    )
                continue
            prefix = "analysis_completion_policy:"
            if reference.startswith(prefix):
                policy_sha256 = reference[len(prefix) :]
                if (
                    self.analysis_completion_policy is None
                    or policy_sha256
                    != self.analysis_completion_policy.policy_sha256
                ):
                    raise ContractError(
                        "scientific decision cites an unknown analysis policy"
                    )
                continue
            prefix = "result_functional_resolution:"
            if reference.startswith(prefix):
                receipt_sha256 = reference[len(prefix) :]
                if receipt_sha256 not in self.result_functional_evidence:
                    raise ContractError(
                        "scientific decision cites unknown result functional evidence"
                    )
                functional_resolution_refs.add(receipt_sha256)
        convention_narrative = " ".join(
            (
                *values["assumptions"],
                values["method_rationale"],
                *values["uncertainties"],
                *values["diagnostics"],
            )
        )
        if (
            re.search(
                r"(?i)(?<![a-z0-9])(?:vwn\s*[35]|b3lypg|b3lyp5)(?![a-z0-9])",
                convention_narrative,
            )
            and not functional_resolution_refs
        ):
            raise ContractError(
                "functional-convention claims require a host resolution receipt"
            )
        record = build_scientific_decision_record(
            decision_id=values["decision_id"],
            task_spec_sha256=task_spec_sha256,
            stage_order=values["stage_order"],
            assumptions=values["assumptions"],
            method_rationale=values["method_rationale"],
            alternatives=values["alternatives"],
            uncertainties=values["uncertainties"],
            diagnostics=values["diagnostics"],
            evidence_refs=evidence_refs,
        )
        self.scientific_decisions[record.record_sha256] = record
        record_body = canonical_data(record)
        record_body.pop("record_sha256")
        self._emit(
            turn_id,
            EventKind.SCIENTIFIC_DECISION_RECORDED,
            record.record_sha256,
            status="recorded",
            decision_id=record.decision_id,
            record=record_body,
        )
        return record

    def _plan_command_workflow(
        self,
        turn_id: str,
        values: dict,
        *,
        node_annotations: Mapping[str, Mapping[str, Any]] | None = None,
    ) -> Any:
        """Record a broad DAG before execution-grade evidence is available."""

        nodes = []
        findings: list[dict[str, str]] = []
        declared_programs = {
            item.program: item for item in self.registry.programs
        }
        for raw_node in values["nodes"]:
            inputs = tuple(
                sorted(
                    (
                        ArtifactInputIntentV1(
                            binding_id=item["binding_id"],
                            artifact_class=item["artifact_class"],
                            artifact_id=item.get("artifact_id", ""),
                            producer_node_id=item["producer_node_id"],
                            producer_output_id=item["producer_output_id"],
                        )
                        for item in raw_node["inputs"]
                    ),
                    key=lambda item: item.binding_id,
                )
            )
            outputs = tuple(
                sorted(
                    (
                        ArtifactOutputIntentV1(
                            output_id=item["output_id"],
                            artifact_class=item["artifact_class"],
                        )
                        for item in raw_node["expected_outputs"]
                    ),
                    key=lambda item: item.output_id,
                )
            )
            node = CommandNodeIntentV1(
                node_id=raw_node["node_id"],
                program=raw_node["program"],
                jobtype=raw_node["jobtype"],
                project_role=raw_node["project_role"],
                dependencies=tuple(sorted(set(raw_node["dependencies"]))),
                inputs=inputs,
                expected_outputs=outputs,
                unresolved_fields=tuple(
                    sorted(set(raw_node["unresolved_fields"]))
                ),
                node_kind=raw_node.get("node_kind", "program_call"),
                charge=raw_node.get("charge"),
                multiplicity=raw_node.get("multiplicity"),
            )
            if node.node_kind == "aggregate":
                # ChemSmart performs the arithmetic, so there is no program
                # capability to check. The contract already restricted the
                # stage, and the operations live in the expression itself.
                nodes.append(node)
                continue
            capability = declared_programs.get(node.program)
            if capability is None:
                findings.append(
                    {
                        "node_id": node.node_id,
                        "rule_id": "workflow.program.not_declared",
                    }
                )
            elif node.jobtype not in capability.jobtypes:
                findings.append(
                    {
                        "node_id": node.node_id,
                        "rule_id": "workflow.job.not_declared",
                    }
                )
            for item in node.inputs:
                if not item.producer_node_id and (
                    not item.artifact_id
                    or item.artifact_id not in self.artifacts
                ):
                    findings.append(
                        {
                            "node_id": node.node_id,
                            "rule_id": "workflow.input.unresolved",
                        }
                    )
            nodes.append(node)
        task_spec_id = self._resolve_task_spec_reference(
            values, "task_spec_id"
        )
        draft = build_command_workflow_draft(
            workflow_id=values["workflow_id"],
            task_spec_id=task_spec_id,
            nodes=tuple(nodes),
        )
        scientific_plan = (
            self._scientific_plan_from_draft(
                draft,
                findings=findings,
                node_annotations=node_annotations,
            )
            if draft.nodes
            else None
        )
        if scientific_plan is not None:
            self.scientific_workflow_plans[scientific_plan.plan_sha256] = (
                scientific_plan
            )
        context = self._workflow_context(
            draft,
            scientific_plan_sha256=(
                scientific_plan.plan_sha256 if scientific_plan else ""
            ),
        )
        finding_nodes = {item["node_id"] for item in findings} | {
            node.node_id for node in draft.nodes if node.unresolved_fields
        }
        actionable = tuple(
            node_id
            for node_id in context.ready_node_ids
            if node_id not in finding_nodes
        )
        completed = context.completed_node_ids
        unresolved_nodes = (
            set(context.waiting_node_ids)
            | set(context.blocked_node_ids)
            | finding_nodes
        )
        unresolved = tuple(
            node.node_id
            for node in draft.nodes
            if node.node_id in unresolved_nodes
        )
        self._emit(
            turn_id,
            EventKind.WORKFLOW_PLANNED,
            draft.draft_sha256,
            status="planned",
            actionable_node_ids=actionable,
            unresolved_node_ids=unresolved,
            scientific_plan_sha256=(
                scientific_plan.plan_sha256 if scientific_plan else ""
            ),
            scientific_plan_record=(
                canonical_data(scientific_plan) if scientific_plan else {}
            ),
        )
        # A draft becomes host state only after every scientific constraint
        # above accepts it and the corresponding event is durable.  In
        # particular, a replan that differs from a frozen approval raises in
        # ``_scientific_plan_from_draft``.  Registering that rejected draft
        # early made it look like the latest observed workflow even though no
        # WORKFLOW_PLANNED event existed for it, so planned termination could
        # not bind its required receipt to the event stream.
        self.workflow_drafts[draft.draft_sha256] = draft
        result = {
            "workflow_draft": draft,
            "scientific_workflow_plan": scientific_plan,
            "actionable_node_ids": actionable,
            "completed_node_ids": completed,
            "unresolved_node_ids": unresolved,
            "findings": tuple(findings),
        }
        if context is not None:
            result["workflow_context"] = context
        if scientific_plan is not None:
            result["approval_readiness"] = self._approval_readiness(
                scientific_plan
            )
        self._latest_program_workflows[draft.workflow_id] = (
            _ResolvedProgramWorkflow(
                draft=draft,
                scientific_plan=scientific_plan,
                command_result=result,
            )
        )
        return result

    def _plan_scientific_workflow(self, turn_id: str, values: dict) -> Any:
        """Plan calculations and their downstream scientific analysis together.

        The existing command planner remains the authority for calculation
        nodes.  Analysis intent is layered on its producer outputs without
        asking the model to invent future artifact or receipt hashes.
        """

        # Analysis-only workflows have no program invocations.  Treat an
        # omitted calculation list as the natural empty list instead of
        # making the model restate a documentary placeholder.
        calculation_nodes = values.get("calculation_nodes", [])
        node_annotations = {
            item["node_id"]: {
                "produces_observables": tuple(item["produces_observables"]),
                "support_state": item["support_state"],
                "blocked_reason": item["blocked_reason"],
            }
            for item in calculation_nodes
        }
        command_result = self._plan_command_workflow(
            turn_id,
            {
                "workflow_id": values["workflow_id"],
                "nodes": calculation_nodes,
                **(
                    {"task_spec_id": values["task_spec_id"]}
                    if "task_spec_id" in values
                    else {}
                ),
            },
            node_annotations=node_annotations,
        )
        draft = command_result["workflow_draft"]
        analysis_nodes = []
        for raw_node in values["analysis_nodes"]:
            analysis_kind = str(raw_node["analysis_kind"])
            artifact_id = str(raw_node.get("artifact_id", "")).strip()
            raw_inputs = tuple(raw_node["inputs"])
            if artifact_id and raw_inputs:
                raise ContractError(
                    "an analysis node must choose a registered result or a "
                    "future producer output, not both"
                )
            if artifact_id:
                artifact = self._artifact(artifact_id)
                result_program = self._analysis_result_program_for_kind(
                    artifact.kind
                )
                if (
                    analysis_kind == "thermochemistry"
                    and result_program == "xyz"
                ):
                    raise ContractError(
                        "thermochemistry requires a complete typed program "
                        "result, not a geometry-only registered artifact"
                    )
                analysis_inputs = (
                    RegisteredResultInputIntentV1(
                        input_id="registered-result",
                        artifact_id=artifact.artifact_id,
                    ),
                )
            else:
                analysis_inputs = tuple(
                    sorted(
                        (
                            AnalysisInputIntentV1(
                                input_id=item["input_id"],
                                source_kind=item["source_kind"],
                                producer_node_id=item["producer_node_id"],
                                producer_output_id=item["producer_output_id"],
                            )
                            for item in raw_inputs
                        ),
                        key=lambda item: item.input_id,
                    )
                )
            analysis_nodes.append(
                AnalysisNodeIntentV1(
                    node_id=raw_node["node_id"],
                    analysis_kind=analysis_kind,
                    dependencies=tuple(sorted(set(raw_node["dependencies"]))),
                    inputs=analysis_inputs,
                    selectors=tuple(
                        sorted(
                            (
                                AnalysisSelectorIntentV1(
                                    quantity_id=item["quantity_id"],
                                    selector=item["selector"],
                                )
                                for item in raw_node["selectors"]
                            ),
                            key=lambda item: item.quantity_id,
                        )
                    ),
                    outputs=tuple(
                        sorted(
                            (
                                AnalysisOutputIntentV1(
                                    output_id=item["output_id"],
                                    quantity_kind=item["quantity_kind"],
                                    unit=item["unit"],
                                )
                                for item in raw_node["outputs"]
                            ),
                            key=lambda item: item.output_id,
                        )
                    ),
                    expression_nodes=tuple(raw_node["expression_nodes"]),
                    expression_output_node_ids=tuple(
                        raw_node["expression_output_node_ids"]
                    ),
                    temperature_k=raw_node.get("temperature_k"),
                    pressure_atm=raw_node.get("pressure_atm"),
                    support_state=raw_node["support_state"],
                    blocked_reason=raw_node["blocked_reason"],
                    concentration_mol_l=raw_node.get("concentration_mol_l"),
                    entropy_method=raw_node.get("entropy_method", "rrho"),
                    entropy_cutoff_cm1=raw_node.get("entropy_cutoff_cm1"),
                    enthalpy_cutoff_cm1=raw_node.get("enthalpy_cutoff_cm1"),
                    alpha=raw_node.get("alpha", 4),
                    use_weighted_mass=raw_node.get(
                        "use_weighted_mass", False
                    ),
                    frequency_scale_factor=raw_node.get(
                        "frequency_scale_factor", 1.0
                    ),
                    validation_rules=tuple(
                        sorted(
                            (
                                AnalysisValidationRuleIntentV1(
                                    rule_id=item["rule_id"],
                                    predicate=item["predicate"],
                                    input_ids=tuple(
                                        sorted(set(item["input_ids"]))
                                    ),
                                    threshold=item.get("threshold"),
                                    expected_count=item.get("expected_count"),
                                    unit=item.get("unit", ""),
                                )
                                for item in raw_node.get(
                                    "validation_rules", ()
                                )
                            ),
                            key=lambda item: item.rule_id,
                        )
                    ),
                )
            )
        observables = {
            item["node_id"]: tuple(item["produces_observables"])
            for item in calculation_nodes
        }
        plan = build_scientific_toolchain_plan(
            plan_id=values["plan_id"],
            workflow_id=values["workflow_id"],
            command_workflow_draft_sha256=draft.draft_sha256,
            calculation_nodes=draft.nodes,
            calculation_observables=observables,
            analysis_nodes=analysis_nodes,
            required_output_ids=values["required_output_ids"],
        )
        self.scientific_toolchain_plans[plan.plan_sha256] = plan
        self._scientific_toolchain_command_results[plan.plan_sha256] = (
            command_result
        )
        self._bind_program_toolchain(plan, command_result)
        frontier = project_scientific_toolchain_frontier(
            plan,
            actionable_calculation_node_ids=command_result[
                "actionable_node_ids"
            ],
            unresolved_calculation_node_ids=command_result[
                "unresolved_node_ids"
            ],
            completed_calculation_node_ids=command_result.get(
                "completed_node_ids", ()
            ),
        )
        return {
            "scientific_toolchain_plan": plan,
            "calculation_plan": command_result,
            "workflow_frontier": frontier,
        }

    def _bind_program_toolchain(
        self,
        plan: ScientificToolchainPlanV1,
        command_result: Mapping[str, Any],
    ) -> None:
        """Bind analysis only to the exact command workflow it extends."""

        resolved = self._latest_program_workflows.get(plan.workflow_id)
        if resolved is None or resolved.command_result is not command_result:
            raise ContractError(
                "scientific toolchain has no exact command workflow binding"
            )
        if plan.command_workflow_draft_sha256 != resolved.draft.draft_sha256:
            raise ContractError(
                "scientific toolchain belongs to another command workflow"
            )
        scientific_plan = command_result.get("scientific_workflow_plan")
        if resolved.draft.nodes:
            if (
                not isinstance(scientific_plan, ScientificWorkflowPlanV2)
                or scientific_plan is not resolved.scientific_plan
            ):
                raise ContractError(
                    "scientific toolchain lacks its task-bound scientific plan"
                )
        elif (
            scientific_plan is not None or resolved.scientific_plan is not None
        ):
            raise ContractError(
                "analysis-only toolchain must not invent a scientific plan"
            )
        self._latest_program_workflows[plan.workflow_id] = replace(
            resolved, scientific_toolchain_plan=plan
        )

    @staticmethod
    def _analysis_result_program_for_kind(artifact_kind: str) -> str:
        """Resolve a registered result kind through the existing readers."""

        from chemsmart.analysis.result_readers import RESULT_READERS

        programs = {
            reader.program
            for reader in RESULT_READERS.values()
            if reader.artifact_kind == artifact_kind
        }
        if artifact_kind == "pyscf_hdf5":
            programs.add("pyscf")
        if len(programs) != 1:
            raise ContractError(
                "registered analysis input must be a supported result "
                f"artifact; kind {artifact_kind!r} maps to {sorted(programs)}"
            )
        return next(iter(programs))

    def _rebind_scientific_workflow_projects(
        self, turn_id: str, values: dict
    ) -> Any:
        """Clone the latest scientific DAG while changing project roles only."""

        workflow_id = values["workflow_id"]
        candidates = tuple(
            plan
            for plan in self.scientific_toolchain_plans.values()
            if plan.workflow_id == workflow_id
        )
        if not candidates:
            raise ContractError("unknown scientific workflow ID")
        current_plan = candidates[-1]
        current_result = self._scientific_toolchain_command_results[
            current_plan.plan_sha256
        ]
        current_draft = current_result["workflow_draft"]
        nodes_by_id = {node.node_id: node for node in current_draft.nodes}
        replacements = {
            item["node_id"]: item["project_role"]
            for item in values["replacements"]
        }
        if len(replacements) != len(values["replacements"]):
            raise ContractError("project rebindings must name each node once")
        unknown = sorted(set(replacements).difference(nodes_by_id))
        if unknown:
            raise ContractError(
                f"project rebindings reference unknown nodes {unknown}"
            )

        for node_id, project_role in replacements.items():
            node = nodes_by_id[node_id]
            project = self._artifact(project_role)
            if project.kind != "project_yaml":
                raise ContractError(
                    f"replacement for node {node_id!r} is not project YAML"
                )
            valid = tuple(
                receipt
                for receipt in self.project_validations.values()
                if receipt.project_artifact_id == project.artifact_id
                and receipt.project_sha256 == project.sha256
                and receipt.program == node.program
                and receipt.jobtype == node.jobtype
                and receipt.status == "valid"
            )
            if not valid:
                raise ContractError(
                    f"replacement project {project_role!r} is not validated "
                    f"for {node.program}/{node.jobtype}"
                )

        revised_nodes = tuple(
            replace(
                node,
                project_role=replacements.get(node.node_id, node.project_role),
            )
            for node in current_draft.nodes
        )
        scientific_v2 = current_result.get("scientific_workflow_plan")
        annotations = {
            node.node_id: {
                "produces_observables": node.produces_observables,
                "support_state": node.support_state,
                "blocked_reason": node.blocked_reason,
            }
            for node in (scientific_v2.nodes if scientific_v2 else ())
        }
        command_result = self._plan_command_workflow(
            turn_id,
            {
                "workflow_id": current_draft.workflow_id,
                "task_spec_id": current_draft.task_spec_id,
                "nodes": tuple(
                    {
                        "node_id": node.node_id,
                        "node_kind": node.node_kind,
                        "program": node.program,
                        "jobtype": node.jobtype,
                        "project_role": node.project_role,
                        **(
                            {
                                "charge": node.charge,
                                "multiplicity": node.multiplicity,
                            }
                            if node.charge is not None
                            else {}
                        ),
                        "dependencies": node.dependencies,
                        "inputs": tuple(
                            {
                                "binding_id": item.binding_id,
                                "artifact_id": item.artifact_id,
                                "artifact_class": item.artifact_class,
                                "producer_node_id": item.producer_node_id,
                                "producer_output_id": item.producer_output_id,
                            }
                            for item in node.inputs
                        ),
                        "expected_outputs": tuple(
                            {
                                "output_id": item.output_id,
                                "artifact_class": item.artifact_class,
                            }
                            for item in node.expected_outputs
                        ),
                        "unresolved_fields": node.unresolved_fields,
                    }
                    for node in revised_nodes
                ),
            },
            node_annotations=annotations,
        )
        draft = command_result["workflow_draft"]
        observables = dict(current_plan.calculation_observables)
        revised_plan = build_scientific_toolchain_plan(
            plan_id=current_plan.plan_id,
            workflow_id=current_plan.workflow_id,
            command_workflow_draft_sha256=draft.draft_sha256,
            calculation_nodes=draft.nodes,
            calculation_observables=observables,
            analysis_nodes=current_plan.analysis_nodes,
            required_output_ids=current_plan.required_output_ids,
        )
        self.scientific_toolchain_plans[revised_plan.plan_sha256] = (
            revised_plan
        )
        self._scientific_toolchain_command_results[
            revised_plan.plan_sha256
        ] = command_result
        self._bind_program_toolchain(revised_plan, command_result)
        return {
            "scientific_toolchain_plan": revised_plan,
            "calculation_plan": command_result,
            "workflow_frontier": project_scientific_toolchain_frontier(
                revised_plan,
                actionable_calculation_node_ids=command_result[
                    "actionable_node_ids"
                ],
                unresolved_calculation_node_ids=command_result[
                    "unresolved_node_ids"
                ],
                completed_calculation_node_ids=command_result.get(
                    "completed_node_ids", ()
                ),
            ),
            "replacements": tuple(
                {
                    "node_id": node_id,
                    "project_role": project_role,
                }
                for node_id, project_role in sorted(replacements.items())
            ),
        }

    def _resolve_program_workflow(
        self, workflow_id: str
    ) -> _ResolvedProgramWorkflow:
        """Resolve the latest workflow surface with exact host bindings."""

        resolved = self._latest_program_workflows.get(workflow_id)
        if resolved is None:
            raise ContractError("unknown scientific workflow ID")
        if (
            self.workflow_drafts.get(resolved.draft.draft_sha256)
            is not resolved.draft
            or resolved.draft.workflow_id != workflow_id
        ):
            raise ContractError(
                "workflow resolver lost its exact command draft binding"
            )
        scientific_plan = resolved.scientific_plan
        if scientific_plan is not None and (
            scientific_plan.workflow_id != workflow_id
            or self.scientific_workflow_plans.get(scientific_plan.plan_sha256)
            is not scientific_plan
        ):
            raise ContractError(
                "workflow resolver lost its exact scientific plan binding"
            )
        toolchain = resolved.scientific_toolchain_plan
        if toolchain is not None and (
            toolchain.workflow_id != workflow_id
            or toolchain.command_workflow_draft_sha256
            != resolved.draft.draft_sha256
            or self.scientific_toolchain_plans.get(toolchain.plan_sha256)
            is not toolchain
            or self._scientific_toolchain_command_results.get(
                toolchain.plan_sha256
            )
            is not resolved.command_result
        ):
            raise ContractError(
                "workflow resolver lost its exact scientific toolchain binding"
            )
        return resolved

    def _inspect_workflow_frontier(self, turn_id: str, values: dict) -> Any:
        """Return the latest connected frontier for a named workflow."""

        del turn_id
        resolved = self._resolve_program_workflow(values["workflow_id"])
        command_result = resolved.command_result
        draft = resolved.draft
        scientific_v2 = resolved.scientific_plan
        context = self._workflow_context(
            draft,
            scientific_plan_sha256=(
                scientific_v2.plan_sha256 if scientific_v2 else ""
            ),
        )
        materialized_inputs, _completed = self._observed_workflow_state(
            draft,
            scientific_plan_sha256=(
                scientific_v2.plan_sha256 if scientific_v2 else ""
            ),
        )
        finding_nodes = {
            item["node_id"] for item in command_result.get("findings", ())
        } | {
            node.node_id
            for node in draft.nodes
            if _remaining_node_unresolved_fields(node, materialized_inputs)
        }
        actionable = tuple(
            node_id
            for node_id in context.ready_node_ids
            if node_id not in finding_nodes
        )
        unresolved = tuple(
            node.node_id
            for node in draft.nodes
            if node.node_id
            in (
                set(context.waiting_node_ids)
                | set(context.blocked_node_ids)
                | finding_nodes
            )
        )
        plan = resolved.scientific_toolchain_plan
        if plan is None:
            return {
                "scientific_workflow_plan": scientific_v2,
                "workflow_context": context,
            }
        analysis_receipts = self._scientific_toolchain_analysis_receipts(
            plan,
            task_spec_sha256=draft.task_spec_id,
        )
        return {
            "scientific_toolchain_plan": plan,
            "workflow_context": context,
            "workflow_frontier": project_scientific_toolchain_frontier(
                plan,
                actionable_calculation_node_ids=actionable,
                unresolved_calculation_node_ids=unresolved,
                completed_calculation_node_ids=context.completed_node_ids,
                completed_analysis_node_ids=analysis_receipts,
            ),
        }

    def _scientific_toolchain_analysis_receipts(
        self,
        plan: ScientificToolchainPlanV1,
        *,
        task_spec_sha256: str,
    ) -> dict[str, tuple[str, ...]]:
        """Relate planned analysis nodes to already observed typed evidence.

        Analysis tools intentionally do not accept a plan-node handle.  That
        keeps the normal extraction, thermochemistry, dimensional-expression,
        claims, and decision surface useful for creative DAGs and combined
        expressions.  The host can still recover the relation from semantics
        already present on both sides: registered result identity, selectors,
        typed output IDs, and receipt-level input dependencies.

        A scientific-validation node is complete only when the host has
        evaluated the exact predicates sealed into that plan node and recorded
        a typed validation receipt. A failed predicate is still an evaluated
        scientific determination; completion does not imply a positive verdict.
        """

        nodes = {node.node_id: node for node in plan.analysis_nodes}
        matched: dict[str, tuple[str, ...]] = {}

        def _dependency_receipts(node: AnalysisNodeIntentV1) -> dict[str, set[str]]:
            return {
                dependency: set(matched.get(dependency, ()))
                for dependency in node.dependencies
            }

        def _decision_evidence(decision: ScientificDecisionRecordV1) -> set[str]:
            evidence: set[str] = set()
            for reference in decision.evidence_refs:
                parsed = _postprocessing_evidence_reference(reference)
                if parsed is not None:
                    evidence.add(parsed[1])
            return evidence

        def _unit_dimension(unit: str) -> tuple[int, ...] | None:
            try:
                _value, _canonical_unit, dimension = normalize_numeric_value(
                    0.0, unit
                )
            except (ContractError, ValueError):
                return None
            return tuple(dimension)

        task_decisions = tuple(
            decision
            for decision in self.scientific_decisions.values()
            if decision.task_spec_sha256 == task_spec_sha256
        )
        task_claims = tuple(
            record
            for record in self.analysis_claim_records.values()
            if record.task_spec_sha256 == task_spec_sha256
        )

        for node_id in plan.node_order:
            node = nodes.get(node_id)
            if node is None or node.support_state != "planned":
                continue
            dependencies = _dependency_receipts(node)
            if any(not receipts for receipts in dependencies.values()):
                continue

            if node.analysis_kind == "result_extraction":
                registered = tuple(
                    item
                    for item in node.inputs
                    if isinstance(item, RegisteredResultInputIntentV1)
                )
                if len(registered) != 1:
                    continue
                # Typed extraction evidence is the triple (registered result,
                # selector, value).  ``quantity_id`` is the model's own label
                # for that value, and a plan node and the extraction call that
                # satisfies it are authored independently, so the same
                # scientist may reasonably name the same selector differently
                # in each.  Match on the selector set, exactly as the
                # policy-driven completion gate already does; requiring the
                # two labels to coincide would gate completion on naming
                # luck rather than on evidence.
                selectors = frozenset(
                    selector.selector for selector in node.selectors
                )
                exact_candidates = tuple(
                    receipt.receipt_sha256
                    for receipt in self.quantity_extractions.values()
                    if receipt.status == "extracted"
                    and receipt.artifact_id == registered[0].artifact_id
                    and selectors.issubset(
                        self.quantity_extraction_selectors.get(
                            receipt.receipt_sha256, ()
                        )
                    )
                )
                candidates = exact_candidates
                if not candidates and task_claims and task_decisions:
                    # A parser may expose only the supported subset of a
                    # broader analysis intent.  Once the final task-bound
                    # claims and decision explicitly cite that typed subset,
                    # the extraction stage was performed; the senior
                    # scientist, not this matcher, judges whether the stated
                    # limitation is acceptable.  Unknown or mismatched
                    # selectors still cannot enter through this fallback.
                    subset_candidates = {
                        receipt.receipt_sha256
                        for receipt in self.quantity_extractions.values()
                        if receipt.status == "extracted"
                        and receipt.artifact_id == registered[0].artifact_id
                        and (
                            observed := self.quantity_extraction_selectors.get(
                                receipt.receipt_sha256, ()
                            )
                        )
                        and selectors.issuperset(observed)
                    }
                    claim_digests = {
                        record.receipt_sha256 for record in task_claims
                    }
                    cited = {
                        digest
                        for decision in task_decisions
                        if _decision_evidence(decision).intersection(
                            claim_digests
                        )
                        for digest in _decision_evidence(decision)
                    }
                    candidates = tuple(sorted(subset_candidates & cited))
                if candidates:
                    matched[node_id] = tuple(sorted(set(candidates)))
                continue

            if node.analysis_kind == "thermochemistry":
                registered = tuple(
                    item
                    for item in node.inputs
                    if isinstance(item, RegisteredResultInputIntentV1)
                )
                if len(registered) == 1:
                    source_artifact_ids = {registered[0].artifact_id}
                else:
                    source_artifact_ids = {
                        receipt.artifact_id
                        for dependency_receipts in dependencies.values()
                        for digest in dependency_receipts
                        if (
                            receipt := self.quantity_extractions.get(digest)
                        )
                        is not None
                    }
                required_kinds = {
                    canonical_thermochemistry_quantity(
                        output.quantity_kind
                    ): output.unit
                    for output in node.outputs
                }
                candidates: list[str] = []
                for receipt in self.thermochemistry_receipts.values():
                    if (
                        receipt.status != "derived"
                        or receipt.artifact_id not in source_artifact_ids
                        or not math.isclose(
                            receipt.temperature_k,
                            float(node.temperature_k),
                            rel_tol=0.0,
                            abs_tol=1.0e-12,
                        )
                        or not math.isclose(
                            receipt.pressure_atm,
                            float(node.pressure_atm),
                            rel_tol=0.0,
                            abs_tol=1.0e-12,
                        )
                        or receipt.concentration_mol_l
                        != node.concentration_mol_l
                        or receipt.entropy_method != node.entropy_method
                        or receipt.entropy_cutoff_cm1
                        != node.entropy_cutoff_cm1
                        or receipt.enthalpy_cutoff_cm1
                        != node.enthalpy_cutoff_cm1
                        or receipt.alpha != node.alpha
                        or receipt.use_weighted_mass
                        is not node.use_weighted_mass
                        or not math.isclose(
                            receipt.frequency_scale_factor,
                            node.frequency_scale_factor,
                            rel_tol=0.0,
                            abs_tol=1.0e-12,
                        )
                    ):
                        continue
                    quantities = {
                        quantity.quantity_id: quantity
                        for quantity in receipt.quantities
                    }
                    if not set(required_kinds).issubset(quantities):
                        continue
                    compatible = True
                    for quantity_kind, unit in required_kinds.items():
                        dimension = _unit_dimension(unit)
                        if dimension is None:
                            compatible = False
                            break
                        if tuple(quantities[quantity_kind].dimension) != dimension:
                            compatible = False
                            break
                    if compatible:
                        candidates.append(receipt.receipt_sha256)
                if candidates:
                    matched[node_id] = tuple(sorted(set(candidates)))
                continue

            if node.analysis_kind == "quantity_expression":
                required_outputs = set(node.expression_output_node_ids)
                planned_output_operations = {
                    str(expression.get("operation") or "")
                    for expression in node.expression_nodes
                    if str(expression.get("node_id") or "")
                    in required_outputs
                }
                planned_dimensions: set[tuple[int, ...]] = set()
                for output in node.outputs:
                    dimension = _unit_dimension(output.unit)
                    if dimension is not None:
                        planned_dimensions.add(dimension)
                candidates: list[str] = []
                for receipt in self.quantity_expression_receipts.values():
                    if receipt.status != "derived":
                        continue
                    observed_outputs = {
                        quantity.quantity_id: quantity
                        for quantity in receipt.outputs
                    }
                    selected_outputs = set(required_outputs)
                    if not required_outputs.issubset(observed_outputs):
                        # A direct typed expression may be a scientifically
                        # stronger expansion of a planned abstract output (for
                        # example, several N-H...O distances instead of one
                        # placeholder contact).  Preserve that flexibility only
                        # when the host-observed operation and dimension still
                        # match the planned typed result.
                        selected_outputs = {
                            dependency.output_id
                            for dependency in receipt.output_dependencies
                            if planned_output_operations.intersection(
                                dependency.convention_operations
                            )
                            and dependency.output_id in observed_outputs
                            and (
                                not planned_dimensions
                                or tuple(
                                    observed_outputs[
                                        dependency.output_id
                                    ].dimension
                                )
                                in planned_dimensions
                            )
                        }
                        observed_operations = {
                            operation
                            for dependency in receipt.output_dependencies
                            if dependency.output_id in selected_outputs
                            for operation in dependency.convention_operations
                        }
                        if (
                            not selected_outputs
                            or not planned_output_operations.issubset(
                                observed_operations
                            )
                        ):
                            continue
                    sources = {
                        source
                        for dependency in receipt.output_dependencies
                        if dependency.output_id in selected_outputs
                        for source in dependency.source_receipt_sha256s
                    }
                    if all(
                        bool(sources.intersection(dependency_receipts))
                        for dependency_receipts in dependencies.values()
                    ):
                        candidates.append(receipt.receipt_sha256)
                if candidates:
                    matched[node_id] = tuple(sorted(set(candidates)))
                continue

            if node.analysis_kind == "scientific_validation":
                candidates = [
                    receipt.receipt_sha256
                    for receipt in self.scientific_validation_receipts.values()
                    if receipt.status == "evaluated"
                    and receipt.workflow_id == plan.workflow_id
                    and receipt.plan_sha256 == plan.plan_sha256
                    and receipt.node_id == node.node_id
                    and all(
                        bool(
                            set(receipt.source_receipt_sha256s).intersection(
                                dependency_receipts
                            )
                        )
                        for dependency_receipts in dependencies.values()
                    )
                ]
                if candidates:
                    matched[node_id] = tuple(sorted(set(candidates)))
                continue

            if node.analysis_kind == "claim_rendering":
                candidates: list[tuple[str, str]] = []
                for claim_record in task_claims:
                    claim_sources = {
                        claim.source_receipt_sha256
                        for claim in claim_record.claims
                    }
                    validation_inputs = tuple(
                        item
                        for item in node.inputs
                        if isinstance(item, AnalysisInputIntentV1)
                        and (
                            producer := nodes.get(item.producer_node_id)
                        )
                        is not None
                        and producer.analysis_kind
                        == "scientific_validation"
                    )
                    if any(
                        not any(
                            claim.source_receipt_sha256
                            in dependencies.get(item.producer_node_id, set())
                            and claim.quantity_id == item.producer_output_id
                            for claim in claim_record.claims
                        )
                        for item in validation_inputs
                    ):
                        continue
                    for decision in task_decisions:
                        evidence = _decision_evidence(decision)
                        if claim_record.receipt_sha256 not in evidence:
                            continue
                        if all(
                            bool(
                                (
                                    claim_sources
                                    | evidence
                                    | {decision.record_sha256}
                                ).intersection(
                                    dependency_receipts
                                )
                            )
                            for dependency_receipts in dependencies.values()
                        ):
                            candidates.append(
                                (
                                    claim_record.receipt_sha256,
                                    decision.record_sha256,
                                )
                            )
                if candidates:
                    matched[node_id] = tuple(
                        sorted({digest for pair in candidates for digest in pair})
                    )
        return matched

    def _prepare_program_node(self, turn_id: str, values: dict) -> Any:
        """Resolve and safe-preview a planned node without model-carried hashes.

        The model chooses the scientific workflow and node.  The host performs
        the relational join across the already observed capability,
        environment, project, identity, and artifact records.  Ambiguity is a
        semantic next action, never an invitation to guess a receipt digest.
        """

        workflow_id = values["workflow_id"]
        node_id = values["node_id"]
        resolved = self._resolve_program_workflow(workflow_id)
        draft = resolved.draft
        nodes = tuple(node for node in draft.nodes if node.node_id == node_id)
        if len(nodes) != 1:
            raise ContractError("workflow has no unique calculation node ID")
        node = nodes[0]
        if not node.inputs:
            return {
                "status": "needs_clarification",
                "workflow_id": workflow_id,
                "node_id": node_id,
                "finding": "program node declares no molecular input",
                "next_action": "repair the workflow input binding",
            }

        scientific_v2 = resolved.scientific_plan
        if scientific_v2 is None:
            raise ContractError("workflow has no task-bound scientific plan")
        plan_sha256 = scientific_v2.plan_sha256
        materialized_inputs, _completed = self._observed_workflow_state(
            draft,
            scientific_plan_sha256=plan_sha256,
        )
        workflow_context = self._workflow_context(
            draft,
            scientific_plan_sha256=plan_sha256,
        )
        node_context = workflow_context.node(node_id)
        if node_context is not None and node_context.state == "completed":
            return {
                "status": "completed",
                "workflow_id": workflow_id,
                "node_id": node_id,
                "next_action": "use the validated result already recorded",
            }
        if node_context is not None and node_context.state == "blocked":
            return {
                "status": "blocked",
                "workflow_id": workflow_id,
                "node_id": node_id,
                "finding": node_context.reason,
                "next_action": "continue an independent ready branch or close the failed branch",
            }
        if node_context is not None and node_context.state == "waiting":
            return {
                "status": "waiting_for_artifact",
                "workflow_id": workflow_id,
                "node_id": node_id,
                "producer_inputs": tuple(
                    canonical_data(item)
                    for item in node_context.unsatisfied_inputs
                ),
                "waiting_on": node_context.waiting_on,
                "next_action": "execute and validate the named producer first",
            }

        remaining_unresolved = _remaining_node_unresolved_fields(
            node, materialized_inputs
        )
        if remaining_unresolved:
            return {
                "status": "needs_clarification",
                "workflow_id": workflow_id,
                "node_id": node_id,
                "unresolved_fields": remaining_unresolved,
                "next_action": "resolve the node's scientific settings",
            }
        resolved_inputs: dict[str, TrustedArtifactRefV1] = {}
        waiting_producers = []
        waiting_artifacts = []
        for item in node.inputs:
            if item.producer_node_id:
                artifact = materialized_inputs.get(
                    (
                        node_id,
                        item.binding_id,
                        item.producer_node_id,
                        item.producer_output_id,
                    )
                )
                if artifact is None:
                    waiting_producers.append(
                        {
                            "binding_id": item.binding_id,
                            "producer_node_id": item.producer_node_id,
                            "producer_output_id": item.producer_output_id,
                        }
                    )
                    continue
            else:
                artifact = self.artifacts.get(item.artifact_id)
                if artifact is None:
                    waiting_artifacts.append(
                        {
                            "binding_id": item.binding_id,
                            "artifact_id": item.artifact_id,
                        }
                    )
                    continue
            resolved_inputs[item.binding_id] = artifact
        if waiting_producers or waiting_artifacts:
            return {
                "status": "waiting_for_artifact",
                "workflow_id": workflow_id,
                "node_id": node_id,
                "producer_inputs": tuple(waiting_producers),
                "external_inputs": tuple(waiting_artifacts),
                "next_action": "materialize the declared workflow inputs",
            }
        if not resolved_inputs:
            return {
                "status": "needs_clarification",
                "workflow_id": workflow_id,
                "node_id": node_id,
                "finding": "program node declares no molecular input",
                "next_action": "repair the workflow input binding",
            }
        if len(resolved_inputs) == 1:
            input_artifact = next(iter(resolved_inputs.values()))
            job_artifact_options: tuple[
                tuple[str, TrustedArtifactRefV1], ...
            ] = ()
        else:
            input_artifact = resolved_inputs.get("filename")
            if input_artifact is None:
                return {
                    "status": "needs_clarification",
                    "workflow_id": workflow_id,
                    "node_id": node_id,
                    "finding": (
                        "a multi-file ChemSmart job needs binding_id "
                        "'filename' for its primary geometry; bind each "
                        "additional artifact by its live job-option name"
                    ),
                    "next_action": "repair the workflow input binding roles",
                }
            job_artifact_options = tuple(
                sorted(
                    (
                        (binding_id, artifact)
                        for binding_id, artifact in resolved_inputs.items()
                        if binding_id != "filename"
                    ),
                    key=lambda item: item[0],
                )
            )

        identity_candidates = tuple(
            identity
            for identity in self.scientific_identities.values()
            if identity.task_spec_sha256 == draft.task_spec_id
            and identity.geometry_artifact_id == input_artifact.artifact_id
            and identity.geometry_artifact_sha256 == input_artifact.sha256
        )
        identities = (
            tuple(
                identity
                for identity in identity_candidates
                if (identity.charge, identity.multiplicity)
                == (node.charge, node.multiplicity)
            )
            if node.charge is not None
            else identity_candidates
        )
        if len(identities) != 1:
            return {
                "status": "needs_clarification",
                "workflow_id": workflow_id,
                "node_id": node_id,
                "finding": (
                    "input has no task-bound identity for the node's explicit "
                    "charge and multiplicity"
                    if node.charge is not None and not identities
                    else "input has no unique task-bound electronic state"
                ),
                "candidate_states": tuple(
                    {
                        "charge": item.charge,
                        "multiplicity": item.multiplicity,
                    }
                    for item in identity_candidates
                ),
                "next_action": "bind one charge and multiplicity to this input",
            }
        identity = identities[0]
        for parameter_name, artifact in job_artifact_options:
            if artifact.kind not in {"geometry_xyz", "xyz"}:
                continue
            secondary_candidates = tuple(
                candidate
                for candidate in self.scientific_identities.values()
                if candidate.task_spec_sha256 == draft.task_spec_id
                and candidate.geometry_artifact_id == artifact.artifact_id
                and candidate.geometry_artifact_sha256 == artifact.sha256
            )
            secondary_identities = tuple(
                candidate
                for candidate in secondary_candidates
                if (candidate.charge, candidate.multiplicity)
                == (identity.charge, identity.multiplicity)
            )
            if len(secondary_identities) != 1:
                return {
                    "status": "needs_clarification",
                    "workflow_id": workflow_id,
                    "node_id": node_id,
                    "finding": (
                        f"{parameter_name} has no unique task-bound "
                        "electronic state"
                    ),
                    "next_action": (
                        "bind charge and multiplicity to every geometry input"
                    ),
                }
            secondary = secondary_identities[0]
            if (secondary.charge, secondary.multiplicity) != (
                identity.charge,
                identity.multiplicity,
            ):
                return {
                    "status": "needs_clarification",
                    "workflow_id": workflow_id,
                    "node_id": node_id,
                    "finding": (
                        f"{parameter_name} is on a different charge or "
                        "multiplicity surface"
                    ),
                    "next_action": "repair the scientific endpoint states",
                }

        capabilities = tuple(
            receipt
            for receipt in self.capabilities.values()
            if receipt.query.program == node.program
            and receipt.query.jobtype == node.jobtype
            and str(receipt.status.value) in {"supported", "preview_only"}
        )
        if len(capabilities) != 1:
            return {
                "status": "needs_capability_selection",
                "workflow_id": workflow_id,
                "node_id": node_id,
                "program": node.program,
                "jobtype": node.jobtype,
                "candidate_engines": tuple(
                    sorted({item.query.engine for item in capabilities})
                ),
                "next_action": "inspect one explicit program engine",
            }
        capability = capabilities[0]
        engine_bindings = tuple(
            binding
            for binding in self.engine_bindings.values()
            if binding.program == node.program
            and binding.capability_receipt_sha256 == capability.receipt_sha256
            and binding.state != "blocked"
        )
        if len(engine_bindings) != 1:
            return {
                "status": "needs_engine_selection",
                "workflow_id": workflow_id,
                "node_id": node_id,
                "candidate_engines": tuple(
                    sorted({item.engine for item in engine_bindings})
                ),
                "next_action": "resolve one observed program engine",
            }
        engine_binding = engine_bindings[0]
        program_binding = self.program_bindings.get(
            engine_binding.program_binding_sha256
        )
        if program_binding is None:
            raise ContractError("engine binding lacks its program binding")

        project = self.artifacts.get(node.project_role)
        if project is None:
            return {
                "status": "needs_project",
                "workflow_id": workflow_id,
                "node_id": node_id,
                "project_role": node.project_role,
                "next_action": "render, promote, and validate this project role",
            }
        validations = tuple(
            receipt
            for receipt in self.project_validations.values()
            if receipt.project_artifact_id == project.artifact_id
            and receipt.project_sha256 == project.sha256
            and receipt.capability_receipt_sha256 == capability.receipt_sha256
            and receipt.program == node.program
            and receipt.jobtype == node.jobtype
            and receipt.status == "valid"
        )
        if len(validations) != 1:
            return {
                "status": "needs_project_validation",
                "workflow_id": workflow_id,
                "node_id": node_id,
                "project_artifact_id": project.artifact_id,
                "next_action": "validate the project for this program stage",
            }
        validation = validations[0]
        execution_target = (
            self.execution_resources.execution_target
            if self.surface.profile == "command_compiled_approved_execution"
            and self.execution_resources is not None
            else "run"
        )
        proposal = CommandProposalV1(
            node_id=node.node_id,
            execution_target=execution_target,
            program=node.program,
            jobtype=node.jobtype,
            project_artifact_id=project.artifact_id,
            input_artifact_id=input_artifact.artifact_id,
            scientific_identity_sha256=identity.binding_sha256,
            charge=identity.charge,
            multiplicity=identity.multiplicity,
        )
        invocation = compile_command(
            proposal,
            capability=capability,
            binding=engine_binding,
            project=project,
            project_validation=validation,
            input_artifact=input_artifact,
            scientific_identity=identity,
            job_artifact_options=dict(job_artifact_options),
            live_schema=self.live_schema,
            server=(
                self.execution_server
                if self.surface.profile
                == "command_compiled_approved_execution"
                else self.preview_server
            ),
        )
        context = _CommandContext(
            proposal=proposal,
            capability=capability,
            program_binding=program_binding,
            engine_binding=engine_binding,
            project_artifact=project,
            project_validation=validation,
            input_artifact=input_artifact,
            scientific_identity=identity,
            job_artifact_options=job_artifact_options,
        )
        compiled = self._record_compiled_command(turn_id, invocation, context)
        if plan_sha256:
            self._invocation_workflow_plan_sha256s[
                invocation.invocation_sha256
            ] = plan_sha256
        if compiled["inspection"].status != "valid":
            return {
                "status": "blocked",
                "workflow_id": workflow_id,
                "node_id": node_id,
                "command": compiled,
                "next_action": "repair the command inspection finding",
            }
        preview = self._preview_command(
            turn_id, {"invocation_sha256": invocation.invocation_sha256}
        )
        preview_status = preview["safe_preview"].status
        preflight = None
        if preview_status == "previewed" and not preview["critical_findings"]:
            # Preparation already holds every host-owned input needed by the
            # deterministic preflight.  Completing it here makes scientific
            # readiness independent of whether a model memorizes a second
            # bookkeeping tool name; it still grants no execution authority.
            preflight = self._preflight_program_node(
                turn_id,
                {
                    "node_id": node.node_id,
                    "capability_receipt_sha256": capability.receipt_sha256,
                    "program_binding_sha256": program_binding.binding_sha256,
                    "engine_binding_sha256": engine_binding.binding_sha256,
                    "geometry_artifact_sha256": input_artifact.sha256,
                    "scientific_identity_sha256": identity.binding_sha256,
                    "charge": identity.charge,
                    "multiplicity": identity.multiplicity,
                    "project_validation_receipt_sha256": validation.receipt_sha256,
                    "invocation_sha256": invocation.invocation_sha256,
                    "command_inspection_receipt_sha256": compiled[
                        "inspection"
                    ].receipt_sha256,
                    "safe_preview_receipt_sha256": preview[
                        "safe_preview"
                    ].receipt_sha256,
                },
            )
        return {
            "status": (
                "previewed"
                if preview_status == "previewed"
                else "preview_failed"
            ),
            "workflow_id": workflow_id,
            "node_id": node_id,
            "program": node.program,
            "jobtype": node.jobtype,
            "project_artifact_id": project.artifact_id,
            "input_artifact_id": input_artifact.artifact_id,
            "job_artifact_inputs": tuple(
                {
                    "parameter_name": parameter_name,
                    "artifact_id": artifact.artifact_id,
                }
                for parameter_name, artifact in job_artifact_options
            ),
            "electronic_state": {
                "charge": identity.charge,
                "multiplicity": identity.multiplicity,
            },
            "command": compiled,
            "preview": preview,
            "preflight": preflight,
            "next_action": (
                "inspect the workflow frontier"
                if preview_status == "previewed"
                else "inspect the generated-input validation findings"
            ),
        }

    def _observed_workflow_state(
        self,
        draft: CommandWorkflowDraftV1,
        *,
        scientific_plan_sha256: str,
    ) -> tuple[
        dict[tuple[str, str, str, str], TrustedArtifactRefV1],
        tuple[str, ...],
    ]:
        """Join this exact approved plan to its executed nodes and handoffs."""

        frozen = self.frozen_workflow_approval
        if (
            frozen is None
            or not scientific_plan_sha256
            or frozen.workflow_id != draft.workflow_id
            or frozen.plan_sha256 != scientific_plan_sha256
        ):
            return {}, ()

        durable_run = self._durable_workflow_run_state(
            workflow_id=draft.workflow_id,
            plan_sha256=scientific_plan_sha256,
        )
        durable_states = {
            node.node_id: node.state
            for node in (durable_run.nodes if durable_run is not None else ())
        }
        completed = tuple(
            sorted(
                node.node_id
                for node in draft.nodes
                if durable_states.get(node.node_id) == "validated"
                or (
                    (receipt := self.execution_receipts.get(node.node_id))
                    is not None
                    and receipt.validated
                )
            )
        )
        resolved: dict[tuple[str, str, str, str], TrustedArtifactRefV1] = {}
        for node in draft.nodes:
            for item in node.inputs:
                if item.artifact_class == "geometry_xyz":
                    handoff = self.handoffs.get(node.node_id)
                elif item.artifact_class == "orca_hessian":
                    handoff = self.hessian_handoffs.get(node.node_id)
                else:
                    handoff = None
                if handoff is None or handoff.status != "validated_handoff":
                    continue
                artifact = self.artifacts.get(handoff.selected_artifact_id)
                if (
                    artifact is None
                    or artifact.sha256 != handoff.selected_artifact_sha256
                ):
                    continue
                if (
                    item.producer_node_id != handoff.producer_node_id
                    or node.node_id != handoff.consumer_node_id
                ):
                    continue
                matching_rules = tuple(
                    rule
                    for rule in frozen.producer_edge_rules
                    if rule.source_node_id == item.producer_node_id
                    and rule.target_node_id == node.node_id
                    and rule.consumer_input_id == item.binding_id
                    and rule.producer_output_id == item.producer_output_id
                    and rule.artifact_class == item.artifact_class
                )
                if len(matching_rules) != 1:
                    continue
                resolved[
                    (
                        node.node_id,
                        item.binding_id,
                        item.producer_node_id,
                        item.producer_output_id,
                    )
                ] = artifact
        return resolved, completed

    def _workflow_context(
        self,
        draft: CommandWorkflowDraftV1,
        *,
        scientific_plan_sha256: str = "",
    ) -> Any:
        """Derive the dependency context the model would otherwise reconstruct.

        Host-derived and read-only: the model is told which nodes are runnable
        and what each waiting node is waiting for, and can never assert it.
        """

        materialized_inputs, completed = self._observed_workflow_state(
            draft,
            scientific_plan_sha256=scientific_plan_sha256,
        )
        durable_run = self._durable_workflow_run_state(
            workflow_id=draft.workflow_id,
            plan_sha256=scientific_plan_sha256,
        )
        blocked_reasons = {}
        if durable_run is not None:
            for node in durable_run.nodes:
                if node.state in {"failed", "blocked", "ambiguous"}:
                    blocked_reasons[node.node_id] = (
                        node.state
                        + " in the durable workflow run"
                        + (
                            ": " + ", ".join(node.failure_rule_ids)
                            if node.failure_rule_ids
                            else ""
                        )
                    )
                elif node.state in {"running", "engine_complete"}:
                    blocked_reasons[node.node_id] = (
                        "execution already started; reconcile its durable receipt"
                    )
        return project_workflow_context(
            workflow_id=draft.workflow_id,
            nodes=draft.nodes,
            materialized_artifact_ids=self.artifacts,
            materialized_producer_inputs=materialized_inputs,
            completed_node_ids=completed,
            blocked_reasons=blocked_reasons,
        )

    def _durable_workflow_run_state(
        self, *, workflow_id: str, plan_sha256: str
    ) -> WorkflowRunStateV1 | None:
        """Return the replayed run only when it owns this exact plan."""

        frozen = self.frozen_workflow_approval
        if (
            frozen is None
            or not plan_sha256
            or frozen.workflow_id != workflow_id
            or frozen.plan_sha256 != plan_sha256
        ):
            return None
        frontier = self.event_store.workflow_frontier(
            workflow_id=workflow_id,
            run_id="run." + frozen.approval_id,
        )
        run = frontier.run_state
        if run is not None and run.plan_sha256 != plan_sha256:
            raise ContractError("durable workflow run belongs to another plan")
        return run

    def _scientific_plan_from_draft(
        self,
        draft: CommandWorkflowDraftV1,
        *,
        findings: list[dict[str, str]],
        node_annotations: Mapping[str, Mapping[str, Any]] | None = None,
    ) -> ScientificWorkflowPlanV2 | None:
        """Project a V1 model draft into the host-owned scientific DAG."""

        if draft.task_spec_id not in self.task_spec_sha256s:
            findings.append(
                {
                    "node_id": draft.workflow_id,
                    "rule_id": "workflow.task_spec.unbound",
                }
            )
            return None
        external_artifact_ids = {
            item.artifact_id
            for node in draft.nodes
            for item in node.inputs
            if not item.producer_node_id and item.artifact_id
        }
        identities = tuple(
            sorted(
                {
                    identity.binding_sha256
                    for identity in self.scientific_identities.values()
                    if identity.task_spec_sha256 == draft.task_spec_id
                    and identity.geometry_artifact_id in external_artifact_ids
                }
            )
        )
        if not identities:
            findings.append(
                {
                    "node_id": draft.workflow_id,
                    "rule_id": "workflow.scientific_identity.unbound",
                }
            )
            return None
        scientific_identity_sha256 = (
            identities[0]
            if len(identities) == 1
            else canonical_sha256({"scientific_identity_sha256s": identities})
        )
        annotations = dict(node_annotations or {})
        scientific_nodes = []
        for node in draft.nodes:
            matching_capabilities = tuple(
                receipt
                for receipt in self.capabilities.values()
                if receipt.query.program == node.program
                and receipt.query.jobtype == node.jobtype
            )
            engines = tuple(
                sorted(
                    {receipt.query.engine for receipt in matching_capabilities}
                )
            )
            engine = engines[0] if len(engines) == 1 else "unresolved"
            unresolved = set(node.unresolved_fields)
            if node.node_kind == "aggregate":
                # The host is the engine and it is always present, so this is
                # resolved by construction rather than left for discovery.
                engine = AGGREGATE_NODE_PROGRAM
                unresolved.discard("engine")
            elif engine == "unresolved":
                unresolved.add("engine")
            requested_programs = {
                binding.requested_program
                for binding in self.program_bindings.values()
                if binding.selected_program == node.program
                and any(
                    receipt.receipt_sha256 == binding.capability_receipt_sha256
                    for receipt in matching_capabilities
                )
            }
            requested_program = (
                next(iter(requested_programs))
                if len(requested_programs) == 1
                else node.program
            )
            annotation = dict(annotations.get(node.node_id) or {})
            declared_support = str(
                annotation.get("support_state") or "planned"
            )
            support_state = (
                "blocked_unsupported"
                if declared_support == "blocked_unsupported"
                else "unresolved_future"
                if unresolved
                else "resolvable"
            )
            blocked_reason = (
                str(annotation.get("blocked_reason") or "")
                if support_state == "blocked_unsupported"
                else ""
            )
            scientific_nodes.append(
                ScientificWorkflowNodeV2(
                    node_id=node.node_id,
                    stage=node.jobtype,
                    requested_program=requested_program,
                    program=node.program,
                    engine=engine,
                    project_role=node.project_role,
                    unresolved_fields=tuple(sorted(unresolved)),
                    produces_observables=tuple(
                        sorted(
                            set(annotation.get("produces_observables") or ())
                        )
                    ),
                    support_state=support_state,
                    blocked_reason=blocked_reason,
                    charge=node.charge,
                    multiplicity=node.multiplicity,
                )
            )
        edges = []
        for node in draft.nodes:
            data_sources = set()
            for item in node.inputs:
                if not item.producer_node_id:
                    continue
                data_sources.add(item.producer_node_id)
                edges.append(
                    ScientificWorkflowEdgeV2(
                        edge_id=(
                            "data."
                            + item.producer_node_id
                            + "."
                            + node.node_id
                            + "."
                            + item.binding_id
                        ),
                        source_node_id=item.producer_node_id,
                        target_node_id=node.node_id,
                        edge_kind="data",
                        artifact_class=item.artifact_class,
                        producer_output_id=item.producer_output_id,
                        consumer_input_id=item.binding_id,
                    )
                )
            for dependency in node.dependencies:
                if dependency in data_sources:
                    continue
                edges.append(
                    ScientificWorkflowEdgeV2(
                        edge_id=("control." + dependency + "." + node.node_id),
                        source_node_id=dependency,
                        target_node_id=node.node_id,
                        edge_kind="control",
                    )
                )
        plan = build_scientific_workflow_plan(
            workflow_id=draft.workflow_id,
            task_spec_sha256=draft.task_spec_id,
            scientific_identity_sha256=scientific_identity_sha256,
            nodes=tuple(scientific_nodes),
            edges=tuple(sorted(edges, key=lambda edge: edge.edge_id)),
            required_observables=tuple(
                sorted(
                    {
                        observable
                        for node in scientific_nodes
                        for observable in node.produces_observables
                    }
                )
            ),
        )
        if (
            self.frozen_workflow_approval is not None
            and getattr(
                self.frozen_workflow_approval,
                "workflow_id",
                plan.workflow_id,
            )
            == plan.workflow_id
            and self.frozen_workflow_approval.plan_sha256 != plan.plan_sha256
        ):
            raise ContractError(
                "planned workflow differs from frozen execution approval"
            )
        self._refuse_observable_regression(plan)
        self.scientific_plans[plan.workflow_id] = plan
        return plan

    def _refuse_observable_regression(self, plan) -> None:
        """Refuse a replan that drops a stage the previous plan carried.

        Repair is scored on whether findings clear, and deleting the node that
        carries the findings is the cheapest way to clear them -- which silently
        discards the stage the task asked for.  A stage that cannot be
        materialized has to stay in the plan as ``blocked_unsupported`` with a
        reason, so an honest plan and a complete plan are the same plan.
        """

        previous = self.scientific_plans.get(plan.workflow_id)
        if previous is None:
            return
        current_ids = {node.node_id for node in plan.nodes}
        dropped = sorted(
            node.node_id
            for node in previous.nodes
            if node.node_id not in current_ids
        )
        if dropped:
            raise ContractError(
                f"replanning removed workflow stage(s) {dropped} that the "
                "previous plan carried; a stage that cannot be materialized "
                "must be kept with support_state='blocked_unsupported' and a "
                "blocked_reason instead of being deleted"
            )

    def _resolve_project_validation(
        self, *, project, capability, program: str, jobtype: str
    ):
        """Find the validation receipt the caller already earned.

        ``project_validation_receipt_sha256`` is the one optional argument on
        ``synthesize_command``, and it shipped undescribed.  Omitting it made
        the compiler refuse with

            project 'X' is bound but has no validation receipt; call
            validate_project_yaml on it first

        which names the wrong repair: the caller had called it, and the receipt
        was sitting in the host registry.  The bookkeeping is mechanical once
        the project, capability, program and stage are fixed -- ``prepare_
        program_node`` already resolves it with exactly this predicate -- so
        the host does it rather than making the model thread a digest it
        cannot get wrong in any interesting way.

        Returns ``None`` when nothing matches, which leaves the existing
        message to say, correctly this time, that validation is missing.
        """

        if project is None:
            return None
        matches = tuple(
            receipt
            for receipt in self.project_validations.values()
            if receipt.project_artifact_id == project.artifact_id
            and receipt.project_sha256 == project.sha256
            and receipt.capability_receipt_sha256 == capability.receipt_sha256
            and receipt.program == program
            and receipt.jobtype == jobtype
            and receipt.status == "valid"
        )
        if len(matches) > 1:
            raise ContractError(
                f"project {project.artifact_id!r} has "
                f"{len(matches)} valid {program} {jobtype} validation "
                "receipts; pass project_validation_receipt_sha256 to say "
                "which one this command was compiled against"
            )
        return matches[0] if matches else None

    def _node_is_previewed(
        self, node_id: str, *, plan_sha256: str = ""
    ) -> bool:
        """Whether a node holds the green preview an approval will demand."""

        preflight = self._preflight_by_node.get(node_id)
        safe_preview = (
            self.safe_previews.get(preflight.safe_preview_receipt_sha256)
            if preflight is not None
            else None
        )
        if preflight is None or safe_preview is None:
            return False
        if (
            plan_sha256
            and self._invocation_workflow_plan_sha256s.get(
                safe_preview.invocation_sha256
            )
            != plan_sha256
        ):
            return False
        try:
            invocation, _context = self._latest_invocation_for_node(
                node_id, plan_sha256=plan_sha256
            )
        except ContractError:
            return False
        return (
            safe_preview.invocation_sha256 == invocation.invocation_sha256
            and preflight.plan_state == "previewed"
            and not preflight.critical_finding_sha256s
        )

    def _bounded_deferred_target_ids(self, plan: Any) -> set[str]:
        """Return causal future nodes that bounded execution can defer.

        A consumer of an optimized geometry cannot be compiled or previewed
        before its producer runs.  The bounded execution contract admits
        that exact dependency without weakening the preview requirement for
        any initially runnable node.  Keep the predicate here aligned with
        the producer-edge subset accepted by ``_admit_bounded_workflow`` so
        planning feedback does not tell the model to delete required science.
        """

        envelope = self.bounded_execution_envelope
        if envelope is None:
            return set()
        nodes = {node.node_id: node for node in getattr(plan, "nodes", ())}
        data_edges = tuple(
            edge
            for edge in getattr(plan, "edges", ())
            if edge.edge_kind == "data"
        )
        incoming_counts: dict[str, int] = {}
        for edge in data_edges:
            incoming_counts[edge.target_node_id] = (
                incoming_counts.get(edge.target_node_id, 0) + 1
            )
        deferred = set()
        for edge in data_edges:
            producer = nodes.get(edge.source_node_id)
            target = nodes.get(edge.target_node_id)
            if (
                producer is None
                or target is None
                or incoming_counts[edge.target_node_id] != 1
                or not is_validated_optimized_geometry_edge(plan, edge)
                or producer.program not in {"gaussian", "orca", "pyscf", "xtb"}
                or target.support_state
                not in {"resolvable", "unresolved_future"}
                or not envelope.allows(target.program, target.engine)
            ):
                continue
            deferred.add(target.node_id)
        return deferred

    def _approval_readiness(self, plan: Any) -> dict[str, Any]:
        """Say which nodes still stand between this plan and execution.

        Exact approvals require every materialized node to hold a green
        preview.  Bounded local execution additionally permits an exact
        producer-data target to remain deferred until its optimized geometry
        exists.  That causal future node is not a preview blocker and must not
        be presented as a stage the model should delete.
        """

        nodes = []
        blocking = []
        deferred_ids = self._bounded_deferred_target_ids(plan)
        non_executable_ids = (
            self._release_non_executable_node_ids(plan)
            if self.bounded_execution_envelope is not None
            else frozenset()
        )
        planned_ids = {
            node.node_id for node in getattr(plan, "nodes", ())
        }
        executable_ids = planned_ids - non_executable_ids
        for node in getattr(plan, "nodes", ()):
            node_id = node.node_id
            previewed = self._node_is_previewed(
                node_id, plan_sha256=plan.plan_sha256
            )
            non_executable = node_id in non_executable_ids
            deferred = (
                not previewed
                and not non_executable
                and node_id in deferred_ids
            )
            blocks_approval = (
                not previewed and not deferred and not non_executable
            )
            if blocks_approval:
                blocking.append(node_id)
            nodes.append(
                {
                    "node_id": node_id,
                    "program": getattr(node, "program", ""),
                    "previewed": previewed,
                    "deferred_admissible": deferred,
                    "non_executable": non_executable,
                    "approval_state": (
                        "non_executable"
                        if non_executable
                        else (
                            "previewed"
                            if previewed
                            else (
                                "deferred_admissible"
                                if deferred
                                else "preview_required"
                            )
                        )
                    ),
                    "blocks_approval": blocks_approval,
                }
            )
        return {
            "approvable": not blocking and bool(executable_ids),
            "authorization_mode": (
                "bounded_local"
                if getattr(self, "bounded_execution_envelope", None)
                is not None
                else "exact_preview"
            ),
            "blocking_node_ids": tuple(blocking),
            "deferred_node_ids": tuple(
                node["node_id"]
                for node in nodes
                if node["deferred_admissible"]
            ),
            "non_executable_node_ids": tuple(
                node["node_id"]
                for node in nodes
                if node["non_executable"]
            ),
            "workflow_blocked_reason": (
                "the workflow has no release-executable stage to review"
                if planned_ids and not executable_ids
                else ""
            ),
            "nodes": tuple(nodes),
            "rule": (
                (
                    "every initially runnable node needs a green preview "
                    "before execution. Under bounded local execution, "
                    "an exact producer-data target is deferred_admissible "
                    "until its producer materializes the optimized geometry. "
                    "A release-unsupported stage marked non_executable is "
                    "retained as scientific intent but is not approved or "
                    "launched and does not require a green preview. "
                    "At least one release-executable stage is required for "
                    "human execution review. "
                    "Repair a preview_required node using the findings "
                    "returned by preview_command; do not delete a "
                    "scientifically required causal stage merely because "
                    "its producer output does not exist yet."
                )
                if getattr(self, "bounded_execution_envelope", None)
                is not None
                else (
                    "every materialized node needs a green preview before "
                    "exact approval. Repair a preview_required node using "
                    "the findings returned by preview_command."
                )
            ),
        }

    def _invocation_identity(
        self, node_id: str, *, plan_sha256: str = ""
    ) -> str:
        """Path-independent identity of a node's latest compiled command.

        Computed from one place so the digest a materialization freezes and
        the digest execution presents cannot be assembled differently.
        """

        try:
            invocation, context = self._latest_invocation_for_node(
                node_id, plan_sha256=plan_sha256
            )
        except ContractError:
            return ""
        project_artifact = context.project_artifact
        identity = context.scientific_identity
        if project_artifact is None or identity is None:
            return ""
        return invocation_identity_sha256(
            program=invocation.command_path[1],
            engine=context.engine_binding.engine,
            jobtype=invocation.command_path[-1],
            project_sha256=project_artifact.sha256,
            input_sha256=context.input_artifact.sha256,
            scientific_identity_sha256=identity.binding_sha256,
            argv=invocation.argv,
            auxiliary_input_bindings=invocation.auxiliary_input_bindings,
        )

    def _environment_identity_for(self, receipt_sha256: str) -> str:
        """Identity of an observed environment, or "" when unresolvable."""

        observed = self.environments.get(receipt_sha256)
        return (
            "" if observed is None else environment_identity_sha256(observed)
        )

    def _environment_identity_is_approved(
        self, observed_sha256: str, approved_sha256s
    ) -> bool:
        """Accept the approved *machine* even when the receipt digest moved.

        An environment receipt's digest folds in its capability receipt, and a
        capability receipt changes with the active overlay, so the receipt a
        plan session records is never the one an execution session computes --
        even on the same interpreter, with identical versions. Pinning the
        digest therefore rejected the machine the approval named, and no
        reviewer could supply the right digest because it does not exist until
        execution is already authorised.

        Comparing environment identity keeps the property the check exists for:
        a different interpreter, a different dependency set, or a different
        accelerator still fails. Only the authorisation flavour is ignored.
        Falls back to refusing when either side's receipt body is unavailable.
        """

        approved = {digest for digest in approved_sha256s if digest}
        if not approved:
            return False
        observed = self.environments.get(observed_sha256)
        if observed is None:
            return False
        identity = environment_identity_sha256(observed)
        for digest in approved:
            candidate = self.environments.get(digest)
            if candidate is None:
                continue
            if environment_identity_sha256(candidate) == identity:
                return True
        approved_identities = getattr(
            self, "approved_environment_identities", ()
        )
        return identity in set(approved_identities or ())

    def _resolve_safe_preview(self, invocation_sha256: str):
        """Find the safe preview this invocation already produced.

        Keyed by invocation, so there is nothing for a caller to choose. When
        an invocation has been previewed more than once the newest receipt is
        the one the preflight should carry, because it reflects the current
        project bytes.
        """

        matches = [
            receipt
            for receipt in self.safe_previews.values()
            if receipt.invocation_sha256 == invocation_sha256
        ]
        return matches[-1] if matches else None

    def _synthesize_command(self, turn_id: str, values: dict) -> Any:
        capability = self._get(
            self.capabilities,
            values["capability_receipt_sha256"],
            "capability receipt",
        )
        binding = self._get(
            self.engine_bindings,
            values["engine_binding_sha256"],
            "engine binding",
        )
        input_artifact = self._artifact(values["input_artifact_id"])
        identity = self._get(
            self.scientific_identities,
            values["scientific_identity_sha256"],
            "scientific identity",
        )
        project = (
            self._artifact(values["project_artifact_id"])
            if values["project_artifact_id"]
            else None
        )
        validation_digest = values.get("project_validation_receipt_sha256", "")
        validation = (
            self._get(
                self.project_validations,
                validation_digest,
                "project validation receipt",
            )
            if validation_digest
            else self._resolve_project_validation(
                project=project,
                capability=capability,
                program=values["program"],
                jobtype=values["jobtype"],
            )
        )
        execution_target = (
            self.execution_resources.execution_target
            if self.surface.profile == "command_compiled_approved_execution"
            and self.execution_resources is not None
            else "run"
        )
        supplied_target = str(values.get("execution_target") or "").strip()
        if supplied_target and supplied_target != execution_target:
            raise ContractError(
                "execution target is host-owned and differs from the active profile"
            )
        proposal = CommandProposalV1(
            node_id=values["node_id"],
            execution_target=execution_target,
            program=values["program"],
            jobtype=values["jobtype"],
            project_artifact_id=values["project_artifact_id"],
            input_artifact_id=values["input_artifact_id"],
            scientific_identity_sha256=values["scientific_identity_sha256"],
            charge=values["charge"],
            multiplicity=values["multiplicity"],
        )
        invocation = compile_command(
            proposal,
            capability=capability,
            binding=binding,
            project=project,
            project_validation=validation,
            input_artifact=input_artifact,
            scientific_identity=identity,
            live_schema=self.live_schema,
            server=(
                self.execution_server
                if self.surface.profile
                == "command_compiled_approved_execution"
                else self.preview_server
            ),
        )
        program_binding = self._get(
            self.program_bindings,
            binding.program_binding_sha256,
            "program binding",
        )
        context = _CommandContext(
            proposal=proposal,
            capability=capability,
            program_binding=program_binding,
            engine_binding=binding,
            project_artifact=project,
            project_validation=validation,
            input_artifact=input_artifact,
            scientific_identity=identity,
        )
        return self._record_compiled_command(turn_id, invocation, context)

    def _repair_command(self, turn_id: str, values: dict) -> Any:
        parent = self._get(
            self.invocations,
            values["invocation_sha256"],
            "canonical invocation",
        )
        context = self._get(
            self._command_contexts,
            parent.invocation_sha256,
            "command context",
        )
        counterexample = self._get(
            self.counterexamples,
            values["counterexample_id"],
            "counterexample",
        )
        if counterexample.invocation_sha256 != parent.invocation_sha256:
            raise ContractError("counterexample targets another invocation")
        if (
            counterexample.task_spec_sha256
            != context.scientific_identity.task_spec_sha256
        ):
            raise ContractError("counterexample targets another task spec")
        invocation = compile_command(
            context.proposal,
            capability=context.capability,
            binding=context.engine_binding,
            project=context.project_artifact,
            project_validation=context.project_validation,
            input_artifact=context.input_artifact,
            scientific_identity=context.scientific_identity,
            job_artifact_options=dict(context.job_artifact_options),
            live_schema=self.live_schema,
            server=(
                self.execution_server
                if self.surface.profile
                == "command_compiled_approved_execution"
                else self.preview_server
            ),
            repair_parent_sha256=parent.invocation_sha256,
            counterexample_sha256=counterexample.counterexample_sha256,
            repair_attempt=parent.repair_attempt + 1,
        )
        return self._record_compiled_command(turn_id, invocation, context)

    def _preview_command(self, turn_id: str, values: dict) -> Any:
        invocation = self._get(
            self.invocations,
            values["invocation_sha256"],
            "canonical invocation",
        )
        context = self._get(
            self._command_contexts,
            invocation.invocation_sha256,
            "command context",
        )
        expectation = build_preview_expectation(
            program=context.proposal.program,
            jobtype=context.proposal.jobtype,
            input_artifact=context.input_artifact,
            project=context.project_validation,
            charge=context.scientific_identity.charge,
            multiplicity=context.scientific_identity.multiplicity,
        )
        receipt = execute_safe_preview(
            invocation,
            input_artifact=context.input_artifact,
            project_artifact=context.project_artifact,
            expectation=expectation,
            auxiliary_input_artifacts=dict(context.job_artifact_options),
        )
        validator = validator_receipt_from_safe_preview(
            node_id=context.proposal.node_id,
            program=context.proposal.program,
            scientific_identity_sha256=(
                context.scientific_identity.binding_sha256
            ),
            safe_preview=receipt,
        )
        self.safe_previews[receipt.receipt_sha256] = receipt
        self.validators[validator.receipt_sha256] = validator
        self._emit(
            turn_id,
            EventKind.SAFE_PREVIEWED,
            receipt.receipt_sha256,
            status=receipt.status,
            program_validation_status=receipt.program_validation_status,
            critical_finding_count=len(receipt.critical_finding_sha256s),
            invocation_sha256=receipt.invocation_sha256,
        )
        self._emit(
            turn_id,
            EventKind.VALIDATOR_OBSERVED,
            validator.receipt_sha256,
            status=validator.status,
            critical_finding_count=len(validator.critical_finding_sha256s),
            # Record what the validator objected to, not only how many
            # objections there were.  A count cannot be reviewed: when a live
            # session was told twice that a correct ORCA project produced
            # invalid input, nothing durable said which field disagreed, so
            # the defect had to be reproduced by hand to be found at all.
            findings=_public_validator_findings(validator),
            source_receipt_sha256=validator.source_receipt_sha256,
        )
        # The model is the one that has to act on these. Disclosing them
        # only to the event log is what left a live session recompiling
        # against hashes it could not resolve.
        return {
            "safe_preview": receipt,
            "validator": validator,
            "critical_findings": _public_validator_findings(validator),
        }

    def _preflight_program_node(self, turn_id: str, values: dict) -> Any:
        capability = self._get(
            self.capabilities,
            values["capability_receipt_sha256"],
            "capability receipt",
        )
        program_binding = self._get(
            self.program_bindings,
            values["program_binding_sha256"],
            "program binding",
        )
        engine_binding = self._get(
            self.engine_bindings,
            values["engine_binding_sha256"],
            "engine binding",
        )
        invocation = self._get(
            self.invocations,
            values["invocation_sha256"],
            "canonical invocation",
        )
        invocation_plan_sha256 = self._invocation_workflow_plan_sha256s.get(
            invocation.invocation_sha256, ""
        )
        if not invocation_plan_sha256:
            planned = any(
                any(node.node_id == values["node_id"] for node in plan.nodes)
                for plan in self.scientific_workflow_plans.values()
            ) or any(
                any(node.node_id == values["node_id"] for node in draft.nodes)
                for draft in getattr(self, "workflow_drafts", {}).values()
            )
            if planned:
                candidate_plan = self._current_execution_plan_for_node(
                    values["node_id"]
                )
                current_invocation, _context = (
                    self._plan_invocation_for_node(
                        plan=candidate_plan,
                        node_id=values["node_id"],
                    )
                )
                if (
                    current_invocation.invocation_sha256
                    != invocation.invocation_sha256
                ):
                    raise ContractError(
                        "preflight invocation is not the current invocation "
                        "for the selected scientific workflow"
                    )
                invocation_plan_sha256 = candidate_plan.plan_sha256
        if invocation_plan_sha256:
            plan = self.scientific_workflow_plans.get(
                invocation_plan_sha256
            )
            if plan is None:
                raise ContractError(
                    "prepared invocation has no registered scientific plan"
                )
            durable_run = self._durable_workflow_run_state(
                workflow_id=plan.workflow_id,
                plan_sha256=plan.plan_sha256,
            )
            if durable_run is not None:
                frontier = self.event_store.workflow_frontier(
                    workflow_id=plan.workflow_id,
                    run_id=durable_run.run_id,
                )
                observed = next(
                    item
                    for item in durable_run.nodes
                    if item.node_id == values["node_id"]
                )
                ready = derive_ready_node_ids(
                    plan,
                    durable_run,
                    frontier.data_edge_bindings,
                )
                if observed.state != "pending" or values["node_id"] not in ready:
                    raise ContractError(
                        "prepared node is not runnable in the current durable "
                        "workflow frontier"
                    )
        command_inspection = self._get(
            self.command_inspections,
            values["command_inspection_receipt_sha256"],
            "command inspection receipt",
        )
        identity = self._get(
            self.scientific_identities,
            values["scientific_identity_sha256"],
            "scientific identity",
        )
        # Both of these are optional arguments naming receipts the host itself
        # produced for this very invocation, moments earlier.  Omitting either
        # made the preflight ``blocked``; execution then refused three layers
        # away with "node requires a green safe-preview preflight", which names
        # neither omission.  The model's decisions are which project and which
        # node -- carrying the digests back is bookkeeping it cannot get
        # interestingly wrong, so the host does it and an explicit value still
        # wins.
        project_digest = values.get("project_validation_receipt_sha256", "")
        project = (
            self._get(
                self.project_validations,
                project_digest,
                "project validation receipt",
            )
            if project_digest
            # The invocation records the exact validation receipt it was
            # compiled against, so there is nothing to re-derive here.
            else self.project_validations.get(
                invocation.project_receipt_sha256
            )
        )
        preview_digest = values.get("safe_preview_receipt_sha256", "")
        safe_preview = (
            self._get(
                self.safe_previews, preview_digest, "safe preview receipt"
            )
            if preview_digest
            else self._resolve_safe_preview(invocation.invocation_sha256)
        )
        derived_validators = tuple(
            sorted(
                (
                    receipt
                    for receipt in self.validators.values()
                    if safe_preview is not None
                    and receipt.source_receipt_sha256
                    == safe_preview.receipt_sha256
                ),
                key=lambda receipt: receipt.receipt_sha256,
            )
        )
        supplied_validator_ids = tuple(
            values.get("validator_receipt_sha256s") or ()
        )
        if supplied_validator_ids:
            supplied_validators = tuple(
                self._get(
                    self.validators,
                    digest,
                    "program validator receipt",
                )
                for digest in supplied_validator_ids
            )
            if {item.receipt_sha256 for item in supplied_validators} != {
                item.receipt_sha256 for item in derived_validators
            }:
                raise ContractError(
                    "validator receipts differ from the selected safe preview"
                )
        validators = derived_validators
        request = build_program_node_preflight_request(
            node_id=values["node_id"],
            capability_receipt_sha256=capability.receipt_sha256,
            program_binding_sha256=program_binding.binding_sha256,
            engine_binding_sha256=engine_binding.binding_sha256,
            environment_receipt_sha256=(
                engine_binding.environment_receipt_sha256
            ),
            geometry_artifact_sha256=values["geometry_artifact_sha256"],
            scientific_identity_sha256=identity.binding_sha256,
            charge=values["charge"],
            multiplicity=values["multiplicity"],
            project_receipt_sha256=project.receipt_sha256 if project else "",
            invocation_sha256=invocation.invocation_sha256,
            command_inspection_sha256=command_inspection.receipt_sha256,
            validator_receipts=validators,
        )
        receipt = evaluate_program_node_preflight(
            request,
            capability=capability,
            program_binding=program_binding,
            engine_binding=engine_binding,
            project=project,
            invocation=invocation,
            command_inspection=command_inspection,
            scientific_identity=identity,
            validator_receipts=validators,
            safe_preview=safe_preview,
        )
        self.preflights[receipt.receipt_sha256] = receipt
        self._preflight_by_node[values["node_id"]] = receipt
        completion = {
            capability.receipt_sha256,
            program_binding.binding_sha256,
            engine_binding.binding_sha256,
            invocation.invocation_sha256,
            command_inspection.receipt_sha256,
            receipt.receipt_sha256,
            *(item.receipt_sha256 for item in validators),
        }
        if project is not None:
            completion.add(project.receipt_sha256)
        if safe_preview is not None:
            completion.add(safe_preview.receipt_sha256)
        self._completion_sets[receipt.receipt_sha256] = tuple(
            sorted(completion)
        )
        self._emit(
            turn_id,
            EventKind.PROGRAM_PREFLIGHTED,
            receipt.receipt_sha256,
            plan_state=receipt.plan_state,
            critical_finding_count=len(receipt.critical_finding_sha256s),
            safe_preview_receipt_sha256=(receipt.safe_preview_receipt_sha256),
            execution_ready=receipt.execution_ready,
        )
        self._materialize_scientific_workflow(
            turn_id=turn_id, node_id=values["node_id"]
        )
        return receipt

    def _materialize_scientific_workflow(
        self, *, turn_id: str, node_id: str
    ) -> MaterializedWorkflowV1 | None:
        """Ground the latest plan containing ``node_id`` from host receipts."""

        plan = next(
            (
                candidate
                for candidate in reversed(
                    tuple(self.scientific_workflow_plans.values())
                )
                if any(node.node_id == node_id for node in candidate.nodes)
            ),
            None,
        )
        if plan is None:
            return None
        materialized_nodes = []
        unresolved_node_ids = []
        for planned_node in plan.nodes:
            try:
                invocation, context = self._latest_invocation_for_node(
                    planned_node.node_id
                )
            except ContractError:
                unresolved_node_ids.append(planned_node.node_id)
                continue
            bound_plan_sha256 = self._invocation_workflow_plan_sha256s.get(
                invocation.invocation_sha256, ""
            )
            if bound_plan_sha256 and bound_plan_sha256 != plan.plan_sha256:
                unresolved_node_ids.append(planned_node.node_id)
                continue
            self._invocation_workflow_plan_sha256s[
                invocation.invocation_sha256
            ] = plan.plan_sha256
            project = context.project_validation
            project_artifact = context.project_artifact
            environment_sha256 = (
                context.engine_binding.environment_receipt_sha256
            )
            if (
                project is None
                or project.status != "valid"
                or project_artifact is None
                or not environment_sha256
            ):
                unresolved_node_ids.append(planned_node.node_id)
                continue
            preflight = self._preflight_by_node.get(planned_node.node_id)
            previewed = self._node_is_previewed(
                planned_node.node_id,
                plan_sha256=plan.plan_sha256,
            )
            materialized_nodes.append(
                MaterializedNodeV1(
                    node_id=planned_node.node_id,
                    input_artifact_sha256=context.input_artifact.sha256,
                    project_artifact_sha256=project_artifact.sha256,
                    project_validation_receipt_sha256=project.receipt_sha256,
                    environment_receipt_sha256=environment_sha256,
                    invocation_sha256=invocation.invocation_sha256,
                    invocation_identity_sha256=(
                        self._invocation_identity(
                            planned_node.node_id,
                            plan_sha256=plan.plan_sha256,
                        )
                    ),
                    preflight_receipt_sha256=(
                        preflight.receipt_sha256 if previewed else ""
                    ),
                    state="previewed" if previewed else "compiled",
                    auxiliary_input_bindings=(
                        invocation.auxiliary_input_bindings
                    ),
                )
            )
        if unresolved_node_ids:
            status = "partial"
        elif materialized_nodes and all(
            node.state == "previewed" for node in materialized_nodes
        ):
            status = "previewed"
        else:
            status = "materialized"
        resource_sha256 = (
            self.execution_resources.resource_sha256
            if self.execution_resources is not None
            else canonical_sha256(
                {
                    "schema_version": "chemsmart.preview-resource.v1",
                    "chemistry_engine_calls": 0,
                }
            )
        )
        workflow = build_materialized_workflow(
            plan=plan,
            live_cli_schema_sha256=self.live_schema.schema_sha256,
            resource_sha256=resource_sha256,
            nodes=tuple(materialized_nodes),
            unresolved_node_ids=tuple(unresolved_node_ids),
            status=status,
        )
        self.materialized_workflows[workflow.materialized_sha256] = workflow
        self.event_store.record_materialized_workflow(
            turn_id=turn_id, workflow=workflow
        )
        return workflow

    def _completion_receipts_for_latest_analysis_toolchain(
        self,
    ) -> tuple[str, ...]:
        """Aggregate a completed registered-result toolchain, if one exists.

        The scientific toolchain itself is the transparent completion policy:
        required output producers must have typed evidence, and a planned
        claim-rendering stage must have a task-bound claim record cited by a
        task-bound scientific decision.  This observes execution of the
        model's declared analysis DAG; it does not grade the numerical result
        or force one scientifically preferred DAG.
        """

        if not self.scientific_toolchain_plans:
            return ()
        plan = next(reversed(self.scientific_toolchain_plans.values()))
        if (
            not self.workflow_drafts
            or next(reversed(self.workflow_drafts))
            != plan.command_workflow_draft_sha256
        ):
            return ()
        command_result = self._scientific_toolchain_command_results.get(
            plan.plan_sha256
        )
        if command_result is None:
            return ()
        draft = command_result.get("workflow_draft")
        if (
            not isinstance(draft, CommandWorkflowDraftV1)
            or plan.calculation_node_ids
            or draft.nodes
        ):
            return ()
        matched = self._scientific_toolchain_analysis_receipts(
            plan,
            task_spec_sha256=draft.task_spec_id,
        )
        required_validation_nodes = {
            node.node_id
            for node in plan.analysis_nodes
            if node.analysis_kind == "scientific_validation"
            and node.support_state == "planned"
        }
        if required_validation_nodes.difference(matched):
            return ()
        output_producers = {
            output.output_id: node.node_id
            for node in plan.analysis_nodes
            for output in node.outputs
        }
        required_producers = {
            output_producers[output_id]
            for output_id in plan.required_output_ids
        }
        analysis_nodes = {
            node.node_id: node for node in plan.analysis_nodes
        }
        missing_required = {
            node_id for node_id in required_producers if node_id not in matched
        }
        unresolved_required = {
            node_id
            for node_id in missing_required
            if analysis_nodes[node_id].support_state != "blocked_unsupported"
        }
        if unresolved_required:
            return ()
        claim_nodes = tuple(
            node.node_id
            for node in plan.analysis_nodes
            if node.analysis_kind == "claim_rendering"
            and node.support_state == "planned"
        )
        rendered_claim_nodes = {
            node_id for node_id in claim_nodes if node_id in matched
        }
        if claim_nodes and not rendered_claim_nodes:
            return ()
        # A required blocked_unsupported output is a completed limitation only
        # when the same task reached its recorded claims/decision stage.  Keep
        # the node visibly blocked in the frontier; do not pretend a missing
        # parser or Hessian produced a value.
        if missing_required and not rendered_claim_nodes:
            return ()
        source_receipts = tuple(
            sorted(
                {
                    digest
                    for receipts in matched.values()
                    for digest in receipts
                }
            )
        )
        if not source_receipts:
            return ()
        body = {
            "schema_version": "chemsmart.analysis-completion-receipt.v1",
            # A registered-result scientific toolchain is already a visible,
            # typed output contract.  Its digest fills the existing aggregate
            # gate's policy identity without inventing a parallel policy file.
            "policy_sha256": plan.plan_sha256,
            "task_spec_sha256": draft.task_spec_id,
            "source_receipt_sha256s": source_receipts,
            "status": "passed",
            "findings": (),
        }
        completion = AnalysisCompletionReceiptV1(
            **body, receipt_sha256=canonical_sha256(body)
        )
        self.analysis_completion_receipts[completion.receipt_sha256] = (
            completion
        )
        completion_record = canonical_data(completion)
        completion_record.pop("receipt_sha256")
        turn_id = self.event_store.state().turn_id or "analysis-toolchain"
        self.event_store.append(
            turn_id=turn_id,
            kind=EventKind.ANALYSIS_COMPLETION_EVALUATED.value,
            payload={
                "receipt_sha256": completion.receipt_sha256,
                "policy_sha256": completion.policy_sha256,
                "task_spec_sha256": completion.task_spec_sha256,
                "source_receipt_sha256s": completion.source_receipt_sha256s,
                "status": completion.status,
                "critical_finding_count": 0,
                "completion_kind": "scientific_toolchain",
                "record": completion_record,
            },
            idempotency_key=(
                "scientific-toolchain-completion:" + completion.receipt_sha256
            ),
        )
        return (completion.receipt_sha256,)

    def completion_receipts_for_latest_preflight(self) -> tuple[str, ...]:
        """Return a green analysis-toolchain or command-preflight gate."""

        analysis_completion = (
            self._completion_receipts_for_latest_analysis_toolchain()
        )
        if analysis_completion:
            return analysis_completion

        if not self.preflights:
            raise ContractError("no node preflight has been observed")
        latest = next(reversed(self.preflights.values()))
        if latest.plan_state != "previewed" or latest.critical_finding_sha256s:
            raise ContractError("latest node preflight is not green")
        completion = list(self._completion_sets[latest.receipt_sha256])
        if self.surface.profile == "command_compiled_approved_execution":
            if not self.execution_receipts:
                raise ContractError(
                    "approved workflow has not executed a node"
                )
            approval = self.workflow_execution_approval
            if approval is None:
                raise ContractError("bounded workflow has not been admitted")
            required_nodes = tuple(
                item.node_id for item in approval.node_bindings
            )
            missing = [
                node_id
                for node_id in required_nodes
                if node_id not in self.execution_receipts
                or not self.execution_receipts[node_id].validated
            ]
            if missing:
                raise ContractError(
                    "approved workflow execution remains incomplete: "
                    + ", ".join(missing)
                )
            completion.extend(
                self.execution_receipts[node_id].receipt_sha256
                for node_id in required_nodes
            )
        return tuple(sorted(set(completion)))

    def completion_receipts_for_analysis(
        self, *, turn_id: str
    ) -> tuple[str, ...]:
        """Return receipts only when the task-owned analysis policy is met."""

        policy = self.analysis_completion_policy
        if policy is None:
            raise ContractError("no analysis completion policy is active")
        selected: list[str] = []
        evidence_receipts: list[str] = []
        downstream_source_receipts = {
            claim.source_receipt_sha256
            for record in self.analysis_claim_records.values()
            for claim in record.claims
        }
        downstream_source_receipts.update(
            digest
            for receipt in self.quantity_expression_receipts.values()
            for dependency in receipt.output_dependencies
            for digest in dependency.source_receipt_sha256s
        )
        downstream_expression_receipts = {
            claim.source_receipt_sha256
            for record in self.analysis_claim_records.values()
            for claim in record.claims
            if claim.source_kind == "quantity_expression"
        }
        decision_receipts = {
            parsed[1]
            for decision in self.scientific_decisions.values()
            for reference in decision.evidence_refs
            if (parsed := _postprocessing_evidence_reference(reference))
            is not None
        }
        target_artifacts = set(policy.target_artifact_sha256s)
        stage_receipts: dict[str, dict[str, str]] = {
            "quantity_extraction": {},
            "thermochemistry": {},
        }

        if "quantity_extraction" in policy.required_stages:
            for artifact_sha256 in sorted(target_artifacts):
                matches = tuple(
                    receipt
                    for receipt in self.quantity_extractions.values()
                    if (
                        receipt.status == "extracted"
                        and receipt.artifact_sha256 == artifact_sha256
                        and set(policy.required_extraction_selectors).issubset(
                            self.quantity_extraction_selectors.get(
                                receipt.receipt_sha256, ()
                            )
                        )
                    )
                )
                if not matches:
                    raise AnalysisIncompleteError(
                        "required quantity extraction has not been observed "
                        f"for artifact {artifact_sha256}"
                    )
                preferred = (
                    tuple(
                        receipt
                        for receipt in matches
                        if receipt.receipt_sha256 in downstream_source_receipts
                    )
                    or matches
                )
                receipt = sorted(
                    preferred, key=lambda item: item.receipt_sha256
                )[-1]
                stage_receipts["quantity_extraction"][artifact_sha256] = (
                    receipt.receipt_sha256
                )
                selected.append(receipt.receipt_sha256)
                evidence_receipts.append(receipt.receipt_sha256)

        if "thermochemistry" in policy.required_stages:
            required_ids = set(policy.required_thermochemistry_quantity_ids)
            for artifact_sha256 in sorted(target_artifacts):
                matches = []
                for receipt in self.thermochemistry_receipts.values():
                    observed_ids = {
                        quantity.quantity_id for quantity in receipt.quantities
                    }
                    if (
                        receipt.status == "derived"
                        and receipt.artifact_sha256 == artifact_sha256
                        and required_ids.issubset(observed_ids)
                        and (
                            policy.temperature_k is None
                            or math.isclose(
                                receipt.temperature_k,
                                policy.temperature_k,
                                rel_tol=0.0,
                                abs_tol=1.0e-12,
                            )
                        )
                        and (
                            policy.pressure_atm is None
                            or math.isclose(
                                receipt.pressure_atm,
                                policy.pressure_atm,
                                rel_tol=0.0,
                                abs_tol=1.0e-12,
                            )
                        )
                    ):
                        matches.append(receipt)
                if not matches:
                    raise AnalysisIncompleteError(
                        "required thermochemistry receipt has not been observed "
                        f"for artifact {artifact_sha256}"
                    )
                preferred = tuple(
                    receipt
                    for receipt in matches
                    if receipt.receipt_sha256 in downstream_source_receipts
                ) or tuple(matches)
                receipt = sorted(
                    preferred, key=lambda item: item.receipt_sha256
                )[-1]
                stage_receipts["thermochemistry"][artifact_sha256] = (
                    receipt.receipt_sha256
                )
                selected.append(receipt.receipt_sha256)
                evidence_receipts.append(receipt.receipt_sha256)

        if "quantity_expression" in policy.required_stages:
            matches: list[Any] = []
            expression_stage_receipts: dict[str, str] = {}
            selected_source_refs = tuple(
                f"receipt:{receipt_sha256}"
                for receipt_sha256 in evidence_receipts
            )
            for requirement in policy.required_expressions:
                requirement_matches = []
                for receipt in self.quantity_expression_receipts.values():
                    observed_output_ids = {
                        quantity.quantity_id for quantity in receipt.outputs
                    }
                    output_dependencies = {
                        dependency.output_id: set(
                            dependency.source_receipt_sha256s
                        )
                        for dependency in receipt.output_dependencies
                    }
                    required_sources_by_output: dict[str, set[str]] = {
                        output_id: set()
                        for output_id in requirement.required_output_ids
                    }
                    for source in requirement.required_sources:
                        source_artifacts = (
                            set(source.artifact_sha256s)
                            if source.artifact_sha256s
                            else target_artifacts
                        )
                        source_receipts = {
                            stage_receipts[source.stage].get(
                                artifact_sha256, ""
                            )
                            for artifact_sha256 in source_artifacts
                        }
                        source_receipts.discard("")
                        target_outputs = (
                            source.output_ids
                            or requirement.required_output_ids
                        )
                        for output_id in target_outputs:
                            required_sources_by_output[output_id].update(
                                source_receipts
                            )
                    if any(required_sources_by_output.values()):
                        dependency_ok = all(
                            required_sources_by_output[output_id].issubset(
                                output_dependencies.get(output_id, set())
                            )
                            for output_id in requirement.required_output_ids
                        )
                    elif selected_source_refs:
                        selected_source_receipts = {
                            reference.removeprefix("receipt:")
                            for reference in selected_source_refs
                        }
                        dependency_ok = all(
                            bool(
                                selected_source_receipts.intersection(
                                    output_dependencies.get(output_id, set())
                                )
                            )
                            for output_id in requirement.required_output_ids
                        )
                    else:
                        dependency_ok = True
                    if (
                        receipt.status == "derived"
                        and receipt.expression_id == requirement.expression_id
                        and set(requirement.required_output_ids).issubset(
                            observed_output_ids
                        )
                        and dependency_ok
                        and (
                            not requirement.semantic_signature_sha256
                            or receipt.semantic_signature_sha256
                            == requirement.semantic_signature_sha256
                        )
                    ):
                        requirement_matches.append(receipt)
                if not requirement_matches:
                    raise AnalysisIncompleteError(
                        "required named quantity expression has not been observed: "
                        + requirement.expression_id
                    )
                preferred = tuple(
                    receipt
                    for receipt in requirement_matches
                    if receipt.receipt_sha256
                    in downstream_expression_receipts | decision_receipts
                ) or tuple(requirement_matches)
                chosen_expression = sorted(
                    preferred, key=lambda item: item.receipt_sha256
                )[-1]
                matches.append(chosen_expression)
                expression_stage_receipts[requirement.expression_id] = (
                    chosen_expression.receipt_sha256
                )
            if len(matches) < policy.minimum_expression_receipts:
                raise AnalysisIncompleteError(
                    "required quantity expressions have not been observed"
                )
            selected.extend(receipt.receipt_sha256 for receipt in matches)
            evidence_receipts.extend(
                receipt.receipt_sha256 for receipt in matches
            )

        if "analysis_claims" in policy.required_stages:
            matching_records = []
            for record in self.analysis_claim_records.values():
                if record.task_spec_sha256 != policy.task_spec_sha256:
                    continue
                claims = {claim.claim_id: claim for claim in record.claims}
                if all(
                    requirement.claim_id in claims
                    and claims[requirement.claim_id].source_kind
                    == requirement.source_kind
                    and claims[requirement.claim_id].quantity_id
                    == requirement.quantity_id
                    and claims[requirement.claim_id].display_unit
                    == requirement.display_unit
                    and claims[requirement.claim_id].source_receipt_sha256
                    in evidence_receipts
                    and (
                        not requirement.source_artifact_sha256s
                        or claims[requirement.claim_id].source_receipt_sha256
                        in {
                            stage_receipts[requirement.source_kind][artifact]
                            for artifact in requirement.source_artifact_sha256s
                        }
                    )
                    and (
                        not requirement.source_expression_id
                        or claims[requirement.claim_id].source_receipt_sha256
                        == expression_stage_receipts.get(
                            requirement.source_expression_id
                        )
                    )
                    and (
                        not requirement.source_selector
                        or self.quantity_extraction_bindings.get(
                            claims[requirement.claim_id].source_receipt_sha256,
                            {},
                        ).get(requirement.quantity_id)
                        == requirement.source_selector
                    )
                    for requirement in policy.required_claims
                ):
                    matching_records.append(record)
            if not matching_records:
                raise AnalysisIncompleteError(
                    "required host-rendered analysis claims have not been observed"
                )
            preferred_records = tuple(
                record
                for record in matching_records
                if record.receipt_sha256 in decision_receipts
            ) or tuple(matching_records)
            record = sorted(
                preferred_records, key=lambda item: item.receipt_sha256
            )[-1]
            selected.append(record.receipt_sha256)
            evidence_receipts.append(record.receipt_sha256)

        if "scientific_decision" in policy.required_stages:
            matches = tuple(
                decision
                for decision in self.scientific_decisions.values()
                if decision.task_spec_sha256 == policy.task_spec_sha256
                and (
                    not policy.require_decision_evidence_binding
                    or set(evidence_receipts).issubset(
                        {
                            parsed[1]
                            for reference in decision.evidence_refs
                            if (
                                parsed := _postprocessing_evidence_reference(
                                    reference
                                )
                            )
                            is not None
                        }
                    )
                )
            )
            if not matches:
                raise AnalysisIncompleteError(
                    "required evidence-bound scientific decision is absent"
                )
            selected.append(
                sorted(matches, key=lambda item: item.record_sha256)[
                    -1
                ].record_sha256
            )

        completion = build_analysis_completion_receipt(
            policy=policy,
            source_receipt_sha256s=tuple(sorted(set(selected))),
        )
        self.analysis_completion_receipts[completion.receipt_sha256] = (
            completion
        )
        completion_record = canonical_data(completion)
        completion_record.pop("receipt_sha256")
        self.event_store.append(
            turn_id=turn_id,
            kind=EventKind.ANALYSIS_COMPLETION_EVALUATED.value,
            payload={
                "receipt_sha256": completion.receipt_sha256,
                "policy_sha256": completion.policy_sha256,
                "task_spec_sha256": completion.task_spec_sha256,
                "source_receipt_sha256s": (completion.source_receipt_sha256s),
                "status": completion.status,
                "critical_finding_count": 0,
                "completion_kind": "numerical_analysis",
                "record": completion_record,
                "policy_record": policy.public_record(),
            },
            idempotency_key=(
                "analysis-completion:" + completion.receipt_sha256
            ),
        )
        return (completion.receipt_sha256,)

    def render_completed_analysis_report(
        self, completion_receipt_sha256: str
    ) -> str:
        """Render the authoritative result from host-owned typed records.

        Provider prose remains visible in the public transcript, but it is not
        allowed to become the numerical answer.  This renderer copies values
        from the exact claim record admitted by the completion receipt and
        uses only the corresponding structured scientific decision for
        interpretation.
        """

        try:
            completion = self.analysis_completion_receipts[
                completion_receipt_sha256
            ]
        except KeyError as exc:
            raise ContractError("unknown analysis completion receipt") from exc
        source_receipts = set(completion.source_receipt_sha256s)
        claim_records = tuple(
            record
            for digest, record in self.analysis_claim_records.items()
            if digest in source_receipts
        )
        decisions = tuple(
            decision
            for decision in self.scientific_decisions.values()
            if decision.record_sha256 in source_receipts
        )
        if len(claim_records) != 1 or len(decisions) != 1:
            raise ContractError(
                "completed analysis must bind one claim record and one decision"
            )
        claims = claim_records[0]
        decision = decisions[0]
        policy = self.analysis_completion_policy
        if policy is None or policy.policy_sha256 != completion.policy_sha256:
            raise ContractError(
                "completed analysis policy is absent or differs from the receipt"
            )
        lines = [
            "# Host-validated structured analysis",
            "",
            f"Completion receipt: `{completion.receipt_sha256}`",
            f"Claim record: `{claims.receipt_sha256}`",
        ]
        if policy.required_conditions:
            lines.extend(
                (
                    "",
                    "## Task-owned conditions",
                    "",
                    "| Condition | Value | Unit | Origin | Evidence |",
                    "|---|---:|---|---|---|",
                )
            )
            for condition in policy.required_conditions:
                lines.append(
                    f"| {condition.condition_id} | `{condition.value}` | "
                    f"{condition.unit} | {condition.origin} | "
                    f"`{condition.evidence_ref}` |"
                )
        lines.extend(
            (
                "",
                "## Host-rendered numerical claims",
                "",
                "| Claim | Value | Unit | Source receipt |",
                "|---|---:|---|---|",
            )
        )
        for claim in claims.claims:
            display = json.dumps(
                canonical_data(claim.display_value),
                ensure_ascii=False,
                separators=(",", ":"),
            )
            lines.append(
                f"| {claim.claim_id} | `{display}` | {claim.display_unit} | "
                f"`{claim.source_receipt_sha256}` |"
            )
        sections = (
            ("Method rationale", (decision.method_rationale,)),
            ("Assumptions", decision.assumptions),
            ("Diagnostics", decision.diagnostics),
            ("Uncertainties", decision.uncertainties),
            ("Alternatives", decision.alternatives),
        )
        for title, values in sections:
            entries = tuple(value for value in values if value)
            if not entries:
                continue
            lines.extend(("", f"## {title}", ""))
            lines.extend(f"- {value}" for value in entries)
        return "\n".join(lines)

    def _execute_approved_program_node(
        self, turn_id: str, values: dict
    ) -> Any:
        """Resolve and execute one approved node without model-authored argv."""

        if self.surface.profile != "command_compiled_approved_execution":
            raise ContractError("execution tool is absent from this profile")
        node_id = values["node_id"]
        frozen_approval = getattr(self, "frozen_workflow_approval", None)
        if (
            self.workflow_execution_approval is not None
            and frozen_approval is None
        ):
            raise ContractError(
                "legacy V1 approval is preview-only; Runtime V2 frozen approval "
                "is required for execution"
            )
        current_plan = self._current_execution_plan_for_node(node_id)
        invocation, context = self._plan_invocation_for_node(
            plan=current_plan,
            node_id=node_id,
        )
        if self.workflow_execution_approval is None:
            raise ContractError(
                "execution requires a human-approved one-shot workflow bundle; "
                "a resource envelope alone never grants authority"
            )
        approval = self.workflow_execution_approval
        if approval is None:  # pragma: no cover - admission narrows this
            raise ContractError("execution approval was not created")
        scientific_plan = self._execution_scientific_plan()
        if scientific_plan.plan_sha256 != current_plan.plan_sha256:
            raise ContractError(
                "a separate workflow owns the immutable execution approval; "
                "start a new independent attempt for the current workflow"
            )
        if frozen_approval is None:
            raise ContractError(
                "legacy V1 approval is preview-only; Runtime V2 frozen approval "
                "is required for execution"
            )
        v2_run_id = "run." + approval.approval_id
        replayed = self.event_store.replayed_execution_receipt(
            workflow_id=scientific_plan.workflow_id,
            run_id=v2_run_id,
            node_id=node_id,
        )
        if replayed is not None:
            self.execution_receipts[node_id] = replayed
            result_validation = self.result_validation_receipts.get(
                replayed.result_validation_receipt_sha256
            )
            return {
                "execution": replayed,
                "idempotent_replay": True,
                "handoff": self.handoffs.get(node_id),
                "result_validation": result_validation,
            }
        existing = self.execution_receipts.get(node_id)
        if existing is not None:
            raise ContractError(
                "process-local execution receipt lacks replay evidence"
            )
        approved_node = approval.node(node_id)
        if (
            context.project_artifact is None
            or context.project_validation is None
        ):
            raise ContractError(
                "execution requires a validated project artifact"
            )
        if context.project_validation.status != "valid":
            raise ContractError("execution project loader gate is red")
        if (
            context.project_validation.settings_sha256
            != approved_node.settings_sha256
        ):
            raise ContractError(
                "effective project settings differ from approval"
            )
        if (
            context.project_artifact.sha256
            != approved_node.project_artifact_sha256
        ):
            raise ContractError("project bytes differ from workflow approval")
        if not self._invocation_has_green_preflight(
            invocation.invocation_sha256
        ):
            raise ContractError("node requires a green safe-preview preflight")

        handoff = self.handoffs.get(node_id)
        real_argv = self._real_execution_argv(invocation)
        execution_invocation = build_program_execution_invocation(
            node_id=node_id,
            approval=approval,
            project_artifact=context.project_artifact,
            input_artifact=context.input_artifact,
            scientific_identity_sha256=(
                context.scientific_identity.binding_sha256
            ),
            environment_receipt_sha256=(
                context.engine_binding.environment_receipt_sha256
            ),
            resources=self.execution_resources,
            argv=real_argv,
            handoff=handoff,
            environment_identity_sha256=self._environment_identity_for(
                context.engine_binding.environment_receipt_sha256
            ),
            # Identity of the *compiled* command, not of ``real_argv``: the
            # preview froze the compiled one, and the host rewrite that adds
            # --no-fake, successful-run scratch cleanup, and the resource
            # flags is deterministic from it.
            invocation_identity_sha256=self._invocation_identity(
                node_id,
                plan_sha256=scientific_plan.plan_sha256,
            ),
            auxiliary_input_artifacts=dict(context.job_artifact_options),
            auxiliary_handoffs=(
                {"hess_filename": self.hessian_handoffs[node_id]}
                if "hess_filename" in dict(context.job_artifact_options)
                and node_id in self.hessian_handoffs
                else {}
            ),
        )
        frozen_preview = frozen_approval.preview_binding(node_id)
        if (
            frozen_preview is not None
            and frozen_preview.auxiliary_input_bindings
            != execution_invocation.auxiliary_input_bindings
        ):
            raise ContractError(
                "auxiliary inputs differ from the frozen safe preview"
            )
        if self.frozen_workflow_approval is not None:
            future_rules = self.frozen_workflow_approval.producer_rules_for(
                node_id
            )
            if future_rules:
                # A future producer rule still pins a receipt digest, so its
                # identity has to be resolved through the observed receipts.
                approved = self._environment_identity_is_approved(
                    execution_invocation.environment_receipt_sha256,
                    {item.environment_receipt_sha256 for item in future_rules},
                )
            else:
                approved = (
                    execution_invocation.environment_identity_sha256
                    in (
                        set(
                            self.frozen_workflow_approval.environment_identity_sha256s
                        )
                    )
                )
            if not approved:
                raise ContractError(
                    "execution environment differs from the exact frozen "
                    "node approval"
                )
        node_workspace = self.approved_workspace / "nodes" / node_id
        _prepare_execution_node_workspace(node_workspace)
        started = datetime.now(timezone.utc).isoformat()
        if frozen_approval.plan_sha256 != scientific_plan.plan_sha256:
            raise ContractError(
                "execution plan differs from frozen workflow approval"
            )
        frontier = self.event_store.workflow_frontier(
            workflow_id=scientific_plan.workflow_id,
            run_id=v2_run_id,
        )
        materialized = frontier.materialized_workflow
        if materialized is None:
            materialized = self.materialized_workflows.get(
                frozen_approval.materialized_workflow_sha256
            )
        if materialized is None:
            raise ContractError(
                "frozen workflow approval has no canonical materialization"
            )
        effective_timeout_seconds = self._require_bounded_launch_budget()
        fence = self.event_store.reserve_workflow_node_launch(
            turn_id=turn_id,
            plan=scientific_plan,
            materialized_workflow=materialized,
            approval=frozen_approval,
            invocation=execution_invocation,
            run_id=v2_run_id,
            timestamp=started,
        )
        if fence.status == "terminal_replay":
            replayed = fence.execution_receipt
            if replayed is None:  # pragma: no cover - contract narrows this
                raise ContractError(
                    "terminal replay lacks an execution receipt"
                )
            self.execution_receipts[node_id] = replayed
            return {
                "execution": replayed,
                "idempotent_replay": True,
                "handoff": self.handoffs.get(node_id),
            }
        command = list(real_argv)
        environment = _program_process_environment(
            overrides=self.execution_environment,
            remove=self.execution_environment_remove,
        )
        source_root = str(Path(__file__).resolve().parents[2])
        current_pythonpath = environment.get("PYTHONPATH", "")
        environment["PYTHONPATH"] = (
            source_root
            if not current_pythonpath
            else source_root + os.pathsep + current_pythonpath
        )
        _require_current_auxiliary_inputs(context.job_artifact_options)
        launch_ambiguous = False
        with ProcessSignalGuard() as signal_guard:
            try:
                process = subprocess.Popen(
                    command,
                    cwd=node_workspace,
                    env=environment,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                    start_new_session=True,
                )
            except OSError as exc:
                process_observation = launch_failure_observation(
                    timeout_seconds=effective_timeout_seconds,
                    memory_limit_mb=(
                        self.execution_resources.memory_gb * 1024.0
                    ),
                    error_type=type(exc).__name__,
                )
                stdout_text = ""
                stderr_text = type(exc).__name__
            else:
                process_result = observe_process(
                    process,
                    timeout_seconds=effective_timeout_seconds,
                    memory_limit_mb=(
                        self.execution_resources.memory_gb * 1024.0
                    ),
                    signal_guard=signal_guard,
                )
                process_observation = process_result.observation
                stdout_text = _public_process_stream(process_result.stdout)
                stderr_text = _public_process_stream(process_result.stderr)
        _write_host_execution_artifact(
            node_workspace / "controller.stdout", stdout_text
        )
        _write_host_execution_artifact(
            node_workspace / "controller.stderr", stderr_text
        )
        _write_host_execution_artifact(
            node_workspace / "execution-resource.receipt.json",
            json.dumps(
                process_observation.as_dict(),
                sort_keys=True,
                separators=(",", ":"),
            )
            + "\n",
        )
        wrapper_exit_status = process_observation.returncode
        launch_ambiguous = process_observation.state.endswith("_ambiguous")
        outputs = self._execution_output_artifacts(
            node_id,
            node_workspace,
            program=context.proposal.program,
            auxiliary_inputs=tuple(
                artifact for _name, artifact in context.job_artifact_options
            ),
        )
        pyscf_engine = (
            _inspect_pyscf_engine_observation(
                outputs, launch_ambiguous=launch_ambiguous
            )
            if context.proposal.program == "pyscf"
            else None
        )
        capability_environment = self.environments.get(
            execution_invocation.environment_receipt_sha256
        )
        if capability_environment is None:
            raise ContractError(
                "execution environment capability receipt is unavailable"
            )
        evaluation = self._evaluate_execution_outputs(
            program=context.proposal.program,
            jobtype=context.proposal.jobtype,
            charge=context.scientific_identity.charge,
            multiplicity=context.scientific_identity.multiplicity,
            expected_settings=dict(context.project_validation.settings),
            expected_input_artifact=context.input_artifact,
            expected_project_artifact=context.project_artifact,
            output_artifacts=outputs,
            exit_status=wrapper_exit_status,
            expected_environment_receipt_sha256=(
                execution_invocation.environment_receipt_sha256
            ),
            capability_environment_receipt=capability_environment,
            pyscf_engine_observation=pyscf_engine,
            stationary_point_policy=(
                self.stationary_point_policy
                if self.stationary_point_policy is not None
                and self.stationary_point_policy.hessian_node_id == node_id
                else None
            ),
            approved_stationary_point_policy_sha256=(
                frozen_approval.stationary_point_policy_sha256
                if self.stationary_point_policy is not None
                and self.stationary_point_policy.hessian_node_id == node_id
                else ""
            ),
            approved_hessian_node_id=(
                node_id if approved_node.jobtype == "hess" else ""
            ),
            process_observation=process_observation,
        )
        staged_auxiliary_findings = _staged_auxiliary_input_findings(
            node_workspace=node_workspace,
            job_artifact_options=context.job_artifact_options,
        )
        if staged_auxiliary_findings:
            evaluation = replace(
                evaluation,
                findings=tuple(
                    sorted(
                        {
                            *evaluation.findings,
                            *staged_auxiliary_findings,
                        }
                    )
                ),
            )
        result_validation_receipt = build_program_result_validation_receipt(
            invocation=execution_invocation,
            validator_id=evaluation.validator_id,
            validator_schema_version=(evaluation.validator_schema_version),
            validator_version=evaluation.validator_version,
            input_artifact=context.input_artifact,
            project_artifact=context.project_artifact,
            capability_environment_receipt_sha256=(
                capability_environment.receipt_sha256
            ),
            output_artifacts=outputs,
            observations=evaluation.observations,
            findings=evaluation.findings,
            run_environment_receipt_sha256=(
                evaluation.run_environment_receipt_sha256
            ),
            environment_validation_sha256=(
                evaluation.environment_validation_sha256
            ),
            stationary_point_policy_sha256=(
                self.stationary_point_policy.policy_sha256
                if approved_node.jobtype == "hess"
                and self.stationary_point_policy is not None
                else ""
            ),
        )
        self.result_validation_receipts[
            result_validation_receipt.receipt_sha256
        ] = result_validation_receipt
        self._emit(
            turn_id,
            EventKind.RESULT_VERIFIED,
            result_validation_receipt.receipt_sha256,
            status=result_validation_receipt.state,
            critical_finding_count=len(result_validation_receipt.findings),
            node_id=node_id,
            record=canonical_data(result_validation_receipt),
        )
        validated = result_validation_receipt.state == "valid"
        validator_sha256s = (result_validation_receipt.receipt_sha256,)
        findings = result_validation_receipt.findings
        if launch_ambiguous:
            execution_state = "ambiguous"
            engine_complete = False
            child_exit_status = None
            engine_receipt_sha256 = ""
        elif pyscf_engine is not None:
            engine_complete = pyscf_engine.engine_complete
            child_exit_status = pyscf_engine.child_exit_status
            engine_receipt_sha256 = pyscf_engine.run_receipt_sha256
            execution_state = (
                "failed"
                if process_observation.state != "exited"
                else "engine_complete"
                if engine_complete
                else "failed"
            )
        else:
            engine_complete = wrapper_exit_status == 0
            if context.proposal.program == "orca":
                engine_complete = bool(
                    engine_complete
                    and (evaluation.observations.get("orca") or {}).get(
                        "normal_termination"
                    )
                )
            child_exit_status = wrapper_exit_status
            engine_receipt_sha256 = ""
            execution_state = (
                "engine_complete"
                if engine_complete and process_observation.state == "exited"
                else "failed"
            )
        if validated:
            if execution_state != "engine_complete" or not engine_complete:
                validated = False
                findings = tuple(
                    sorted(
                        {
                            *findings,
                            "execution.validation_without_engine_completion",
                        }
                    )
                )
            else:
                execution_state = "validated"
        finished = datetime.now(timezone.utc).isoformat()
        receipt = build_program_execution_receipt(
            execution_invocation,
            execution_state=execution_state,
            exit_status=wrapper_exit_status,
            child_exit_status=child_exit_status,
            engine_complete=engine_complete,
            validated=validated,
            engine_receipt_sha256=engine_receipt_sha256,
            result_validation_receipt_sha256=(
                result_validation_receipt.receipt_sha256
            ),
            output_artifacts=outputs,
            validator_receipt_sha256s=validator_sha256s,
            findings=findings,
            started_at=started,
            finished_at=finished,
        )
        self.event_store.record_program_execution_receipt(
            turn_id=turn_id,
            workflow_id=scientific_plan.workflow_id,
            run_id=v2_run_id,
            receipt=receipt,
        )
        self.execution_receipts[node_id] = receipt
        produced_handoffs = []
        pending_data_edges = []
        if receipt.validated and context.proposal.program in {
            "gaussian",
            "orca",
            "pyscf",
            "xtb",
        }:
            outgoing_edges = tuple(
                sorted(
                    (
                        edge
                        for edge in approval.producer_edges
                        if edge.producer_node_id == node_id
                    ),
                    key=lambda item: (
                        item.selection_rule
                        != "validated_optimized_geometry",
                        item.consumer_node_id,
                    ),
                )
            )
            for edge in outgoing_edges:
                if edge.producer_node_id != node_id:
                    continue
                consumer_binding = approval.node(edge.consumer_node_id)
                if edge.selection_rule == "validated_final_orca_ts_hessian":
                    if (
                        context.proposal.program != "orca"
                        or context.proposal.jobtype != "ts"
                    ):
                        raise ContractError(
                            "final ORCA Hessian handoff requires an ORCA TS"
                        )
                    result_candidates = tuple(
                        item for item in outputs if item.kind == "orca_output"
                    )
                    hessian_candidates = tuple(
                        item for item in outputs if item.kind == "orca_hessian"
                    )
                    if len(result_candidates) != 1 or not hessian_candidates:
                        raise ContractError(
                            "validated ORCA TS requires one output and at least "
                            "one native Hessian candidate"
                        )
                    artifact, observed = handoff_final_orca_ts_hessian(
                        producer_receipt=receipt,
                        result_artifact=result_candidates[0],
                        hessian_candidates=hessian_candidates,
                        producer_edge=edge,
                        approved_workspace=self.approved_workspace,
                        hessian_artifact_id=(
                            f"hessian.{edge.producer_node_id}-to-"
                            f"{edge.consumer_node_id}"
                        ),
                        expected_charge=context.scientific_identity.charge,
                        expected_multiplicity=(
                            context.scientific_identity.multiplicity
                        ),
                    )
                    geometry_handoff = self.handoffs.get(edge.consumer_node_id)
                    if geometry_handoff is None:
                        raise ContractError(
                            "ORCA TS Hessian handoff requires its final geometry"
                        )
                    if observed.consumer_state != (
                        consumer_binding.charge,
                        consumer_binding.multiplicity,
                    ):
                        raise ContractError(
                            "ORCA IRC must remain on the transition-state "
                            "charge and multiplicity surface"
                        )
                    geometry = self.artifacts.get(
                        geometry_handoff.geometry_artifact_id
                    )
                    if geometry is None:
                        raise ContractError(
                            "ORCA TS Hessian lacks its selected geometry"
                        )
                    identity = build_scientific_identity_binding(
                        task_spec_sha256=approval.task_spec_sha256,
                        geometry_artifact=geometry,
                        charge=consumer_binding.charge,
                        multiplicity=consumer_binding.multiplicity,
                    )
                    self.artifacts[artifact.artifact_id] = artifact
                    self.hessian_handoffs[edge.consumer_node_id] = observed
                elif context.proposal.program == "pyscf":
                    candidates = tuple(
                        item for item in outputs if item.kind == "pyscf_hdf5"
                    )
                    if len(candidates) != 1:
                        raise ContractError(
                            "validated PySCF OPT requires exactly one HDF5 result"
                        )
                    artifact, observed = handoff_optimized_pyscf_geometry(
                        producer_receipt=receipt,
                        result_artifact=candidates[0],
                        producer_edge=edge,
                        approved_workspace=self.approved_workspace,
                        geometry_artifact_id=(
                            f"geometry.{edge.producer_node_id}-to-"
                            f"{edge.consumer_node_id}"
                        ),
                        expected_charge=context.scientific_identity.charge,
                        expected_multiplicity=(
                            context.scientific_identity.multiplicity
                        ),
                        consumer_charge=consumer_binding.charge,
                        consumer_multiplicity=(consumer_binding.multiplicity),
                    )
                elif context.proposal.program == "xtb":
                    candidates = tuple(
                        item
                        for item in outputs
                        if item.kind == "geometry_xyz"
                        and Path(item.path).name == "xtbopt.xyz"
                    )
                    if len(candidates) != 1:
                        raise ContractError(
                            "validated xTB OPT requires exactly one xtbopt.xyz"
                        )
                    artifact, observed = handoff_optimized_xtb_geometry(
                        producer_receipt=receipt,
                        result_artifact=candidates[0],
                        input_artifact=context.input_artifact,
                        producer_edge=edge,
                        approved_workspace=self.approved_workspace,
                        geometry_artifact_id=(
                            f"geometry.{edge.producer_node_id}-to-"
                            f"{edge.consumer_node_id}"
                        ),
                        expected_charge=context.scientific_identity.charge,
                        expected_multiplicity=(
                            context.scientific_identity.multiplicity
                        ),
                        consumer_charge=consumer_binding.charge,
                        consumer_multiplicity=(consumer_binding.multiplicity),
                    )
                else:
                    output_kind = f"{context.proposal.program}_output"
                    candidates = tuple(
                        item for item in outputs if item.kind == output_kind
                    )
                    if len(candidates) != 1:
                        raise ContractError(
                            f"validated {context.proposal.program} OPT/TS "
                            f"requires exactly one {output_kind}"
                        )
                    artifact, observed = handoff_optimized_native_geometry(
                        program=context.proposal.program,
                        producer_receipt=receipt,
                        result_artifact=candidates[0],
                        input_artifact=context.input_artifact,
                        producer_edge=edge,
                        approved_workspace=self.approved_workspace,
                        geometry_artifact_id=(
                            f"geometry.{edge.producer_node_id}-to-"
                            f"{edge.consumer_node_id}"
                        ),
                        expected_charge=context.scientific_identity.charge,
                        expected_multiplicity=(
                            context.scientific_identity.multiplicity
                        ),
                        consumer_charge=consumer_binding.charge,
                        consumer_multiplicity=(consumer_binding.multiplicity),
                    )
                if edge.selection_rule == "validated_optimized_geometry":
                    consumer_charge, consumer_multiplicity = (
                        observed.consumer_state
                    )
                    identity = build_scientific_identity_binding(
                        task_spec_sha256=approval.task_spec_sha256,
                        geometry_artifact=artifact,
                        charge=consumer_charge,
                        multiplicity=consumer_multiplicity,
                    )
                    self.artifacts[artifact.artifact_id] = artifact
                    self.scientific_identities[identity.binding_sha256] = identity
                    self.handoffs[edge.consumer_node_id] = observed
                scientific_edge = next(
                    (
                        item
                        for item in scientific_plan.edges
                        if item.edge_kind == "data"
                        and item.source_node_id == edge.producer_node_id
                        and item.target_node_id == edge.consumer_node_id
                        and item.artifact_class == edge.artifact_kind
                    ),
                    None,
                )
                if scientific_edge is None:
                    raise ContractError(
                        "producer handoff lacks an exact scientific data edge"
                    )
                produced_handoffs.append(
                    {
                        "handoff": observed,
                        "artifact": artifact,
                        "scientific_identity": identity,
                    }
                )
                pending_data_edges.append(
                    (scientific_edge, edge, observed, identity)
                )
                if edge.selection_rule == "validated_optimized_geometry":
                    self._emit(
                        turn_id,
                        EventKind.OPTIMIZED_GEOMETRY_HANDED_OFF,
                        observed.receipt_sha256,
                        status=observed.status,
                        producer_node_id=edge.producer_node_id,
                        consumer_node_id=edge.consumer_node_id,
                    )
        if self.frozen_workflow_approval is not None:
            output_sha256s = tuple(
                sorted(
                    {
                        *(artifact.sha256 for artifact in outputs),
                        *(
                            item["artifact"].sha256
                            for item in produced_handoffs
                        ),
                    }
                )
            )
            if receipt.engine_complete:
                self.event_store.transition_workflow_run_node(
                    turn_id=turn_id,
                    run_id=v2_run_id,
                    node_id=node_id,
                    new_state="engine_complete",
                    execution_receipt_sha256=receipt.receipt_sha256,
                    output_artifact_sha256s=output_sha256s,
                    timestamp=finished,
                )
                if receipt.validated:
                    self.event_store.transition_workflow_run_node(
                        turn_id=turn_id,
                        run_id=v2_run_id,
                        node_id=node_id,
                        new_state="validated",
                        validator_receipt_sha256s=(
                            receipt.validator_receipt_sha256s
                        ),
                        result_validation_receipt=(result_validation_receipt),
                        timestamp=finished,
                    )
                    for (
                        scientific_edge,
                        producer_edge,
                        observed_handoff,
                        consumer_identity,
                    ) in pending_data_edges:
                        binding = build_validated_data_edge_binding(
                            run_id=v2_run_id,
                            plan=scientific_plan,
                            approval=frozen_approval,
                            scientific_edge=scientific_edge,
                            producer_edge=producer_edge,
                            producer_invocation=execution_invocation,
                            producer_receipt=receipt,
                            handoff=observed_handoff,
                            producer_scientific_identity_sha256=(
                                context.scientific_identity.binding_sha256
                            ),
                            consumer_scientific_identity_sha256=(
                                consumer_identity.binding_sha256
                            ),
                        )
                        self.event_store.record_validated_data_edge_binding(
                            turn_id=turn_id,
                            binding=binding,
                        )
                        for item in produced_handoffs:
                            if item["handoff"] == observed_handoff:
                                item["data_edge_binding"] = binding
                else:
                    self.event_store.transition_workflow_run_node(
                        turn_id=turn_id,
                        run_id=v2_run_id,
                        node_id=node_id,
                        new_state="failed",
                        execution_receipt_sha256=receipt.receipt_sha256,
                        output_artifact_sha256s=output_sha256s,
                        failure_rule_ids=(
                            receipt.findings
                            or ("execution.validation.failed",)
                        ),
                        timestamp=finished,
                    )
            else:
                self.event_store.transition_workflow_run_node(
                    turn_id=turn_id,
                    run_id=v2_run_id,
                    node_id=node_id,
                    new_state=receipt.execution_state,
                    execution_receipt_sha256=receipt.receipt_sha256,
                    output_artifact_sha256s=output_sha256s,
                    failure_rule_ids=(
                        receipt.findings
                        or ("execution.state." + receipt.execution_state,)
                    ),
                    timestamp=finished,
                )
        return {
            "execution": receipt,
            "idempotent_replay": False,
            "produced_handoffs": tuple(produced_handoffs),
            "process_observation": process_observation.as_dict(),
            "result_validation": result_validation_receipt,
        }

    def _require_bounded_launch_budget(self) -> float:
        """Return this launch's timeout while preserving analysis time.

        The envelope's node timeout is an upper bound, not a requirement that
        the full duration remain available before every node.  Planning and
        earlier nodes may consume part of the episode; a later launch receives
        the smaller live window left after reserving postprocessing time.
        """

        envelope = self.bounded_execution_envelope
        if envelope is None:
            return float(self.execution_resources.node_timeout_seconds)
        if len(self.execution_receipts) >= envelope.max_engine_calls:
            raise ContractError(
                "bounded execution engine-call budget exhausted"
            )
        remaining = envelope.episode_wall_time_seconds - (
            time.monotonic() - self._bounded_execution_started_at
        )
        execution_window = remaining - envelope.postprocess_reserve_seconds
        effective_timeout = min(
            float(self.execution_resources.node_timeout_seconds),
            execution_window,
        )
        if effective_timeout < 1.0:
            raise ContractError(
                "insufficient episode time for a usable node launch while "
                "protecting the postprocessing reserve"
            )
        return effective_timeout

    def _release_non_executable_node_ids(
        self, plan: ScientificWorkflowPlanV2
    ) -> frozenset[str]:
        """Plan stages this release cannot execute, whatever the plan asked for.

        Release maturity is a host fact, not a planning choice.  A scientific
        workflow is required to keep a stage it cannot materialize rather than
        silently drop it, so such a stage must not also block the human review
        of the stages that *can* run.  It is displayed with the workflow,
        excluded from the approval, and never launched.
        """

        non_executable: set[str] = set()
        for node in plan.nodes:
            capability = self.registry.get(node.program)
            executable = (
                set(capability.execution_engine_job_pairs)
                if capability is not None
                else set()
            )
            if (node.engine, node.stage) in executable:
                continue
            if node.support_state != "blocked_unsupported":
                continue
            non_executable.add(node.node_id)
        return frozenset(non_executable)

    def execution_review_ineligibility_reason(
        self,
        *,
        plan: ScientificWorkflowPlanV2,
        planned_node: ScientificWorkflowNodeV2,
    ) -> str:
        """Explain why one valid plan node cannot enter human execution review.

        A stage this release declares preview-only is deferred rather than
        ineligible: see ``_release_non_executable_node_ids``.
        """

        non_executable = self._release_non_executable_node_ids(plan)
        if planned_node.node_id in non_executable:
            if non_executable == {node.node_id for node in plan.nodes}:
                return (
                    "workflow has no release-executable stage to approve; "
                    "retain this stage as non-executable scientific intent"
                )
            return ""

        envelope = self.bounded_execution_envelope
        if envelope is None:
            return "an explicit execution envelope is required"
        if not envelope.allows(planned_node.program, planned_node.engine):
            return "program/engine is outside the execution envelope"
        program_capability = self.registry.get(planned_node.program)
        executable_pairs = (
            set(program_capability.execution_engine_job_pairs)
            if program_capability is not None
            else set()
        )
        if (planned_node.engine, planned_node.stage) not in executable_pairs:
            return (
                "job is supported for planning or preview, not Agent "
                "execution; keep the stage and declare it blocked_unsupported "
                "so the executable stages can still be reviewed"
            )
        if planned_node.program == "orca" and planned_node.stage == "ts":
            data_target_ids = {
                edge.target_node_id
                for edge in plan.edges
                if edge.edge_kind == "data"
            }
            context = self._bounded_node_context(
                plan=plan,
                planned_node=planned_node,
                data_target_ids=data_target_ids,
            )
            if context.project_validation is None:
                raise ContractError(
                    "ORCA transition-state review lacks validated project settings"
                )
            settings = dict(context.project_validation.settings)
            if not bool(settings.get("freq") or settings.get("numfreq")):
                return (
                    "ORCA transition-state execution requires a requested "
                    "frequency analysis to establish a first-order saddle"
                )
        return ""

    def build_execution_review(
        self,
        *,
        workspace: str | Path,
        request_id: str = "",
    ) -> WorkflowExecutionReviewV1:
        """Freeze review evidence without creating execution authority.

        The resource envelope is advisory input during planning.  This method
        is the only product path from that preview session: it emits an inert
        packet which a separate human-only command may approve later.  It
        deliberately does not populate either approval field on this host.
        """

        envelope = self.bounded_execution_envelope
        resources = self.execution_resources
        if envelope is None or resources is None:
            raise ContractError(
                "execution review requires an explicit resource envelope"
            )
        if resources.resource_sha256 != envelope.resources.resource_sha256:
            raise ContractError("review resources differ from execution envelope")
        plans = tuple(self.scientific_workflow_plans.values())
        if not plans:
            raise ContractError("execution review requires a scientific workflow")
        plan = plans[-1]
        if plan.task_spec_sha256 not in self.task_spec_sha256s:
            raise ContractError("review workflow belongs to another task")
        non_executable_ids = self._release_non_executable_node_ids(plan)
        executable_nodes = tuple(
            item
            for item in plan.nodes
            if item.node_id not in non_executable_ids
        )
        if not executable_nodes:
            raise ContractError(
                "execution review requires at least one executable node"
            )
        if len(executable_nodes) > envelope.max_engine_calls:
            raise ContractError(
                "scientific workflow exceeds bounded engine-call budget: "
                f"{len(executable_nodes)} nodes for "
                f"{envelope.max_engine_calls} calls"
            )
        data_edges = tuple(
            edge
            for edge in plan.edges
            if edge.edge_kind == "data"
            and edge.target_node_id not in non_executable_ids
        )
        if any(
            edge.source_node_id in non_executable_ids for edge in data_edges
        ):
            raise ContractError(
                "an executed node cannot consume the output of a stage this "
                "release cannot execute"
            )
        data_target_ids = {edge.target_node_id for edge in data_edges}
        unsupported = tuple(
            item.node_id
            for item in executable_nodes
            if item.support_state == "blocked_unsupported"
            or (
                item.support_state == "unresolved_future"
                and item.node_id not in data_target_ids
            )
            or item.support_state not in {"resolvable", "unresolved_future"}
        )
        if unsupported:
            raise ContractError(
                "execution review refuses unsupported or non-causal unresolved "
                "nodes: " + ", ".join(unsupported)
            )
        materialized = self._latest_bounded_materialization(plan)
        producer_edges = []
        for edge in data_edges:
            producer = next(
                item
                for item in plan.nodes
                if item.node_id == edge.source_node_id
            )
            if is_validated_optimized_geometry_edge(plan, edge):
                selection_rule = "validated_optimized_geometry"
            elif is_validated_orca_ts_hessian_edge(plan, edge):
                selection_rule = "validated_final_orca_ts_hessian"
            else:
                raise ContractError(
                    "execution review has no exact selection rule for data "
                    f"edge {edge.edge_id!r}; expected optimized geometry or "
                    "an ORCA final-TS Hessian for IRC"
                )
            if (
                selection_rule == "validated_optimized_geometry"
                and producer.program
                not in {"gaussian", "orca", "pyscf", "xtb"}
            ):
                raise ContractError(
                    "execution review has no optimized-geometry handoff for "
                    f"producer program {producer.program!r}"
                )
            producer_edges.append(
                build_producer_edge_rule(
                    producer_node_id=edge.source_node_id,
                    consumer_node_id=edge.target_node_id,
                    artifact_kind=edge.artifact_class,
                    selection_rule=selection_rule,
                )
            )
        geometry_edges = tuple(
            edge
            for edge in producer_edges
            if edge.selection_rule == "validated_optimized_geometry"
        )
        edge_by_target = {
            edge.consumer_node_id: edge
            for edge in geometry_edges
        }
        if (
            len(edge_by_target) != len(geometry_edges)
            or set(edge_by_target) != data_target_ids
        ):
            raise ContractError(
                "every producer-dependent calculation requires exactly one "
                "validated geometry input"
            )
        node_bindings = []
        environment_bindings = []
        node_reviews: list[WorkflowExecutionNodeReviewV1] = []
        unbindable: list[tuple[str, str]] = []
        for planned_node in executable_nodes:
            ineligibility = self.execution_review_ineligibility_reason(
                plan=plan,
                planned_node=planned_node,
            )
            if ineligibility:
                raise ContractError(
                    "workflow node is not eligible for Agent execution: "
                    f"{planned_node.program}/{planned_node.engine}/"
                    f"{planned_node.stage}: {ineligibility}"
                )
            try:
                context = self._bounded_node_context(
                    plan=plan,
                    planned_node=planned_node,
                    data_target_ids=data_target_ids,
                )
            except ContractError as exc:
                # The release can run this stage, but its project, capability
                # or environment evidence does not resolve to exactly one
                # record -- most often because the plan explored more than one
                # program and left ambiguous artifacts behind.  Collect every
                # such node before reporting, so one message names all of them
                # and their reasons instead of stopping at whichever happened
                # to be visited first.
                #
                # The host does not mark them non-executable on the model's
                # behalf: a plan states its own unsupported stages, its digest
                # covers that statement, and rewriting it here would forge the
                # provenance of a plan the model did not author.
                logger.warning(
                    "node %s cannot enter execution review: %s",
                    planned_node.node_id,
                    exc,
                )
                unbindable.append((planned_node.node_id, str(exc)))
                continue
            environment_receipt_sha256 = (
                context.engine_binding.environment_receipt_sha256
            )
            environment_identity = self._environment_identity_for(
                environment_receipt_sha256
            )
            if not environment_identity:
                raise ContractError(
                    f"node {planned_node.node_id!r} lacks environment identity"
                )
            environment_receipt = self.environments.get(
                environment_receipt_sha256
            )
            if environment_receipt is None:
                raise ContractError(
                    f"node {planned_node.node_id!r} lacks environment evidence"
                )
            environment_bindings.append(
                WorkflowEnvironmentBindingV1(
                    node_id=planned_node.node_id,
                    program=planned_node.program,
                    engine=planned_node.engine,
                    environment_receipt_sha256=environment_receipt_sha256,
                    environment_identity_sha256=environment_identity,
                )
            )
            edge = edge_by_target.get(planned_node.node_id)
            project = context.project_artifact
            validation = context.project_validation
            if (
                project is None
                or validation is None
                or validation.status != "valid"
            ):
                raise ContractError(
                    f"node {planned_node.node_id!r} lacks valid project evidence"
                )
            target_charge = (
                planned_node.charge
                if planned_node.charge is not None
                else context.scientific_identity.charge
            )
            target_multiplicity = (
                planned_node.multiplicity
                if planned_node.multiplicity is not None
                else context.scientific_identity.multiplicity
            )
            if edge is None and (target_charge, target_multiplicity) != (
                context.scientific_identity.charge,
                context.scientific_identity.multiplicity,
            ):
                raise ContractError(
                    f"initial node {planned_node.node_id!r} explicit state "
                    "differs from its task-bound molecular input"
                )
            if edge is None:
                invocation, invocation_context = self._plan_invocation_for_node(
                    plan=plan, node_id=planned_node.node_id
                )
                review_input_sha256 = context.input_artifact.sha256
                input_binding = (
                    context.input_artifact.cli_value,
                    ("molecular-input", review_input_sha256),
                )
                coordinate_identity = {
                    "kind": "exact-input-artifact",
                    "geometry_artifact_sha256": review_input_sha256,
                }
            else:
                invocation_context = context
                scientific_identity = replace(
                    context.scientific_identity,
                    charge=target_charge,
                    multiplicity=target_multiplicity,
                    binding_sha256=canonical_sha256(
                        {
                            "schema_version": (
                                "chemsmart.scientific-identity-binding.v1"
                            ),
                            "task_spec_sha256": (
                                context.scientific_identity.task_spec_sha256
                            ),
                            "geometry_artifact_id": (
                                context.scientific_identity.geometry_artifact_id
                            ),
                            "geometry_artifact_sha256": (
                                context.scientific_identity.geometry_artifact_sha256
                            ),
                            "charge": target_charge,
                            "multiplicity": target_multiplicity,
                        }
                    ),
                )
                proposal = replace(
                    context.proposal,
                    scientific_identity_sha256=(
                        scientific_identity.binding_sha256
                    ),
                    charge=target_charge,
                    multiplicity=target_multiplicity,
                )
                invocation = compile_command(
                    proposal,
                    capability=context.capability,
                    binding=context.engine_binding,
                    project=context.project_artifact,
                    project_validation=context.project_validation,
                    input_artifact=context.input_artifact,
                    scientific_identity=scientific_identity,
                    job_artifact_options=dict(context.job_artifact_options),
                    live_schema=self.live_schema,
                    server=self.preview_server,
                )
                review_input_sha256 = edge.edge_sha256
                input_binding = (
                    context.input_artifact.cli_value,
                    ("producer-geometry", review_input_sha256),
                )
                coordinate_identity = {
                    "kind": "validated-producer-output",
                    "producer_edge_sha256": edge.edge_sha256,
                    "selection_rule": edge.selection_rule,
                    "reference_geometry_sha256": context.input_artifact.sha256,
                }
            approved_identity = next(
                (
                    identity
                    for identity in self.approved_molecular_identities.values()
                    if identity.geometry_sha256
                    == context.input_artifact.sha256
                ),
                None,
            )
            if approved_identity is None:
                molecular_identity = {
                    "identity_evidence_status": "task-bound-geometry-only",
                    **_review_molecule_identity(context.input_artifact),
                    "coordinate_identity": coordinate_identity,
                    "input_binding_sha256": review_input_sha256,
                    "charge": target_charge,
                    "multiplicity": target_multiplicity,
                    "electronic_state": "charge-and-multiplicity-specified",
                    "scientific_identity_sha256": (
                        context.scientific_identity.binding_sha256
                        if edge is None
                        else "deferred-until-producer-output"
                    ),
                }
            else:
                molecular_identity = {
                    "identity_evidence_status": "approved-molecular-identity",
                    **approved_identity.public_record(),
                    "atom_count": len(approved_identity.atom_order),
                    "formula": _formula_from_symbols(
                        approved_identity.atom_order
                    ),
                    "coordinate_identity": coordinate_identity,
                    "input_binding_sha256": review_input_sha256,
                    "charge": target_charge,
                    "multiplicity": target_multiplicity,
                    "electronic_state": "charge-and-multiplicity-specified",
                    "scientific_identity_sha256": (
                        context.scientific_identity.binding_sha256
                        if edge is None
                        else "deferred-until-producer-output"
                    ),
                }
            server_profile_sha256 = execution_server_profile_sha256(
                resources=resources,
                scratch_root=envelope.scratch_root,
            )
            server_token = execution_path_placeholder(
                "server-profile", server_profile_sha256
            )
            real_argv = build_real_execution_argv(
                compiled_argv=invocation.argv,
                command_path=invocation.command_path,
                resources=resources,
                server=server_token,
            )
            path_bindings = {
                sys.executable: (
                    "controller-python",
                    file_sha256(Path(sys.executable).resolve()),
                ),
                project.cli_value: ("project-yaml", project.sha256),
                input_binding[0]: input_binding[1],
            }
            for auxiliary in invocation.auxiliary_input_bindings:
                artifact = dict(invocation_context.job_artifact_options).get(
                    auxiliary.parameter_name
                )
                if artifact is None:
                    raise ContractError(
                        "review invocation lacks an auxiliary artifact"
                    )
                path_bindings[artifact.cli_value] = (
                    "auxiliary-" + auxiliary.parameter_name,
                    auxiliary.artifact_sha256,
                )
            node_reviews.append(
                build_workflow_execution_node_review(
                    node_id=planned_node.node_id,
                    stage=planned_node.stage,
                    program=planned_node.program,
                    engine=planned_node.engine,
                    molecular_identity=molecular_identity,
                    project_artifact_sha256=project.sha256,
                    project_settings_sha256=validation.settings_sha256,
                    project_settings=validation.settings,
                    real_execution_argv=project_real_execution_argv(
                        real_argv, path_bindings=path_bindings
                    ),
                    execution_resources=resources,
                    environment_summary=environment_review_summary(
                        environment_receipt
                    ),
                    server_profile_sha256=server_profile_sha256,
                    environment_receipt_sha256=environment_receipt_sha256,
                    environment_identity_sha256=environment_identity,
                )
            )
            node_bindings.append(
                ApprovedNodeBindingV1(
                    node_id=planned_node.node_id,
                    program=planned_node.program,
                    engine=planned_node.engine,
                    jobtype=planned_node.stage,
                    project_artifact_sha256=project.sha256,
                    settings_sha256=validation.settings_sha256,
                    charge=target_charge,
                    multiplicity=target_multiplicity,
                    input_mode="producer" if edge is not None else "initial",
                    initial_artifact_id=(
                        "" if edge is not None else context.input_artifact.artifact_id
                    ),
                    initial_artifact_sha256=(
                        "" if edge is not None else context.input_artifact.sha256
                    ),
                    scientific_identity_sha256=(
                        ""
                        if edge is not None
                        else context.scientific_identity.binding_sha256
                    ),
                    producer_edge_sha256=(
                        edge.edge_sha256 if edge is not None else ""
                    ),
                    auxiliary_input_bindings=(
                        self._latest_invocation_for_node(
                            planned_node.node_id,
                            plan_sha256=plan.plan_sha256,
                        )[0].auxiliary_input_bindings
                        if planned_node.node_id not in data_target_ids
                        else ()
                    ),
                )
            )
        if unbindable:
            detail = "; ".join(
                f"{node_id}: {reason}" for node_id, reason in sorted(unbindable)
            )
            raise ContractError(
                "these workflow nodes have no unique project, capability and "
                f"environment evidence, so they cannot be reviewed: {detail}. "
                "Resolve the ambiguity, or declare the stage "
                "blocked_unsupported in the plan with its reason so the rest "
                "of the workflow can still be reviewed"
            )
        if not node_reviews:
            raise ContractError(
                "no workflow node could be bound to unique project, "
                "capability and environment evidence, so there is nothing to "
                "review"
            )
        review_request_id = (
            str(request_id).strip()
            or "review-" + plan.plan_sha256[:16]
        )
        request = build_workflow_approval_request(
            request_id=review_request_id,
            workflow_id=plan.workflow_id,
            workflow_sha256=plan.plan_sha256,
            task_spec_sha256=plan.task_spec_sha256,
            workspace=workspace,
            resources=resources,
            node_bindings=tuple(node_bindings),
            producer_edges=tuple(producer_edges),
        )
        return build_workflow_execution_review(
            request=request,
            scientific_plan=plan,
            materialized_workflow=materialized,
            execution_resources=resources,
            execution_envelope=canonical_data(envelope),
            environment_bindings=tuple(
                sorted(environment_bindings, key=lambda item: item.node_id)
            ),
            node_reviews=tuple(
                sorted(node_reviews, key=lambda item: item.node_id)
            ),
            stationary_point_policy=self.stationary_point_policy,
            non_executable_node_ids=tuple(sorted(non_executable_ids)),
        )

    def _admit_bounded_workflow(
        self, *, node_id: str, plan_sha256: str = ""
    ) -> None:
        """Freeze current ChemSmart evidence under the operating envelope.

        This is the deferred equivalent of loading an approval file.  It does
        not approve a model-authored command or path: every binding comes from
        the host's current plan, validated project, compiled invocation, green
        preview, observed environment, and exact scientific data edge.
        """

        envelope = self.bounded_execution_envelope
        if envelope is None:
            raise ContractError("workflow has no approval or bounded envelope")
        self._require_bounded_launch_budget()
        plans = tuple(
            plan
            for plan in self.scientific_workflow_plans.values()
            if any(item.node_id == node_id for item in plan.nodes)
            and (not plan_sha256 or plan.plan_sha256 == plan_sha256)
        )
        if not plans:
            raise ContractError(
                "bounded execution requires a current scientific workflow "
                "containing the requested node"
            )
        # Node IDs are workflow-local.  Choose the current exact plan first,
        # then require its own materialization.  Falling back to an older plan
        # with the same node name is what turned a missing identity/preparation
        # step in an edge-free diagnostic into an unrelated future-edge error.
        plan = plans[-1]
        frozen = self.frozen_workflow_approval
        if frozen is not None:
            if frozen.plan_sha256 == plan.plan_sha256:
                raise ContractError("bounded workflow is already admitted")
            raise ContractError(
                "a separate workflow already owns the immutable bounded "
                "approval; preserve that evidence and start a new independent "
                "attempt for this workflow"
            )
        current_materializations = tuple(
            workflow
            for workflow in self.materialized_workflows.values()
            if workflow.plan_sha256 == plan.plan_sha256
        )
        if not current_materializations:
            raise ContractError(
                "current bounded workflow has not been materialized; prepare "
                "and preflight this workflow's node after resolving its exact "
                "molecular identity and project evidence"
            )
        if plan.task_spec_sha256 not in self.task_spec_sha256s:
            raise ContractError("bounded workflow belongs to another task")
        if len(plan.nodes) > envelope.max_engine_calls:
            raise ContractError(
                "scientific workflow exceeds bounded engine-call budget: "
                f"{len(plan.nodes)} nodes for {envelope.max_engine_calls} calls"
            )
        data_edges = tuple(
            edge for edge in plan.edges if edge.edge_kind == "data"
        )
        data_target_ids = {edge.target_node_id for edge in data_edges}
        unsupported = tuple(
            item.node_id
            for item in plan.nodes
            if item.support_state == "blocked_unsupported"
            or (
                item.support_state == "unresolved_future"
                and item.node_id not in data_target_ids
            )
            or item.support_state not in {"resolvable", "unresolved_future"}
        )
        if unsupported:
            raise ContractError(
                "bounded execution refuses unsupported or non-causal unresolved "
                "nodes: " + ", ".join(unsupported)
            )
        materialized = self._latest_bounded_materialization(plan)
        node_bindings = []
        producer_edges = []
        for edge in data_edges:
            producer = next(
                item
                for item in plan.nodes
                if item.node_id == edge.source_node_id
            )
            if is_validated_optimized_geometry_edge(plan, edge):
                selection_rule = "validated_optimized_geometry"
            elif is_validated_orca_ts_hessian_edge(plan, edge):
                selection_rule = "validated_final_orca_ts_hessian"
            else:
                raise ContractError(
                    "bounded execution has no exact selection rule for data "
                    f"edge {edge.edge_id!r}; expected optimized geometry or "
                    "an ORCA final-TS Hessian for IRC"
                )
            if (
                selection_rule == "validated_optimized_geometry"
                and producer.program
                not in {"gaussian", "orca", "pyscf", "xtb"}
            ):
                raise ContractError(
                    "bounded execution has no optimized-geometry handoff for "
                    f"producer program {producer.program!r}"
                )
            producer_edges.append(
                build_producer_edge_rule(
                    producer_node_id=edge.source_node_id,
                    consumer_node_id=edge.target_node_id,
                    artifact_kind=edge.artifact_class,
                    selection_rule=selection_rule,
                )
            )
        geometry_edges = tuple(
            edge
            for edge in producer_edges
            if edge.selection_rule == "validated_optimized_geometry"
        )
        edge_by_target = {
            edge.consumer_node_id: edge
            for edge in geometry_edges
        }
        if (
            len(edge_by_target) != len(geometry_edges)
            or set(edge_by_target) != data_target_ids
        ):
            raise ContractError(
                "every producer-dependent calculation requires exactly one "
                "validated geometry input"
            )
        environment_identities = set()
        future_environments = {}
        for planned_node in plan.nodes:
            if not envelope.allows(planned_node.program, planned_node.engine):
                raise ContractError(
                    "workflow uses program/engine outside bounded allowlist: "
                    f"{planned_node.program}/{planned_node.engine}"
                )
            context = self._bounded_node_context(
                plan=plan,
                planned_node=planned_node,
                data_target_ids=data_target_ids,
            )
            environment_identity = self._environment_identity_for(
                context.engine_binding.environment_receipt_sha256
            )
            if not environment_identity:
                raise ContractError(
                    f"node {planned_node.node_id!r} lacks environment identity"
                )
            environment_identities.add(environment_identity)
            if planned_node.node_id in data_target_ids:
                future_environments[planned_node.node_id] = (
                    context.engine_binding.environment_receipt_sha256
                )
            edge = edge_by_target.get(planned_node.node_id)
            project = context.project_artifact
            validation = context.project_validation
            if (
                project is None
                or validation is None
                or validation.status != "valid"
            ):
                raise ContractError(
                    f"node {planned_node.node_id!r} lacks valid project evidence"
                )
            target_charge = (
                planned_node.charge
                if planned_node.charge is not None
                else context.scientific_identity.charge
            )
            target_multiplicity = (
                planned_node.multiplicity
                if planned_node.multiplicity is not None
                else context.scientific_identity.multiplicity
            )
            if edge is None and (target_charge, target_multiplicity) != (
                context.scientific_identity.charge,
                context.scientific_identity.multiplicity,
            ):
                raise ContractError(
                    f"initial node {planned_node.node_id!r} explicit state "
                    "differs from its task-bound molecular input"
                )
            node_bindings.append(
                ApprovedNodeBindingV1(
                    node_id=planned_node.node_id,
                    program=planned_node.program,
                    engine=planned_node.engine,
                    jobtype=planned_node.stage,
                    project_artifact_sha256=project.sha256,
                    settings_sha256=validation.settings_sha256,
                    charge=target_charge,
                    multiplicity=target_multiplicity,
                    input_mode="producer" if edge is not None else "initial",
                    initial_artifact_id=(
                        ""
                        if edge is not None
                        else context.input_artifact.artifact_id
                    ),
                    initial_artifact_sha256=(
                        ""
                        if edge is not None
                        else context.input_artifact.sha256
                    ),
                    scientific_identity_sha256=(
                        ""
                        if edge is not None
                        else context.scientific_identity.binding_sha256
                    ),
                    producer_edge_sha256=(
                        edge.edge_sha256 if edge is not None else ""
                    ),
                    auxiliary_input_bindings=(
                        self._latest_invocation_for_node(
                            planned_node.node_id,
                            plan_sha256=plan.plan_sha256,
                        )[0].auxiliary_input_bindings
                        if planned_node.node_id not in data_target_ids
                        else ()
                    ),
                )
            )
        approval_id = "bounded-" + plan.plan_sha256[:16]
        approval = build_workflow_execution_approval(
            approval_id=approval_id,
            workflow_id=plan.workflow_id,
            workflow_sha256=plan.plan_sha256,
            task_spec_sha256=plan.task_spec_sha256,
            approved_workspace=self.approved_workspace,
            resources=self.execution_resources,
            node_bindings=tuple(node_bindings),
            producer_edges=tuple(producer_edges),
        )
        receipt_identity_map = {
            receipt_sha256: self._environment_identity_for(receipt_sha256)
            for receipt_sha256 in {
                node.environment_receipt_sha256 for node in materialized.nodes
            }.union(future_environments.values())
        }
        frozen = build_frozen_workflow_approval(
            approval_id=approval_id,
            plan=plan,
            materialized_workflow=materialized,
            resources=self.execution_resources,
            environment_identity_sha256s=tuple(sorted(environment_identities)),
            future_node_environment_identity_sha256s=future_environments,
            environment_identity_by_receipt=receipt_identity_map,
            stationary_point_policy=self.stationary_point_policy,
        )
        self.workflow_execution_approval = approval
        self.frozen_workflow_approval = frozen
        self.materialized_workflows[materialized.materialized_sha256] = (
            materialized
        )

    def _latest_bounded_materialization(
        self, plan: ScientificWorkflowPlanV2
    ) -> MaterializedWorkflowV1:
        matches = tuple(
            workflow
            for workflow in self.materialized_workflows.values()
            if workflow.plan_sha256 == plan.plan_sha256
        )
        if not matches:
            raise ContractError("bounded workflow has not been materialized")
        materialized = matches[-1]
        previewed_ids = {
            item.node_id
            for item in materialized.nodes
            if item.state == "previewed"
        }
        non_executable_ids = self._release_non_executable_node_ids(plan)
        data_targets = {
            edge.target_node_id
            for edge in plan.edges
            if edge.edge_kind == "data"
            and edge.source_node_id not in non_executable_ids
            and edge.target_node_id not in non_executable_ids
        }
        initial_ids = (
            {item.node_id for item in plan.nodes}
            - data_targets
            - non_executable_ids
        )
        if not initial_ids or not initial_ids.issubset(previewed_ids):
            raise ContractError(
                "every initial workflow node requires a green preview before "
                "bounded execution"
            )
        unresolved_ids = (
            set(materialized.unresolved_node_ids) - non_executable_ids
        )
        if unresolved_ids != data_targets:
            raise ContractError(
                "only exact producer-dependent nodes may remain unresolved"
            )
        return materialized

    def _bounded_node_context(
        self,
        *,
        plan: ScientificWorkflowPlanV2,
        planned_node: ScientificWorkflowNodeV2,
        data_target_ids: set[str],
    ) -> _CommandContext:
        """Resolve current or future context without model-carried receipts."""

        if planned_node.node_id not in data_target_ids:
            _invocation, context = self._plan_invocation_for_node(
                plan=plan, node_id=planned_node.node_id
            )
            return context
        draft = next(
            (
                draft
                for draft in reversed(tuple(self.workflow_drafts.values()))
                if draft.workflow_id == plan.workflow_id
            ),
            None,
        )
        if draft is None:
            raise ContractError(
                "future bounded node lacks command workflow draft"
            )
        node = next(
            (
                item
                for item in draft.nodes
                if item.node_id == planned_node.node_id
            ),
            None,
        )
        producer_inputs = tuple(
            item for item in (node.inputs if node is not None else ())
            if item.producer_node_id
        )
        geometry_inputs = tuple(
            item
            for item in producer_inputs
            if item.artifact_class == "geometry_xyz"
        )
        auxiliary_inputs = tuple(
            item for item in producer_inputs if item not in geometry_inputs
        )
        valid_orca_irc_auxiliary = bool(
            planned_node.program == "orca"
            and planned_node.stage == "irc"
            and len(auxiliary_inputs) == 1
            and len(geometry_inputs) == 1
            and geometry_inputs[0].binding_id == "filename"
            and auxiliary_inputs[0].binding_id == "hess_filename"
            and auxiliary_inputs[0].artifact_class == "orca_hessian"
        )
        if (
            node is None
            or len(geometry_inputs) != 1
            or (auxiliary_inputs and not valid_orca_irc_auxiliary)
        ):
            raise ContractError(
                "future bounded node requires one filename/geometry_xyz input; "
                "ORCA IRC may additionally consume one "
                "hess_filename/orca_hessian input"
            )
        project = self._artifact(node.project_role)
        capabilities = tuple(
            receipt
            for receipt in self.capabilities.values()
            if receipt.query.program == planned_node.program
            and receipt.query.jobtype == planned_node.stage
            and receipt.query.engine == planned_node.engine
            and str(receipt.status.value) in {"supported", "preview_only"}
        )
        if len(capabilities) != 1:
            raise ContractError(
                f"future node {planned_node.node_id!r} lacks one exact "
                "program/stage/engine capability"
            )
        capability = capabilities[0]
        bindings = tuple(
            binding
            for binding in self.engine_bindings.values()
            if binding.capability_receipt_sha256 == capability.receipt_sha256
            and binding.program == planned_node.program
            and binding.selected_engine == planned_node.engine
            # Provider planning deliberately omits the execution tool, so a
            # correctly observed program environment is resolved but not
            # execution-ready in this host.  Review freezes that observed
            # environment; the provider-free executor rechecks readiness
            # before any launch.
            and binding.state == "resolved"
            and bool(binding.environment_receipt_sha256)
        )
        validation = self._resolve_project_validation(
            project=project,
            capability=capability,
            program=planned_node.program,
            jobtype=planned_node.stage,
        )
        if validation is None or len(bindings) != 1:
            raise ContractError(
                f"future node {planned_node.node_id!r} lacks unique project/environment evidence"
            )
        producer = next(
            (
                edge.source_node_id
                for edge in plan.edges
                if edge.edge_kind == "data"
                and edge.target_node_id == planned_node.node_id
                and edge.artifact_class == "geometry_xyz"
                and edge.consumer_input_id == "filename"
            ),
            None,
        )
        if producer is None:
            # This node was classified as consuming a producer's geometry, but
            # the plan carries no such edge for it.  Say so: a bare
            # StopIteration escapes every caller that expects a contract
            # failure and leaves the operator with no statement of what is
            # wrong.
            raise ContractError(
                f"future node {planned_node.node_id!r} consumes a produced "
                "geometry but the plan declares no geometry_xyz data edge "
                "into its 'filename' input"
            )
        _producer_invocation, producer_context = (
            self._plan_invocation_for_node(plan=plan, node_id=producer)
        )
        return _CommandContext(
            proposal=CommandProposalV1(
                node_id=planned_node.node_id,
                execution_target="run",
                program=planned_node.program,
                jobtype=planned_node.stage,
                project_artifact_id=project.artifact_id,
                input_artifact_id=producer_context.input_artifact.artifact_id,
                scientific_identity_sha256=(
                    producer_context.scientific_identity.binding_sha256
                ),
                charge=producer_context.scientific_identity.charge,
                multiplicity=producer_context.scientific_identity.multiplicity,
            ),
            capability=capability,
            program_binding=self.program_bindings[
                bindings[0].program_binding_sha256
            ],
            engine_binding=bindings[0],
            project_artifact=project,
            project_validation=validation,
            input_artifact=producer_context.input_artifact,
            scientific_identity=producer_context.scientific_identity,
        )

    def _execution_scientific_plan(self) -> ScientificWorkflowPlanV2:
        """Resolve the exact execution DAG without inferring tuple order."""

        approval = self.workflow_execution_approval
        if approval is None:
            raise ContractError("workflow has not been admitted for execution")
        frozen_approval = getattr(self, "frozen_workflow_approval", None)
        plans = getattr(self, "scientific_workflow_plans", None)
        if plans is None:
            plans = {}
            self.scientific_workflow_plans = plans
        if frozen_approval is not None:
            plan = plans.get(frozen_approval.plan_sha256)
            if plan is None:
                raise ContractError(
                    "frozen workflow approval has no registered scientific plan"
                )
            return plan
        plan = _scientific_plan_from_v1_approval(approval)
        plans.setdefault(plan.plan_sha256, plan)
        return plan

    def _current_execution_plan_for_node(
        self, node_id: str
    ) -> ScientificWorkflowPlanV2:
        """Resolve the current workflow before considering same-name history."""

        plans = tuple(
            plan
            for plan in getattr(self, "scientific_workflow_plans", {}).values()
            if any(node.node_id == node_id for node in plan.nodes)
        )
        for draft in reversed(
            tuple(getattr(self, "workflow_drafts", {}).values())
        ):
            if not any(node.node_id == node_id for node in draft.nodes):
                continue
            current = tuple(
                plan for plan in plans if plan.workflow_id == draft.workflow_id
            )
            if not current:
                raise ContractError(
                    "current workflow node has no task-bound scientific "
                    "identity; prepare it before requesting execution"
                )
            return current[-1]
        if not plans:
            raise ContractError(
                "execution requires a current scientific workflow containing "
                "the requested node"
            )
        return plans[-1]

    def _latest_invocation_for_node(
        self, node_id: str, *, plan_sha256: str = ""
    ) -> tuple[CanonicalCommandInvocationV1, _CommandContext]:
        invocations = getattr(self, "invocations", {})
        contexts = getattr(self, "_command_contexts", {})
        plan_bindings = getattr(
            self, "_invocation_workflow_plan_sha256s", {}
        )
        for invocation in reversed(tuple(invocations.values())):
            context = contexts[invocation.invocation_sha256]
            if (
                context.proposal.node_id == node_id
                and (
                    not plan_sha256
                    or plan_bindings.get(invocation.invocation_sha256)
                    == plan_sha256
                )
            ):
                return invocation, context
        message = "node has no compiled command invocation"
        if plan_sha256:
            message += " for the selected scientific workflow"
        raise ContractError(message)

    def _plan_invocation_for_node(
        self,
        *,
        plan: ScientificWorkflowPlanV2,
        node_id: str,
    ) -> tuple[CanonicalCommandInvocationV1, _CommandContext]:
        """Resolve and bind one invocation without cross-plan fallback."""

        try:
            return self._latest_invocation_for_node(
                node_id, plan_sha256=plan.plan_sha256
            )
        except ContractError:
            invocation, context = self._latest_invocation_for_node(node_id)
        bindings = self._invocation_workflow_plan_sha256s
        observed = bindings.get(invocation.invocation_sha256, "")
        if observed and observed != plan.plan_sha256:
            raise ContractError(
                "node invocation belongs to another scientific workflow"
            )
        if not any(node.node_id == node_id for node in plan.nodes):
            raise ContractError("scientific workflow has no such node")
        bindings[invocation.invocation_sha256] = plan.plan_sha256
        return invocation, context

    def _invocation_has_green_preflight(self, invocation_sha256: str) -> bool:
        return any(
            invocation_sha256
            in self._completion_sets.get(receipt.receipt_sha256, ())
            and receipt.plan_state == "previewed"
            and receipt.execution_ready
            and not receipt.critical_finding_sha256s
            for receipt in self.preflights.values()
        )

    def _real_execution_argv(
        self, invocation: CanonicalCommandInvocationV1
    ) -> tuple[str, ...]:
        return build_real_execution_argv(
            compiled_argv=invocation.argv,
            command_path=invocation.command_path,
            resources=self.execution_resources,
            server=self.execution_server,
        )

    def verify_reviewed_real_execution_argv(
        self,
        *,
        node_id: str,
        invocation_sha256: str,
        review: WorkflowExecutionNodeReviewV1,
    ) -> tuple[str, ...]:
        """Recheck the reviewed scientific meaning before a real launch.

        The live ChemSmart compiler is authoritative for the launch command.
        Human approval covers molecular state, program/stage, effective project
        settings, resources and the displayed causal DAG; an argv digest is not
        a second execution authority.
        """

        if review.node_id != node_id:
            raise ContractError("execution node differs from reviewed command")
        invocation = self._get(
            self.invocations, invocation_sha256, "canonical invocation"
        )
        context = self._get(
            self._command_contexts,
            invocation.invocation_sha256,
            "command context",
        )
        if (review.program, review.engine, review.stage) != (
            context.proposal.program,
            context.engine_binding.engine,
            context.proposal.jobtype,
        ):
            raise ContractError("program command differs from human review")
        if self.execution_resources is None or canonical_data(
            self.execution_resources
        ) != review.execution_resources:
            raise ContractError("execution resources differ from human review")
        if context.project_artifact is None or context.project_validation is None:
            raise ContractError("reviewed execution requires a validated project")
        reviewed_settings = json.loads(review.project_settings_text)
        effective_settings = dict(context.project_validation.settings)
        if canonical_data(effective_settings) != canonical_data(reviewed_settings):
            raise ContractError("effective project settings differ from human review")
        current_environment = self.environments.get(
            context.engine_binding.environment_receipt_sha256
        )
        if (
            current_environment is None
            or environment_review_summary(current_environment)
            != review.environment_summary
        ):
            raise ContractError("execution environment facts differ from human review")
        molecular = review.molecular_identity
        if (
            context.scientific_identity.charge != molecular.get("charge")
            or context.scientific_identity.multiplicity
            != molecular.get("multiplicity")
        ):
            raise ContractError("molecular electronic state differs from human review")
        coordinate = molecular.get("coordinate_identity")
        if not isinstance(coordinate, Mapping):
            raise ContractError("human review lacks coordinate identity")
        path_bindings = {
            sys.executable: (
                "controller-python",
                file_sha256(Path(sys.executable).resolve()),
            ),
            context.project_artifact.cli_value: (
                "project-yaml",
                context.project_artifact.sha256,
            )
        }
        if coordinate.get("kind") == "exact-input-artifact":
            expected_input = str(coordinate.get("geometry_artifact_sha256", ""))
            if context.input_artifact.sha256 != expected_input:
                raise ContractError("molecular input bytes differ from human review")
            if (
                molecular.get("scientific_identity_sha256")
                != context.scientific_identity.binding_sha256
            ):
                raise ContractError("molecular identity differs from human review")
            input_role = "molecular-input"
            input_digest = expected_input
        elif coordinate.get("kind") == "validated-producer-output":
            binding = self.workflow_execution_approval.node(node_id)
            input_digest = str(coordinate.get("producer_edge_sha256", ""))
            if (
                binding.input_mode != "producer"
                or binding.producer_edge_sha256 != input_digest
            ):
                raise ContractError("producer geometry edge differs from human review")
            if self.handoffs.get(node_id) is None:
                raise ContractError("producer geometry lacks validated handoff")
            input_role = "producer-geometry"
        else:
            raise ContractError("unknown reviewed coordinate identity")
        path_bindings[context.input_artifact.cli_value] = (
            input_role,
            input_digest,
        )
        if self.execution_server:
            server_path = Path(self.execution_server)
            if (
                not server_path.is_file()
                or server_path.is_symlink()
            ):
                raise ContractError("execution server profile is unavailable")
            path_bindings[self.execution_server] = (
                "server-profile",
                review.server_profile_sha256,
            )
        auxiliary_by_name = dict(context.job_artifact_options)
        for binding in invocation.auxiliary_input_bindings:
            artifact = auxiliary_by_name.get(binding.parameter_name)
            if artifact is None or artifact.sha256 != binding.artifact_sha256:
                raise ContractError("auxiliary input differs from human review")
            path_bindings[artifact.cli_value] = (
                "auxiliary-" + binding.parameter_name,
                binding.artifact_sha256,
            )
        projected = project_real_execution_argv(
            self._real_execution_argv(invocation),
            path_bindings=path_bindings,
        )
        if projected != review.real_execution_argv:
            raise ContractError(
                "recompiled ChemSmart CLI operation differs from human review"
            )
        return projected

    def _execution_output_artifacts(
        self,
        node_id: str,
        workspace: Path,
        *,
        program: str = "",
        auxiliary_inputs: Sequence[TrustedArtifactRefV1] = (),
    ) -> tuple[TrustedArtifactRefV1, ...]:
        auxiliary_by_basename = {
            (Path(item.path).name, item.size_bytes, item.sha256)
            for item in auxiliary_inputs
        }
        artifacts = []
        ordinal = 0
        for path in sorted(workspace.rglob("*")):
            if path.is_symlink():
                raise ContractError("execution emitted a symbolic link")
            if not path.is_file() or path.name.startswith("controller."):
                continue
            size_bytes = path.stat().st_size
            if any(
                basename == path.name
                and expected_size == size_bytes
                and expected_sha256 == file_sha256(path)
                for basename, expected_size, expected_sha256 in (
                    auxiliary_by_basename
                )
            ):
                # A staged multi-file input is launch evidence, not a newly
                # produced calculation result.  It is checked separately by
                # _staged_auxiliary_input_findings below.
                continue
            ordinal += 1
            suffix = path.suffix.lower()
            kind = "program_output"
            if suffix == ".h5":
                kind = "pyscf_hdf5"
            elif suffix == ".xyz":
                kind = "geometry_xyz"
            elif program == "orca" and suffix == ".hess":
                kind = "orca_hessian"
            elif suffix == ".json":
                kind = "json"
            elif program == "xtb" and suffix == ".out":
                kind = "xtb_output"
            elif program == "orca" and suffix == ".out":
                kind = "orca_output"
            elif program == "gaussian" and suffix in {".log", ".out"}:
                kind = "gaussian_output"
            before = path.stat()
            observed_sha256 = file_sha256(path)
            after = path.stat()
            if (
                before.st_size != after.st_size
                or before.st_mtime_ns != after.st_mtime_ns
            ):
                raise ContractError(
                    "execution output changed while it was being bound"
                )
            artifact = TrustedArtifactRefV1(
                artifact_id=f"result.{node_id}.{ordinal}",
                kind=kind,
                sha256=observed_sha256,
                size_bytes=after.st_size,
                path=str(path.resolve()),
                cli_value=str(path.resolve()),
            )
            artifacts.append(artifact)
            self.artifacts[artifact.artifact_id] = artifact
        return tuple(artifacts)

    @staticmethod
    def _evaluate_execution_outputs(
        *,
        program: str,
        jobtype: str,
        charge: int,
        multiplicity: int,
        expected_settings: Mapping[str, Any] | None = None,
        expected_input_artifact: TrustedArtifactRefV1 | None = None,
        expected_project_artifact: TrustedArtifactRefV1 | None = None,
        output_artifacts: tuple[TrustedArtifactRefV1, ...],
        exit_status: int | None,
        expected_environment_receipt_sha256: str = "",
        capability_environment_receipt: (
            EnvironmentCapabilityReceiptV1 | None
        ) = None,
        pyscf_engine_observation: _PySCFEngineObservation | None = None,
        stationary_point_policy: (
            StationaryPointValidationPolicyV1 | None
        ) = None,
        approved_stationary_point_policy_sha256: str = "",
        approved_hessian_node_id: str = "",
        process_observation: ProcessObservationV1 | None = None,
    ) -> _ExecutionValidationEvaluation:
        findings: list[str] = []
        observation: dict[str, Any] = {
            "program": program,
            "jobtype": jobtype,
            "wrapper_exit_status": exit_status,
        }
        if process_observation is not None:
            observation["process_observation"] = canonical_data(
                process_observation.as_dict()
            )
            findings.extend(_process_observation_findings(process_observation))
        for artifact in output_artifacts:
            try:
                _current_artifact_path(
                    artifact, field_name="execution output artifact"
                )
            except ContractError:
                findings.append("execution.output.artifact_binding_mismatch")
        if program == "pyscf":
            engine = pyscf_engine_observation or (
                _inspect_pyscf_engine_observation(
                    output_artifacts, launch_ambiguous=exit_status is None
                )
            )
            findings.extend(engine.findings)
            run_receipt = engine.run_receipt
            if capability_environment_receipt is not None:
                if (
                    expected_environment_receipt_sha256
                    != capability_environment_receipt.receipt_sha256
                ):
                    findings.append(
                        "pyscf.environment.capability_binding_mismatch"
                    )
                if capability_environment_receipt.program != program:
                    findings.append(
                        "pyscf.environment.program_binding_mismatch"
                    )
            environment_observation, environment_findings = (
                _pyscf_environment_evidence(
                    output_artifacts=output_artifacts,
                    run_receipt=run_receipt,
                    capability_environment=capability_environment_receipt,
                )
                if capability_environment_receipt is not None
                else ({}, ())
            )
            findings.extend(environment_findings)
            observation.update(
                {
                    "child_exit_status": engine.child_exit_status,
                    "engine_complete": engine.engine_complete,
                    "engine_receipt_sha256": engine.run_receipt_sha256,
                    "runner_scientifically_validated": (
                        run_receipt.get("scientifically_validated")
                        if run_receipt is not None
                        else None
                    ),
                    "environment_validation": environment_observation,
                }
            )
            if run_receipt is not None:
                runner_findings = tuple(
                    sorted(
                        {
                            str(item.get("rule_id") or "unknown")
                            for item in run_receipt.get("findings") or ()
                            if isinstance(item, Mapping)
                        }
                    )
                )
                observation["runner_findings"] = runner_findings
                observation["runner_state"] = run_receipt.get("state")
                deferred_hessian_classification = (
                    _runner_defers_hessian_classification(
                        run_receipt=run_receipt,
                        jobtype=jobtype,
                        hessian_node_id=approved_hessian_node_id,
                        engine_complete=engine.engine_complete,
                        stationary_point_policy=stationary_point_policy,
                        approved_stationary_point_policy_sha256=(
                            approved_stationary_point_policy_sha256
                        ),
                    )
                )
                if deferred_hessian_classification:
                    observation["runner_validation_delegation"] = (
                        "approved_stationary_point_policy"
                        if stationary_point_policy is not None
                        else "downstream_scientific_analysis"
                    )
                else:
                    observation["runner_validation_delegation"] = "none"
                if not deferred_hessian_classification:
                    if run_receipt.get("scientifically_validated") is not True:
                        findings.append(
                            "pyscf.run_receipt.scientific_validation_failed"
                        )
                    elif runner_findings:
                        findings.append(
                            "pyscf.run_receipt.validation_state_inconsistent"
                        )
                    if run_receipt.get("state") != "validated":
                        findings.append(
                            "pyscf.run_receipt.state_not_validated"
                        )
                if run_receipt.get("state") == "validated" and (
                    run_receipt.get("scientifically_validated") is not True
                    or runner_findings
                ):
                    findings.append(
                        "pyscf.run_receipt.validation_state_inconsistent"
                    )
                if (
                    expected_project_artifact is not None
                    and run_receipt.get("project_yaml_sha256")
                    != expected_project_artifact.sha256
                ):
                    findings.append(
                        "pyscf.run_receipt.project_digest_mismatch"
                    )
                if expected_input_artifact is not None:
                    expected_geometry_sha256 = _pyscf_input_geometry_sha256(
                        expected_input_artifact,
                        charge=charge,
                        multiplicity=multiplicity,
                    )
                    if (
                        not expected_geometry_sha256
                        or run_receipt.get("input_geometry_sha256")
                        != expected_geometry_sha256
                    ):
                        findings.append(
                            "pyscf.run_receipt.input_geometry_digest_mismatch"
                        )
                    observed_artifact_sha256 = run_receipt.get(
                        "input_artifact_sha256"
                    )
                    if (
                        observed_artifact_sha256
                        and observed_artifact_sha256
                        != expected_input_artifact.sha256
                    ):
                        findings.append(
                            "pyscf.run_receipt.input_digest_mismatch"
                        )

            result = engine.result_artifact
            if result is not None:
                try:
                    from chemsmart.jobs.pyscf.validation import (
                        validate_pyscf_result,
                    )

                    expected_symbols, expected_positions = (
                        _pyscf_input_geometry(expected_input_artifact)
                    )
                    expected_receipt = _pyscf_result_receipt_expectation(
                        run_receipt
                    )
                    result_validation = validate_pyscf_result(
                        _current_artifact_path(
                            result, field_name="PySCF HDF5 result"
                        ),
                        settings=expected_settings or {},
                        expected_jobtype=jobtype,
                        expected_charge=charge,
                        expected_multiplicity=multiplicity,
                        expected_symbols=expected_symbols,
                        expected_positions=(
                            expected_positions
                            if jobtype in {"sp", "hess"}
                            else None
                        ),
                        expected_receipt=expected_receipt,
                        stationary_point_policy=stationary_point_policy,
                    )
                    observation["result_validation"] = canonical_data(
                        result_validation
                    )
                    if stationary_point_policy is not None:
                        observation["stationary_point_policy_sha256"] = (
                            stationary_point_policy.policy_sha256
                        )
                    findings.extend(
                        str(item.rule_id)
                        for item in result_validation["findings"]
                    )
                    result_state = result_validation.get("state")
                    unclassified_for_downstream = bool(
                        deferred_hessian_classification
                        and result_state == "unclassified"
                    )
                    if (
                        result_state != "validated"
                        and not unclassified_for_downstream
                        and not result_validation.get("findings")
                    ):
                        findings.append(
                            "pyscf.result.validation_state_inconsistent"
                        )
                    from chemsmart.io.pyscf.output import read_pyscf_h5

                    result_spec, _provenance, _status, _results = (
                        read_pyscf_h5(
                            _current_artifact_path(
                                result, field_name="PySCF HDF5 program binding"
                            )
                        )
                    )
                    expected_engine = (
                        capability_environment_receipt.engine
                        if capability_environment_receipt is not None
                        else str((expected_settings or {}).get("engine") or "")
                    )
                    if result_spec.get("program") != "pyscf":
                        findings.append(
                            "pyscf.result.program_binding_mismatch"
                        )
                    if expected_engine and result_spec.get("engine") != (
                        expected_engine
                    ):
                        findings.append("pyscf.result.engine_binding_mismatch")
                except Exception as exc:
                    observation["pyscf_result_error_type"] = type(exc).__name__
                    findings.append("pyscf.result.unreadable")
            elif (
                run_receipt is not None
                and observation.get("runner_validation_delegation")
                == "approved_stationary_point_policy"
            ):
                findings.append("pyscf.result.policy_validation_unavailable")
        elif exit_status != 0 and program not in {"orca", "gaussian"}:
            findings.append("execution.process.nonzero_or_unknown")
        elif program == "orca":
            if exit_status != 0:
                findings.append("execution.process.nonzero_or_unknown")
            candidates = tuple(
                artifact
                for artifact in output_artifacts
                if artifact.kind == "orca_output"
            )
            orca_observation: dict[str, Any] = {
                "output_count": len(candidates),
                "normal_termination": False,
                "optimization_converged": False,
                "charge": None,
                "multiplicity": None,
                "energy_hartree": None,
                "vibrational_mode_count": 0,
                "transition_count": 0,
                "consequential_imaginary_mode_count": 0,
            }
            if len(candidates) != 1:
                findings.append("orca.result.output_count")
            else:
                try:
                    from chemsmart.io.orca import (
                        normalize_orca_neb_joboption,
                    )
                    from chemsmart.io.orca.output import (
                        ORCANEBOutput,
                        ORCAOutput,
                    )

                    output_class = (
                        ORCANEBOutput if jobtype == "neb" else ORCAOutput
                    )
                    output = output_class(
                        str(
                            _current_artifact_path(
                                candidates[0], field_name="ORCA result output"
                            )
                        )
                    )
                    if not output.normal_termination:
                        from chemsmart.io.native_failure import (
                            summarize_orca_native_failure,
                        )

                        diagnostic_artifacts = _native_diagnostic_artifacts(
                            output_artifacts
                        )
                        failure_summary = summarize_orca_native_failure(
                            output.contents,
                            diagnostic_lines=_iter_native_diagnostic_lines(
                                diagnostic_artifacts
                            ),
                        )
                        if failure_summary is not None:
                            orca_observation["native_failure"] = (
                                _bound_native_failure_summary(
                                    failure_summary,
                                    artifacts=(
                                        candidates[0],
                                        *diagnostic_artifacts,
                                    ),
                                )
                            )
                            findings.append(
                                "orca.native_failure."
                                + failure_summary.error_class
                            )
                    frequencies = tuple(output.vibrational_frequencies or ())
                    transitions = tuple(output.excited_state_records or ())
                    finite_frequencies = bool(frequencies) and all(
                        math.isfinite(float(value)) for value in frequencies
                    )
                    consequential_imaginary_modes = tuple(
                        float(value)
                        for value in frequencies
                        if math.isfinite(float(value))
                        and float(value) < -20.0
                    )
                    orca_observation.update(
                        {
                            "normal_termination": bool(
                                output.normal_termination
                            ),
                            "optimization_converged": bool(output.converged),
                            "charge": output.charge,
                            "multiplicity": output.multiplicity,
                            "energy_hartree": output.final_energy,
                            "vibrational_mode_count": len(frequencies),
                            "transition_count": len(transitions),
                            "consequential_imaginary_mode_count": len(
                                consequential_imaginary_modes
                            ),
                        }
                    )
                    if jobtype == "neb":
                        expected_joboption = normalize_orca_neb_joboption(
                            (expected_settings or {}).get("joboption")
                        )
                        observed_joboption = output.route_object.neb_joboption
                        neb_converged = bool(output.neb_converged)
                        ts_converged = bool(output.ts_converged)
                        ts_required = bool(
                            (expected_joboption or observed_joboption or "")
                            .upper()
                            .endswith("-TS")
                        )
                        orca_observation.update(
                            {
                                "neb_joboption": observed_joboption,
                                "expected_neb_joboption": expected_joboption,
                                "neb_converged": neb_converged,
                                "ts_converged": ts_converged,
                                "optimization_converged": bool(
                                    neb_converged
                                    and (not ts_required or ts_converged)
                                ),
                            }
                        )
                        if expected_joboption is None:
                            findings.append(
                                "orca.result.neb_joboption_missing"
                            )
                        elif observed_joboption != expected_joboption:
                            findings.append(
                                "orca.result.neb_joboption_mismatch"
                            )
                        if not neb_converged:
                            findings.append("orca.result.neb_not_converged")
                        if ts_required and not ts_converged:
                            findings.append("orca.result.neb_ts_not_converged")
                    if not output.normal_termination:
                        findings.append("orca.result.normal_termination")
                    if output.charge != charge:
                        findings.append("orca.result.charge_mismatch")
                    if output.multiplicity != multiplicity:
                        findings.append("orca.result.multiplicity_mismatch")
                    if (
                        jobtype in {"opt", "ts"}
                        and output.converged is not True
                    ):
                        findings.append(
                            "orca.result.optimization_not_converged"
                        )
                    if output.final_energy is None:
                        findings.append("orca.result.energy_missing")
                    requested_frequency_analysis = any(
                        bool((expected_settings or {}).get(field))
                        for field in ("freq", "numfreq", "vpt2")
                    )
                    if requested_frequency_analysis and not frequencies:
                        findings.append("orca.result.frequencies_missing")
                    elif requested_frequency_analysis and not finite_frequencies:
                        findings.append("orca.result.frequencies_invalid")
                    if jobtype == "ts":
                        if not finite_frequencies:
                            if "orca.result.frequencies_invalid" not in findings:
                                findings.append("orca.result.frequencies_invalid")
                        elif len(consequential_imaginary_modes) != 1:
                            findings.append(
                                "orca.result.ts_imaginary_mode_count"
                            )
                    if jobtype == "td":
                        requested_states = int(
                            (expected_settings or {}).get("nstates") or 1
                        )
                        if len(transitions) < requested_states:
                            findings.append("orca.result.transitions_missing")
                except Exception as exc:
                    orca_observation["parser_error_type"] = type(exc).__name__
                    findings.append("orca.result.unreadable")
            observation["orca"] = orca_observation
        elif program == "xtb":
            receipts: list[Path] = []
            for artifact in output_artifacts:
                if (
                    artifact.kind != "json"
                    or "xtb-result-receipt" not in Path(artifact.path).name
                ):
                    continue
                try:
                    path = _current_artifact_path(
                        artifact, field_name="xTB result receipt"
                    )
                    json.loads(path.read_text(encoding="utf-8"))
                    receipts.append(path)
                except (ContractError, OSError, json.JSONDecodeError):
                    pass
            if len(receipts) != 1:
                findings.append("xtb.result.receipt_count")
            else:
                from chemsmart.jobs.xtb.validation import (
                    audit_xtb_result_receipt,
                )

                xtb_observation, xtb_findings = audit_xtb_result_receipt(
                    receipts[0],
                    expected_jobtype=jobtype,
                    expected_charge=charge,
                    expected_multiplicity=multiplicity,
                    expected_settings=expected_settings,
                    expected_source_sha256=(
                        expected_input_artifact.sha256
                        if expected_input_artifact is not None
                        else ""
                    ),
                    expected_project_sha256=(
                        expected_project_artifact.sha256
                        if expected_project_artifact is not None
                        else ""
                    ),
                )
                observation["xtb"] = xtb_observation
                findings.extend(xtb_findings)
        elif program == "gaussian":
            if exit_status != 0:
                findings.append("execution.process.nonzero_or_unknown")
            candidates = tuple(
                artifact
                for artifact in output_artifacts
                if artifact.kind == "gaussian_output"
            )
            requested = expected_settings or {}
            expected_result_jobtype = (
                str(requested.get("jobtype") or "").strip().lower()
                if jobtype == "link"
                else jobtype
            )
            requested_direction = str(requested.get("direction") or "").lower()
            if jobtype == "irc":
                expected_directions = (
                    {"ircf"}
                    if requested_direction == "forward"
                    else (
                        {"ircr"}
                        if requested_direction == "reverse"
                        else {"ircf", "ircr"}
                    )
                )
            else:
                expected_directions = set()
            required_count = (
                len(expected_directions) if expected_directions else 1
            )
            gaussian_observation: dict[str, Any] = {
                "output_count": len(candidates),
                "required_output_count": required_count,
                "outputs": (),
            }
            if len(candidates) != required_count:
                findings.append("gaussian.result.output_count")
            else:
                from chemsmart.io.gaussian.output import Gaussian16Output

                rows = []
                observed_directions = set()

                def _level_token(value: Any) -> str:
                    return "".join(
                        character
                        for character in str(value or "").casefold()
                        if character.isalnum()
                    )

                expected_method = next(
                    (
                        requested.get(field)
                        for field in (
                            "functional",
                            "ab_initio",
                            "semiempirical",
                        )
                        if requested.get(field) is not None
                    ),
                    None,
                )
                for artifact in candidates:
                    row: dict[str, Any] = {
                        "artifact_sha256": artifact.sha256,
                        "normal_termination": False,
                        "optimization_converged": False,
                        "charge": None,
                        "multiplicity": None,
                        "jobtype": None,
                        "method": None,
                        "basis": None,
                        "energy_hartree": None,
                        "vibrational_mode_count": 0,
                        "transition_count": 0,
                        "wavefunction_stability_history": (),
                    }
                    try:
                        path = _current_artifact_path(
                            artifact, field_name="Gaussian result output"
                        )
                        output = Gaussian16Output(str(path))
                        if not output.normal_termination:
                            from chemsmart.io.native_failure import (
                                summarize_gaussian_native_failure,
                            )

                            failure_summary = (
                                summarize_gaussian_native_failure(
                                    output.contents
                                )
                            )
                            if failure_summary is not None:
                                row["native_failure"] = (
                                    _bound_native_failure_summary(
                                        failure_summary,
                                        artifacts=(artifact,),
                                    )
                                )
                                findings.append(
                                    "gaussian.native_failure."
                                    + failure_summary.error_class
                                )
                        energies = tuple(
                            float(value) for value in output.energies
                        )
                        energy = energies[-1] if energies else None
                        frequencies = tuple(
                            float(value)
                            for value in (output.vibrational_frequencies or ())
                        )
                        transitions = tuple(output.tddft_transitions or ())
                        stability_history = tuple(
                            output.wavefunction_stability_history or ()
                        )
                        observed_jobtype = str(output.jobtype or "").lower()
                        optimization_converged = (
                            any(
                                "Optimization completed." in line
                                for line in output.contents
                            )
                            and not output.convergence_criterion_not_met
                        )
                        row.update(
                            {
                                "normal_termination": bool(
                                    output.normal_termination
                                ),
                                "optimization_converged": (
                                    optimization_converged
                                ),
                                "charge": output.charge,
                                "multiplicity": output.multiplicity,
                                "jobtype": observed_jobtype,
                                "method": output.method,
                                "basis": output.basis,
                                "energy_hartree": energy,
                                "vibrational_mode_count": len(frequencies),
                                "transition_count": len(transitions),
                                "wavefunction_stability_history": (
                                    stability_history
                                ),
                            }
                        )
                        if observed_jobtype in {"ircf", "ircr"}:
                            observed_directions.add(observed_jobtype)
                        if not output.normal_termination:
                            findings.append(
                                "gaussian.result.normal_termination"
                            )
                        if output.charge != charge:
                            findings.append("gaussian.result.charge_mismatch")
                        if output.multiplicity != multiplicity:
                            findings.append(
                                "gaussian.result.multiplicity_mismatch"
                            )
                        if energy is None or not math.isfinite(energy):
                            findings.append("gaussian.result.energy_missing")
                        if (
                            expected_result_jobtype in {"opt", "ts"}
                            and not optimization_converged
                        ):
                            findings.append(
                                "gaussian.result.optimization_not_converged"
                            )
                        if (
                            jobtype not in {"td", "irc"}
                            and observed_jobtype != expected_result_jobtype
                        ):
                            findings.append("gaussian.result.jobtype_mismatch")
                        if jobtype == "link" and (
                            not stability_history
                            or stability_history[-1]
                            != "stable_under_considered_perturbations"
                        ):
                            findings.append(
                                "gaussian.result.wavefunction_not_stable"
                            )
                        if (
                            bool(requested.get("freq")) or jobtype == "freq"
                        ) and not frequencies:
                            findings.append(
                                "gaussian.result.frequencies_missing"
                            )
                        if jobtype == "td":
                            requested_states = int(
                                requested.get("nstates") or 1
                            )
                            if len(transitions) < requested_states:
                                findings.append(
                                    "gaussian.result.transitions_missing"
                                )
                        if expected_method is not None and _level_token(
                            output.method
                        ) != _level_token(expected_method):
                            findings.append("gaussian.result.method_mismatch")
                        expected_basis = requested.get("basis")
                        if expected_basis is not None and _level_token(
                            output.basis
                        ) != _level_token(expected_basis):
                            findings.append("gaussian.result.basis_mismatch")
                    except Exception as exc:
                        row["parser_error_type"] = type(exc).__name__
                        findings.append("gaussian.result.unreadable")
                    rows.append(row)
                if (
                    expected_directions
                    and observed_directions != expected_directions
                ):
                    findings.append("gaussian.result.irc_direction_mismatch")
                gaussian_observation["outputs"] = tuple(rows)
                if expected_directions:
                    gaussian_observation["irc_directions"] = tuple(
                        sorted(observed_directions)
                    )
            observation["gaussian"] = gaussian_observation
        elif program not in {"pyscf", "xtb", "orca", "gaussian"}:
            findings.append("execution.program.validator_unavailable")
        normalized_findings = tuple(sorted(set(findings)))
        validator_schema_version = str(
            (observation.get("result_validation") or {}).get("schema_version")
            or "chemsmart.generic-result-validation.v1"
        )
        environment_validation = (
            observation.get("environment_validation") or {}
        )
        return _ExecutionValidationEvaluation(
            validator_id=(
                "pyscf-result-validator"
                if program == "pyscf"
                else (
                    "xtb-result-validator"
                    if program == "xtb"
                    else (
                        "orca-result-validator"
                        if program == "orca"
                        else (
                            "gaussian-result-validator"
                            if program == "gaussian"
                            else "program-result-validator"
                        )
                    )
                )
            ),
            validator_schema_version=validator_schema_version,
            validator_version="1",
            observations=canonical_data(observation),
            findings=normalized_findings,
            run_environment_receipt_sha256=str(
                environment_validation.get(
                    "run_environment_receipt_sha256", ""
                )
            ),
            environment_validation_sha256=str(
                environment_validation.get("validation_sha256", "")
            ),
        )

    @staticmethod
    def _validate_execution_outputs(
        **values: Any,
    ) -> tuple[bool, tuple[str, ...], tuple[str, ...]]:
        """Compatibility projection used by focused legacy tests.

        Runtime execution stores the complete ProgramResultValidationReceiptV1
        produced from ``_evaluate_execution_outputs`` instead of this digest.
        """

        evaluation = CommandCompiledToolHostV1._evaluate_execution_outputs(
            **values
        )
        digest = canonical_sha256(
            {
                "observations": evaluation.observations,
                "findings": evaluation.findings,
                "validator_schema_version": (
                    evaluation.validator_schema_version
                ),
                "validator_version": evaluation.validator_version,
            }
        )
        return evaluation.validated, (digest,), evaluation.findings

    def latest_workflow_draft_receipt(self) -> str:
        """Return the latest useful plan identity without claiming readiness."""

        if not self.workflow_drafts:
            raise ContractError("no command workflow draft has been observed")
        return next(reversed(self.workflow_drafts))

    def _inspect_calculation_artifact(self, turn_id: str, values: dict) -> Any:
        artifact = self._artifact(values["artifact_id"])
        project = self._artifact(values["project_artifact_id"])
        settings = self._get(
            self.settings_objects, values["settings_id"], "settings object"
        )
        run_receipt = self._get(
            self.run_receipts, values["run_receipt_id"], "run receipt"
        )
        receipt = inspect_generated_artifact(
            program=values["program"],
            settings=settings,
            artifact=artifact,
            project_artifact=project,
            expected_receipt=run_receipt,
        )
        self.result_inspections[receipt.receipt_sha256] = receipt
        self._emit(
            turn_id,
            EventKind.RESULT_VERIFIED,
            receipt.receipt_sha256,
            status=receipt.status,
            artifact_sha256=receipt.artifact_sha256,
            expected_receipt_sha256=receipt.expected_receipt_sha256,
        )
        return receipt

    def _extract_result_quantities(self, turn_id: str, values: dict) -> Any:
        artifact = self._artifact(values["artifact_id"])
        selectors = tuple(
            QuantitySelectorV1(
                quantity_id=str(item["quantity_id"]),
                selector=str(item["selector"]),
            )
            for item in values["selectors"]
        )
        receipt = extract_trusted_result_quantities(
            artifact=artifact,
            program=values["program"],
            selectors=selectors,
        )
        self.quantity_extractions[receipt.receipt_sha256] = receipt
        self.quantity_extraction_selectors[receipt.receipt_sha256] = tuple(
            sorted({selector.selector for selector in selectors})
        )
        self.quantity_extraction_bindings[receipt.receipt_sha256] = {
            selector.quantity_id: selector.selector for selector in selectors
        }
        record = canonical_data(receipt)
        record.pop("receipt_sha256")
        self._emit(
            turn_id,
            EventKind.RESULT_QUANTITIES_EXTRACTED,
            receipt.receipt_sha256,
            status=receipt.status,
            artifact_sha256=receipt.artifact_sha256,
            quantity_ids=tuple(
                item.quantity_id for item in receipt.quantities
            ),
            selector_bindings=self.quantity_extraction_bindings[
                receipt.receipt_sha256
            ],
            record=record,
        )
        return receipt

    def _derive_thermochemistry(self, turn_id: str, values: dict) -> Any:
        artifact = self._artifact(values["artifact_id"])
        receipt = derive_trusted_thermochemistry(
            artifact=artifact,
            program=values["program"],
            temperature_k=float(values["temperature_k"]),
            pressure_atm=float(values["pressure_atm"]),
            concentration_mol_l=(
                float(values["concentration_mol_l"])
                if "concentration_mol_l" in values
                else None
            ),
            entropy_method=str(values.get("entropy_method", "rrho")),
            entropy_cutoff_cm1=(
                float(values["entropy_cutoff_cm1"])
                if "entropy_cutoff_cm1" in values
                else None
            ),
            enthalpy_cutoff_cm1=(
                float(values["enthalpy_cutoff_cm1"])
                if "enthalpy_cutoff_cm1" in values
                else None
            ),
            alpha=int(values.get("alpha", 4)),
            use_weighted_mass=bool(values.get("use_weighted_mass", False)),
            frequency_scale_factor=float(
                values.get("frequency_scale_factor", 1.0)
            ),
        )
        self.thermochemistry_receipts[receipt.receipt_sha256] = receipt
        record = canonical_data(receipt)
        record.pop("receipt_sha256")
        # Default PySCF RRHO receipts retain the original v1 canonical body so
        # existing archives and newly derived values share one identity.  The
        # dataclass also exposes newer optional controls with their defaults;
        # omitting those defaults from the persisted Runtime record keeps the
        # event representation identical to the receipt that was actually
        # derived instead of rejecting valid thermochemistry as a digest
        # mismatch.
        if canonical_sha256(record) != receipt.receipt_sha256:
            legacy_defaults = {
                "concentration_mol_l": None,
                "entropy_method": "rrho",
                "entropy_cutoff_cm1": None,
                "enthalpy_cutoff_cm1": None,
                "alpha": 4,
                "use_weighted_mass": False,
                "frequency_scale_factor": 1.0,
            }
            if receipt.program == "pyscf" and all(
                record.get(key) == value
                for key, value in legacy_defaults.items()
            ):
                for key in legacy_defaults:
                    record.pop(key, None)
        self._emit(
            turn_id,
            EventKind.THERMOCHEMISTRY_DERIVED,
            receipt.receipt_sha256,
            status=receipt.status,
            artifact_sha256=receipt.artifact_sha256,
            temperature_k=receipt.temperature_k,
            pressure_atm=receipt.pressure_atm,
            record=record,
        )
        return receipt

    def _evaluate_quantity_expression(self, turn_id: str, values: dict) -> Any:
        inputs: list[QuantityValueV1] = []
        for item in values["inputs"]:
            input_id = str(item["input_id"])
            receipt_sha256 = str(item["receipt_sha256"])
            quantity_id = str(item["quantity_id"])
            _kind, _receipt, source = self._typed_quantity_from_receipt(
                receipt_sha256=receipt_sha256,
                quantity_id=quantity_id,
                operation="quantity expression",
            )
            semantic_role = str(item.get("semantic_role", "")).strip()
            semantic_role_ref = (
                f";semantic-role:{semantic_role}" if semantic_role else ""
            )
            inputs.append(
                make_quantity_value(
                    quantity_id=input_id,
                    source_value=source.source_value,
                    source_unit=source.source_unit,
                    value=source.value,
                    unit=source.unit,
                    dimension=source.dimension,
                    evidence_ref=(
                        f"{source.evidence_ref};receipt:{receipt_sha256};"
                        f"quantity:{source.quantity_id}{semantic_role_ref}"
                    ),
                    data_kind=source.data_kind,
                )
            )
        nodes = tuple(
            QuantityExpressionNodeV1(
                node_id=str(item["node_id"]),
                operation=str(item["operation"]),
                input_ids=tuple(
                    str(value) for value in item.get("input_ids", ())
                ),
                reference=str(item.get("reference", "")),
                indices=tuple(int(value) for value in item.get("indices", ())),
                literal_value=item.get("literal_value"),
                literal_unit=str(item.get("literal_unit", "1")),
                scale_factor=(
                    float(item["scale_factor"])
                    if "scale_factor" in item
                    else None
                ),
                target_unit=str(item.get("target_unit", "")),
                cardinal_numbers=tuple(
                    int(value) for value in item.get("cardinal_numbers", ())
                ),
                extrapolation_exponent=(
                    float(item["extrapolation_exponent"])
                    if "extrapolation_exponent" in item
                    else None
                ),
            )
            for item in values["nodes"]
        )
        request = QuantityExpressionRequestV1(
            schema_version="chemsmart.quantity-expression-request.v1",
            expression_id=str(values["expression_id"]),
            inputs=tuple(inputs),
            nodes=nodes,
            output_node_ids=tuple(
                str(item) for item in values["output_node_ids"]
            ),
        )
        receipt = evaluate_typed_quantity_expression(request)
        self.quantity_expression_receipts[receipt.receipt_sha256] = receipt
        self.quantity_expression_requests[receipt.receipt_sha256] = request
        record = canonical_data(receipt)
        record.pop("receipt_sha256")
        self._emit(
            turn_id,
            EventKind.QUANTITY_EXPRESSION_EVALUATED,
            receipt.receipt_sha256,
            status=receipt.status,
            output_ids=tuple(item.quantity_id for item in receipt.outputs),
            semantic_signature_sha256=receipt.semantic_signature_sha256,
            record=record,
        )
        return receipt

    def _typed_quantity_from_receipt(
        self,
        *,
        receipt_sha256: str,
        quantity_id: str,
        operation: str,
    ) -> tuple[str, Any, QuantityValueV1]:
        """Resolve one exact typed quantity without accepting a model path."""

        registries = (
            ("quantity_extraction", self.quantity_extractions, "quantities"),
            ("thermochemistry", self.thermochemistry_receipts, "quantities"),
            (
                "quantity_expression",
                self.quantity_expression_receipts,
                "outputs",
            ),
            (
                "scientific_validation",
                self.scientific_validation_receipts,
                "outputs",
            ),
        )
        for source_kind, registry, collection in registries:
            receipt = registry.get(receipt_sha256)
            if receipt is None:
                continue
            matches = tuple(
                quantity
                for quantity in getattr(receipt, collection)
                if quantity.quantity_id == quantity_id
            )
            if len(matches) != 1:
                raise ContractError(
                    f"{operation} input is absent or ambiguous in its receipt"
                )
            return source_kind, receipt, matches[0]
        raise ContractError(f"{operation} references an unknown receipt")

    def _evaluate_scientific_validation(
        self, turn_id: str, values: dict
    ) -> ScientificValidationReceiptV1:
        """Evaluate one exact validation node without model-authored rules."""

        workflow_id = str(values["workflow_id"])
        node_id = str(values["node_id"])
        resolved = self._resolve_program_workflow(workflow_id)
        plan = resolved.scientific_toolchain_plan
        if plan is None:
            raise ContractError("workflow has no scientific analysis toolchain")
        nodes = tuple(
            node for node in plan.analysis_nodes if node.node_id == node_id
        )
        if len(nodes) != 1:
            raise ContractError("workflow has no unique validation node ID")
        node = nodes[0]
        if (
            node.analysis_kind != "scientific_validation"
            or node.support_state != "planned"
        ):
            raise ContractError("requested node is not planned validation")

        matched = self._scientific_toolchain_analysis_receipts(
            plan,
            task_spec_sha256=resolved.draft.task_spec_id,
        )
        input_intents = {
            item.input_id: item
            for item in node.inputs
            if isinstance(item, AnalysisInputIntentV1)
        }
        supplied = tuple(values["inputs"])
        supplied_ids = tuple(sorted(str(item["input_id"]) for item in supplied))
        if (
            len(supplied_ids) != len(set(supplied_ids))
            or supplied_ids != tuple(sorted(input_intents))
        ):
            raise ContractError(
                "scientific validation requires exactly its planned inputs"
            )
        producer_nodes = {
            item.node_id: item for item in plan.analysis_nodes
        }
        bound_inputs: dict[str, tuple[str, QuantityValueV1]] = {}
        for item in supplied:
            input_id = str(item["input_id"])
            intent = input_intents[input_id]
            source_receipt_sha256 = str(item["receipt_sha256"])
            require_sha256(
                source_receipt_sha256,
                "scientific validation source receipt",
            )
            if source_receipt_sha256 not in matched.get(
                intent.producer_node_id, ()
            ):
                raise ContractError(
                    "scientific validation input is not typed evidence from "
                    "its planned producer"
                )
            producer = producer_nodes.get(intent.producer_node_id)
            if producer is None:
                raise ContractError(
                    "scientific validation input producer is not analysis"
                )
            expected_kinds = {
                "result_extraction": "quantity_extraction",
                "thermochemistry": "thermochemistry",
                "quantity_expression": "quantity_expression",
                "scientific_validation": "scientific_validation",
            }
            expected_kind = expected_kinds.get(producer.analysis_kind)
            if expected_kind is None:
                raise ContractError(
                    "scientific validation producer has no typed quantities"
                )
            source_kind, source_receipt, quantity = (
                self._typed_quantity_from_receipt(
                    receipt_sha256=source_receipt_sha256,
                    quantity_id=str(item["quantity_id"]),
                    operation="scientific validation",
                )
            )
            if source_kind != expected_kind:
                raise ContractError(
                    "scientific validation receipt kind differs from producer"
                )
            declared_outputs = tuple(
                output
                for output in producer.outputs
                if output.output_id == intent.producer_output_id
            )
            if len(declared_outputs) != 1:
                raise ContractError(
                    "scientific validation input lacks a unique planned output"
                )
            declared = declared_outputs[0]
            try:
                _value, _unit, dimension = normalize_numeric_value(
                    0.0, declared.unit
                )
            except (ContractError, ValueError) as exc:
                raise ContractError(
                    "scientific validation producer unit is invalid"
                ) from exc
            if tuple(quantity.dimension) != tuple(dimension):
                raise ContractError(
                    "scientific validation quantity dimension differs from "
                    "its planned producer output"
                )
            available_ids = {
                candidate.quantity_id
                for candidate in getattr(
                    source_receipt,
                    "outputs"
                    if source_kind
                    in {"quantity_expression", "scientific_validation"}
                    else "quantities",
                )
            }
            if (
                intent.producer_output_id in available_ids
                and quantity.quantity_id != intent.producer_output_id
            ):
                raise ContractError(
                    "scientific validation selected another quantity despite "
                    "the planned producer output being present"
                )
            if source_kind == "quantity_extraction":
                selector = self.quantity_extraction_bindings.get(
                    source_receipt_sha256, {}
                ).get(quantity.quantity_id)
                if selector not in {
                    planned.selector for planned in producer.selectors
                }:
                    raise ContractError(
                        "scientific validation extraction quantity is outside "
                        "the planned selector set"
                    )
            bound_inputs[input_id] = (
                source_receipt_sha256,
                quantity,
            )

        receipt = evaluate_planned_scientific_validation(
            workflow_id=workflow_id,
            plan_sha256=plan.plan_sha256,
            node=node,
            inputs=bound_inputs,
        )
        self.scientific_validation_receipts[receipt.receipt_sha256] = receipt
        record = canonical_data(receipt)
        record.pop("receipt_sha256")
        self.event_store.append(
            turn_id=turn_id,
            kind=EventKind.SCIENTIFIC_VALIDATION_EVALUATED.value,
            payload={
                "receipt_sha256": receipt.receipt_sha256,
                "workflow_id": receipt.workflow_id,
                "plan_sha256": receipt.plan_sha256,
                "node_id": receipt.node_id,
                "source_receipt_sha256s": (
                    receipt.source_receipt_sha256s
                ),
                "all_rules_passed": receipt.all_rules_passed,
                "status": receipt.status,
                "record": record,
            },
            idempotency_key=(
                "scientific-validation:" + receipt.receipt_sha256
            ),
        )
        return receipt

    def _record_analysis_claims(self, turn_id: str, values: dict) -> Any:
        """Render reportable values from exact typed receipt outputs."""

        task_spec_sha256 = self._resolve_task_spec_reference(
            values, "task_spec_sha256"
        )
        if task_spec_sha256 not in self.task_spec_sha256s:
            raise ContractError("analysis claims target an unknown task spec")
        registries = (
            ("quantity_extraction", self.quantity_extractions, "quantities"),
            ("thermochemistry", self.thermochemistry_receipts, "quantities"),
            (
                "quantity_expression",
                self.quantity_expression_receipts,
                "outputs",
            ),
            (
                "scientific_validation",
                self.scientific_validation_receipts,
                "outputs",
            ),
        )
        claims = []
        for item in values["claims"]:
            receipt_sha256 = str(item["receipt_sha256"])
            require_sha256(receipt_sha256, "analysis claim receipt_sha256")
            source_kind = ""
            source_receipt = None
            quantity_collection = ""
            for kind, registry, collection in registries:
                if receipt_sha256 in registry:
                    source_kind = kind
                    source_receipt = registry[receipt_sha256]
                    quantity_collection = collection
                    break
            if source_receipt is None:
                raise ContractError("analysis claim cites an unknown receipt")
            quantity_id = str(item["quantity_id"])
            quantities = getattr(source_receipt, quantity_collection)
            matches = tuple(
                quantity
                for quantity in quantities
                if quantity.quantity_id == quantity_id
            )
            if len(matches) != 1:
                raise ContractError(
                    "analysis claim quantity is absent or ambiguous"
                )
            quantity = matches[0]
            if quantity.data_kind in {"text", "text_vector"}:
                raise ContractError("analysis claims must be numerical")
            display_unit = str(item["display_unit"])
            display_value = convert_normalized_value(
                quantity.value, quantity.dimension, display_unit
            )
            claims.append(
                AnalysisReportedQuantityV1(
                    claim_id=str(item["claim_id"]),
                    source_kind=source_kind,
                    source_receipt_sha256=receipt_sha256,
                    quantity_id=quantity.quantity_id,
                    quantity_value_sha256=quantity.value_sha256,
                    display_value=display_value,
                    display_unit=display_unit,
                    canonical_value=quantity.value,
                    canonical_unit=quantity.unit,
                    dimension=quantity.dimension,
                    data_kind=quantity.data_kind,
                )
            )
        record = build_analysis_claim_record(
            task_spec_sha256=task_spec_sha256,
            claims=tuple(claims),
        )
        self.analysis_claim_records[record.receipt_sha256] = record
        record_body = canonical_data(record)
        record_body.pop("receipt_sha256")
        self.event_store.append(
            turn_id=turn_id,
            kind=EventKind.ANALYSIS_CLAIMS_RECORDED.value,
            payload={
                "receipt_sha256": record.receipt_sha256,
                "task_spec_sha256": record.task_spec_sha256,
                "status": record.status,
                "source_receipt_sha256s": tuple(
                    sorted(
                        {
                            claim.source_receipt_sha256
                            for claim in record.claims
                        }
                    )
                ),
                "claim_ids": tuple(claim.claim_id for claim in record.claims),
                "critical_finding_count": 0,
                "record": record_body,
            },
            idempotency_key="analysis-claims:" + record.receipt_sha256,
        )
        return record

    def _record_compiled_command(
        self,
        turn_id: str,
        invocation: CanonicalCommandInvocationV1,
        context: _CommandContext,
    ) -> dict[str, Any]:
        inspection = inspect_command(invocation, live_schema=self.live_schema)
        self.invocations[invocation.invocation_sha256] = invocation
        self.command_inspections[inspection.receipt_sha256] = inspection
        self._command_contexts[invocation.invocation_sha256] = context
        # Carry the node's identity, not only its digests.  An approval binds
        # exact nodes, so a record that says a command was compiled without
        # saying which node it was, on which program, at what charge, cannot
        # be reviewed into one -- the reviewer would have to re-plan and hope
        # for a byte-identical workflow.
        self._emit(
            turn_id,
            EventKind.COMMAND_COMPILED,
            invocation.invocation_sha256,
            status=invocation.status,
            node_id=invocation.node_id,
            program=context.proposal.program,
            jobtype=context.proposal.jobtype,
            execution_target=context.proposal.execution_target,
            charge=context.proposal.charge,
            multiplicity=context.proposal.multiplicity,
            display_command=invocation.display_command,
            input_sha256=invocation.input_sha256,
            project_sha256=invocation.project_sha256,
            project_receipt_sha256=invocation.project_receipt_sha256,
            program_engine_binding_sha256=(
                invocation.program_engine_binding_sha256
            ),
            scientific_identity_sha256=(invocation.scientific_identity_sha256),
        )
        self._emit(
            turn_id,
            EventKind.COMMAND_INSPECTED,
            inspection.receipt_sha256,
            status=inspection.status,
            invocation_sha256=inspection.invocation_sha256,
        )
        return {"invocation": invocation, "inspection": inspection}

    def _emit_binding(self, turn_id: str, kind: EventKind, value: Any) -> None:
        self.event_store.append(
            turn_id=turn_id,
            kind=kind.value,
            payload={
                "binding_sha256": value.binding_sha256,
                "state": value.state,
                "program": getattr(value, "program", ""),
            },
            idempotency_key=f"{kind.value}:{value.binding_sha256}",
        )

    def _emit(
        self,
        turn_id: str,
        kind: EventKind,
        receipt_sha256: str,
        **payload: Any,
    ) -> None:
        body = {"receipt_sha256": receipt_sha256, **payload}
        self.event_store.append(
            turn_id=turn_id,
            kind=kind.value,
            payload=body,
            idempotency_key=f"{kind.value}:{receipt_sha256}",
        )

    def _artifact(self, artifact_id: str) -> TrustedArtifactRefV1:
        return self._get(self.artifacts, artifact_id, "trusted artifact")

    @staticmethod
    def _get(values: Mapping[str, Any], key: str, label: str) -> Any:
        try:
            return values[key]
        except KeyError as exc:
            # This is the shared lookup behind every host-bound object, so the
            # message it raises is the one a caller sees for most mistaken IDs.
            # "unknown X ID" names neither what was asked for nor what exists,
            # which leaves retrying blind; listing the bound IDs makes the
            # rejection something the caller can act on.
            known = sorted(values)
            if not known:
                # An empty registry is the case where listing IDs helps least
                # and naming the producer helps most.  Observed live: a session
                # asked twice for a counterexample that no failure had yet
                # produced, because being told the registry was empty did not
                # say what fills it.
                detail = f"no {label} is bound yet"
                producer = REGISTRY_PRODUCERS.get(label)
                if producer:
                    detail += f"; one is bound {producer}"
            elif len(known) <= 8:
                detail = f"bound {label} IDs: {known}"
            else:
                detail = (
                    f"bound {label} IDs include {known[:8]} "
                    f"and {len(known) - 8} more"
                )
            raise ContractError(
                f"unknown {label} ID {key!r}; {detail}"
            ) from exc


#: How many validator findings an event records, and how much of each value.
#: A rejection is for reading, not a channel for arbitrary bytes, so both the
#: count and the rendered values stay bounded.
_RECORDED_FINDINGS = 8
_FINDING_VALUE_CHARS = 120


def _public_validator_findings(validator: Any) -> tuple[dict[str, str], ...]:
    """Render preview-validation findings as bounded, reviewable records."""

    def _text(value: Any) -> str:
        if isinstance(value, (dict, list, tuple)):
            rendered = f"a {type(value).__name__} of {len(value)} entries"
        else:
            rendered = repr(value)
        if len(rendered) > _FINDING_VALUE_CHARS:
            rendered = rendered[: _FINDING_VALUE_CHARS - 3] + "..."
        return rendered

    findings = tuple(getattr(validator, "findings", ()) or ())
    recorded = tuple(
        {
            "rule_id": str(getattr(item, "rule_id", "")),
            "field": str(getattr(item, "field", "")),
            "expected": _text(getattr(item, "expected", None)),
            "observed": _text(getattr(item, "observed", None)),
        }
        for item in findings[:_RECORDED_FINDINGS]
    )
    if len(findings) > _RECORDED_FINDINGS:
        recorded += (
            {
                "rule_id": "record.truncated",
                "field": "",
                "expected": f"{len(findings)} findings",
                "observed": f"first {_RECORDED_FINDINGS} recorded",
            },
        )
    return recorded


def _validate_tool_arguments(
    surface: AgentToolSurfaceV1, tool_name: str, arguments: dict[str, Any]
) -> None:
    definition = next(
        (
            item["function"]
            for item in surface.tool_definitions
            if item["function"]["name"] == tool_name
        ),
        None,
    )
    if definition is None:
        raise ContractError("tool is not exposed by this profile")
    schema = definition["parameters"]
    required = set(schema.get("required", ()))
    properties = schema.get("properties", {})
    missing = sorted(required.difference(arguments))
    if missing:
        # A caller that is told only the field name has to go back to the
        # schema; carrying the field's own description makes the rejection
        # self-contained.
        described = "; ".join(
            (
                f"{name} ({properties[name]['description']})"
                if isinstance(properties.get(name), Mapping)
                and properties[name].get("description")
                else name
            )
            for name in missing
        )
        raise ContractError(f"{tool_name} requires {described}")
    unknown = sorted(set(arguments).difference(properties))
    if unknown:
        raise ContractError(
            f"{tool_name} does not accept {unknown}; it accepts "
            f"{sorted(properties)}"
        )
    for name, value in arguments.items():
        _validate_json_value(name, value, properties[name])


#: How much of an offending value a rejection may quote back.  Long enough to
#: identify a wrong identifier or a malformed path, short enough that a large
#: payload cannot be echoed through a refusal.
_REJECTED_VALUE_CHARS = 120


def _offending_text(value: Any) -> str:
    """Quote a rejected value compactly, never a whole payload."""

    if isinstance(value, (dict, list)):
        text = f"a {type(value).__name__} of {len(value)} entries"
    else:
        text = repr(value)
    if len(text) > _REJECTED_VALUE_CHARS:
        text = text[: _REJECTED_VALUE_CHARS - 3] + "..."
    return text


def _validate_json_value(
    name: str, value: Any, schema: Mapping[str, Any]
) -> None:
    """Reject an argument by naming the path, the value, and the rule.

    A caller that is told only which field is wrong has to guess at the value
    and the constraint, and a model that guesses generally resubmits the same
    argument.  Every message here therefore carries all three.
    """

    alternatives = schema.get("oneOf")
    if isinstance(alternatives, list):
        matches = 0
        for alternative in alternatives:
            if not isinstance(alternative, Mapping):
                continue
            try:
                _validate_json_value(name, value, alternative)
            except ContractError:
                continue
            matches += 1
        if matches != 1:
            shapes = ", ".join(
                str(item.get("type", "?"))
                for item in alternatives
                if isinstance(item, Mapping)
            )
            raise ContractError(
                f"tool argument {name} is {_offending_text(value)}, which "
                f"matches {matches} of the allowed shapes; exactly one must "
                f"match. Allowed shapes: {shapes}"
            )
        return
    expected = schema.get("type")

    def _matches(kind: str) -> bool:
        return {
            "string": isinstance(value, str),
            "integer": isinstance(value, int) and not isinstance(value, bool),
            "number": (
                isinstance(value, (int, float))
                and not isinstance(value, bool)
                and math.isfinite(float(value))
            ),
            "boolean": isinstance(value, bool),
            "array": isinstance(value, list),
            "object": isinstance(value, dict),
            "null": value is None,
        }.get(kind, True)

    # A union type such as ["number", "null"] is how a field says the concept
    # may not apply here.  Treating it as an unknown type accepted anything.
    if isinstance(expected, list):
        allowed = [str(item) for item in expected]
        if not any(_matches(kind) for kind in allowed):
            raise ContractError(
                f"tool argument {name} must be one of {allowed}, but got "
                f"{type(value).__name__} {_offending_text(value)}"
            )
        if value is None:
            # No further keyword applies to an explicit null.
            return
    elif not _matches(expected):
        raise ContractError(
            f"tool argument {name} must be {expected}, but got "
            f"{type(value).__name__} {_offending_text(value)}"
        )
    if "enum" in schema and value not in schema["enum"]:
        allowed = list(schema["enum"])
        raise ContractError(
            f"tool argument {name} is {_offending_text(value)}, which is not "
            f"one of {allowed}"
        )
    if isinstance(value, str) and schema.get("pattern"):
        pattern = str(schema["pattern"])
        if re.fullmatch(pattern, value) is None:
            raise ContractError(
                f"tool argument {name} is {_offending_text(value)}, which "
                f"does not match the required pattern {pattern}"
            )
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        if "minimum" in schema and value < float(schema["minimum"]):
            raise ContractError(
                f"tool argument {name} is {value}, below its minimum "
                f"{schema['minimum']}"
            )
        if "maximum" in schema and value > float(schema["maximum"]):
            raise ContractError(
                f"tool argument {name} is {value}, above its maximum "
                f"{schema['maximum']}"
            )
        if "exclusiveMinimum" in schema and value <= float(
            schema["exclusiveMinimum"]
        ):
            raise ContractError(
                f"tool argument {name} is {value}, which must be greater than "
                f"{schema['exclusiveMinimum']}"
            )
        if "exclusiveMaximum" in schema and value >= float(
            schema["exclusiveMaximum"]
        ):
            raise ContractError(
                f"tool argument {name} is {value}, which must be less than "
                f"{schema['exclusiveMaximum']}"
            )
    if isinstance(value, list):
        if "minItems" in schema and len(value) < int(schema["minItems"]):
            raise ContractError(
                f"tool argument {name} has {len(value)} items, fewer than the "
                f"required {schema['minItems']}"
            )
        if "maxItems" in schema and len(value) > int(schema["maxItems"]):
            raise ContractError(
                f"tool argument {name} has {len(value)} items, more than the "
                f"allowed {schema['maxItems']}"
            )
    if isinstance(value, list) and isinstance(schema.get("items"), Mapping):
        for index, item in enumerate(value):
            # Index the path: "outputs[]" tells a caller with eight outputs
            # nothing about which one to change.
            _validate_json_value(f"{name}[{index}]", item, schema["items"])
    if isinstance(value, dict):
        properties = schema.get("properties", {})
        additional = schema.get("additionalProperties")
        required = set(schema.get("required", ()))
        missing = sorted(required.difference(value))
        if missing:
            raise ContractError(
                f"tool argument {name} is missing {missing}; it supplied "
                f"{sorted(value)}"
            )
        if schema.get("additionalProperties") is False:
            unknown = sorted(set(value).difference(properties))
            if unknown:
                raise ContractError(
                    f"tool argument {name} supplied {unknown}, which this "
                    f"object does not accept; it accepts {sorted(properties)}"
                )
        for child_name, child_value in value.items():
            child_schema = properties.get(child_name)
            if child_schema is None and isinstance(additional, Mapping):
                child_schema = additional
            if child_schema is not None:
                _validate_json_value(
                    f"{name}.{child_name}", child_value, child_schema
                )


def _require_registry_keys(
    values: Mapping[str, Any], attribute: str, label: str
) -> None:
    if any(key != getattr(value, attribute) for key, value in values.items()):
        raise ContractError(f"{label} registry key mismatch")


def _program_process_environment(
    *, overrides: Mapping[str, str], remove: tuple[str, ...]
) -> dict[str, str]:
    """Build an engine environment without provider credential labels."""

    environment = os.environ.copy()
    environment.update(
        {str(key): str(value) for key, value in overrides.items()}
    )
    for key in remove:
        environment.pop(key, None)
    return environment


def _public_process_stream(value: str | bytes | None) -> str:
    """Return process output as bounded public text, never a Python repr."""

    if value is None:
        return ""
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    return str(value)


def _write_host_execution_artifact(path: Path, payload: str) -> None:
    """Write one host-owned node artifact without replacing child output."""

    try:
        with path.open("x", encoding="utf-8") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
    except FileExistsError as exc:
        raise ContractError(
            "execution emitted a reserved host artifact name"
        ) from exc


def _process_observation_findings(
    observation: ProcessObservationV1,
) -> tuple[str, ...]:
    """Map a digest-valid process observation to deterministic findings."""

    findings: set[str] = set()
    if observation.state == "launch_failed":
        findings.add("execution.process.launch_failed")
    if observation.timed_out:
        findings.add("execution.process.timeout")
    if observation.memory_limit_exceeded:
        findings.add("execution.process.memory_limit_exceeded")
    if observation.termination_requested and not (
        observation.termination_confirmed
    ):
        findings.add("execution.process.termination_ambiguous")
    if observation.pid is not None and not observation.process_group_owned:
        findings.add("execution.process.group_not_owned")
    if (
        observation.state == "exited"
        and observation.memory_observation_state != "observed"
    ):
        findings.add("execution.process.memory_observation_unavailable")
    return tuple(sorted(findings))


def _prepare_execution_node_workspace(path: Path) -> Path:
    """Create one empty node directory without deleting prior evidence."""

    if path.is_symlink():
        raise ContractError("execution node workspace cannot be a symlink")
    if path.exists():
        if not path.is_dir():
            raise ContractError("execution node workspace is not a directory")
        if any(path.iterdir()):
            raise ContractError(
                "execution workspace already contains outputs; import or "
                "inspect them before authorizing another launch"
            )
    else:
        path.mkdir(parents=True)
    return path


def _require_current_auxiliary_inputs(
    job_artifact_options: tuple[tuple[str, TrustedArtifactRefV1], ...],
) -> None:
    """Recheck every additional command input immediately before launch."""

    for parameter_name, artifact in job_artifact_options:
        source = Path(artifact.path)
        if not source.is_file():
            raise ContractError(
                f"auxiliary input {parameter_name} is unavailable"
            )
        before = source.stat()
        observed_sha256 = file_sha256(source)
        after = source.stat()
        if (
            before.st_size != after.st_size
            or before.st_mtime_ns != after.st_mtime_ns
        ):
            raise ContractError(
                f"auxiliary input {parameter_name} changed while inspected"
            )
        if (
            after.st_size != artifact.size_bytes
            or observed_sha256 != artifact.sha256
        ):
            raise ContractError(
                f"auxiliary input {parameter_name} differs from approval"
            )


def _remaining_node_unresolved_fields(
    node: Any,
    materialized_producer_inputs: Mapping[
        tuple[str, str, str, str], TrustedArtifactRefV1
    ],
) -> tuple[str, ...]:
    """Clear only stale input markers backed by an exact validated handoff.

    Planning legitimately marks a future optimized geometry unresolved.  The
    draft stays immutable after approval, so that marker must be interpreted
    against current host evidence rather than copied forever.  Scientific
    fields such as a method, solvent, state, or IRC direction are never cleared
    here; only names derived from the exact resolved producer input are.
    """

    resolved_input_markers: set[str] = set()
    for item in getattr(node, "inputs", ()) or ():
        producer_node_id = str(getattr(item, "producer_node_id", "") or "")
        producer_output_id = str(
            getattr(item, "producer_output_id", "") or ""
        )
        if not producer_node_id or not producer_output_id:
            continue
        key = (
            str(node.node_id),
            str(item.binding_id),
            producer_node_id,
            producer_output_id,
        )
        if key not in materialized_producer_inputs:
            continue
        binding_id = str(item.binding_id)
        artifact_class = str(getattr(item, "artifact_class", "") or "")
        resolved_input_markers.update(
            {
                binding_id,
                producer_output_id,
                artifact_class,
                "input." + binding_id,
                "input." + artifact_class,
                "input_artifact",
            }
        )
        if artifact_class == "geometry_xyz":
            resolved_input_markers.add("geometry")

    normalized_markers = {
        marker.replace("-", "_") for marker in resolved_input_markers
    }
    return tuple(
        field
        for field in getattr(node, "unresolved_fields", ()) or ()
        if str(field).replace("-", "_") not in normalized_markers
    )


def _staged_auxiliary_input_findings(
    *,
    node_workspace: Path,
    job_artifact_options: tuple[tuple[str, TrustedArtifactRefV1], ...],
) -> tuple[str, ...]:
    """Confirm a multi-file writer staged the approved bytes it references."""

    findings = []
    for parameter_name, artifact in job_artifact_options:
        basename = Path(artifact.path).name
        candidates = tuple(
            path
            for path in node_workspace.rglob(basename)
            if path.is_file() and not path.is_symlink()
        )
        matching = tuple(
            path
            for path in candidates
            if path.stat().st_size == artifact.size_bytes
            and file_sha256(path) == artifact.sha256
        )
        if len(matching) != 1:
            findings.append(
                f"execution.auxiliary_input_not_staged.{parameter_name}"
            )
    return tuple(sorted(findings))


def _scientific_plan_from_v1_approval(
    approval: WorkflowExecutionApprovalV1,
) -> ScientificWorkflowPlanV2:
    """Compatibility projection that does not turn tuple order into edges."""

    identity_sha256s = tuple(
        sorted(
            {
                node.scientific_identity_sha256
                for node in approval.node_bindings
                if node.scientific_identity_sha256
            }
        )
    )
    if not identity_sha256s:
        raise ContractError("V1 workflow approval lacks scientific identity")
    scientific_identity_sha256 = (
        identity_sha256s[0]
        if len(identity_sha256s) == 1
        else canonical_sha256(
            {"scientific_identity_sha256s": identity_sha256s}
        )
    )
    nodes = tuple(
        ScientificWorkflowNodeV2(
            node_id=node.node_id,
            stage=node.jobtype,
            requested_program=node.program,
            program=node.program,
            engine=node.engine,
            project_role="approved." + node.program,
            unresolved_fields=(),
        )
        for node in approval.node_bindings
    )
    edges = tuple(
        sorted(
            (
                ScientificWorkflowEdgeV2(
                    edge_id=(
                        "data."
                        + edge.producer_node_id
                        + "."
                        + edge.consumer_node_id
                        + "."
                        + edge.artifact_kind
                    ),
                    source_node_id=edge.producer_node_id,
                    target_node_id=edge.consumer_node_id,
                    edge_kind="data",
                    artifact_class=edge.artifact_kind,
                    producer_output_id=edge.selection_rule,
                    consumer_input_id="geometry",
                )
                for edge in approval.producer_edges
            ),
            key=lambda edge: edge.edge_id,
        )
    )
    return build_scientific_workflow_plan(
        workflow_id=approval.workflow_id,
        task_spec_sha256=approval.task_spec_sha256,
        scientific_identity_sha256=scientific_identity_sha256,
        nodes=nodes,
        edges=edges,
    )


def _project_v1_execution_run_state(
    plan: ScientificWorkflowPlanV2,
    approval: WorkflowExecutionApprovalV1,
    receipts: Mapping[str, ProgramExecutionReceiptV1],
) -> WorkflowRunStateV1:
    """Project legacy receipts into V2 solely for deterministic readiness."""

    node_states = []
    for planned_node in sorted(plan.nodes, key=lambda node: node.node_id):
        receipt = receipts.get(planned_node.node_id)
        if receipt is None:
            node_states.append(
                WorkflowNodeRunStateV1(
                    node_id=planned_node.node_id,
                    state="pending",
                    invocation_sha256="",
                    execution_receipt_sha256="",
                    validator_receipt_sha256s=(),
                    output_artifact_sha256s=(),
                    failure_rule_ids=(),
                )
            )
            continue
        if receipt.validated:
            state = "validated"
        elif receipt.execution_state == "engine_complete":
            state = "engine_complete"
        elif receipt.execution_state in {"failed", "ambiguous", "running"}:
            state = receipt.execution_state
        else:
            state = "blocked"
        failure_rule_ids = ()
        if state in {"failed", "ambiguous", "blocked"}:
            failure_rule_ids = tuple(
                sorted(
                    set(receipt.findings)
                    or {"execution.state." + receipt.execution_state}
                )
            )
        node_states.append(
            WorkflowNodeRunStateV1(
                node_id=planned_node.node_id,
                state=state,
                invocation_sha256=receipt.invocation_sha256,
                execution_receipt_sha256=receipt.receipt_sha256,
                validator_receipt_sha256s=(receipt.validator_receipt_sha256s),
                output_artifact_sha256s=tuple(
                    sorted(
                        artifact.sha256
                        for artifact in receipt.output_artifacts
                    )
                ),
                failure_rule_ids=failure_rule_ids,
            )
        )
    observed_states = {node.state for node in node_states}
    if "ambiguous" in observed_states:
        workflow_state = "ambiguous"
    elif "failed" in observed_states:
        workflow_state = "failed"
    elif "blocked" in observed_states:
        workflow_state = "blocked"
    elif observed_states == {"validated"}:
        workflow_state = "validated"
    elif observed_states == {"pending"}:
        workflow_state = "approved"
    else:
        workflow_state = "running"
    started_values = tuple(
        receipt.started_at
        for receipt in receipts.values()
        if receipt.started_at
    )
    finished_values = tuple(
        receipt.finished_at
        for receipt in receipts.values()
        if receipt.finished_at
    )
    body = {
        "schema_version": "chemsmart.workflow-run-state.v1",
        "run_id": "compat." + approval.approval_id,
        "workflow_id": plan.workflow_id,
        "plan_sha256": plan.plan_sha256,
        "approval_id": approval.approval_id,
        "approval_sha256": approval.approval_sha256,
        "approval_consumed": True,
        "state": workflow_state,
        "nodes": tuple(node_states),
        "started_at": min(started_values) if started_values else "",
        "finished_at": (
            max(finished_values)
            if workflow_state in {"validated", "failed"} and finished_values
            else ""
        ),
    }
    return WorkflowRunStateV1(**body, run_state_sha256=canonical_sha256(body))


def _native_diagnostic_artifacts(
    artifacts: tuple[TrustedArtifactRefV1, ...],
) -> tuple[TrustedArtifactRefV1, ...]:
    """Return at most two bound native stderr sidecars for parser diagnosis."""

    return tuple(
        artifact
        for artifact in artifacts
        if artifact.kind == "program_output"
        and Path(artifact.path).suffix.lower() == ".err"
    )[:2]


def _iter_native_diagnostic_lines(
    artifacts: tuple[TrustedArtifactRefV1, ...],
):
    """Stream bound sidecar lines without retaining or exposing the log."""

    for artifact in artifacts:
        try:
            path = _current_artifact_path(
                artifact, field_name="native diagnostic sidecar"
            )
            with path.open(encoding="utf-8", errors="replace") as handle:
                yield from handle
        except (ContractError, OSError):
            continue


def _bound_native_failure_summary(
    summary: Any,
    *,
    artifacts: tuple[TrustedArtifactRefV1, ...],
) -> dict[str, Any]:
    """Attach path-free trusted artifact pointers to a parser summary."""

    record = dict(summary.as_dict())
    record["artifact_refs"] = tuple(
        {
            "artifact_id": artifact.artifact_id,
            "kind": artifact.kind,
            "sha256": artifact.sha256,
            "size_bytes": artifact.size_bytes,
        }
        for artifact in artifacts
    )
    return record


def _model_visible_data(value: Any) -> Any:
    """Remove host-local resolution details from model-visible tool results."""

    if isinstance(value, dict):
        return {
            key: _model_visible_data(item)
            for key, item in value.items()
            if key not in {"path", "cli_value", "workspace"}
        }
    if isinstance(value, (list, tuple)):
        return [_model_visible_data(item) for item in value]
    return value


__all__ = ["CommandCompiledToolHostV1"]
