"""Host-owned dispatcher for the command-compiled model tool surface."""

from __future__ import annotations

import hashlib
import json
import logging
import math
import os
import re
import shlex
import socket
import subprocess
import sys
import time
from dataclasses import dataclass, replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from chemsmart.agent._contracts import (
    ContractError,
    RoutedContractError,
    TrustedArtifactRefV1,
    canonical_data,
    canonical_sha256,
    file_sha256,
    require_identifier,
    require_sha256,
)
from chemsmart.agent.analysis_claims import (
    AnalysisReportedQuantityV1,
    analysis_claim_record_from_record,
    build_analysis_claim_record,
)
from chemsmart.agent.analysis_completion import (
    AnalysisCompletionPolicyV1,
    AnalysisCompletionReceiptV1,
    AnalysisIncompleteError,
    build_analysis_completion_receipt,
)
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
    CommandInspectionReceiptV1,
    CommandProposalV1,
    ScientificIdentityBindingV1,
    build_scientific_identity_binding,
    compile_command,
    inspect_command,
    native_coordinate_options,
)
from chemsmart.agent.delivery import _OPEN_STATES as _OPEN_SUFFICIENCY_STATES
from chemsmart.agent.delivery import (
    judge_sufficiency,
)
from chemsmart.agent.execution import (
    DEFERRABLE_GEOMETRY_PRODUCER_STAGES,
    HESSIAN_CONSUMER_ROLES,
    AnomalyObservationV1,
    ApprovedNodeBindingV1,
    AtomAppendReceiptV1,
    DatabaseRecordExtractionReceiptV1,
    ExecutionResourceSpecV1,
    FrozenWorkflowApprovalV1,
    GeometryEditReceiptV1,
    MolecularCompositionReceiptV1,
    MolecularDerivationReceiptV1,
    OptimizedGeometryHandoffV1,
    ORCAHessianHandoffV1,
    ProgramExecutionReceiptV1,
    ProgramResultValidationReceiptV1,
    ProjectArtifactPromotionV1,
    PubchemGeometryReceiptV1,
    ScientificDecisionRecordV1,
    SymmetryBreakReceiptV1,
    WorkflowEnvironmentBindingV1,
    WorkflowExecutionApprovalV1,
    WorkflowExecutionNodeReviewV1,
    WorkflowExecutionReviewV1,
    WorkflowNodeRunStateV1,
    WorkflowRunStateV1,
    anomaly_standing,
    append_trusted_molecular_atom,
    bind_project_promotion_validation,
    break_trusted_molecular_symmetry,
    build_anomaly_observation,
    build_frozen_workflow_approval,
    build_producer_edge_rule,
    build_program_execution_invocation,
    build_program_execution_receipt,
    build_program_result_validation_receipt,
    build_reached_geometry,
    build_real_execution_argv,
    build_scientific_decision_record,
    build_stationary_point_characterisation,
    build_validated_data_edge_binding,
    build_workflow_approval_request,
    build_workflow_execution_approval,
    build_workflow_execution_node_review,
    build_workflow_execution_review,
    compose_trusted_molecular_arrangement,
    derive_ready_node_ids,
    derive_trusted_molecular_species,
    displace_trusted_geometry_along_mode,
    environment_review_summary,
    execution_path_placeholder,
    execution_server_profile_sha256,
    existing_node_branches,
    extract_trusted_database_record_geometry,
    fetch_trusted_pubchem_geometry,
    handoff_final_orca_ts_hessian,
    handoff_optimized_native_geometry,
    handoff_optimized_pyscf_geometry,
    handoff_optimized_xtb_geometry,
    handoff_scan_minimum_geometry,
    handoff_validated_orca_producer_hessian,
    hessian_role_for_rule,
    invocation_identity_sha256,
    is_validated_optimized_geometry_edge,
    is_validated_orca_ts_hessian_edge,
    is_validated_producer_orca_hessian_edge,
    is_validated_scan_minimum_geometry_edge,
    node_branch_directory,
    project_real_execution_argv,
    promote_project_candidate,
    transform_trusted_molecular_geometry,
)
from chemsmart.agent.execution_envelope import BoundedExecutionEnvelopeV1
from chemsmart.agent.guides import guide_for_tool
from chemsmart.agent.identity import (
    ApprovedMolecularIdentityV1,
    refuse_impossible_electronic_state,
)
from chemsmart.agent.inspection import (
    GeneratedArtifactInspectionReceiptV1,
    inspect_generated_artifact,
)
from chemsmart.agent.knowledge import (
    FunctionalEquivalenceReceiptV1,
    ProgramSubstitutionApprovalV1,
    ProgramSubstitutionReceiptV1,
    ScientificClaimEvidenceV1,
    assess_typed_program_substitution,
    build_program_substitution_request,
)
from chemsmart.agent.postprocessing import (
    derive_trusted_thermochemistry,
    evaluate_typed_quantity_expression,
    extract_trusted_result_quantities,
)
from chemsmart.agent.preflight import (
    ProgramNodePreflightReceiptV1,
    ProgramValidatorReceiptV1,
    build_program_node_preflight_request,
    evaluate_program_node_preflight,
    validator_receipt_from_safe_preview,
)
from chemsmart.agent.preview import SafePreviewReceiptV1, execute_safe_preview
from chemsmart.agent.program_verifiers import build_preview_expectation
from chemsmart.agent.projects import (
    ProjectDocumentV1,
    ProjectRenderReceiptV1,
    ProjectValidationReceiptV1,
    PySCFFunctionalResolutionReceiptV1,
    project_document,
    project_scientific_materializations,
    project_section_application_observation,
    read_project_yaml,
    render_project_yaml,
    validate_project_yaml,
)
from chemsmart.agent.report_format import (
    ABSENCES_HEADING,
    CLAIM_RECORD_LABEL,
    CLAIMS_HEADING,
    COMPLETION_RECEIPT_LABEL,
    CONDITIONS_HEADING,
    DECISION_SECTIONS,
    EVIDENCE_COLUMN,
    FINDINGS_HEADING,
    HOST_REPORT_TITLE,
    LITERATURE_CONSTANTS_HEADING,
    NO_DECISION_PREFIX,
    PARTIAL_STATUS_PREFIX,
    PREDICTIONS_HEADING,
    RECOVERY_PREFIX,
    SOURCE_RECEIPT_COLUMN,
    SURVIVING_HEADING,
    THERMO_CONDITIONS_HEADING,
    TOOLCHAIN_PLAN_LABEL,
    VERDICTS_HEADING,
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
from chemsmart.agent.skills import resolve_skill
from chemsmart.agent.terminal_states import (
    GEOMETRY_SEARCH_JOBTYPES,
    STATIONARY_POINT_PROMISES,
    consequential_imaginary_mode_count,
    expected_imaginary_mode_count,
    stationary_point_order_finding,
)
from chemsmart.agent.tool_specs import (
    REGISTRY_PRODUCERS,
    AgentToolSurfaceV1,
    build_approved_execution_tool_surface,
    build_command_compiled_tool_surface,
)
from chemsmart.agent.workflow_context import (
    project_workflow_context,
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
    build_command_workflow_draft,
    build_materialized_workflow,
    build_scientific_workflow_plan,
)
from chemsmart.analysis.literature_constants import (
    UnknownLiteratureConstantError,
    literature_constant,
)
from chemsmart.analysis.quantity_expressions import (
    QuantityExpressionError,
    QuantityExpressionRequestV1,
    canonical_unit_for_dimension,
    convert_normalized_value,
    expression_node_from_plan,
    normalize_numeric_value,
    quantity_expression_receipt_from_record,
    unit_dimension,
)
from chemsmart.analysis.result_quantities import (
    QuantityExtractionReceiptV1,
    QuantitySelectorV1,
    QuantityValueV1,
    ThermochemistryReceiptV1,
    canonical_extraction_receipt_body,
    canonical_thermochemistry_quantity,
    make_quantity_value,
    quantity_extraction_receipt_from_record,
    thermochemistry_receipt_from_record,
)
from chemsmart.analysis.result_readers import (
    reader_for,
    registered_reader_programs,
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
    #: The driven or held internal coordinates of this node, carried so the
    #: approved binding can record them. Without that the executor has nothing
    #: to rebuild a scan's range from.
    internal_coordinates: Mapping[str, Any] | None = None
    #: The anomaly this node investigates, when it is an excursion.
    excursion: str = ""
    #: The goal's original bound geometry behind a reached or displaced
    #: input; None when the input is the goal's own start.
    root_artifact: TrustedArtifactRefV1 | None = None


def _undeferrable_producer_finding(
    waiting_producers: Sequence[Mapping[str, Any]],
) -> dict[str, str]:
    """Explain a wait that no amount of model effort can end.

    A consumer waiting on an ``opt`` or ``ts`` geometry is deferrable: that
    producer ends at one stationary structure, so the stage can sit inside the
    same approval and take it when it exists. A consumer waiting on a relaxed
    scan cannot, because a scan ends at a surface and which point to carry
    forward is a scientific judgement the surface has to inform.

    Told only to "materialize the declared workflow inputs", a session tries to
    do the impossible and its node blocks approval for ever with no reason
    given. Naming the real constraint is the difference between a dead end and
    a decision: run the scan, read its surface, choose a point, and plan the
    next stage against that geometry -- which is a changed molecular input and
    therefore a new workflow with its own review.
    """

    undeferrable = tuple(
        item
        for item in waiting_producers
        if not item.get("deferrable_within_one_approval", True)
    )
    if not undeferrable:
        return {"next_action": "materialize the declared workflow inputs"}
    names = ", ".join(
        f"{item.get('producer_node_id')} ({item.get('producer_stage')})"
        for item in undeferrable
    )
    return {
        "finding": (
            f"this node waits on {names}, which does not end at a single "
            "stationary geometry; which structure to carry forward is a "
            "scientific choice the computed result has to inform, so it "
            "cannot be settled inside this approval"
        ),
        "next_action": (
            "either retain this stage as declared non-executable intent with "
            "its reason, or drop it from this workflow and plan it once the "
            "producer's result exists and a structure has been chosen from it"
        ),
    }


def _node_coordinates(node) -> dict[str, str]:
    """Render a planning node's internal coordinates into program options."""

    return native_coordinate_options(
        node.program, getattr(node, "internal_coordinates", None)
    )


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
    #: Host-detected surprises, each a signal id with the numbers that
    #: tripped it. They ride beside the verdict and never enter it: the
    #: validation receipt is built from observations and findings alone.
    anomalies: tuple[Mapping[str, Any], ...] = ()

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
    candidates: list[tuple[TrustedArtifactRefV1, Mapping[str, Any] | None]] = (
        []
    )
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


def _spin_square_observation(output: Any, multiplicity: Any) -> dict[str, Any]:
    """<S^2> as the program printed it beside what the bound state implies.

    An observation, deliberately not a finding: contamination is a fact
    about the method a scientist weighs, and the analysis chain's own
    predicates are where a threshold belongs.
    """

    # The readers expose the table as a property; calling it raised a
    # swallowed TypeError, so this observation had been empty on every
    # ORCA and Gaussian result since it was written (audit, 2026-09-03).
    try:
        history = getattr(output, "spin_square_history", None)
        if callable(history):
            history = history()
        history = tuple(history or ())
        if not history:
            return {}
        last = history[-1]
        if isinstance(last, Mapping):
            # Gaussian records the SCF value before and after
            # annihilation; the wavefunction's own value is the
            # diagnostic, the annihilated one is a projected estimate.
            last = last.get("before_annihilation")
        observed = float(last)
    except Exception:  # noqa: BLE001 - a reader without the table
        return {}
    record: dict[str, Any] = {"spin_square_observed": observed}
    try:
        spin = (int(multiplicity) - 1) / 2.0
    except (TypeError, ValueError):
        return record
    expected = spin * (spin + 1.0)
    record["spin_square_expected"] = expected
    record["spin_square_deviation"] = observed - expected
    return record


def _output_artifact_kind(program: str, path: Path) -> str:
    """The typed kind of one file an engine left behind.

    xTB writes g98.out beside its own log -- a Gaussian-98-style
    frequency table for other tools -- and registering it as a second
    xtb_output made a live hess node's analysis refuse to choose
    between two logs. It is a sidecar, and reads as one.
    """

    suffix = path.suffix.lower()
    if suffix == ".h5":
        return "pyscf_hdf5"
    if suffix == ".xyz":
        return "geometry_xyz"
    if program == "orca" and suffix == ".hess":
        return "orca_hessian"
    if suffix == ".json":
        return "json"
    if program == "xtb" and suffix == ".out":
        if path.name.lower() == "g98.out":
            return "program_output"
        return "xtb_output"
    if program == "orca" and suffix == ".out":
        return "orca_output"
    if program == "gaussian" and suffix in {".log", ".out"}:
        return "gaussian_output"
    return "program_output"


def _xtb_log_frequencies(logs: tuple[Any, ...]) -> tuple[float, ...]:
    """The frequencies an xTB log printed, read by the same parser the
    typed layer uses; empty when the run printed none or the log is
    unreadable, which the rule treats as no claim. Accepts registered
    artifacts (rehashed before reading) or bare paths."""

    from chemsmart.io.xtb.file import XTBMainOut

    for artifact in logs:
        try:
            path = (
                Path(artifact)
                if isinstance(artifact, (str, Path))
                else _current_artifact_path(artifact, field_name="xTB log")
            )
            values = XTBMainOut(str(path)).vibrational_frequencies
        except Exception:  # noqa: BLE001 - an unreadable log makes no claim
            continue
        if values:
            try:
                return tuple(float(item) for item in values)
            except (TypeError, ValueError):
                return ()
    return ()


def _scan_boundary_sensor(
    profile: Sequence[Mapping[str, float]],
) -> dict[str, Any] | None:
    """An extremum of a scanned surface that sits on the grid's edge.

    Four of four scan extrema in one window sat on a boundary: two
    product-side ring-opening scans rose monotonically to their last
    point and two hydrogen-transfer scans reached their edge, and every
    "barrier position" the extremum operation returned was the end of
    the grid. One session found it by hand four hours later; another
    delivered the gap between two boundary points as an answer (NOVEL-3
    po3 and po1, 2026-09-05). An observation with standing, never a
    refusal: a boundary minimum is the right answer to a dissociation
    curve, and only the scientist knows which end was asked for.
    """

    points = tuple(profile)
    if len(points) < 3:
        return None
    energies = [float(point["energy"]) for point in points]
    coordinates = [float(point["coordinate"]) for point in points]
    last = len(points) - 1
    index_max = max(range(len(points)), key=lambda i: energies[i])
    index_min = min(range(len(points)), key=lambda i: energies[i])
    if index_max not in (0, last) and index_min not in (0, last):
        return None
    steps = [b - a for a, b in zip(energies, energies[1:])]
    monotone = all(step > 0 for step in steps) or all(
        step < 0 for step in steps
    )
    kcal = 627.5094740631
    return {
        "signal_id": "scan.extremum_at_grid_boundary",
        "points": len(points),
        "grid_start": coordinates[0],
        "grid_end": coordinates[-1],
        "maximum_index": index_max,
        "maximum_at_boundary": index_max in (0, last),
        "maximum_coordinate": coordinates[index_max],
        "minimum_index": index_min,
        "minimum_at_boundary": index_min in (0, last),
        "minimum_coordinate": coordinates[index_min],
        "monotone": monotone,
        "span_kcal_mol": round((max(energies) - min(energies)) * kcal, 3),
        "last_step_kcal_mol": round(steps[-1] * kcal, 3),
        "first_step_kcal_mol": round(steps[0] * kcal, 3),
    }


def _basin_sensor_inputs(
    input_artifact: TrustedArtifactRefV1 | None, output: Any, jobtype: str
) -> dict[str, Any]:
    """How far an optimisation walked from the structure it was given.

    A converged minimum is a valid minimum whatever basin it lies in,
    and the stationary-point rule types the order of a point, never
    which structure it is: a planar cyclohexane relaxed to the
    twist-boat in two live arms and both delivered its energy as "the
    minimum" with budget unspent (E2 window, 2026-09-03). The host
    measures the walk -- heavy-atom RMSD after Kabsch alignment and
    any bond made or broken -- and records the numbers; what the walk
    means is the scientist's.
    """

    if jobtype not in GEOMETRY_SEARCH_JOBTYPES:
        return {}
    symbols, positions = _pyscf_input_geometry(input_artifact)
    if not symbols:
        return {}
    try:
        import numpy as np

        from chemsmart.agent.execution import _molecule_graph
        from chemsmart.io.molecules.structure import Molecule
        from chemsmart.utils.utils import kabsch_align

        molecule = output.molecule
        out_symbols = tuple(str(item) for item in molecule.chemical_symbols)
        out_positions = np.asarray(molecule.positions, dtype=float)
    except Exception:
        return {}
    if out_symbols != tuple(symbols) or len(out_positions) != len(symbols):
        return {}
    heavy = [index for index, symbol in enumerate(symbols) if symbol != "H"]
    inputs: dict[str, Any] = {}
    if len(heavy) < SENSOR_HEAVY_ATOM_FLOOR:
        inputs["heavy_atom_rmsd_floor_applied"] = True
        inputs["heavy_atom_count"] = len(heavy)
    if len(heavy) >= SENSOR_HEAVY_ATOM_FLOOR:
        before = np.asarray(positions, dtype=float)[heavy]
        after = out_positions[heavy]
        *_rest, rmsd = kabsch_align(before, after)
        inputs["heavy_atom_rmsd_angstrom"] = float(f"{float(rmsd):.4f}")
    try:
        before_graph = _molecule_graph(
            Molecule(
                symbols=list(symbols),
                positions=[list(row) for row in positions],
            )
        )
        after_graph = _molecule_graph(molecule)
        before_edges = {frozenset(edge) for edge in before_graph.edges}
        after_edges = {frozenset(edge) for edge in after_graph.edges}
        inputs["bonds_made"] = sorted(
            sorted(edge) for edge in after_edges - before_edges
        )
        inputs["bonds_broken"] = sorted(
            sorted(edge) for edge in before_edges - after_edges
        )
    except Exception:
        return inputs
    return inputs


#: A Kabsch heavy-atom RMSD needs three heavy atoms to mean anything:
#: below that, two structures align exactly by construction and the
#: number is zero whatever the molecule did. The basin sensor and the
#: same-structure sensor both stop there, and both used to stop in
#: silence -- which is how a digest join that missed every handoff
#: consumer stayed hidden until the first molecule with three heavy
#: atoms arrived (PySCF round 2, 2026-09-13). The floor is declared, and
#: reaching it is a recorded fact rather than an absence.
SENSOR_HEAVY_ATOM_FLOOR = 3


def _same_structure_observations(
    receipts: Mapping[str, Any],
    node_id: str,
    output: Any,
    input_sha256: str = "",
    output_sha256s: Sequence[str] = (),
    handoffs: Mapping[str, Any] | None = None,
) -> tuple[dict[str, Any], ...]:
    """Whether this result is the same structure as one already validated.

    A session that optimises twice and gets one answer reads the
    agreement as convergence. Two live deliveries called two
    twist-boats "degenerate chairs" and one dismissed a contradicting
    1.85 kcal/mol as an artifact (E4', 2026-09-03), and nothing in the
    host had compared the two structures to each other: the basin
    sensor measures a node against its own input, never against a
    sibling result. This measures it, and says nothing about what it
    means -- two indistinguishable results may be an honest repeat, a
    lost perturbation, or a defect, and only the scientist decides.

    It fires only for INDEPENDENT starts that converged on one
    structure. A single point on an optimisation's own geometry, and
    two siblings launched from the same geometry, are one structure by
    construction; recording those would bury the real observation
    under the ordinary shape of a workflow.

    That exclusion is read from the handoff records the host keeps
    (``handoffs``, keyed by consumer node), never inferred from digests:
    a handoff writes a fresh geometry file the producer's receipt does
    not list, so the digest join missed every fixed-geometry consumer of
    a validated handoff and the sensor's own three-heavy-atom floor hid
    it until the first molecule large enough arrived -- acrolein's td on
    its own optimisation and formic acid's three single points on their
    CCSD geometry were each recorded as a surprise (PySCF round 2,
    2026-09-13).
    """

    if output is None:
        return ()
    try:
        import numpy as np

        from chemsmart.analysis.result_readers import reader_for
        from chemsmart.utils.utils import kabsch_align

        molecule = output.molecule
        symbols = tuple(str(item) for item in molecule.chemical_symbols)
        positions = np.asarray(molecule.positions, dtype=float)
        energy = output.final_energy
    except Exception:  # noqa: BLE001 - a reader without a geometry
        return ()
    heavy = [index for index, symbol in enumerate(symbols) if symbol != "H"]
    if len(heavy) < SENSOR_HEAVY_ATOM_FLOOR:
        # Not silence: the block says the comparison was not made and
        # why, so a reader can tell "no sibling matched" from "no
        # comparison was possible".
        return (
            # ``signal_id`` and not some other word for the same thing:
            # every observation this function returns is read by one
            # consumer, which builds an anomaly observation from it and
            # asks for that key by name. Naming the floor block
            # differently from its own sibling below left the two halves
            # disagreeing while the suite stayed green, because the test
            # asserted the producer's spelling instead of driving the
            # consumer -- and the first goal to validate a node with
            # fewer than three heavy atoms died of a KeyError after its
            # engine had already run, so a finished calculation was
            # typed interrupted_mid_engine and the next node never
            # launched. Water has one heavy atom.
            {
                "signal_id": "geometry.same_structure_comparison_not_made",
                "heavy_atom_rmsd_floor_applied": True,
                "heavy_atom_count": len(heavy),
                "node_id": str(node_id),
            },
        )
    # The nodes this one is one structure with by construction: the
    # producer whose geometry it consumed, every sibling that consumed
    # the same producer's geometry, and every consumer of its own.
    joined: set[str] = set()
    for consumer, handoff in (handoffs or {}).items():
        producer = str(getattr(handoff, "producer_node_id", "") or "")
        if str(consumer) == node_id and producer:
            joined.add(producer)
        if producer == node_id:
            joined.add(str(consumer))
    own_producer = str(
        getattr((handoffs or {}).get(node_id), "producer_node_id", "") or ""
    )
    if own_producer:
        for consumer, handoff in (handoffs or {}).items():
            if str(getattr(handoff, "producer_node_id", "")) == own_producer:
                joined.add(str(consumer))
    found: list[dict[str, Any]] = []
    for receipt in receipts.values():
        other_id = str(getattr(receipt, "node_id", "") or "")
        if not other_id or other_id == node_id or other_id in joined:
            continue
        if str(getattr(receipt, "state", "")) != "valid":
            continue
        other_input = str(getattr(receipt, "input_artifact_sha256", "") or "")
        other_outputs = {
            str(getattr(item, "sha256", "") or "")
            for item in getattr(receipt, "output_artifacts", ()) or ()
        }
        if input_sha256 and (
            input_sha256 in other_outputs or input_sha256 == other_input
        ):
            continue
        if other_input and other_input in set(output_sha256s):
            continue
        for artifact in getattr(receipt, "output_artifacts", ()) or ():
            try:
                other = reader_for(
                    str(getattr(receipt, "program", ""))
                ).open_output(str(artifact.path))
                other_molecule = other.molecule
                other_symbols = tuple(
                    str(item) for item in other_molecule.chemical_symbols
                )
                if other_symbols != symbols:
                    continue
                other_positions = np.asarray(
                    other_molecule.positions, dtype=float
                )
                *_rest, rmsd = kabsch_align(
                    other_positions[heavy], positions[heavy]
                )
            except Exception:  # noqa: BLE001 - not a readable geometry
                continue
            if float(rmsd) >= 0.10:
                break
            record = {
                "signal_id": "geometry.results_indistinguishable",
                "other_node_id": other_id,
                "heavy_atom_rmsd_angstrom": float(f"{float(rmsd):.4f}"),
            }
            try:
                gap = (float(energy) - float(other.final_energy)) * 627.5095
                record["energy_difference_kcal_mol"] = float(f"{gap:.4f}")
            except (TypeError, ValueError):
                pass
            found.append(record)
            break
    return tuple(found)


def _imaginary_mode_sensor_inputs(
    output: Any, frequencies: tuple[float, ...]
) -> dict[str, Any]:
    """What the lowest imaginary mode is, for the anomaly sensor.

    The magnitude and the heavy atoms sharing the mode's motion separate
    a structural saddle (an inversion, a symmetry breaking, a hidden
    reaction coordinate) from a rotor: an inversion mode spreads over
    several heavy atoms, a methyl rotor lives on three hydrogens. Absent
    modes leave the spread unknown, never guessed.
    """

    imaginary = [
        (index, float(value))
        for index, value in enumerate(frequencies)
        if math.isfinite(float(value)) and float(value) < -20.0
    ]
    if not imaginary:
        return {}
    index, lowest = min(imaginary, key=lambda item: item[1])
    inputs: dict[str, Any] = {"lowest_imaginary_cm1": float(f"{lowest:.2f}")}
    try:
        from chemsmart.analysis.result_readers import (
            _vibrational_mode_atom_participation,
        )

        shares = _vibrational_mode_atom_participation(output)[index]
        symbols = list(getattr(output.molecule, "chemical_symbols", ()))
    except Exception:
        return inputs
    heavy_shares = [
        (atom + 1, float(share))
        for atom, share in enumerate(shares)
        if atom < len(symbols) and symbols[atom] != "H"
    ]
    inputs["participating_heavy_atoms"] = [
        atom for atom, share in heavy_shares if share >= 0.05
    ]
    inputs["heavy_atom_share"] = float(
        f"{sum(share for _atom, share in heavy_shares):.3f}"
    )
    return inputs


def _observed_spin_deviation(
    observation: Mapping[str, Any], program: str
) -> float | None:
    """The ⟨S²⟩ deviation a program's observation carries, or None."""

    block = observation.get(program)
    if not isinstance(block, Mapping):
        return None
    if program == "gaussian":
        rows = tuple(block.get("outputs") or ())
        if len(rows) != 1 or not isinstance(rows[0], Mapping):
            return None
        block = rows[0]
    value = block.get("spin_square_deviation")
    try:
        return float(value) if value is not None else None
    except (TypeError, ValueError):
        return None


def _observed_imaginary_mode_count(
    observation: Mapping[str, Any], program: str
) -> int | None:
    """The consequential imaginary-mode count a program's observation
    carries, or None when the run printed no frequencies."""

    block = observation.get(program)
    if not isinstance(block, Mapping):
        return None
    if program == "gaussian":
        rows = tuple(block.get("outputs") or ())
        if len(rows) != 1 or not isinstance(rows[0], Mapping):
            return None
        block = rows[0]
    value = block.get("consequential_imaginary_mode_count")
    return int(value) if isinstance(value, int) else None


#: A bond-forming saddle's imaginary mode is hundreds of wavenumbers; a
#: mode inside this band is an intermolecular or torsional motion that
#: happens to fall past the 20 cm-1 noise convention.
SOFT_IMAGINARY_MODE_BAND_CM1 = 50.0

#: geomeTRIC's ``convergence_gmax`` (Eh/Bohr), the optimiser's own word for
#: "the gradient is zero".  A Hessian computed at a geometry whose gradient
#: exceeds it is an observation with standing: PySCF's harmonic analysis
#: projects rotations out, so the projected spectrum at a non-stationary
#: point can be entirely real (HF/STO-3G water at 1.10 A / 120 deg: 2015,
#: 2868, 3225 cm-1; the archived stretched-water Hessian: three real modes
#: at max|g| = 0.0185 Eh/Bohr), and zero imaginary modes then proves
#: nothing about stationarity.  Never a refusal: a Hessian off a stationary
#: point is a legitimate thing to ask for.  Owner ruling, 2026-09-12.
HESS_STATIONARITY_GRADIENT_EH_PER_BOHR = 4.5e-4


def _gradient_anomaly(gradient: float) -> dict[str, Any]:
    """The projected spectrum can be all-real at a geometry that is not
    stationary; the gradient the Hessian stage recorded says how far from
    stationary it was, against the optimiser's own criterion.  An
    observation with standing, never a verdict."""

    return {
        "signal_id": "stationary_point.gradient_above_optimizer_criterion",
        "max_abs_gradient_eh_per_bohr": float(f"{gradient:.6g}"),
        "optimizer_criterion_eh_per_bohr": (
            HESS_STATIONARITY_GRADIENT_EH_PER_BOHR
        ),
        "policy_id": "hess_stationarity_gradient",
    }


def _neutral_sensor_facts(
    *,
    program: str,
    jobtype: str,
    multiplicity: Any,
    output_artifacts: Sequence[TrustedArtifactRefV1],
    expected_input_artifact: TrustedArtifactRefV1 | None,
    expected_root_artifact: TrustedArtifactRefV1 | None,
) -> tuple[dict[str, Any], dict[str, Any]]:
    """The facts every host sensor consumes, read once through the reader.

    The stationary-point rule, the spin observation, the basin walk and
    the imaginary-mode inputs are program-neutral where they are computed
    and were fed per program: the ORCA branch wrote all of them, xTB and
    Gaussian wrote the count, PySCF wrote none, so the neutral word never
    reached PySCF and the basin sensor never reached anyone but ORCA.
    This reads the same facts from whatever the program's reader opened,
    and a branch that already wrote a key keeps its own value.

    Returns ``(block_facts, sensor_inputs)``; both empty when the primary
    result cannot be opened, which is an absence the findings already
    name.
    """

    from chemsmart.analysis.result_readers import reader_for

    reader = reader_for(program)
    if reader is None:
        return {}, {}
    primary = next(
        (
            artifact
            for artifact in output_artifacts
            if artifact.kind == reader.artifact_kind
        ),
        None,
    )
    if primary is None:
        return {}, {}
    try:
        output = reader.open_output(
            _current_artifact_path(primary, field_name="primary result")
        )
    except Exception:  # noqa: BLE001 - unreadable is an absence here
        return {}, {}
    block: dict[str, Any] = {}
    inputs: dict[str, Any] = {}
    # The electronic surface this result is on, where its reader can say
    # it. Every organ that asks whether two results describe the same
    # thing -- the characterisation join, the same-structure sensor --
    # reads this one fact rather than each deciding for itself, and a
    # reader that cannot say leaves it absent rather than guessing.
    try:
        surface = reader.surface_for_output(output)
    except Exception:  # noqa: BLE001 - a reader without an identity
        surface = None
    if surface:
        from chemsmart.analysis.result_readers import surface_token

        block["surface"] = canonical_data(surface)
        block["surface_id"] = surface_token(surface)
    raw_frequencies = getattr(output, "vibrational_frequencies", None)
    frequencies: tuple[float, ...] = ()
    if raw_frequencies is not None:
        try:
            frequencies = tuple(float(item) for item in raw_frequencies)
        except (TypeError, ValueError):
            frequencies = ()
    # How many modes the program printed, for every program: an
    # optimisation that printed none made no claim about its
    # stationary point, and the projections that word a delivered
    # number read this count (the ORCA branch wrote its own; PySCF
    # and xTB rows carried None and read as zero).
    block["vibrational_mode_count"] = len(frequencies)
    count = consequential_imaginary_mode_count(frequencies)
    if count is not None:
        block["consequential_imaginary_mode_count"] = count
        block["imaginary_frequencies_cm1"] = [
            value for value in frequencies if value <= -20.0
        ]
    for name in ("charge", "multiplicity"):
        value = getattr(output, name, None)
        if isinstance(value, int) and not isinstance(value, bool):
            block[name] = value
    try:
        observed = float(reader.read(output, "spin_square")[0])
        target = float(reader.read(output, "spin_square_target")[0])
        block["spin_square_observed"] = observed
        block["spin_square_expected"] = target
        block["spin_square_deviation"] = observed - target
    except Exception:  # noqa: BLE001 - a closed shell or no diagnostic
        pass
    if frequencies:
        inputs.update(_imaginary_mode_sensor_inputs(output, frequencies))
    if jobtype in GEOMETRY_SEARCH_JOBTYPES:
        try:
            inputs["basin"] = _basin_sensor_inputs(
                expected_input_artifact, output, jobtype
            )
            if expected_root_artifact is not None and (
                expected_input_artifact is None
                or expected_root_artifact.sha256
                != expected_input_artifact.sha256
            ):
                inputs["basin_root"] = {
                    **_basin_sensor_inputs(
                        expected_root_artifact, output, jobtype
                    ),
                    "reference": "goal_root",
                    "reference_artifact_id": expected_root_artifact.artifact_id,
                }
        except Exception:  # noqa: BLE001 - no readable input geometry
            pass
    forces = getattr(output, "forces", None)
    if (
        jobtype in {"hess", "freq"}
        and forces is not None
        and getattr(output, "forces_unit", None) == "Eh/Bohr"
    ):
        try:
            import numpy as np

            gradient = float(np.max(np.abs(np.asarray(forces, dtype=float))))
            block["max_abs_gradient_eh_per_bohr"] = gradient
            inputs["stationarity_gradient"] = gradient
        except (TypeError, ValueError):
            pass
    # The response stage's own numbers, as facts with no threshold and no
    # signal: which root an optimisation followed, how far it ended from
    # the ground state and from its neighbour, and how many requested
    # roots the program's positive-eigenvalue filter dropped. A root is
    # an index, never a state identity; a small gap is what a single-
    # reference response cannot describe, and the session reads it here
    # before it calls the delivered geometry "the S1 minimum". A numeric
    # policy is registered only after a sealed case shows one is needed.
    block.update(_excited_state_sensor_facts(output))
    return block, inputs


def _plain_reader_value(value: Any) -> Any:
    """Return a reader's value in a plain Python container.

    Nothing makes two readers agree on the container: the HDF5 path
    serves numpy arrays and numpy scalars where a log parser serves
    lists and floats for the same selector. Both idioms a sensor
    naturally reaches for are shape-sensitive -- ``if values:`` raises
    on an array of more than one element, and an array is not a list or
    a tuple, so a membership branch silently takes the scalar path -- so
    the container is removed once, here, and every fact below is written
    against plain Python.
    """

    tolist = getattr(value, "tolist", None)
    return tolist() if callable(tolist) else value


def _excited_state_sensor_facts(output: Any) -> dict[str, Any]:
    """Gap and filter facts of a response stage, read through the reader."""

    facts: dict[str, Any] = {}
    td_stage = getattr(output, "td_stage", None)
    if isinstance(td_stage, Mapping):
        for key, name in (
            ("nstates_requested", "excited_state_roots_requested"),
            ("nstates_obtained", "excited_state_roots_obtained"),
            ("roots_filtered", "excited_state_roots_filtered"),
            ("unconverged_roots", "excited_state_unconverged_roots"),
        ):
            value = _plain_reader_value(td_stage.get(key))
            if value is not None:
                facts[name] = (
                    [int(item) for item in value]
                    if isinstance(value, (list, tuple))
                    else int(value)
                )
        excitations = _plain_reader_value(
            getattr(output, "excitation_energies", None)
        )
        if excitations is not None and len(excitations):
            try:
                facts["excited_state_lowest_root_ev"] = (
                    float(excitations[0]) * 27.211386245988
                )
            except (TypeError, ValueError, IndexError):
                pass
    record = getattr(output, "excited_state_record", None)
    if isinstance(record, Mapping):
        for key, name in (
            ("root", "excited_state_followed_root"),
            (
                "root_gap_to_ground_end_ev",
                "excited_state_root_gap_to_ground_ev",
            ),
            (
                "root_gap_to_neighbour_end_ev",
                "excited_state_root_gap_to_neighbour_ev",
            ),
            (
                "followed_root_converged_end",
                "excited_state_followed_root_converged",
            ),
            (
                "final_gradient_max_eh_per_bohr",
                "excited_state_final_gradient_max",
            ),
        ):
            value = record.get(key)
            if value is None:
                continue
            if key == "root":
                facts[name] = int(value)
            elif isinstance(value, bool):
                facts[name] = bool(value)
            else:
                facts[name] = float(value)
    return facts


#: Every function that turns a reader's numbers into sensor facts,
#: enumerated by the host rather than by whichever ones a test author
#: remembered. A reader may serve a Python list where another serves a
#: numpy array for the same quantity, and the two must produce the same
#: facts: ``if values:`` raises on an array of more than one element and
#: ``isinstance(value, (list, tuple))`` is false for one, so a sensor
#: that reads a list correctly can be silently wrong or fatal on the
#: same numbers from another program. The shape-invariance test drives
#: this tuple and fails when an entry has no driver.
SENSOR_FACT_FUNCTIONS = (
    _scan_boundary_sensor,
    _basin_sensor_inputs,
    _same_structure_observations,
    _imaginary_mode_sensor_inputs,
    _gradient_anomaly,
    _excited_state_sensor_facts,
)


def _observed_soft_imaginary_mode(
    observation: Mapping[str, Any], program: str, *, jobtype: str
) -> float | None:
    """The one imaginary mode of a validated saddle when it lies inside
    the soft band, else None.

    REACH-1 po3 cycle 4 (2026-09-06): a transition-state search on a
    dual-contact guess relaxed to a van der Waals complex with the
    azide 2.9-3.3 A from an unreacted alkyne and one imaginary mode at
    -22.8 cm-1, 2.8 cm-1 past the noise convention; the order rule
    certified it and the word was `validated`. A rule at a threshold
    certifies noise on the far side of the threshold, so the number is
    recorded as an observation with standing, never as a verdict: the
    session decides what a 22.8 cm-1 mode means.
    """

    if expected_imaginary_mode_count(jobtype) != 1:
        return None
    block = observation.get(program)
    if not isinstance(block, Mapping):
        return None
    if program == "gaussian":
        rows = tuple(block.get("outputs") or ())
        if len(rows) != 1 or not isinstance(rows[0], Mapping):
            return None
        block = rows[0]
    if block.get("consequential_imaginary_mode_count") != 1:
        return None
    modes = tuple(
        float(value)
        for value in (block.get("imaginary_frequencies_cm1") or ())
        if float(value) <= -20.0
    )
    if len(modes) != 1:
        return None
    return modes[0] if abs(modes[0]) < SOFT_IMAGINARY_MODE_BAND_CM1 else None


#: ORCA caps a geometry optimisation at 3N iterations by default. A
#: 21-atom Fe(II) ammine complex hit that cap in two sealed windows while
#: geom_maxiter sat, rendered and unused, on the project tool: the lever
#: was visible and nothing said this system needed it (NOVEL-1/2 ino1,
#: 2026-09-04). The host knows the program, the job type, the settings
#: and the atom count at compile time, so it says so there.
_ORCA_GEOMETRY_CAP_JOBTYPES = GEOMETRY_SEARCH_JOBTYPES


def _artifact_id_taken(
    requested: Any, taken: list[str]
) -> RoutedContractError:
    """The refusal a re-used artifact id meets, stated once for every
    surface that mints one. A live run collided five times in a row on
    the bare message; a session under NOVEL-3 met it on a re-derived
    species and re-promoted project alike."""

    return RoutedContractError(
        gate="artifact.id_is_unused",
        invariant=(
            "an artifact id stands for exactly one set of bytes for the "
            "life of the workspace."
        ),
        diagnosis=(
            f"artifact ID {requested!r} is already registered; taken ids: "
            f"{taken}."
        ),
        route=(
            "choose an id not in the taken list; a re-derived or "
            "re-promoted artifact takes a fresh id and the earlier one "
            "stays evidence."
        ),
    )


def _is_restatement_of(written: float, exact: float) -> bool:
    """Whether a written tolerance is the exact one, to its own digits.

    A supersession corrects how an observable is named or measured and
    never how good the answer has to be, so a replacement's tolerance
    must be the retired one carried across. Exact equality is the wrong
    test: 0.05 eV is 1.1530 kcal/mol and a scientist writes 1.15, which
    is the same requirement written to fewer digits. A fixed percentage
    band would be a magic number nobody could defend.

    So the comparison is made at the precision the session actually
    wrote: round the host's conversion to the written value's own
    significant figures and require equality. 1.15 for 1.1530 is a
    restatement; 200 for 2 is a new requirement, which is how a
    2 kJ/mol obligation was relabelled as 200 kJ/mol under the same
    meaning, the same unit and the same quoted source. Writing fewer
    digits may loosen the number by less than one unit in its last
    place, which is what rounding is; anything tighter is a stricter
    obligation the session took on itself and is never an escape.
    """

    if written <= exact:
        return True
    return _reads_as(written, exact)


def _reads_as(written: float, exact: float) -> bool:
    """Whether a written number is an exact one, to its own digits.

    Rounding the host's number to the digits the session wrote is the
    comparison: 1.15 reads 1.1530 and 200 does not read 2. Unlike a
    restated *tolerance*, a stated *uncertainty* is not free to be
    smaller than its source -- understating what you measured is the
    escape this check exists to close -- so the uncertainty comparison
    uses this in both directions and never the looser rule above.
    """

    if exact == 0.0 or written == 0.0:
        return written == exact
    digits = (
        len(
            f"{abs(written):.12g}".replace(".", "")
            .replace("-", "")
            .lstrip("0")
        )
        or 1
    )
    return float(f"{abs(exact):.{digits}g}") == float(f"{abs(written):.12g}")


def _restate_display_value(
    value: float, from_unit: str, to_unit: str
) -> float | None:
    """A displayed number in another unit the host can reach.

    Within one dimension through the unit table; between a wavenumber
    and a molar energy through h*c*N_A. None where no conversion the
    host owns applies.
    """

    from chemsmart.analysis.quantity_expressions import (
        ENERGY,
        FREQUENCY,
        QuantityExpressionError,
        _unit_spec,
    )

    try:
        from_dimension, _from_canonical, from_scale = _unit_spec(from_unit)
        to_dimension, _to_canonical, to_scale = _unit_spec(to_unit)
    except QuantityExpressionError:
        return None
    hartree_in_cm1 = 219474.6313705
    canonical = float(value) * from_scale
    if from_dimension == to_dimension:
        return canonical / to_scale
    if from_dimension == FREQUENCY and to_dimension == ENERGY:
        return (canonical / hartree_in_cm1) / to_scale
    if from_dimension == ENERGY and to_dimension == FREQUENCY:
        return (canonical * hartree_in_cm1) / to_scale
    return None


def compile_time_observations(
    *,
    program: str,
    jobtype: str,
    settings: Mapping[str, Any] | Sequence[tuple[str, Any]],
    atom_count: int,
    geometry: Any = None,
) -> tuple[str, ...]:
    """Facts the host can state about a compiled node before it runs.

    Observations, never refusals: each names a default the program will
    apply and the project field that changes it, computed from the plan
    alone; with the input geometry in hand, its symmetry estimate too.
    """

    resolved = (
        dict(settings) if not isinstance(settings, Mapping) else settings
    )
    observations: list[str] = []
    # Every job type that promises a stationary point: the two that
    # search for one and the two that evaluate a Hessian on one. The
    # union this replaced added the geometry-cap set to a hand-written
    # copy of itself plus hess and freq, and named the same four.
    if geometry is not None and jobtype in STATIONARY_POINT_PROMISES:
        try:
            from chemsmart.agent.symmetry import symmetry_observation

            observations.append(
                symmetry_observation(
                    tuple(geometry.chemical_symbols), geometry.positions
                )
            )
        except Exception:
            pass
    if program == "orca" and jobtype in _ORCA_GEOMETRY_CAP_JOBTYPES:
        if resolved.get("geom_maxiter") in (None, "", 0):
            cap = 3 * int(atom_count)
            observations.append(
                "ORCA's default geometry-iteration cap for this "
                f"{int(atom_count)}-atom input is 3N = {cap}; a floppy or "
                "metal-centred system may need more, and geom_maxiter on "
                "the project raises it (opt_convergence changes the "
                "criterion, not the count)"
            )
    return tuple(observations)


def _periodic_degrees(value: float) -> float:
    """A dihedral is periodic: 285 is -75, and a refusal on the range
    costs a call for nothing (REACH-1 ino3, one refusal)."""

    return float(((value + 180.0) % 360.0) - 180.0)


def _geometry_for_observation(artifact: Any) -> Any:
    """The molecule behind a geometry_xyz artifact, or None."""

    if getattr(artifact, "kind", "") != "geometry_xyz":
        return None
    try:
        from chemsmart.io.molecules.structure import Molecule

        return Molecule.from_filepath(str(artifact.path))
    except Exception:
        return None


def promotion_field_observations(
    render: Any, earlier: Sequence[tuple[str, Any]]
) -> tuple[str, ...]:
    """Fields an earlier promoted project of the same program set, in a
    section this project also carries, that this project leaves unset.

    A promoted project carries no field it does not state. A session under
    NOVEL-3 promoted a quartet project with ``geom_maxiter: 150``, was
    refused on another field, re-promoted a ``-v2`` without the cap, and
    two nodes then stopped at ORCA's 3N default (2026-09-05). Nothing
    said the lever had been lost. An observation, never a refusal: a
    dropped field can be deliberate.
    """

    import yaml

    def sections(text: str) -> dict[str, dict[str, Any]]:
        try:
            document = yaml.safe_load(text) or {}
        except yaml.YAMLError:
            return {}
        if not isinstance(document, Mapping):
            return {}
        return {
            str(name): dict(fields)
            for name, fields in document.items()
            if isinstance(fields, Mapping)
        }

    current = sections(str(getattr(render, "rendered_yaml", "") or ""))
    observations: list[str] = []
    for artifact_id, receipt in earlier:
        if getattr(receipt, "program", None) != getattr(
            render, "program", None
        ):
            continue
        previous = sections(str(getattr(receipt, "rendered_yaml", "") or ""))
        dropped = [
            f"{section}.{key} ({value!r})"
            for section, fields in sorted(previous.items())
            if section in current
            for key, value in sorted(fields.items())
            if key not in current[section]
        ]
        if dropped:
            observations.append(
                f"against {artifact_id!r} ({render.program}): "
                f"{', '.join(dropped)} set there and not here; a promoted "
                "project carries no field it does not state, so a "
                "re-promotion that means to keep one restates it"
            )
    return tuple(observations)


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
        run_evidence_root: str | Path | None = None,
        cycle_label: str | None = None,
        execution_resources: ExecutionResourceSpecV1 | None = None,
        workflow_execution_approval: WorkflowExecutionApprovalV1 | None = None,
        frozen_workflow_approval: FrozenWorkflowApprovalV1 | None = None,
        bounded_execution_envelope: BoundedExecutionEnvelopeV1 | None = None,
        engine_calls_remaining: int | None = None,
        excursion_calls_remaining: int | None = None,
        offered_repair_routes: Sequence[str] = (),
        preview_retention_root: Path | None = None,
        input_check_executable: Path | None = None,
        input_check_env: Mapping[str, str] | None = None,
        input_check_cap_seconds: float = 20.0,
        wall_seconds_remaining: float | None = None,
        revisions_remaining: int | None = None,
        goal_delivered_declared_ids: Sequence[str] = (),
        prior_anomaly_observations: Sequence[Mapping[str, Any]] = (),
        approved_environment_identities: tuple[str, ...] = (),
        materialized_workflow: MaterializedWorkflowV1 | None = None,
        approved_requested_observable_declarations: Sequence[
            Mapping[str, Any]
        ] = (),
        approved_scientific_toolchain_plan: (
            ScientificToolchainPlanV1 | None
        ) = None,
        scientific_workflow_plan: ScientificWorkflowPlanV2 | None = None,
        preview_server: str = "",
        execution_server: str = "",
        execution_server_file_sha256: str = "",
        execution_environment: Mapping[str, str] = {},
        execution_environment_remove: tuple[str, ...] = (),
        active_guides: Iterable[str] = (),
        execute_analysis_only_plans: bool = False,
        analysis_only_run_directory: str | Path = "",
        analysis_only_workspace: str | Path = "",
    ) -> None:
        self.event_store = event_store
        #: Under a goal wake with a previous run, a plan with no
        #: calculation node is walked by the host the moment it is
        #: planned (executor.execute_analysis_only_toolchain); outside a
        #: goal a planned chain stays a plan, as before.
        self.execute_analysis_only_plans = bool(execute_analysis_only_plans)
        self.analysis_only_run_directory = (
            Path(analysis_only_run_directory)
            if analysis_only_run_directory
            else None
        )
        self.analysis_only_workspace = (
            Path(analysis_only_workspace) if analysis_only_workspace else None
        )
        self.registry = registry or load_program_capabilities()
        self.live_schema = live_schema or build_live_click_schema()
        preview_overlay = build_command_compiled_preview_overlay(
            self.registry,
            conformance_receipts=component_conformance_receipts,
            live_schema=self.live_schema,
        )
        self.active_guides: set[str] = set(active_guides)
        command_surface = build_command_compiled_tool_surface(
            self.registry, guides=tuple(sorted(self.active_guides))
        )
        execution_surface = build_approved_execution_tool_surface(
            self.registry
        )
        if tool_surface is not None and tool_surface.profile not in {
            command_surface.profile,
            execution_surface.profile,
        }:
            raise ContractError(
                "injected tool surface is not a canonical profile"
            )
        self.surface = tool_surface or command_surface
        self.artifacts = dict(artifacts)
        self.task_spec_sha256s = frozenset(task_spec_sha256s)
        # The session's own restatement of what the task asks for:
        # observable_id -> {unit, dimension, meaning}.  A commitment the
        # completion gate checks by kind and unit, never value.
        #
        # The declarations are made on the planning session's host and the
        # gate runs on the provider-free executor's, which is a different
        # host object built from the approval bundle rather than handed
        # down from the session.  Left unseeded there, the
        # gate and the expectation rows took their empty branch on every
        # execution this product has ever performed, so a declaration was
        # advertised as a commitment and never once checked.  The approved
        # bundle carries them across.
        self.requested_observable_declarations: dict[str, dict[str, Any]] = {
            str(item["observable_id"]): dict(item)
            for item in approved_requested_observable_declarations
            if item.get("observable_id")
        }
        self.approved_environment_identities = tuple(
            approved_environment_identities
        )
        self.approved_workspace = (
            Path(approved_workspace).resolve()
            if approved_workspace is not None
            else None
        )
        # approved_workspace means a different root per surface: the
        # planning host binds its private preview root, the executor
        # binds the user workspace. Recorded runs live only under the
        # user workspace, so run inspection resolves against this
        # field and never against approved_workspace.
        self.run_evidence_root = (
            Path(run_evidence_root).resolve()
            if run_evidence_root is not None
            else None
        )
        #: The cycle this host is executing, as a folder name. A branch
        #: lives under it, so one calculation's whole record sits in one
        #: place a chemist can open. ``None`` keeps the legacy layout, so
        #: every recorded run keeps the path it was written at.
        self.cycle_label = str(cycle_label) if cycle_label else None
        #: The wave this session selected, in the order it selected it.
        #: Declared here rather than created on first use, so the session
        #: that carries it to the dispatcher reads an attribute that
        #: always exists: a `getattr` default would turn a renamed
        #: attribute into an empty wave, which dispatches as a single job
        #: and looks exactly like a session that chose not to select one.
        self.selected_execution_wave: tuple[str, ...] = ()
        self.execution_resources = execution_resources
        self.workflow_execution_approval = workflow_execution_approval
        self.frozen_workflow_approval = frozen_workflow_approval
        self.bounded_execution_envelope = bounded_execution_envelope
        #: What the goal has left of its engine-call budget after earlier
        #: cycles, when a wake carries it; the envelope alone bounds a
        #: first cycle.
        self.excursion_calls_remaining = (
            None
            if excursion_calls_remaining is None
            else int(excursion_calls_remaining)
        )
        #: The repair-menu routes the wake offered this session, keyed
        #: as the menu keys them: the terminal states of the previous
        #: run. A disposition names one of these or is refused.
        self.offered_repair_routes = tuple(
            str(route) for route in offered_repair_routes
        )
        #: Where a safe preview's emitted files are kept by digest, and
        #: the program's own checker the host may run on them. Both come
        #: from the host-owned server profile; neither is model-visible.
        self.preview_retention_root = (
            None
            if preview_retention_root is None
            else Path(preview_retention_root)
        )
        self.input_check_executable = (
            None
            if input_check_executable is None
            else Path(input_check_executable)
        )
        self.input_check_env = (
            None if input_check_env is None else dict(input_check_env)
        )
        self.input_check_cap_seconds = float(input_check_cap_seconds)
        self._input_check_by_node: dict[str, Any] = {}
        self.wall_seconds_remaining = (
            None
            if wall_seconds_remaining is None
            else float(wall_seconds_remaining)
        )
        self.revisions_remaining = (
            None if revisions_remaining is None else int(revisions_remaining)
        )
        self.goal_delivered_declared_ids = frozenset(
            str(item) for item in goal_delivered_declared_ids
        )
        self.verified_unreachable_ids: set[str] = set()
        self.engine_calls_remaining = (
            None
            if engine_calls_remaining is None
            else int(engine_calls_remaining)
        )
        self._bounded_execution_started_at = time.monotonic()
        self.preview_server = str(preview_server)
        self.execution_server = str(execution_server)
        self.execution_server_file_sha256 = str(execution_server_file_sha256)
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
            # The approved executor carries BOTH by design: the approval
            # names what may run and the envelope bounds how much.  The
            # old mutual-exclusion encoded the deleted execution-enabled
            # session plane and fired the moment the envelope started
            # riding the bundle into the executor; the coherence that
            # matters is the resource check below.
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
        #: Durable provenance of those consultations -- id, version, and
        #: digests only, rehydrated from events so a review built later can
        #: display what the session read without re-reading transcripts.
        #: Pure provenance: nothing gates on it and it grants no authority.
        self.consulted_skill_records = {}
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
        # The approved analysis chain, in the executor host. The validation
        # evaluator and the completion gate resolve their sealed plan from
        # here when no planning session ever ran in this process.
        self.approved_scientific_toolchain_plan = (
            approved_scientific_toolchain_plan
        )
        if (
            self.analysis_completion_policy is not None
            and self.analysis_completion_policy.task_spec_sha256
            not in self.task_spec_sha256s
        ):
            raise ContractError(
                "analysis completion policy targets another task spec"
            )
        self.molecular_compositions: dict[
            str, MolecularCompositionReceiptV1
        ] = {}
        self.molecular_derivations: dict[str, MolecularDerivationReceiptV1] = (
            {}
        )
        self.pubchem_geometries: dict[str, PubchemGeometryReceiptV1] = {}
        self.database_extractions: dict[
            str, DatabaseRecordExtractionReceiptV1
        ] = {}
        self.geometry_edits: dict[str, GeometryEditReceiptV1] = {}
        self.mode_displacements: dict[str, Any] = {}
        self.symmetry_breaks: dict[str, SymmetryBreakReceiptV1] = {}
        #: The review the loop built at "host readiness gates passed", or
        #: the refusal it recorded, so the session runner reuses one and
        #: never appends the other after the stream is sealed.
        self.prepared_execution_review: WorkflowExecutionReviewV1 | None = None
        self.execution_review_refusal: dict[str, str] = {}
        self.atom_appends: dict[str, AtomAppendReceiptV1] = {}
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
        self._declared_observable_join_fields = {}
        self._reply_observations: tuple[dict[str, Any], ...] = ()
        #: The current sufficiency assessment of each declared
        #: requirement, by observable id.
        self.requirement_assessments: dict[str, dict[str, Any]] = {}
        self.analysis_completion_receipts: dict[str, Any] = {}
        self.workflow_drafts: dict[str, CommandWorkflowDraftV1] = {}
        self.scientific_toolchain_plans: dict[
            str, ScientificToolchainPlanV1
        ] = {}
        if approved_scientific_toolchain_plan is not None:
            self.scientific_toolchain_plans[
                approved_scientific_toolchain_plan.plan_sha256
            ] = approved_scientific_toolchain_plan
        self._scientific_toolchain_command_results: dict[
            str, dict[str, Any]
        ] = {}
        self._latest_program_workflows: dict[str, _ResolvedProgramWorkflow] = (
            {}
        )
        #: Last accepted scientific plan per workflow, so a later repair can be
        #: checked against the question and not only against its own findings.
        self.scientific_plans: dict[str, Any] = {}
        self.scientific_workflow_plans: dict[str, ScientificWorkflowPlanV2] = (
            {}
        )
        self.materialized_workflows: dict[str, MaterializedWorkflowV1] = {}
        if materialized_workflow is not None:
            self.materialized_workflows[
                materialized_workflow.materialized_sha256
            ] = materialized_workflow
        #: Which registered plan is the one the session is standing on.
        #: The registry is keyed by digest, and a dict re-insertion keeps
        #: a key in its original position -- so after plan A, plan B,
        #: and an amendment back to A, ``values()[-1]`` is still B, the
        #: plan the session abandoned. po3-r17 (2026-09-11) reverted a
        #: change exactly that way: its readiness frontier reported the
        #: restored nodes previewed and approvable, and the execution
        #: review, reading the last-inserted plan, refused the workflow
        #: over the superseded plan's red previews. Reverting a change
        #: is self-correction, which this harness exists to support, so
        #: the current plan is named rather than inferred from ordering.
        self.current_scientific_plan_sha256: str = ""
        if scientific_workflow_plan is not None:
            self.scientific_workflow_plans[
                scientific_workflow_plan.plan_sha256
            ] = scientific_workflow_plan
            self.current_scientific_plan_sha256 = (
                scientific_workflow_plan.plan_sha256
            )
        self.project_promotions: dict[str, ProjectArtifactPromotionV1] = {}
        self.scientific_decisions: dict[str, ScientificDecisionRecordV1] = {}
        self.execution_receipts: dict[str, ProgramExecutionReceiptV1] = {}
        self.result_validation_receipts: dict[
            str, ProgramResultValidationReceiptV1
        ] = {}
        self.anomaly_observations: dict[str, AnomalyObservationV1] = {}
        self.stationary_point_characterisations: dict[str, Any] = {}
        self.reached_geometries: dict[str, Any] = {}
        #: Anomalies earlier cycles recorded, seeded from the wake so a
        #: claim-only cycle can cite and carry them.
        self.prior_anomaly_observations: tuple[Mapping[str, Any], ...] = tuple(
            dict(item) for item in prior_anomaly_observations
        )
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
                values["limitation_output_ids"] = tuple(
                    values.get("limitation_output_ids") or ()
                )
                # The fourth field, added after its three siblings and
                # never threaded here. JSON returns a list and the
                # dataclass compares against a tuple, so a host rebuilt
                # over a stream that recorded any anomaly -- the exact
                # streams worth reconstructing -- raised "anomaly output
                # ids must be non-empty, sorted, unique" before doing
                # anything. A hand-listed field set beside a dataclass
                # that grew.
                values["anomaly_output_ids"] = tuple(
                    values.get("anomaly_output_ids") or ()
                )
                receipt = AnalysisCompletionReceiptV1(
                    **values, receipt_sha256=receipt_sha256
                )
                self.analysis_completion_receipts[receipt_sha256] = receipt
            elif event.kind == EventKind.DOMAIN_SKILL_CONSULTED.value:
                values = {
                    key: str(record.get(key) or "")
                    for key in (
                        "skill_id",
                        "skill_version",
                        "origin",
                        "body_sha256",
                        "document_sha256",
                    )
                }
                require_sha256(
                    values["document_sha256"], "consulted skill digest"
                )
                self.consulted_skill_records[values["document_sha256"]] = (
                    values
                )
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

    #: Every tool the host handles, by the method that handles it. The
    #: capability registry reads this to say a tool is wired.
    TOOL_HANDLERS: Mapping[str, str] = {
        "inspect_program": "_inspect_program",
        "project_yaml": "_project_yaml",
        "compile_command": "_prepare_program_node",
        "inspect_run": "_inspect_run",
        "inspect_program_capability": "_inspect_program_capability",
        "inspect_program_environment": "_inspect_program_environment",
        "assess_program_candidate": "_assess_program_candidate",
        "render_project_yaml": "_render_project_yaml",
        "promote_project_yaml": "_promote_project_yaml",
        "establish_project": "_establish_project",
        "bind_scientific_identity": "_bind_scientific_identity",
        "bind_scan_point_geometry": "_bind_scan_point_geometry",
        "compose_molecular_arrangement": "_compose_molecular_arrangement",
        "derive_molecular_species": "_derive_molecular_species",
        "fetch_pubchem_geometry": "_fetch_pubchem_geometry",
        "edit_molecular_geometry": "_edit_molecular_geometry",
        "break_symmetry": "_break_symmetry",
        "append_molecular_atom": "_append_molecular_atom",
        "displace_along_vibrational_mode": "_displace_along_vibrational_mode",
        "bind_reached_geometry": "_bind_reached_geometry",
        "characterise_stationary_point": "_characterise_stationary_point",
        "inspect_database_records": "_inspect_database_records",
        "extract_database_record_geometry": "_extract_database_record_geometry",
        "read_project_yaml": "_read_project_yaml",
        "validate_project_yaml": "_validate_project_yaml",
        "plan_scientific_workflow": "_plan_scientific_workflow",
        "amend_scientific_workflow": "_amend_scientific_workflow",
        "inspect_workflow_frontier": "_inspect_workflow_frontier",
        "select_execution_wave": "_select_execution_wave",
        "prepare_program_node": "_prepare_program_node",
        "synthesize_command": "_synthesize_command",
        "preview_command": "_preview_command",
        "preflight_program_node": "_preflight_program_node",
        "inspect_calculation_artifact": "_inspect_calculation_artifact",
        "inspect_result_selectors": "_inspect_result_selectors",
        "inspect_run_outcome": "_inspect_run_outcome",
        "extract_result_quantities": "_extract_result_quantities",
        "derive_thermochemistry": "_derive_thermochemistry",
        "evaluate_quantity_expression": "_evaluate_quantity_expression",
        "evaluate_scientific_validation": "_evaluate_scientific_validation",
        "record_analysis_claims": "_record_analysis_claims",
        "record_scientific_decision": "_record_scientific_decision",
        "declare_requested_observable": "_declare_requested_observable",
        "execute_approved_program_node": "_execute_approved_program_node",
        "consult_domain_skill": "_consult_domain_skill",
        "open_guide": "_open_guide",
    }

    def dispatch(
        self, *, turn_id: str, tool_name: str, arguments: Mapping[str, Any]
    ) -> dict[str, Any]:
        """Validate a call and invoke exactly one approved host operation."""

        values = dict(arguments)
        # A leaf tool called by name before its guide opened: the model
        # asked, so the guide opens on that signal and the call proceeds.
        owner = guide_for_tool(tool_name)
        opened_by_call: tuple[Any, ...] = ()
        if owner and owner not in self.active_guides:
            opened_by_call = self.activate_guides(
                turn_id, (owner,), signal="model_call"
            )
        _validate_tool_arguments(self.surface, tool_name, values)
        handlers = {
            name: getattr(self, method)
            for name, method in self.TOOL_HANDLERS.items()
        }
        handler = handlers.get(tool_name)
        if handler is None:
            raise ContractError("tool is absent from command-compiled profile")
        self._reply_observations = ()
        result = handler(turn_id, values)
        reply = {
            "schema_version": "chemsmart.tool-result.v1",
            "tool": tool_name,
            "status": "ok",
            "result": _model_visible_data(canonical_data(result)),
        }
        if self._reply_observations:
            # A handler's observations about its own inputs ride beside
            # the receipt rather than inside it, so the receipt's digest
            # stays what the arithmetic makes it.
            reply["observations"] = tuple(self._reply_observations)
            self._reply_observations = ()
        # Every guide this call opened travels back with its body, whether
        # the call itself asked (a leaf tool by name) or the plan did.
        opened = tuple(opened_by_call) + tuple(
            self._guides_from_planning(turn_id, tool_name, result)
        )
        if opened:
            reply["guides_opened"] = tuple(
                self._guide_record(guide) for guide in opened
            )
        return reply

    # -- guides: the leaves of the surface ---------------------------------

    def activate_guides(
        self, turn_id: str, guide_ids: Iterable[str], *, signal: str
    ) -> tuple[Any, ...]:
        """Open guides not yet open; rebuild the surface; record each."""

        from chemsmart.agent.guides import GUIDES_BY_ID

        opened = []
        for guide_id in guide_ids:
            guide = GUIDES_BY_ID.get(str(guide_id))
            if guide is None or guide.guide_id in self.active_guides:
                continue
            self.active_guides.add(guide.guide_id)
            opened.append(guide)
        if not opened:
            return ()
        if self.surface.profile == "command_compiled_preview":
            self.surface = build_command_compiled_tool_surface(
                self.registry, guides=tuple(sorted(self.active_guides))
            )
        for guide in opened:
            self.event_store.append(
                turn_id=turn_id,
                kind=EventKind.GUIDE_ACTIVATED.value,
                payload={
                    "guide_id": guide.guide_id,
                    "signal": signal,
                    "tools": list(guide.tools),
                    "operations": list(guide.operations),
                    "tool_schema_sha256": self.surface.tool_schema_sha256,
                },
                idempotency_key=f"guide:{turn_id}:{guide.guide_id}:{signal}",
            )
        return tuple(opened)

    @staticmethod
    def _guide_record(guide: Any) -> dict[str, Any]:
        # The guide's own rules render inside its body, once, when it
        # opens; they used to render in the stem for every session.
        from chemsmart.agent.rules import render_rules

        leaf_rules = render_rules(f"leaf:{guide.guide_id}")
        return {
            "guide_id": guide.guide_id,
            "title": guide.title,
            "body": (
                guide.body + " " + leaf_rules if leaf_rules else guide.body
            ),
            "tools_now_available": list(guide.tools),
            "operations_now_available": list(guide.operations),
        }

    def _guides_from_planning(
        self, turn_id: str, tool_name: str, result: Any
    ) -> tuple[Any, ...]:
        """The plan-derived signal: what the DAG the model just planned
        needs. Read from the typed plan, never from prose."""

        from chemsmart.agent.guides import guides_from_plan

        if tool_name not in {
            "plan_scientific_workflow",
            "amend_scientific_workflow",
            "declare_requested_observable",
        }:
            return ()
        jobtypes: set[str] = set()
        operations: set[str] = set()
        constants: set[str] = set()
        programs: set[str] = set()
        for plan in self.scientific_workflow_plans.values():
            for node in getattr(plan, "nodes", ()):
                jobtypes.add(str(getattr(node, "jobtype", "")))
                programs.add(str(getattr(node, "program", "")))
        for toolchain in self.scientific_toolchain_plans.values():
            for node in getattr(toolchain, "analysis_nodes", ()):
                for item in getattr(node, "expression_nodes", ()):
                    operations.add(str(item.get("operation", "")))
                    if str(item.get("operation", "")) == "constant":
                        constants.add(str(item.get("constant_name", "")))
        wanted = guides_from_plan(
            jobtypes=jobtypes,
            operations=operations,
            constants=constants,
            programs=programs,
        )
        return self.activate_guides(turn_id, wanted, signal="plan")

    def _open_guide(self, turn_id: str, values: dict) -> Any:
        """The model-pull path: a guide, or an advisory skill."""

        from chemsmart.agent.guides import GUIDES_BY_ID
        from chemsmart.agent.skills import available_skill_ids

        guide_id = require_identifier(values["guide_id"], "guide_id")
        guide = GUIDES_BY_ID.get(guide_id)
        if guide is not None:
            self.activate_guides(turn_id, (guide_id,), signal="model")
            return {
                "schema_version": "chemsmart.guide-opened.v1",
                **self._guide_record(guide),
                "guidance": (
                    "You opened this guide, so work under it from here. Its "
                    "tools and operations are on the surface now. It settles "
                    "no scientific status: readiness, approval, terminal "
                    "state, validity and accuracy come only from typed host "
                    "receipts."
                ),
            }
        if guide_id in available_skill_ids():
            return self._consult_domain_skill(turn_id, {"skill_id": guide_id})
        raise ContractError(
            f"unknown guide {guide_id!r}; guides: {sorted(GUIDES_BY_ID)}; "
            f"skills: {list(available_skill_ids())}"
        )

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

    def _evidence_already_in_hand(self) -> bool:
        """Whether this GOAL already holds numbers for this task.

        A pre-registration predates the physics. OPEN-1 ino3 extracted
        its spin populations, then declared six observables with bands
        drawn around the values it had just read, and the completion
        printed twelve rows saying `agreed` with nothing to separate
        them from the two written before any number existed
        (2026-09-07). The verdict does not change -- being right after
        the fact is still being right -- but the reader is told.

        The promise is task-wide and the check was process-wide, which
        are the same thing only for a goal that never woke. A woken
        cycle is a fresh host with every registry empty, so at cycle 3
        of ino3-r15 a session declared a new tolerance-bearing
        observable as its first act -- quoting the previous cycle's
        delivered numbers verbatim in its own basis, and choosing the
        tolerance after seeing the uncertainty it had to clear -- and
        the host recorded it as a pre-registration. The goal's own
        physics was on this object at that moment. Nothing is refused
        and no verdict moves; the reader is told, which is the whole
        point of the flag (SUFFICIENCY-5, 2026-09-10).

        Deliberately NOT evidence: prior *declarations*. A cycle that
        declared and delivered nothing produced no physics, and a
        genuinely new prediction after a first molecule is ordinary
        science that must stay possible.
        """

        return bool(
            self.quantity_extractions
            or self.thermochemistry_receipts
            or self.quantity_expression_receipts
            or self.analysis_claim_records
            # Physics this cycle holds that the four above miss: a
            # session carrying only a validation verdict, a
            # characterisation or an anomaly has numbers in hand too.
            or self.scientific_validation_receipts
            or self.stationary_point_characterisations
            or self.anomaly_observations
            or self.reached_geometries
            # And physics the GOAL holds from an earlier cycle: a
            # delivered observable is a claim that was made.
            or self.goal_delivered_declared_ids
            or self._recorded_run_analysis()
        )

    def _recorded_run_analysis(self) -> bool:
        """Whether a recorded run of this goal already produced numbers.

        Read from the durable streams the host wrote, so a woken cycle
        inherits the chronology instead of starting it again. Cached:
        the declaration path may ask more than once per session.
        """

        cached = getattr(self, "_recorded_analysis_seen", None)
        if cached is not None:
            return bool(cached)
        root = getattr(self, "run_evidence_root", None)
        seen = False
        if root:
            marks = (
                '"kind": "result_quantities_extracted"',
                '"kind": "thermochemistry_derived"',
                '"kind": "quantity_expression_evaluated"',
                '"kind": "analysis_claims_recorded"',
            )
            for stream in sorted(
                Path(root).glob(".chemsmart-agent/goals/*/runs/*/events.jsonl")
            ):
                try:
                    text = stream.read_text(encoding="utf-8")
                except OSError:
                    continue
                if any(mark in text for mark in marks):
                    seen = True
                    break
        self._recorded_analysis_seen = seen
        return seen

    def _declare_requested_observable(self, turn_id: str, values: dict) -> Any:
        """Bind the session's restatement of what the task asks for.

        A declaration is a commitment, not a plan: the completion gate
        later requires a delivered claim of each declared dimension --
        kind and unit, never value -- and an undelivered declared
        observable is named in the completion receipt exactly like a
        plan output the chain could not fulfil.
        """

        declared = []
        kept_prior: list[str] = []
        declared_after_evidence = self._evidence_already_in_hand()
        for item in values["observables"]:
            observable_id = str(item["observable_id"])
            require_identifier(observable_id, "observable_id")
            unit = str(item["unit"]).strip()
            meaning = str(item["meaning"]).strip()
            if not meaning:
                raise ContractError(
                    "a declared observable requires one sentence of meaning"
                )
            try:
                dimension = unit_dimension(unit)
            except QuantityExpressionError as exc:
                raise RoutedContractError(
                    gate="declaration.unit_is_in_the_typed_vocabulary",
                    invariant=(
                        "a declared observable is joined to its claim by "
                        "id and judged in its dimension, so its unit is "
                        "one the typed plane measures."
                    ),
                    diagnosis=(
                        f"declared unit {unit!r} is not in the typed unit "
                        f"vocabulary: {exc}."
                    ),
                    route=(
                        "declare the unit the answer will be reported in, "
                        "from the typed vocabulary, e.g. 'kcal/mol', "
                        "'kJ/mol', 'eV', 'cm^-1', 'angstrom', 'degree', "
                        "'1' for a count."
                    ),
                ) from None
            # Both the committed declarations and the ones this call has
            # already accepted: the commit is deferred to the end, so a
            # duplicate inside one call is invisible to the dict.
            existing = self.requested_observable_declarations.get(
                observable_id
            ) or next(
                (
                    pending
                    for pending in declared
                    if pending["observable_id"] == observable_id
                ),
                None,
            )
            if existing is not None:
                if tuple(existing["dimension"]) != tuple(
                    int(value) for value in dimension
                ):
                    raise ContractError(
                        f"declared observable {observable_id!r} is already "
                        f"bound to unit {existing['unit']!r}; a changed "
                        "observable is a new declaration under a new "
                        "identifier"
                    )
                # The first declaration stands, and the reply says so: a
                # woken session re-declared its expectations with a
                # flipped sign convention and wider bands, and the
                # completion row printed agreed over a falsified first
                # prior (live, 2026-09-02).
                kept_prior.append(observable_id)
                continue
            expected_sign = str(item.get("expected_sign", "")).strip().lower()
            basis = str(item.get("expectation_basis", "")).strip()
            low = item.get("expected_low")
            high = item.get("expected_high")
            if (low is None) != (high is None):
                raise ContractError(
                    "an expected range needs both ends: expected_low and "
                    "expected_high, in the observable's own unit"
                )
            # A point expectation -- a count of exactly one imaginary mode,
            # a yes-or-no observable -- is a band whose ends coincide.
            # Observed live (W1c, R2c): both sessions declared one and were
            # refused for it, one refusal each, on an otherwise clean turn.
            if low is not None and float(low) > float(high):
                raise ContractError(
                    "expected_low must not exceed expected_high; a point "
                    "expectation is written with both ends equal"
                )
            if expected_sign and expected_sign not in {"positive", "negative"}:
                raise ContractError(
                    "expected_sign is 'positive' or 'negative'; the "
                    "vocabulary grows only when a loss class earns a new "
                    "comparator"
                )
            # A sign the band excludes can never agree: five correct zero
            # imaginary-mode counts printed "diverged" because their
            # expectation carried expected_sign positive with a 0..0 band
            # (live, 2026-09-02). A zero has no sign.
            if low is not None and (
                (expected_sign == "positive" and float(high) <= 0.0)
                or (expected_sign == "negative" and float(low) >= 0.0)
            ):
                raise ContractError(
                    f"expected_sign {expected_sign!r} contradicts the band "
                    f"[{float(low):g}, {float(high):g}]: a zero has no sign, "
                    "so omit expected_sign for a zero-valued expectation and "
                    "let the band speak"
                )
            if (expected_sign or low is not None) and not basis:
                raise ContractError(
                    "an expectation requires expectation_basis: what it "
                    "rests on. An expectation without a reason is a coin "
                    "flip, and the record would carry it as though it were "
                    "reasoning"
                )
            # A diagnostic is the session's own prediction about the
            # route, given standing (owner ruling R3, 2026-09-06): it is
            # joined and scored like a requested expectation and it is
            # never a deliverable, so it needs a prediction to score and
            # a rule saying what its failure changes.
            role = str(item.get("role") or "requested").strip().lower()
            if role not in {"requested", "diagnostic"}:
                raise ContractError(
                    "role is 'requested' or 'diagnostic'; nothing else "
                    "has standing"
                )
            update_rule = str(item.get("failure_update_rule", "")).strip()
            resolution = item.get("method_resolution")
            if role == "diagnostic" and not (expected_sign or low is not None):
                raise ContractError(
                    f"diagnostic {observable_id!r} predicts nothing: a "
                    "diagnostic carries expected_sign or a band, or it "
                    "cannot be scored"
                )
            if role == "diagnostic" and not update_rule:
                raise ContractError(
                    f"diagnostic {observable_id!r} needs "
                    "failure_update_rule: what its falsification changes "
                    "about the route. A prediction whose failure changes "
                    "nothing is not a diagnostic"
                )
            if resolution is not None and float(resolution) < 0.0:
                raise ContractError(
                    "method_resolution is a magnitude in the observable's "
                    "own unit and cannot be negative"
                )
            # The precision the task asks for. It is the one number in
            # this record the model does not choose: it restates it, and
            # the host holds it against the uncertainty a claim states.
            # Nothing here held it before, so a session could deliver a
            # number its own limitations said missed the requester's
            # tolerance and the goal settled achieved over it (OPEN-2
            # ino3-qwen, 2026-09-07).
            tolerance = item.get("required_tolerance")
            tolerance_basis = str(item.get("tolerance_basis", "")).strip()
            tolerance_origin = str(item.get("tolerance_origin", "")).strip()
            if tolerance is not None:
                if float(tolerance) < 0.0:
                    raise ContractError(
                        "required_tolerance is a magnitude in the "
                        "observable's own unit and cannot be negative"
                    )
                if not tolerance_basis:
                    raise ContractError(
                        f"observable {observable_id!r} states a "
                        "required_tolerance without tolerance_basis: what "
                        "in the task fixes it, or that the task fixes "
                        "none and this is your reading. A tolerance with "
                        "no source is a number nobody asked for"
                    )
                # Who set the precision is a fact about the
                # obligation, not a judgement about it. A live session
                # declared a margin with a self-chosen 0.4 V tolerance
                # after seeing that its uncertainty was 0.25 V, said so
                # honestly in prose, and the record could not join on
                # it: `tolerance_basis` is free text and every reader
                # that decides anything reads the number beside it. The
                # route itself is one the host offers and stays legal;
                # what changes is that a reader can tell the two apart
                # (SUFFICIENCY-5, 2026-09-10).
                if tolerance_origin not in {"task", "session"}:
                    # Recorded as unstated, never refused. Requiring the
                    # field would be a gate, and the ruling this repair
                    # implements is to record the origin and gate
                    # nothing: declaring the margin a decision turns on
                    # is a route the host itself offers, and a
                    # legitimate exploratory tolerance must not cost a
                    # refusal. A reader sees the omission instead.
                    tolerance_origin = "unstated"
                if role == "diagnostic":
                    raise ContractError(
                        f"diagnostic {observable_id!r} carries a "
                        "required_tolerance: a diagnostic is your own "
                        "prediction about the route and is never owed, so "
                        "no tolerance is owed on it either"
                    )
            record = {
                "observable_id": observable_id,
                "unit": unit,
                "dimension": tuple(int(value) for value in dimension),
                "meaning": meaning,
            }
            supersedes = str(item.get("supersedes_observable_id", "")).strip()
            if supersedes:
                if supersedes not in self.requested_observable_declarations:
                    raise ContractError(
                        f"supersedes_observable_id {supersedes!r} names no "
                        "observable this session declared"
                    )
                if supersedes == observable_id:
                    raise ContractError(
                        "an observable cannot supersede itself"
                    )
                # Supersession corrects the observable -- its
                # identifier, its unit -- and never retires the
                # obligation the task set. A replacement needed only to
                # name an existing id, so a tolerance-free declaration,
                # or a diagnostic, could make a requested precision
                # requirement vanish by relabelling. The host does not
                # carry the number across, because supersession exists
                # for a wrong unit and a tolerance is stated in the
                # observable's own unit: the session restates it in the
                # corrected unit, which is a judgement only it can make.
                retired_record = self.requested_observable_declarations[
                    supersedes
                ]
                if retired_record.get("required_tolerance") is not None:
                    if role == "diagnostic":
                        raise ContractError(
                            f"observable {supersedes!r} carries a required "
                            "tolerance the task asked for; a diagnostic is "
                            "never owed and cannot retire it. Supersede it "
                            "with a requested observable, or deliver it"
                        )
                    if tolerance is None:
                        raise ContractError(
                            f"observable {supersedes!r} carries a required "
                            f"tolerance of "
                            f"{retired_record['required_tolerance']} "
                            f"{retired_record.get('unit')!r}; the "
                            "declaration replacing it states none. Restate "
                            "the tolerance in this observable's own unit -- "
                            "the host will not convert it, because a "
                            "supersession is often a unit correction"
                        )
                    # Requiring a tolerance preserved its presence and
                    # not the obligation: 2 kJ/mol was replaced by
                    # 200 kJ/mol under the same meaning, the same unit
                    # and the same quoted source, and the host retired
                    # the original. The model authored its finish line
                    # one declaration later. Where the host owns the
                    # conversion the physical tolerance is preserved by
                    # arithmetic; where it does not -- the cross-
                    # dimension correction this route exists for -- the
                    # restatement stands and the record carries what it
                    # replaced, so the change is visible rather than
                    # silent.
                    retired_tolerance = float(
                        retired_record["required_tolerance"]
                    )
                    retired_unit = str(retired_record.get("unit") or "")
                    converted = (
                        _restate_display_value(
                            retired_tolerance, retired_unit, unit
                        )
                        if retired_unit and unit
                        else None
                    )
                    if converted is not None and not _is_restatement_of(
                        float(tolerance), converted
                    ):
                        raise ContractError(
                            f"observable {supersedes!r} carries a required "
                            f"tolerance of {retired_tolerance} "
                            f"{retired_unit!r}, which is {converted:g} "
                            f"{unit!r}; this declaration states "
                            f"{float(tolerance):g} {unit!r}. A supersession "
                            "corrects how an observable is named or "
                            "measured and never how good the answer has to "
                            "be -- the task set that. State the same "
                            "precision, or deliver the observable you "
                            "declared"
                        )
                    record["superseded_required_tolerance"] = retired_tolerance
                    record["superseded_tolerance_unit"] = retired_unit
                record["supersedes_observable_id"] = supersedes
            if role == "diagnostic":
                record["role"] = role
            if declared_after_evidence:
                record["declared_after_evidence"] = True
            if update_rule:
                record["failure_update_rule"] = update_rule
            if resolution is not None:
                record["method_resolution"] = float(resolution)
            if tolerance is not None:
                record["required_tolerance"] = float(tolerance)
                record["tolerance_basis"] = tolerance_basis
                record["tolerance_origin"] = tolerance_origin
            if expected_sign or low is not None:
                record["expectation_basis"] = basis
            if expected_sign:
                record["expected_sign"] = expected_sign
                # A sign the band already implies is not a prediction: a
                # gap declared as next-minus-ground with band 0..60 and
                # sign positive cannot fail its sign (NOVEL-2 ino1,
                # 2026-09-04). Said on the record, so the row never
                # prints agreed on it.
                if low is not None and (
                    (expected_sign == "positive" and float(low) >= 0.0)
                    or (expected_sign == "negative" and float(high) <= 0.0)
                ):
                    record["sign_implied_by_band"] = True
            if low is not None:
                record["expected_low"] = float(low)
                record["expected_high"] = float(high)
            declared.append(record)
        # Commit nothing until every item has passed. The loop above used
        # to write each record as it went, so a call that raised on a
        # later item left the earlier ones in the host while the model was
        # told the call was rejected -- and the event below, which is the
        # only provenance this field has, was never appended. One live run
        # declared seven observables, raised on the third, and the two
        # carrying its entire headline reached the review packet, the
        # approval bundle and the completion receipt while appearing in no
        # declaration event anywhere.
        for record in declared:
            self.requested_observable_declarations[
                str(record["observable_id"])
            ] = record
        if declared:
            self.event_store.append(
                turn_id=turn_id,
                kind=EventKind.REQUESTED_OBSERVABLE_DECLARED.value,
                payload={
                    "observables": tuple(declared),
                    "declared_total": len(
                        self.requested_observable_declarations
                    ),
                },
                idempotency_key=(
                    "requested-observable:"
                    + canonical_sha256(
                        tuple(item["observable_id"] for item in declared)
                    )
                ),
            )
        reach_warnings = {
            str(record["observable_id"]): warning
            for record in declared
            for warning in (
                self._declaration_reach_warning(str(record["meaning"])),
            )
            if warning
        }
        return {
            "declared": tuple(declared),
            "declared_total": len(self.requested_observable_declarations),
            "kept_prior": tuple(kept_prior),
            **({"reach_warnings": reach_warnings} if reach_warnings else {}),
            "kept_prior_meaning": (
                "the goal's first declaration of an observable stands; an "
                "expectation written after the physics exists is not a "
                "prediction"
                if kept_prior
                else ""
            ),
            "completion_consequence": (
                "the completion gate requires, for every declared "
                "observable, a delivered claim carrying its id -- as the "
                "claim's claim_id, or as the quantity id of the receipt "
                "it stands on -- in its dimension; values are never "
                "checked, and a declared observable you cannot deliver "
                "is stated as a limitation"
            ),
        }

    def termination_notice(self) -> dict[str, Any] | None:
        """What the host says, once, when a goal session is about to end
        with declared observables undelivered while budget remains.

        Informational and never a demand (owner ruling, 2026-09-06). None
        outside a goal, when nothing declared is undelivered, or when no
        budget line remains -- ending is then the only thing left.
        """

        if self.engine_calls_remaining is None:
            return None
        undelivered: list[str] = []
        for task_spec_sha256 in sorted(self.task_spec_sha256s):
            _misses, limitations = self._declared_observable_completion(
                task_spec_sha256=task_spec_sha256
            )
            for item in limitations:
                observable_id = str(item).split(":", 1)[-1]
                if (
                    observable_id in self.goal_delivered_declared_ids
                    or observable_id in self.verified_unreachable_ids
                    or observable_id in undelivered
                ):
                    continue
                undelivered.append(observable_id)
        if not undelivered:
            return None
        budgets = {
            "engine_calls_remaining": int(self.engine_calls_remaining),
            "excursion_calls_remaining": int(
                self.excursion_calls_remaining or 0
            ),
            "wall_seconds_remaining": float(
                self.wall_seconds_remaining or 0.0
            ),
            "revisions_remaining": int(self.revisions_remaining or 0),
        }
        if not any(budgets.values()):
            return None
        from chemsmart.agent.rules import rules_by_id

        text = (
            rules_by_id()["wake.termination_notice"].text
            + " Undelivered declared observables: "
            + ", ".join(undelivered)
            + ". Remaining: engine calls "
            + str(budgets["engine_calls_remaining"])
            + ", excursion calls "
            + str(budgets["excursion_calls_remaining"])
            + ", revisions "
            + str(budgets["revisions_remaining"])
            + ", engine wall "
            + f"{budgets['wall_seconds_remaining']:.0f} s."
        )
        return {
            "undelivered_declared_observable_ids": tuple(undelivered),
            "budgets": budgets,
            "text": text,
        }

    def _declaration_reach_warning(self, meaning: str) -> str:
        """Say, at declaration, when a meaning names a quantity kind no
        envelope program declares -- and the two routes.

        The readers' own selector vocabulary defines the kinds: a
        selector whose distinctive words appear in the meaning is what
        the meaning names. A session declared the cation's spin
        distribution as an observable, ran four cycles, and learned only
        at the end that no route to it existed in the envelope (NOVEL-3
        ino3, 2026-09-05). Warn and route, never refuse: the session may
        compose the quantity or refuse it in a form the host verifies.
        """

        import re

        from chemsmart.analysis.result_readers import (
            registered_reader_selectors,
        )

        envelope = self.bounded_execution_envelope
        if envelope is None:
            return ""
        programs = tuple(
            str(program)
            for program, _engines in envelope.allowed_program_engines
        )
        generic = {
            "atomic",
            "total",
            "energy",
            "energies",
            "final",
            "value",
            "per",
            "atom",
            "of",
            "the",
            "and",
            "in",
            "on",
            "last",
            "first",
            "count",
            "number",
            "list",
            "table",
            "block",
            "printed",
        }
        words = {
            token.rstrip("s")
            for token in re.findall(r"[a-z0-9]+", meaning.lower())
        }
        matched: dict[str, set[str]] = {}
        for program, selectors in registered_reader_selectors().items():
            for selector in selectors:
                parts = {
                    part.rstrip("s") for part in selector.split("_")
                } - generic
                # Two distinctive words, always: a selector with one
                # ("dipole") cannot be told from a passing mention, and a
                # meaning that says "Gibbs" names no selector by it.
                if len(parts) < 2 or len(parts & words) < 2:
                    continue
                matched.setdefault(selector, set()).add(str(program))
        if not matched:
            return ""
        served = {
            selector
            for selector, owners in matched.items()
            if owners & set(programs)
        }
        if served:
            return ""
        unserved = sorted(matched)
        elsewhere = sorted({p for owners in matched.values() for p in owners})
        return (
            "the meaning names a quantity kind "
            f"({', '.join(unserved)}) that no envelope program "
            f"({', '.join(programs)}) declares"
            + (
                f"; declared by {', '.join(elsewhere)}, outside the "
                "envelope"
                if elsewhere
                else ""
            )
            + ". Routes: compose it from selectors the envelope's programs "
            "declare, or refuse it in record_scientific_decision's "
            "unreachable_observable_ids with the receipts that show the "
            "gap -- the host verifies a selector absence and settles "
            "unreachable_from_evidence. Nothing is refused here."
        )

    def _declared_observable_completion(
        self, *, task_spec_sha256: str
    ) -> tuple[tuple[str, ...], tuple[str, ...]]:
        """Match declared observables against delivered claim dimensions.

        Kind and unit, never value: a declaration is satisfied by any
        recorded claim on this task whose dimension equals the declared
        unit's.  An unmatched declaration is a limitation on a green
        receipt -- the delivery stands and states what it delivered
        without -- never a finding, because findings mean the chain
        itself broke.  Dimension vectors from different eras differ
        only by trailing bases, so both sides are compared zero-padded.
        A diagnostic -- the session's own prediction about the route --
        is joined when delivered and is never a miss: it is no
        deliverable, so its absence is no limitation and never reaches
        the notice or the settlement (owner ruling R3, 2026-09-06).
        Returns (human-readable misses for the event, limitation ids).
        """

        if not self.requested_observable_declarations:
            return (), ()

        def _padded(dimension: tuple[int, ...]) -> tuple[int, ...]:
            values = tuple(int(value) for value in dimension)
            return values + (0,) * (9 - len(values))

        # A dimension is not an identity. Matching by dimension certified
        # a delivery whose declared endo:exo ratio had no claim, because
        # two imaginary-mode counts share its dimension, while the session
        # had deliberately relabelled its number as a conformer difference
        # (live, 2026-09-03). A claim answers a declaration by carrying
        # its id; the dimension is then checked, never used to guess.
        # The id a claim carries lives in two fields: claim_id, which the
        # session or the executor names, and quantity_id, which the
        # receipt the claim stands on was computed under. A live chain
        # rendered six correct numbers with the declared id in
        # quantity_id and the plan's short input label in claim_id, and
        # this gate, reading claim_id alone, called all six undelivered
        # three fields away from the id it wanted (NOVEL-2 po2,
        # 2026-09-04). The host's word must be true of the record it
        # holds: either field answers, claim_id first, and the join
        # says which field it read.
        claims_by_id: dict[str, tuple[int, ...]] = {}
        joined_on: dict[str, str] = {}
        for record in self.analysis_claim_records.values():
            if getattr(record, "task_spec_sha256", "") != task_spec_sha256:
                continue
            for claim in getattr(record, "claims", ()):
                for field in ("claim_id", "quantity_id"):
                    key = str(getattr(claim, field, "") or "")
                    if key:
                        # The newest claim answers. This gate was
                        # first-wins while the workspace record and the
                        # settlement are last-wins, so a session that
                        # corrected a wrong-dimension claim -- which is
                        # exactly what this gate's own refusal tells it
                        # to do -- was told "undelivered" for ever while
                        # the settlement called it delivered. Two
                        # readers, one record, opposite answers, in the
                        # direction that hides the miss.
                        claims_by_id[key] = _padded(claim.dimension)
                        joined_on[key] = field
        from chemsmart.agent.delivery import superseded_observable_ids

        retired = superseded_observable_ids(
            tuple(self.requested_observable_declarations.values())
        )
        misses = []
        limitations = []
        for observable_id, record in sorted(
            self.requested_observable_declarations.items()
        ):
            delivered = claims_by_id.get(observable_id)
            if delivered is not None and delivered == _padded(
                record["dimension"]
            ):
                self._declared_observable_join_fields[observable_id] = (
                    joined_on[observable_id]
                )
                continue
            if record.get("role") == "diagnostic":
                continue
            if observable_id in retired:
                # A later declaration named this one as its correction:
                # the meaning moved to the new id and this one is not
                # owed. Both stay on the record (OPEN-1 ino3 declared
                # six spin observables in 'e' and corrected them).
                continue
            misses.append(
                f"declared observable {observable_id!r} "
                f"({record['unit']}) has no delivered claim named "
                f"{observable_id!r}"
                + (
                    " of matching dimension: the claim under that id is "
                    f"in {canonical_unit_for_dimension(delivered)!r} and "
                    f"the declaration in {record['unit']!r}; restate it "
                    "with wavenumber_to_energy / energy_to_wavenumber "
                    "(frequency and energy) or convert (one dimension) "
                    "and claim it again under the id"
                    if delivered is not None
                    else "; a claim answers a declaration by carrying its "
                    "id as its claim_id, or by standing on a receipt "
                    "quantity of that id"
                )
            )
            limitations.append(f"declared_observable:{observable_id}")
        return tuple(misses), tuple(limitations)

    #: Which field answered each declared observable at the last
    #: completion: ``claim_id`` or ``quantity_id``. Read by the completion
    #: event so a reader can see a join that would have failed under the
    #: old single-field rule.
    _declared_observable_join_fields: dict[str, str]

    def _declared_observable_predictions(
        self, *, task_spec_sha256: str
    ) -> tuple[dict[str, Any], ...]:
        """Restate each recorded expectation beside what was delivered.

        A session that writes down what it expects before the evidence
        exists has done the honest thing, and the host owes that
        expectation the same visibility as the number that answers it:
        a live session predicted one isomer lower on a steric argument,
        the physics returned the other, and the record carried the
        falsified premise beside the correct value with nothing joining
        them.  What is compared is only the sign of a scalar claim of
        the declared dimension, by arithmetic the host already owns.

        Divergence is never a finding and never a limitation.  A wrong
        prediction is a scientific result -- often the interesting one
        -- while a finding means the chain itself broke; conflating
        them would make a correct delivery look defective and teach
        sessions to predict nothing.  When the dimension does not
        identify one scalar claim the row says so rather than guessing
        which number the expectation was about.
        """

        predicted = {
            observable_id: record
            for observable_id, record in (
                self.requested_observable_declarations.items()
            )
            if record.get("expected_sign") or "expected_low" in record
        }
        if not predicted:
            return ()

        def _padded(dimension: tuple[int, ...]) -> tuple[int, ...]:
            values = tuple(int(value) for value in dimension)
            return values + (0,) * (9 - len(values))

        claims_by_id: dict[str, Any] = {}
        for claim_record in self.analysis_claim_records.values():
            if (
                getattr(claim_record, "task_spec_sha256", "")
                != task_spec_sha256
            ):
                continue
            for claim in getattr(claim_record, "claims", ()):
                # Same two-field join as the completion gate above.
                for field in ("claim_id", "quantity_id"):
                    key = str(getattr(claim, field, "") or "")
                    if key:
                        claims_by_id.setdefault(key, claim)

        rows = []
        for observable_id, record in sorted(predicted.items()):
            row = {
                "observable_id": observable_id,
                "expected_sign": record.get("expected_sign", ""),
                "expected_low": record.get("expected_low", ""),
                "expected_high": record.get("expected_high", ""),
                "expectation_basis": record["expectation_basis"],
                "delivered_claim_id": "",
                "delivered_value": "",
                "delivered_unit": "",
                "agreement": "not_comparable",
                "role": str(record.get("role") or "requested"),
                **(
                    {"declared_after_evidence": True}
                    if record.get("declared_after_evidence")
                    else {}
                ),
            }
            if "failure_update_rule" in record:
                row["failure_update_rule"] = record["failure_update_rule"]
            if "method_resolution" in record:
                row["method_resolution"] = record["method_resolution"]

            def _scalar(claim: Any) -> bool:
                return isinstance(
                    claim.display_value, (int, float)
                ) and not isinstance(claim.display_value, bool)

            # The identifier first, the dimension only as a fallback.  A
            # session names its claims after the observables it declared,
            # and three potentials in volts are three different questions
            # that a dimension cannot tell apart: the first live use
            # declared a sign for each of three volt-valued observables
            # and every row came back "not comparable" while the claims
            # sat there carrying the very same identifiers.
            # By identifier only. The dimension fallback joined a
            # declared endo:exo difference to a number the session had
            # relabelled a conformer difference and printed "agreed"
            # (live, 2026-09-03); a dimension is not an identity.
            named = claims_by_id.get(observable_id)
            matches = [named] if named is not None and _scalar(named) else []
            if len(matches) == 1:
                claim = matches[0]
                value = float(claim.display_value)
                row["delivered_claim_id"] = claim.claim_id
                row["delivered_value"] = claim.display_value
                row["delivered_unit"] = claim.display_unit
                # Every stated expectation is tested and all must hold.
                # A range is tested only in the unit it was declared in;
                # inventing a conversion here would be arithmetic the
                # expectation never authorised.
                verdicts = []
                sign = record.get("expected_sign", "")
                if record.get("sign_implied_by_band"):
                    # The band decides; the sign adds no prediction.
                    row["sign_implied_by_band"] = True
                elif sign:
                    verdicts.append(
                        value != 0.0
                        and (
                            ("positive" if value > 0.0 else "negative") == sign
                        )
                    )
                if "expected_low" in record:
                    if claim.display_unit == record["unit"]:
                        verdicts.append(
                            record["expected_low"]
                            <= value
                            <= record["expected_high"]
                        )
                    else:
                        # A number in another unit is compared after the
                        # host's own conversion -- within one dimension,
                        # or across frequency and energy through h*c*N_A,
                        # which is a definition and not a choice. A
                        # sixfold-consistent ferromagnetic sign against a
                        # declared antiferromagnetic band printed
                        # not_comparable because the claim was in kcal/mol
                        # and the declaration in cm^-1 (NOVEL-3 ino2).
                        restated = _restate_display_value(
                            value, claim.display_unit, record["unit"]
                        )
                        if restated is None:
                            row["band_untestable"] = (
                                f"delivered in {claim.display_unit!r}, "
                                f"declared in {record['unit']!r}, no "
                                "conversion the host owns"
                            )
                            verdicts.append(None)
                        else:
                            row["delivered_value_in_declared_unit"] = restated
                            verdicts.append(
                                record["expected_low"]
                                <= restated
                                <= record["expected_high"]
                            )
                # A falsified sign diverges whether or not the band could
                # be tested; agreed needs every armed verdict.
                # Within the method's own resolution the number carries
                # no sign and no band: the row says indeterminate rather
                # than grading noise either way.
                comparable = (
                    value
                    if claim.display_unit == record["unit"]
                    else row.get("delivered_value_in_declared_unit")
                )
                if comparable is None and claim.display_unit != record["unit"]:
                    comparable = _restate_display_value(
                        value, claim.display_unit, record["unit"]
                    )
                resolution = record.get("method_resolution")
                if (
                    resolution is not None
                    and comparable is not None
                    and abs(float(comparable)) <= float(resolution)
                ):
                    row["agreement"] = "indeterminate"
                    row["within_method_resolution"] = True
                elif any(verdict is False for verdict in verdicts):
                    row["agreement"] = "diverged"
                elif verdicts and None not in verdicts:
                    # An approximation the session declared answers the
                    # question conditionally and must not be reported as
                    # having met it outright. A live claim delivered a
                    # zero-point-corrected electronic difference,
                    # evaluated on a structure the host had typed a
                    # first-order saddle, against a declaration about a
                    # Gibbs difference between minima -- same id, same
                    # dimension, different quantity -- and this row said
                    # `agreed`. The number stays delivered and the word
                    # carries the qualification (SUFFICIENCY-5).
                    stands_in = getattr(claim, "approximates", None) or {}
                    if str(stands_in.get("observable_id") or "") == str(
                        record.get("observable_id") or ""
                    ) and stands_in.get("relationship"):
                        row["agreement"] = "agreed_as_approximation"
                        row["approximation_relationship"] = str(
                            stands_in.get("relationship")
                        )
                        if stands_in.get("basis"):
                            row["approximation_basis"] = str(
                                stands_in.get("basis")
                            )
                    else:
                        row["agreement"] = "agreed"
            elif len(matches) > 1:
                row["delivered_claim_id"] = ",".join(
                    sorted(claim.claim_id for claim in matches)
                )
            rows.append(row)
        return tuple(rows)

    def _consult_domain_skill(self, turn_id: str, values: dict) -> Any:
        """Return one skill body, as instructions the session now works under.

        The body is knowledge, not evidence: it carries no readiness, no
        approval, no terminal state, and no accuracy claim.  The digests make
        the exact text the model read reconstructible from replay.

        Those status disclaimers used to be the only thing the payload said
        about itself, and they were quietly doing a second job.  "This
        establishes no scientific status" and "you need not follow this" are
        different statements, and a payload stamped only ``advisory_only`` says
        both.  A consulted skill is guidance the session has chosen to adopt --
        the reason to fetch it is to be governed by it -- so the payload now
        says so plainly while every status flag stays exactly as it was.
        """

        skill_id = require_identifier(values["skill_id"], "skill_id")
        document = resolve_skill(skill_id)
        if document is None:
            raise ContractError(f"unknown domain skill: {skill_id}")
        self.consulted_skills[document.document_sha256] = document
        record = {
            "skill_id": document.skill_id,
            "skill_version": document.skill_version,
            "origin": document.origin,
            "body_sha256": document.body_sha256,
            "document_sha256": document.document_sha256,
        }
        self.consulted_skill_records[document.document_sha256] = record
        self.event_store.append(
            turn_id=turn_id,
            kind=EventKind.DOMAIN_SKILL_CONSULTED.value,
            payload={"record": record},
            idempotency_key=(
                "domain-skill:" + turn_id + ":" + document.document_sha256
            ),
        )
        return {
            "schema_version": "chemsmart.domain-skill-consultation.v1",
            "skill_id": document.skill_id,
            "skill_version": document.skill_version,
            "origin": document.origin,
            "description": document.description,
            "body": document.body,
            "body_sha256": document.body_sha256,
            "document_sha256": document.document_sha256,
            # Status axis, unchanged: a skill never establishes scientific
            # standing, and nothing here may be cited as evidence.
            "advisory_only": True,
            "readiness_authority": False,
            "accuracy_authority": False,
            # Adoption axis: what the model should now do about it.
            "applies_to": "the rest of this session",
            "guidance": (
                "You consulted this skill, so work under it from here. Follow "
                "its principles when they bear on what you are doing, and say "
                "when a stated fact came from it. It settles no scientific "
                "status: readiness, approval, terminal state, validity and "
                "accuracy come only from typed host receipts, and where a "
                "receipt and this text disagree the receipt is what happened."
            ),
        }

    def _bind_scan_point_geometry(self, turn_id: str, values: dict) -> Any:
        """Register one chosen point of a completed scan as a geometry input.

        A relaxed scan ends at a surface. ChemSmart writes the structure it
        converged at every point, but those files were opaque
        ``program_output`` with nothing tying them to a coordinate or an
        energy, and nothing could carry one into a later calculation. So a
        session could compute a torsional profile, see exactly where the well
        was, and still have no way to optimise the structure sitting in it.

        Which point matters is a scientific judgement and stays with the
        scientist: this binds the index it is given and ranks nothing. The
        surface is readable through the ordinary selectors --
        ``scan_point_indices``, ``scan_coordinate_values``, ``scan_energies``
        -- so the choice is made from evidence rather than from a rule the
        host imposed.

        What the host owns is the lineage: which result, which point, at what
        coordinate and energy. The geometry that comes back is an ordinary
        trusted input, so using it is a changed molecular input and therefore a
        new workflow with its own review, exactly as the charter requires.
        """

        source = self._artifact(values["artifact_id"])
        reader = reader_for(values.get("program") or "orca")
        if source.kind != reader.artifact_kind:
            raise ContractError(
                f"scan point binding needs a {reader.artifact_kind} result; "
                f"{source.artifact_id!r} is {source.kind}"
            )
        parsed = reader.open_output(Path(source.path))
        records = tuple(getattr(parsed, "scan_point_records", ()) or ())
        if not records:
            raise ContractError(
                f"{source.artifact_id!r} records no scan surface, so it has "
                "no points to choose between"
            )
        requested = int(values["point_index"])
        chosen = next(
            (item for item in records if item["index"] == requested), None
        )
        if chosen is None:
            raise ContractError(
                f"this scan has points 1 to {len(records)}; there is no "
                f"point {requested}"
            )
        if not chosen["geometry_file"]:
            raise ContractError(
                f"point {requested} converged but ChemSmart kept no geometry "
                "file for it, so it cannot be carried forward"
            )
        path = Path(chosen["geometry_file"]).resolve()
        if not path.is_file():
            raise ContractError("the chosen scan point geometry is missing")
        artifact = TrustedArtifactRefV1(
            artifact_id=values["artifact_id"] + f".point.{requested:03d}",
            kind="geometry_xyz",
            sha256=file_sha256(path),
            size_bytes=path.stat().st_size,
            path=str(path),
            cli_value=str(path),
        )
        self.artifacts[artifact.artifact_id] = artifact
        return {
            "schema_version": "chemsmart.scan-point-geometry.v1",
            "artifact": artifact,
            "source_result_artifact_id": source.artifact_id,
            "source_result_sha256": source.sha256,
            "point_index": requested,
            "point_count": len(records),
            "coordinate": chosen["coordinate"],
            "energy_hartree": chosen["energy"],
            "selection_owner": "model",
            "next_action": (
                "bind this geometry's charge and multiplicity, then plan the "
                "stage that uses it as a new workflow for review"
            ),
        }

    def _compose_molecular_arrangement(
        self, turn_id: str, values: dict
    ) -> Any:
        """Place two identity-bound fragments into one host-owned arrangement.

        Without this affordance, sessions that needed a bimolecular
        observable ended uncomputable because nothing could join two
        approved monomers -- observed repeatedly, including a session
        probing invented option names looking for it. The host owns the
        placement mathematics and the bytes; the model owns the scientific
        choices (fragments, contact atoms, distance) and must bind the
        arrangement's charge and multiplicity explicitly afterwards --
        composition never infers an electronic state, and the consuming
        stage is a new workflow.
        """

        if self.approved_workspace is None:
            raise ContractError(
                "composition requires an approved workspace to write into"
            )
        composed_artifact_id = str(values["composed_artifact_id"])
        if composed_artifact_id in self.artifacts:
            taken = sorted(self.artifacts)
            raise ContractError(
                f"artifact ID {composed_artifact_id!r} is already "
                f"registered; choose one not in {taken}"
            )
        fragment_a = self._artifact(values["fragment_a_artifact_id"])
        fragment_b = self._artifact(values["fragment_b_artifact_id"])
        identities = {
            binding.geometry_artifact_sha256: binding
            for binding in self.scientific_identities.values()
        }
        for label, fragment in (("A", fragment_a), ("B", fragment_b)):
            if fragment.sha256 not in identities:
                raise ContractError(
                    f"fragment {label} ({fragment.artifact_id!r}) carries "
                    "no scientific identity; composition requires "
                    "identity-bound parents -- call bind_scientific_identity "
                    "for it first"
                )
        artifact, receipt = compose_trusted_molecular_arrangement(
            approved_workspace=self.approved_workspace,
            composed_artifact_id=composed_artifact_id,
            fragment_a=fragment_a,
            fragment_a_identity_sha256=(
                identities[fragment_a.sha256].binding_sha256
            ),
            fragment_b=fragment_b,
            fragment_b_identity_sha256=(
                identities[fragment_b.sha256].binding_sha256
            ),
            fragment_a_atom=int(values["fragment_a_atom"]),
            fragment_b_atom=int(values["fragment_b_atom"]),
            distance_angstrom=float(values["distance_angstrom"]),
            fragment_b_atoms=(
                [int(item) for item in values["fragment_b_atoms"]]
                if values.get("fragment_b_atoms")
                else None
            ),
            fragment_a_atom_2=(
                int(values["fragment_a_atom_2"])
                if values.get("fragment_a_atom_2") is not None
                else None
            ),
            fragment_b_atom_2=(
                int(values["fragment_b_atom_2"])
                if values.get("fragment_b_atom_2") is not None
                else None
            ),
            distance_angstrom_2=(
                float(values["distance_angstrom_2"])
                if values.get("distance_angstrom_2") is not None
                else None
            ),
        )
        self.artifacts[artifact.artifact_id] = artifact
        self.molecular_compositions[artifact.sha256] = receipt
        self.event_store.append(
            turn_id=turn_id,
            kind=EventKind.MOLECULAR_ARRANGEMENT_COMPOSED.value,
            payload={
                "receipt_sha256": receipt.receipt_sha256,
                "composed_artifact_id": artifact.artifact_id,
                "composed_artifact_sha256": artifact.sha256,
                "fragment_a_sha256": fragment_a.sha256,
                "fragment_b_sha256": fragment_b.sha256,
                "placement": receipt.placement,
            },
            idempotency_key=(
                "molecular-composition:" + receipt.receipt_sha256
            ),
        )
        return {
            "composition": receipt,
            "artifact": artifact,
            "next_action": (
                "bind charge and multiplicity explicitly with "
                "bind_scientific_identity -- composition does not infer "
                "electronic state; the stage that consumes this geometry "
                "is a new workflow needing its own review"
            ),
        }

    def _derive_molecular_species(self, turn_id: str, values: dict) -> Any:
        """Take an ordered subset of one identity-bound parent's atoms.

        Composition could join two fragments and nothing could derive one, so
        every route that makes a new species from an old one -- homolysis,
        deprotonation, extracting a fragment -- ended uncomputable from a
        single parent geometry. Observed three times in one campaign: a
        methanol bond-dissociation session planned all four species correctly,
        previewed the parent green, and then declined because no tool could
        remove the hydrogen and it is not permitted to edit the geometry
        itself.

        The separation is composition's, mirrored. The host owns the selection
        arithmetic, the bytes and the lineage; the model owns which atoms and
        why, and must bind the new charge and multiplicity explicitly --
        removing a hydrogen makes a radical or an anion depending on where its
        electron went, and only the model's chemistry says which. The stage
        that consumes the derived geometry is a new workflow for review.
        """

        if self.approved_workspace is None:
            raise ContractError(
                "derivation requires an approved workspace to write into"
            )
        derived_artifact_id = str(values["derived_artifact_id"])
        if derived_artifact_id in self.artifacts:
            taken = sorted(self.artifacts)
            raise ContractError(
                f"artifact ID {derived_artifact_id!r} is already "
                f"registered; choose one not in {taken}"
            )
        parent = self._artifact(values["parent_artifact_id"])
        identities = {
            binding.geometry_artifact_sha256: binding
            for binding in self.scientific_identities.values()
        }
        if parent.sha256 not in identities:
            raise RoutedContractError(
                gate="derivation.parent_is_identity_bound",
                invariant=(
                    "a derived geometry inherits what its parent is, so "
                    "the parent carries a bound identity first."
                ),
                diagnosis=(
                    f"parent ({parent.artifact_id!r}) carries no "
                    "scientific identity."
                ),
                route=(
                    "bind_scientific_identity on "
                    f"{parent.artifact_id!r} (charge, multiplicity), then "
                    "derive_molecular_species again."
                ),
            )
        kept = values.get("kept_atoms")
        removed = values.get("removed_atoms")
        artifact, receipt = derive_trusted_molecular_species(
            approved_workspace=self.approved_workspace,
            derived_artifact_id=derived_artifact_id,
            parent=parent,
            parent_identity_sha256=(identities[parent.sha256].binding_sha256),
            kept_atoms=None if kept is None else [int(x) for x in kept],
            removed_atoms=(
                None if removed is None else [int(x) for x in removed]
            ),
        )
        self.artifacts[artifact.artifact_id] = artifact
        self.molecular_derivations[artifact.sha256] = receipt
        self.event_store.append(
            turn_id=turn_id,
            kind=EventKind.MOLECULAR_SPECIES_DERIVED.value,
            payload={
                "receipt_sha256": receipt.receipt_sha256,
                "derived_artifact_id": artifact.artifact_id,
                "derived_artifact_sha256": artifact.sha256,
                "parent_sha256": parent.sha256,
                "selection_mode": receipt.selection_mode,
                "kept_atoms": list(receipt.kept_atoms),
                "removed_atoms": list(receipt.removed_atoms),
                "formula": receipt.formula,
                "fragment_count": receipt.fragment_count,
            },
            idempotency_key=("molecular-derivation:" + receipt.receipt_sha256),
        )
        return {
            "derivation": receipt,
            "artifact": artifact,
            "next_action": (
                "bind charge and multiplicity explicitly with "
                "bind_scientific_identity -- derivation does not infer "
                "electronic state, and removing an atom decides neither; the "
                "stage that consumes this geometry is a new workflow needing "
                "its own review"
            ),
        }

    def _fetch_pubchem_geometry(self, turn_id: str, values: dict) -> Any:
        """Bring in a molecule the workspace never supplied.

        Every other origin needs the workspace to already hold the
        molecule: a supplied file, a database record, a previous result,
        or a derivation of one of those. So a session that needed a
        reference couple, a calibration standard or a literature
        comparison had no route, and the campaign read that as the model
        declining on cost (OPEN-2 ino3-qwen, 2026-09-07). The human CLI
        has carried -p/--pubchem for the same programs all along.

        The hub invariant is untouched: the model names an identifier,
        the host fetches through the same library call the CLI uses, and
        the host owns the bytes. What arrives is a database conformer
        with no electronic state, so the next act is an explicit
        bind_scientific_identity and the consuming stage is a new
        workflow for review.
        """

        if self.approved_workspace is None:
            raise ContractError(
                "a pubchem geometry requires an approved workspace to "
                "write into"
            )
        artifact_id = str(values["artifact_id"])
        if artifact_id in self.artifacts:
            taken = sorted(self.artifacts)
            raise ContractError(
                f"artifact ID {artifact_id!r} is already registered; "
                f"choose one not in {taken}"
            )
        artifact, receipt = fetch_trusted_pubchem_geometry(
            approved_workspace=self.approved_workspace,
            artifact_id=artifact_id,
            identifier=values["identifier"],
        )
        self.artifacts[artifact.artifact_id] = artifact
        self.pubchem_geometries[artifact.sha256] = receipt
        self.event_store.append(
            turn_id=turn_id,
            kind=EventKind.PUBCHEM_GEOMETRY_FETCHED.value,
            payload={
                "receipt_sha256": receipt.receipt_sha256,
                "artifact_id": artifact.artifact_id,
                "artifact_sha256": artifact.sha256,
                "identifier": receipt.identifier,
                "identifier_kind": receipt.identifier_kind,
                "atom_count": receipt.atom_count,
                "formula": receipt.formula,
                "fragment_count": receipt.fragment_count,
            },
            idempotency_key=("pubchem-geometry:" + receipt.receipt_sha256),
        )
        molecule = None
        try:
            from chemsmart.io.molecules.structure import Molecule

            molecule = Molecule.from_filepath(artifact.path)
        except Exception:
            molecule = None
        # What arrived, in the channel the session actually reads. The
        # receipt carried `formula`, `atom_count` and `fragment_count`
        # all along and the observations channel one line away returned
        # only a point-group estimate -- so a 102-atom record in 27
        # pieces was described to the session as "C1 within 0.1 A". A
        # live session named two numeric CIDs from prior knowledge and
        # got two unrelated molecules; the fact that refutes the label
        # is the formula, and it belongs where the model looks.
        #
        # Measurements, never refusals: an ion pair, a salt, a solvate
        # and a metallocene as deposited are all legitimately more than
        # one piece, and `compose_molecular_arrangement` exists to
        # consume fragments (SUFFICIENCY-5, 2026-09-10).
        arrival = [
            f"identifier {receipt.identifier!r} "
            f"({receipt.identifier_kind}) returned {receipt.formula} "
            f"with {receipt.atom_count} atoms in "
            f"{receipt.fragment_count} connected piece(s); check the "
            "formula against the molecule you meant"
        ]
        if int(getattr(receipt, "fragment_count", 1) or 1) > 1:
            arrival.append(
                f"this record converted to {receipt.fragment_count} "
                "disconnected pieces -- an observation, not a verdict: "
                "compose_molecular_arrangement consumes fragments, and "
                "an ion pair or a solvate is legitimately more than one"
            )
        return {
            "pubchem_geometry": receipt,
            "artifact": artifact,
            "observations": tuple(arrival)
            + tuple(self._symmetry_observations(artifact, molecule)),
            "next_action": (
                "bind charge and multiplicity explicitly with "
                "bind_scientific_identity -- a database record carries no "
                "electronic state; this is a depositor's conformer and not "
                "a relaxed structure, so measure it before assuming it, and "
                "the stage that consumes it is a new workflow needing its "
                "own review"
            ),
        }

    def _identity_bound_parent(
        self, artifact_id: Any, operation: str
    ) -> tuple[Any, str]:
        """Resolve a parent geometry that already carries an identity."""

        parent = self._artifact(artifact_id)
        identities = {
            binding.geometry_artifact_sha256: binding
            for binding in self.scientific_identities.values()
        }
        if parent.sha256 not in identities:
            raise RoutedContractError(
                gate="derivation.parent_is_identity_bound",
                invariant=(
                    "an edited or composed geometry inherits what its "
                    "parent is, so the parent carries a bound identity "
                    "first."
                ),
                diagnosis=(
                    f"parent ({parent.artifact_id!r}) carries no "
                    "scientific identity."
                ),
                route=(
                    "bind_scientific_identity on "
                    f"{parent.artifact_id!r} (charge, multiplicity), then "
                    f"{operation} again."
                ),
            )
        return parent, identities[parent.sha256].binding_sha256

    def _unused_artifact_id(self, artifact_id: Any) -> str:
        """Reject an identifier already standing for other bytes."""

        value = str(artifact_id)
        if value in self.artifacts:
            raise _artifact_id_taken(value, sorted(self.artifacts))
        return value

    def _edit_molecular_geometry(self, turn_id: str, values: dict) -> Any:
        """Set one internal coordinate of an identity-bound parent.

        The model names the coordinate, the value it wants, and which side of
        the coordinate moves; the host owns the transformation, measures the
        coordinate before and after, and enumerates the atoms it moved. No
        energy exists here, so nothing judges whether the requested value is
        a good one: what makes the edit scientific evidence is that an
        optimisation later disagrees with it, in public.
        """

        if self.approved_workspace is None:
            raise ContractError(
                "a geometry edit requires an approved workspace to write into"
            )
        edited_artifact_id = self._unused_artifact_id(
            values["edited_artifact_id"]
        )
        parent, identity_sha256 = self._identity_bound_parent(
            values["input_artifact_id"], "a geometry edit"
        )
        artifact, receipt = transform_trusted_molecular_geometry(
            approved_workspace=self.approved_workspace,
            edited_artifact_id=edited_artifact_id,
            parent=parent,
            parent_identity_sha256=identity_sha256,
            operation=str(values["operation"]),
            atoms=[int(item) for item in values.get("atoms", ())],
            moving_side_atom=int(values["moving_side_atom"]),
            target_value=float(values["target_value"]),
        )
        self.artifacts[artifact.artifact_id] = artifact
        self.geometry_edits[artifact.sha256] = receipt
        self.event_store.append(
            turn_id=turn_id,
            kind=EventKind.MOLECULAR_GEOMETRY_EDITED.value,
            payload={
                "receipt_sha256": receipt.receipt_sha256,
                "edited_artifact_id": artifact.artifact_id,
                "edited_artifact_sha256": artifact.sha256,
                "parent_sha256": parent.sha256,
                "operation": receipt.operation,
                "coordinate_atoms": list(receipt.coordinate_atoms),
                "moving_side_atom": receipt.moving_side_atom,
                "moved_atoms": list(receipt.moved_atoms),
                "value_unit": receipt.value_unit,
                "value_before": receipt.value_before,
                "value_requested": receipt.value_requested,
                "value_achieved": receipt.value_achieved,
                "connectivity_changed": receipt.connectivity_changed,
            },
            idempotency_key=("geometry-edit:" + receipt.receipt_sha256),
        )
        return {
            "geometry_edit": receipt,
            "artifact": artifact,
            "next_action": (
                "bind charge and multiplicity explicitly with "
                "bind_scientific_identity -- an edit does not change or infer "
                "electronic state; the edited geometry is a starting "
                "structure, and the stage that optimises it is a new workflow "
                "needing its own review, where the coordinate you requested "
                "can be measured against the relaxed result"
            ),
        }

    def _break_symmetry(self, turn_id: str, values: dict) -> Any:
        """Perturb every atom of an identity-bound geometry by seed.

        A source geometry carries its builder's symmetry (six live
        cases: idealised D4h and threefold-rotor starts optimised onto
        saddles). The model names the seed and amplitude; the host draws
        the step, records the largest displacement it took, and owns the
        bytes. Refusals are structural only.
        """

        if self.approved_workspace is None:
            raise ContractError(
                "a symmetry break requires an approved workspace to write "
                "into"
            )
        perturbed_artifact_id = self._unused_artifact_id(
            values["perturbed_artifact_id"]
        )
        parent, identity_sha256 = self._identity_bound_parent(
            values["input_artifact_id"], "break_symmetry"
        )
        artifact, receipt = break_trusted_molecular_symmetry(
            approved_workspace=self.approved_workspace,
            perturbed_artifact_id=perturbed_artifact_id,
            parent=parent,
            parent_identity_sha256=identity_sha256,
            seed=values["seed"],
            amplitude_angstrom=float(values["amplitude_angstrom"]),
        )
        self.artifacts[artifact.artifact_id] = artifact
        self.symmetry_breaks[artifact.sha256] = receipt
        self.event_store.append(
            turn_id=turn_id,
            kind=EventKind.SYMMETRY_BROKEN.value,
            payload={
                "receipt_sha256": receipt.receipt_sha256,
                "perturbed_artifact_id": artifact.artifact_id,
                "perturbed_artifact_sha256": artifact.sha256,
                "parent_artifact_id": parent.artifact_id,
                "seed": receipt.seed,
                "amplitude_angstrom": receipt.amplitude_angstrom,
                "max_displacement_angstrom": (
                    receipt.max_displacement_angstrom
                ),
                "point_group_before": receipt.point_group_before,
                "point_group_after": receipt.point_group_after,
                "connectivity_changed": receipt.connectivity_changed,
            },
            idempotency_key=("symmetry-break:" + receipt.receipt_sha256),
        )
        return {
            "symmetry_break": receipt,
            "artifact": artifact,
            "next_action": (
                "bind charge and multiplicity on the perturbed artifact, "
                "then plan the optimisation that decides whether the step "
                "escaped the saddle"
            ),
        }

    def _bind_reached_geometry(self, turn_id: str, values: dict) -> Any:
        """Carry a run's reached structure forward as a starting geometry.

        The repair menu named this route for the whole campaign and no
        tool walked it, so a session facing an optimisation that ran out
        of iterations had to rebuild a start by hand and pay engine calls
        for the difference. The host owns the bytes and the lineage; the
        session owns whether restarting from there is the right science.
        """

        if self.approved_workspace is None:
            raise ContractError(
                "carrying a reached geometry requires an approved "
                "workspace to write into"
            )
        reached_artifact_id = self._unused_artifact_id(
            values["reached_artifact_id"]
        )
        artifact, receipt = build_reached_geometry(
            approved_workspace=self.approved_workspace,
            reached_artifact_id=reached_artifact_id,
            result_artifact=self._artifact(values["artifact_id"]),
            program=str(values["program"]),
            run_evidence_root=self.run_evidence_root,
        )
        self.artifacts[artifact.artifact_id] = artifact
        self.reached_geometries[receipt.receipt_sha256] = receipt
        self._emit(
            turn_id,
            EventKind.REACHED_GEOMETRY_BOUND,
            receipt.receipt_sha256,
            reached_artifact_id=artifact.artifact_id,
            source_result_artifact_id=receipt.source_result_artifact_id,
            recorded_terminal_state=receipt.recorded_terminal_state,
            record=canonical_data(receipt),
        )
        return {
            "schema_version": "chemsmart.reached-geometry.v1",
            "artifact": artifact,
            "reached_geometry": canonical_data(receipt),
            "next_action": (
                "bind this geometry's charge and multiplicity with "
                "bind_scientific_identity -- a reached structure carries "
                "no electronic state -- then plan the stage that "
                "optimises it as a new workflow for review"
            ),
        }

    def _characterise_stationary_point(
        self, turn_id: str, values: dict
    ) -> Any:
        """Say what a result is, checked against its own printed modes.

        A search that missed its promise keeps the failure; this records
        the other true statement about the same bytes, so a number taken
        from them can be delivered as the stationary point it belongs to
        rather than as the one that was asked for. The host owns the
        convention and refuses a claim the frequencies do not support;
        the session owns the claim and everything it means.
        """

        receipt = build_stationary_point_characterisation(
            result_artifact=self._artifact(values["result_artifact_id"]),
            program=str(values["program"]),
            order_claimed=int(values["order_claimed"]),
            node_id=str(values.get("node_id") or ""),
            anomaly_sha256=str(values.get("anomaly_receipt_sha256") or ""),
        )
        self.stationary_point_characterisations[receipt.receipt_sha256] = (
            receipt
        )
        self._emit(
            turn_id,
            EventKind.STATIONARY_POINT_CHARACTERISED,
            receipt.receipt_sha256,
            node_id=receipt.node_id,
            order_claimed=receipt.order_claimed,
            record=canonical_data(receipt),
        )
        return {"stationary_point_characterisation": canonical_data(receipt)}

    def _displace_along_vibrational_mode(
        self, turn_id: str, values: dict
    ) -> Any:
        """Step a completed result's geometry along one of its own modes.

        This is the move a chemist makes when an optimisation converges onto
        a saddle: read the offending mode, step along it, relax again. The
        displacement vectors come from the program's own output and the host
        owns the arithmetic; the model owns which mode and how far.

        Nothing here decides whether the step was a good one. A displaced
        geometry is a starting structure, and the optimisation that consumes
        it is what grades the choice -- the same discipline every other
        geometry-producing operation is held to.
        """

        if self.approved_workspace is None:
            raise ContractError(
                "a mode displacement requires an approved workspace to "
                "write into"
            )
        displaced_artifact_id = self._unused_artifact_id(
            values["displaced_artifact_id"]
        )
        result_artifact = self._artifact(values["result_artifact_id"])
        artifact, receipt = displace_trusted_geometry_along_mode(
            approved_workspace=self.approved_workspace,
            displaced_artifact_id=displaced_artifact_id,
            result_artifact=result_artifact,
            program=str(values["program"]),
            mode_index=int(values["mode_index"]),
            amplitude_angstrom=float(values["amplitude_angstrom"]),
        )
        self.artifacts[artifact.artifact_id] = artifact
        self.mode_displacements[artifact.sha256] = receipt
        self.event_store.append(
            turn_id=turn_id,
            kind=EventKind.GEOMETRY_DISPLACED_ALONG_MODE.value,
            payload={
                "receipt_sha256": receipt.receipt_sha256,
                "displaced_artifact_id": artifact.artifact_id,
                "displaced_artifact_sha256": artifact.sha256,
                "result_sha256": receipt.result_sha256,
                "mode_index": receipt.mode_index,
                "mode_frequency_cm_1": receipt.mode_frequency_cm_1,
                "mode_is_imaginary": receipt.mode_is_imaginary,
                "amplitude_angstrom": receipt.amplitude_angstrom,
                "achieved_max_displacement_angstrom": (
                    receipt.achieved_max_displacement_angstrom
                ),
                "leading_atoms": list(receipt.leading_atoms),
                "connectivity_changed": receipt.connectivity_changed,
            },
            idempotency_key=("mode-displacement:" + receipt.receipt_sha256),
        )
        return {
            "mode_displacement": receipt,
            "artifact": artifact,
            "next_action": (
                "bind charge and multiplicity explicitly with "
                "bind_scientific_identity -- a displacement does not change "
                "or infer electronic state; the displaced geometry is a "
                "starting structure, and the optimisation that relaxes it is "
                "a new workflow needing its own review, where whether the "
                "step escaped the saddle becomes an observation"
            ),
        }

    def _append_molecular_atom(self, turn_id: str, values: dict) -> Any:
        """Add one atom to an identity-bound parent by its coordinates.

        Derivation's mirror. Taking a hydrogen off gives a radical or an
        anion depending on where its electron went; putting one on gives a
        cation or a radical depending on whether it brought one, so the
        appended species binds its charge and multiplicity explicitly too.
        """

        if self.approved_workspace is None:
            raise ContractError(
                "an atom append requires an approved workspace to write into"
            )
        appended_artifact_id = self._unused_artifact_id(
            values["appended_artifact_id"]
        )
        parent, identity_sha256 = self._identity_bound_parent(
            values["input_artifact_id"], "an atom append"
        )
        artifact, receipt = append_trusted_molecular_atom(
            approved_workspace=self.approved_workspace,
            appended_artifact_id=appended_artifact_id,
            parent=parent,
            parent_identity_sha256=identity_sha256,
            element=str(values["element"]),
            anchor_atom=int(values["anchor_atom"]),
            angle_atom=int(values["angle_atom"]),
            dihedral_atom=int(values["dihedral_atom"]),
            bond_length_angstrom=float(values["bond_length_angstrom"]),
            angle_degrees=float(values["angle_degrees"]),
            dihedral_degrees=_periodic_degrees(
                float(values["dihedral_degrees"])
            ),
        )
        self.artifacts[artifact.artifact_id] = artifact
        self.atom_appends[artifact.sha256] = receipt
        self.event_store.append(
            turn_id=turn_id,
            kind=EventKind.MOLECULAR_ATOM_APPENDED.value,
            payload={
                "receipt_sha256": receipt.receipt_sha256,
                "appended_artifact_id": artifact.artifact_id,
                "appended_artifact_sha256": artifact.sha256,
                "parent_sha256": parent.sha256,
                "element": receipt.element,
                "appended_atom_index": receipt.appended_atom_index,
                "anchor_atoms": list(receipt.anchor_atoms),
                "bond_length_angstrom": receipt.bond_length_angstrom,
                "angle_degrees": receipt.angle_degrees,
                "dihedral_degrees": receipt.dihedral_degrees,
                "formula": receipt.formula,
                "fragment_count": receipt.fragment_count,
            },
            idempotency_key=("atom-append:" + receipt.receipt_sha256),
        )
        return {
            "atom_append": receipt,
            "artifact": artifact,
            "next_action": (
                "bind charge and multiplicity explicitly with "
                "bind_scientific_identity -- appending does not infer "
                "electronic state, and whether the added atom brought an "
                "electron decides it; the appended geometry is a starting "
                "structure, and the stage that consumes it is a new workflow "
                "needing its own review"
            ),
        }

    def _inspect_database_records(self, turn_id: str, values: dict) -> Any:
        """Enumerate a workspace database's records as observations.

        A .db row is provenance, not authority: it may store charge,
        multiplicity, and energy, or store no electronic state at all
        (the ASE-kin row guarantees only geometry and metadata).  This
        tool therefore reports stored fields under an explicit
        observation role and never turns them into bindings; the
        recorded, evidence-bearing act is the extraction that follows.
        """

        from chemsmart.database.database import Database
        from chemsmart.database.query import DatabaseQuery

        artifact = self._artifact(values["database_artifact_id"])
        if artifact.kind != "chemsmart_db":
            raise ContractError(
                f"artifact {artifact.artifact_id!r} is {artifact.kind!r}, "
                "not a chemsmart_db database"
            )
        database = Database(artifact.path)
        total = int(database.count_records())
        limit = int(values.get("limit") or 50)
        query_string = str(values.get("query") or "").strip()
        if query_string:
            try:
                rows = DatabaseQuery(
                    artifact.path, query_string, target="records"
                ).query_summaries()
            except ValueError as exc:
                raise ContractError(str(exc)) from exc
            detail_records = [
                database.get_record(record_id=row.get("record_id"))
                for row in rows[:limit]
            ]
            matched_count = len(rows)
        else:
            detail_records = database.get_all_records()[:limit]
            matched_count = total
        records_payload = []
        for record in detail_records:
            if record is None:
                continue
            structures = list(record.get("molecules") or ())
            meta = record.get("meta") or {}
            records_payload.append(
                {
                    "record_index": record.get("record_index"),
                    "record_id": record.get("record_id"),
                    "method": meta.get("method"),
                    "basis": meta.get("basis"),
                    "structure_count": len(structures),
                    "structures": [
                        {
                            "structure_index": position,
                            "formula": entry.get("chemical_formula"),
                            "atom_count": entry.get("number_of_atoms"),
                            "stored_charge": entry.get("charge"),
                            "stored_multiplicity": entry.get("multiplicity"),
                            "stored_energy": entry.get("energy"),
                            "stored_is_optimized": bool(
                                entry.get("is_optimized_structure")
                            ),
                        }
                        for position, entry in enumerate(structures, start=1)
                    ],
                }
            )
        return {
            "schema_version": "chemsmart.database-records-inspection.v1",
            "database_artifact_id": artifact.artifact_id,
            "record_count": total,
            "matched_count": matched_count,
            "returned_records": len(records_payload),
            "query": query_string,
            "records": records_payload,
            "stored_fields_role": (
                "stored_charge, stored_multiplicity, and stored_energy are "
                "observations from the database's own provenance, never "
                "identity bindings; a record may store no electronic state "
                "at all"
            ),
            "next_action": (
                "admit a record's geometry with "
                "extract_database_record_geometry, then bind charge and "
                "multiplicity explicitly with bind_scientific_identity"
            ),
        }

    def _extract_database_record_geometry(
        self, turn_id: str, values: dict
    ) -> Any:
        """Copy one record's coordinates into a host-owned geometry."""

        if self.approved_workspace is None:
            raise ContractError(
                "extraction requires an approved workspace to write into"
            )
        extracted_artifact_id = str(values["extracted_artifact_id"])
        if extracted_artifact_id in self.artifacts:
            taken = sorted(self.artifacts)
            raise ContractError(
                f"artifact ID {extracted_artifact_id!r} is already "
                f"registered; choose one not in {taken}"
            )
        database = self._artifact(values["database_artifact_id"])
        record_index = values.get("record_index")
        record_id = values.get("record_id")
        structure_index = values.get("structure_index")
        artifact, receipt = extract_trusted_database_record_geometry(
            approved_workspace=self.approved_workspace,
            extracted_artifact_id=extracted_artifact_id,
            database_artifact=database,
            record_index=(None if record_index is None else int(record_index)),
            record_id=None if record_id is None else str(record_id),
            structure_index=(
                None if structure_index is None else int(structure_index)
            ),
        )
        self.artifacts[artifact.artifact_id] = artifact
        self.database_extractions[artifact.sha256] = receipt
        self.event_store.append(
            turn_id=turn_id,
            kind=EventKind.DATABASE_RECORD_EXTRACTED.value,
            payload={
                "receipt_sha256": receipt.receipt_sha256,
                "extracted_artifact_id": artifact.artifact_id,
                "extracted_artifact_sha256": artifact.sha256,
                "database_sha256": receipt.database_sha256,
                "record_id": receipt.record_id,
                "record_index": receipt.record_index,
                "structure_index": receipt.structure_index,
                "formula": receipt.formula,
            },
            idempotency_key=("database-extraction:" + receipt.receipt_sha256),
        )
        return {
            "extraction": receipt,
            "artifact": artifact,
            "stored_state_observation": {
                "charge": receipt.stored_charge,
                "multiplicity": receipt.stored_multiplicity,
                "energy": receipt.stored_energy,
                "is_optimized": receipt.stored_is_optimized,
            },
            "next_action": (
                "bind charge and multiplicity explicitly with "
                "bind_scientific_identity -- the stored values above are "
                "observations from the database's provenance, not "
                "bindings; state the electronic state you intend for this "
                "calculation"
            ),
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
        # Read the geometry before binding anything.  The state about to be
        # bound is checked against this molecule's electron count, and the
        # same read supplies the facts returned below -- so the refusal
        # cannot be swallowed by the best-effort guard around those facts.
        try:
            from chemsmart.io.molecules.structure import Molecule

            molecule = Molecule.from_filepath(artifact.path)
            symbols = tuple(molecule.chemical_symbols)
        except Exception:
            molecule = None
            symbols = ()
        if symbols:
            refuse_impossible_electronic_state(
                symbols,
                values["charge"],
                values["multiplicity"],
                context="scientific identity",
            )
        binding = build_scientific_identity_binding(
            task_spec_sha256=task_spec_sha256,
            geometry_artifact=artifact,
            charge=values["charge"],
            multiplicity=values["multiplicity"],
        )
        self.scientific_identities[binding.binding_sha256] = binding
        # The bare binding told a session nothing about the molecule it had
        # just bound, and nothing pointed geometry artifacts at the
        # measurement operations that already exist -- in a four-session
        # observation exactly one session discovered that a bound geometry
        # is measurable and verified its own conformer labels; the others
        # assumed and disclosed. Surface the facts and the route here, where
        # every session already looks. The binding itself stays digest-frozen
        # inside the wrapper.
        try:
            geometry_facts = {
                "atom_count": len(symbols),
                "formula": molecule.get_chemical_formula(),
                "symbols": symbols,
            }
        except Exception:
            geometry_facts = {}
        return {
            "schema_version": "chemsmart.scientific-identity-bound.v1",
            "binding_sha256": binding.binding_sha256,
            "binding": binding,
            "geometry": geometry_facts,
            "observations": self._symmetry_observations(artifact, molecule),
            "measurement_route": (
                "this geometry_xyz artifact is readable without any engine: "
                "extract_result_quantities with program 'xyz' yields "
                "positions and symbols, and evaluate_quantity_expression "
                "offers distance, angle, and dihedral operations over those "
                "positions -- measure a coordinate before assuming its value"
            ),
        }

    def _root_artifact_for(
        self, node_id: str, input_artifact: TrustedArtifactRefV1
    ) -> TrustedArtifactRefV1 | None:
        """The goal's original bound geometry behind this node's input.

        Under a frozen approval the binding says which bytes; the
        executor located them beside the initial artifacts. In a planning
        host the lineage is walked: structural hops to their parents, and
        a reached or displaced geometry to the input of the recorded run
        that produced its source. None when the input is the goal's own
        start or the root bytes are not in hand.
        """

        approval = getattr(self, "workflow_execution_approval", None)
        if approval is not None:
            for binding in approval.node_bindings:
                if binding.node_id != node_id:
                    continue
                root_sha256 = str(
                    getattr(binding, "root_artifact_sha256", "") or ""
                )
                if not root_sha256:
                    return None
                known = self.artifacts.get(
                    str(getattr(binding, "root_artifact_id", "") or "")
                )
                if known is not None and known.sha256 == root_sha256:
                    return known
                return None
            return None
        return self._goal_root_of(input_artifact)

    def _goal_root_of(
        self, artifact: TrustedArtifactRefV1
    ) -> TrustedArtifactRefV1 | None:
        reached_by_sha = {
            str(getattr(item, "reached_artifact_sha256", "")): item
            for item in self.reached_geometries.values()
        }
        cursor = artifact.sha256
        seen: set[str] = set()
        root_sha256 = ""
        while cursor and cursor not in seen and len(seen) < 64:
            seen.add(cursor)
            # The by-parent hops come from `_GEOMETRY_HOPS` rather than
            # a second hand-list: this walk carried three of the four
            # and the review walk carried all four, which is how they
            # drifted. What stays local is the *purpose* -- this walk
            # resolves a result-predecessor through the run that
            # produced it, because it is looking for the goal's
            # original bound geometry, not for a chain to display.
            for kind, registry_name, predecessor in self._GEOMETRY_HOPS:
                if predecessor != "parent_sha256":
                    continue
                registry = getattr(self, registry_name, None)
                receipt = (
                    registry.get(cursor)
                    if isinstance(registry, dict)
                    else None
                )
                if receipt is not None:
                    cursor = str(getattr(receipt, predecessor, "") or "")
                    break
            else:
                displaced = self.mode_displacements.get(cursor)
                reached = reached_by_sha.get(cursor)
                source = (
                    str(getattr(displaced, "result_sha256", "") or "")
                    if displaced is not None
                    else (
                        str(getattr(reached, "source_result_sha256", "") or "")
                        if reached is not None
                        else ""
                    )
                )
                if not source:
                    break
                run_input = self._recorded_run_input_sha256(source)
                if not run_input:
                    break
                cursor = run_input
                root_sha256 = cursor
        if not root_sha256 or root_sha256 == artifact.sha256:
            return None
        for known in self.artifacts.values():
            if known.sha256 == root_sha256 and known.kind == "geometry_xyz":
                return known
        located = self._locate_workspace_bytes(root_sha256)
        if located is None:
            return None
        return TrustedArtifactRefV1(
            artifact_id=f"root-{root_sha256[:8]}",
            kind="geometry_xyz",
            sha256=root_sha256,
            size_bytes=located.stat().st_size,
            path=str(located),
            cli_value=str(located),
        )

    def _recorded_run_input_sha256(self, result_sha256: str) -> str:
        """The input digest of the recorded run whose outputs include the
        given result, read from this workspace's own run streams."""

        root = self.run_evidence_root
        if root is None:
            return ""
        records = root / ".chemsmart-agent"
        for pattern in (
            "goals/*/runs/*/events.jsonl",
            "runs/*/events.jsonl",
            "executions/*/events.jsonl",
        ):
            for path in sorted(records.glob(pattern)):
                if path.is_symlink() or not path.is_file():
                    continue
                try:
                    lines = path.read_text(encoding="utf-8").splitlines()
                except OSError:
                    continue
                for line in lines:
                    if result_sha256 not in line:
                        continue
                    try:
                        event = json.loads(line)
                    except ValueError:
                        continue
                    if event.get("kind") != EventKind.PROGRAM_EXECUTED.value:
                        continue
                    record = (event.get("payload") or {}).get("record") or {}
                    outputs = record.get("output_artifacts") or ()
                    if any(
                        str((item or {}).get("sha256") or "") == result_sha256
                        for item in outputs
                        if isinstance(item, Mapping)
                    ):
                        return str(record.get("input_artifact_sha256") or "")
        return ""

    def _locate_workspace_bytes(self, sha256: str) -> Path | None:
        roots = [
            root
            for root in (self.approved_workspace, self.run_evidence_root)
            if root is not None
        ]
        for root in roots:
            for candidate in sorted(Path(root).rglob("*.xyz")):
                if candidate.is_symlink() or not candidate.is_file():
                    continue
                try:
                    if file_sha256(candidate) == sha256:
                        return candidate
                except OSError:
                    continue
        return None

    def _symmetry_observations(
        self, artifact: Any, molecule: Any
    ) -> tuple[str, ...]:
        """What the host can say about a geometry's symmetry and its
        builder's idealised coordinates, at the moment it is bound."""

        from chemsmart.agent.symmetry import (
            idealised_coordinate_observation,
            idealised_internal_coordinate_count,
            symmetry_observation,
        )

        observations: list[str] = []
        if molecule is None:
            return ()
        try:
            observations.append(
                symmetry_observation(
                    tuple(molecule.chemical_symbols), molecule.positions
                )
            )
        except Exception:
            return ()
        appends = []
        edits = []
        cursor = artifact.sha256
        seen: set[str] = set()
        while cursor and cursor not in seen:
            seen.add(cursor)
            edit_receipt = self.geometry_edits.get(cursor)
            if edit_receipt is not None:
                # Collected, not stepped past: a coordinate *edited* onto
                # the 60° lattice starts on the same saddle as one an
                # append placed there, and the walk saw every edit on its
                # way to the appends and kept none of them.
                edits.append(edit_receipt)
                cursor = edit_receipt.parent_sha256
                continue
            append_receipt = self.atom_appends.get(cursor)
            if append_receipt is not None:
                appends.append(append_receipt)
                cursor = append_receipt.parent_sha256
                continue
            symmetry_receipt = self.symmetry_breaks.get(cursor)
            if symmetry_receipt is not None:
                cursor = symmetry_receipt.parent_sha256
                continue
            break
        if appends or edits:
            observations.append(
                idealised_coordinate_observation(
                    idealised_internal_coordinate_count(
                        appends, edit_receipts=edits
                    )
                )
            )
        return tuple(observations)

    def _inspect_program(self, turn_id: str, values: dict) -> Any:
        """Capability and environment in one call; every receipt returned."""

        capability = self._inspect_program_capability(turn_id, values)
        environment = self._inspect_program_environment(
            turn_id, {"capability_receipt_sha256": capability.receipt_sha256}
        )
        record = (
            dict(environment)
            if isinstance(environment, Mapping)
            else {"environment": environment}
        )
        return {"capability": capability, **record}

    def _project_yaml(self, turn_id: str, values: dict) -> Any:
        """One project tool; the action names which of the five steps."""

        from chemsmart.agent.tool_specs import PROJECT_YAML_ACTION_ARGUMENTS

        action = str(values.get("action") or "").strip()
        if action not in PROJECT_YAML_ACTION_ARGUMENTS:
            raise ContractError(
                f"project_yaml action must be one of "
                f"{sorted(PROJECT_YAML_ACTION_ARGUMENTS)}, not {action!r}"
            )
        needed = PROJECT_YAML_ACTION_ARGUMENTS[action]
        missing = [name for name in needed if name not in values]
        if missing:
            raise ContractError(
                f"project_yaml(action={action}) needs {list(needed)}; "
                f"missing {missing}"
            )
        arguments = {name: values[name] for name in needed}
        handler = {
            "establish": self._establish_project,
            "render": self._render_project_yaml,
            "promote": self._promote_project_yaml,
            "read": self._read_project_yaml,
            "validate": self._validate_project_yaml,
        }[action]
        return handler(turn_id, arguments)

    def _inspect_run(self, turn_id: str, values: dict) -> Any:
        """A run's typed outcome, or one finished result's selectors."""

        if values.get("artifact_id"):
            if not values.get("program"):
                raise RoutedContractError(
                    gate="inspect.result_needs_its_reader",
                    invariant=(
                        "a result's selectors are read by the reader of the "
                        "program that wrote it."
                    ),
                    diagnosis="inspect_run named an artifact_id and no program.",
                    route=(
                        "pass program beside artifact_id (orca, gaussian, "
                        "xtb or pyscf), or name a run to read its outcome."
                    ),
                )
                raise ContractError(  # pragma: no cover - unreachable
                    "inspect_run with artifact_id also needs program, the "
                    "reader that opens the result"
                )
            return self._inspect_result_selectors(
                turn_id,
                {
                    "program": values["program"],
                    "artifact_id": values["artifact_id"],
                },
            )
        return self._inspect_run_outcome(
            turn_id, {"run": str(values.get("run") or "")}
        )

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
            raise _artifact_id_taken(
                values["artifact_id"], sorted(self.artifacts)
            )
        render = self._get(
            self.project_renders,
            values["render_receipt_sha256"],
            "project render receipt",
        )
        earlier = tuple(
            (artifact_id, self.project_renders[item.render_receipt_sha256])
            for artifact_id, item in self.project_promotions.items()
            if item.render_receipt_sha256 in self.project_renders
        )
        observations = promotion_field_observations(render, earlier)
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
        return {
            "artifact": artifact,
            "promotion": promotion,
            "observations": observations,
        }

    def _establish_project(self, turn_id: str, values: dict) -> Any:
        """Render, promote, and validate one project in a single turn.

        These three are never used apart. A node needs all of them before it
        can be prepared, always in this order, and each node needs its own --
        so the ceremony scales with the size of the graph. Measured across
        three live sessions it took 34%, 56% and 61% of the whole tool budget,
        and the paper task ran out of turns inside it, having reasoned about
        the chemistry correctly and never reached a reviewable workflow.

        Nothing is skipped or weakened: each step is the same handler, called
        in the same order, emitting the same events, and every receipt is
        returned. What changes is that a five-node graph spends five turns here
        instead of fifteen.
        """

        rendered = self._render_project_yaml(
            turn_id,
            {"program": values["program"], "sections": values["sections"]},
        )
        promoted = self._promote_project_yaml(
            turn_id,
            {
                "render_receipt_sha256": rendered.receipt_sha256,
                "artifact_id": values["artifact_id"],
            },
        )
        validated = self._validate_project_yaml(
            turn_id,
            {
                "project_artifact_id": values["artifact_id"],
                "capability_receipt_sha256": values[
                    "capability_receipt_sha256"
                ],
            },
        )
        return {
            "schema_version": "chemsmart.project-establishment.v1",
            "rendered": rendered,
            "promoted": promoted,
            "validated": validated,
        }

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
            command_result = self._scientific_toolchain_command_results.get(
                plan.plan_sha256
            )
            if command_result is None:
                # The executor host seeds the approved toolchain plan without
                # any planning-session command result; the relational
                # observation below exists to keep a *planning* session from
                # binding a stale project, which cannot happen here.
                continue
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
                    "amend_scientific_workflow"
                ),
            }
        return {"status": "no_scientific_workflow_planned"}

    _RUN_RECEIPT_KINDS = frozenset(
        {
            "result_quantities_extracted",
            "thermochemistry_derived",
            "quantity_expression_evaluated",
            "scientific_validation_evaluated",
            "analysis_claims_recorded",
        }
    )

    def _recorded_run_receipt(self, receipt_sha256: str) -> bool:
        """Whether a recorded run of this workspace minted the receipt.

        A woken session reads the previous cycle's executed chain through
        inspect_run and was refused when its decision cited one of that
        chain's receipts: the session host held only its own (NOVEL-3
        po1 and ino2, 2026-09-05). The run streams are the host's own
        durable record; a receipt they carry is one the host minted.
        """

        root = getattr(self, "run_evidence_root", None)
        if not root:
            return False
        needle = f'"receipt_sha256": "{receipt_sha256}"'
        for stream in sorted(
            Path(root).glob(".chemsmart-agent/goals/*/runs/*/events.jsonl")
        ):
            try:
                text = stream.read_text(encoding="utf-8")
            except OSError:
                continue
            if needle not in text:
                continue
            for line in text.splitlines():
                if needle not in line:
                    continue
                try:
                    event = json.loads(line)
                except json.JSONDecodeError:
                    continue
                if event.get("kind") in self._RUN_RECEIPT_KINDS:
                    return True
        return False

    def _record_scientific_decision(self, turn_id: str, values: dict) -> Any:
        task_spec_sha256 = self._resolve_task_spec_reference(
            values, "task_spec_sha256"
        )
        if task_spec_sha256 not in self.task_spec_sha256s:
            raise ContractError(
                "scientific decision targets an unknown task spec"
            )
        # Derived from the one receipt authority rather than hand-listed
        # a sixth time: the five copies had drifted and this gate refused
        # receipts the host itself minted.
        postprocessing_registries = tuple(self._receipt_registries().values())
        postprocessing_receipt_sha256s = tuple(
            str(item)
            for item in values.get("postprocessing_receipt_sha256s", ())
        )
        for receipt_sha256 in postprocessing_receipt_sha256s:
            require_sha256(receipt_sha256, "postprocessing_receipt_sha256")
            if not any(
                receipt_sha256 in receipts
                for receipts in postprocessing_registries
            ) and not self._recorded_run_receipt(receipt_sha256):
                raise RoutedContractError(
                    gate="decision.receipt_is_one_the_host_minted",
                    invariant=(
                        "a decision cites only receipts this session or a "
                        "recorded run of this workspace minted."
                    ),
                    diagnosis=(
                        f"{receipt_sha256[:8]} is "
                        + (
                            self._digest_names(receipt_sha256)
                            or "no digest this host minted"
                        )
                        + ", so it cannot stand as postprocessing "
                        "evidence here."
                    ),
                    route=(
                        "cite the receipt_sha256 a tool returned in this "
                        "session, or one inspect_run shows on a recorded "
                        "run of this goal"
                    ),
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
            prefix = "anomaly:"
            if reference.startswith(prefix):
                # A decision cites only anomalies the host recorded, by
                # digest: nothing enters the record that no receipt
                # backs, and the interpretation stays prose beside it.
                receipt_sha256 = reference[len(prefix) :]
                require_sha256(receipt_sha256, "anomaly_receipt_sha256")
                known = set(self.anomaly_observations) | {
                    str(item.get("receipt_sha256") or "")
                    for item in self.prior_anomaly_observations
                }
                if receipt_sha256 not in known:
                    raise ContractError(
                        "scientific decision cites an unknown anomaly "
                        "observation"
                    )
                continue
            prefix = "doubt:"
            if reference.startswith(prefix):
                # A typed doubt: the session names the exact receipt it
                # distrusts. The completion gate intersects these digests
                # with the rendered claims' supporting receipts, so a
                # doubted number cannot ship under a clean certification.
                receipt_sha256 = reference[len(prefix) :]
                require_sha256(receipt_sha256, "doubted_receipt_sha256")
                if not any(
                    receipt_sha256 in receipts
                    for receipts in postprocessing_registries
                ):
                    raise ContractError(
                        "scientific decision doubts an unknown "
                        "postprocessing receipt"
                    )
                continue
            parsed_postprocessing_ref = _postprocessing_evidence_reference(
                reference
            )
            if parsed_postprocessing_ref is not None:
                receipt_kind, receipt_sha256 = parsed_postprocessing_ref
                # A typed prefix names one role and stays narrow: a
                # reader asking for a thermochemistry receipt means that
                # kind. An untyped citation asks the authenticity
                # question instead, and that one is the whole ledger --
                # a completion receipt and a stationary-point
                # characterisation are the host's own and were refused
                # here (SUFFICIENCY-5, 2026-09-10).
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
                    known = self._resolve_receipt(receipt_sha256) is not None
                else:
                    known = receipt_sha256 in registries[receipt_kind]
                # A receipt an earlier cycle's executed chain minted is
                # the host's own; the run stream carries it.
                if not known and self._recorded_run_receipt(receipt_sha256):
                    known = True
                if not known:
                    raise RoutedContractError(
                        gate="decision.receipt_is_one_the_host_minted",
                        invariant=(
                            "a decision cites only receipts this session "
                            "or a recorded run of this workspace minted."
                        ),
                        diagnosis=(
                            f"{receipt_sha256[:8]} is "
                            + (
                                self._digest_names(receipt_sha256)
                                or "no digest this host minted"
                            )
                            + f", and this citation asks for "
                            f"{receipt_kind!r}."
                        ),
                        route=(
                            "cite the receipt_sha256 a tool returned in "
                            "this session, or one inspect_run shows on a "
                            "recorded run of this goal"
                        ),
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
        unreachable = self._verify_unreachable_observables(
            values.get("unreachable_observable_ids") or ()
        )
        self.verified_unreachable_ids.update(
            str(item.get("observable_id") or "")
            for item in unreachable
            if item.get("verified")
        )
        dispositions = self._verify_menu_route_dispositions(
            values.get("menu_route_dispositions") or ()
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
        extra: dict[str, Any] = {}
        if unreachable:
            extra["unreachable_observables"] = unreachable
        if dispositions:
            extra["menu_route_dispositions"] = dispositions
        self._emit(
            turn_id,
            EventKind.SCIENTIFIC_DECISION_RECORDED,
            record.record_sha256,
            status="recorded",
            decision_id=record.decision_id,
            record=record_body,
            **extra,
        )
        if dispositions and not unreachable:
            return {**canonical_data(record), **extra}
        if unreachable:
            return {
                **canonical_data(record),
                **extra,
                "settlement_consequence": (
                    "a refusal the host verified settles the goal "
                    "unreachable_from_evidence once every other declared "
                    "observable is delivered by id; an unverified one "
                    "returns the goal to the human naming it"
                ),
            }
        return record

    def _verify_unreachable_observables(
        self, entries: Sequence[Mapping[str, Any]]
    ) -> tuple[dict[str, Any], ...]:
        """Check each typed refusal against what the host can know.

        Two sessions wrote honest refusals and neither reached the word:
        po3 said "typed refusal" in prose and recorded nothing typed;
        ino3 built five blocked nodes and the settlement joined
        blocked-ness to required output ids, not observable ids (NOVEL-3,
        2026-09-05). The owner ruled (2026-09-06) the refusal is host-
        verified: the session names the observable, the producer it
        would need and the receipts that show the gap; the host checks
        what it can -- a selector no envelope program declares for the
        job type, or a blocked_unsupported node in this session's plan
        whose output is the observable -- and records the basis. A guard
        against refusing one's way out stays: an unverified refusal
        returns to the human, never settles.
        """

        from chemsmart.analysis.result_readers import (
            reader_for,
            registered_reader_programs,
        )

        declared = set(self.requested_observable_declarations)
        envelope = self.bounded_execution_envelope
        programs = (
            tuple(
                str(program)
                for program, _engines in envelope.allowed_program_engines
            )
            if envelope is not None
            else registered_reader_programs()
        )
        verified: list[dict[str, Any]] = []
        for entry in entries:
            observable_id = str(entry.get("observable_id") or "")
            statement = str(entry.get("statement") or "").strip()
            receipts = tuple(
                str(item) for item in (entry.get("receipt_sha256s") or ())
            )
            if observable_id not in declared:
                raise RoutedContractError(
                    gate="decision.unreachable_id_is_declared",
                    invariant=(
                        "a typed refusal names an observable this goal "
                        "declared."
                    ),
                    diagnosis=(
                        f"{observable_id!r} is not among the declared "
                        f"observables {sorted(declared)}."
                    ),
                    route=(
                        "name the declared observable_id verbatim, or "
                        "declare it first with declare_requested_observable"
                    ),
                )
            if not statement or not receipts:
                raise RoutedContractError(
                    gate="decision.refusal_carries_its_evidence",
                    invariant=(
                        "a typed refusal states the producer it would need "
                        "and cites at least one receipt that shows the gap."
                    ),
                    diagnosis=(
                        f"{observable_id!r} was named with "
                        + ("no statement" if not statement else "no receipt")
                        + "."
                    ),
                    route=(
                        "cite the receipt of the probe that showed the "
                        "producer absent -- an inspect_run or "
                        "extract_result_quantities receipt -- and say what "
                        "producer the observable needs"
                    ),
                )
            for receipt_sha256 in receipts:
                require_sha256(receipt_sha256, "refusal receipt_sha256")
                if not self._receipt_known(receipt_sha256):
                    raise RoutedContractError(
                        gate="decision.receipt_is_one_the_host_minted",
                        invariant=(
                            "a refusal cites only receipts this session or "
                            "a recorded run of this workspace minted."
                        ),
                        diagnosis=(
                            f"{receipt_sha256[:8]} is "
                            + (
                                self._digest_names(receipt_sha256)
                                or "no digest this host minted"
                            )
                            + "."
                        ),
                        route=(
                            "cite the receipt_sha256 a tool returned, or "
                            "one inspect_run shows on a recorded run"
                        ),
                    )
            selector = str(entry.get("selector") or "").strip()
            jobtype = str(entry.get("jobtype") or "").strip().lower()
            blocked_node_id = str(entry.get("blocked_node_id") or "").strip()
            basis = ""
            is_verified = False
            if selector:
                declaring = []
                for program in programs:
                    reader = reader_for(program)
                    if reader is None:
                        continue
                    jobtypes = (
                        (jobtype,)
                        if jobtype
                        else tuple(
                            item[0] for item in reader.jobtype_selectors
                        )
                    )
                    for candidate in jobtypes:
                        declared_here = reader.selectors_for_jobtype(candidate)
                        if (
                            declared_here is not None
                            and selector in declared_here
                        ):
                            declaring.append(f"{program}/{candidate}")
                if declaring:
                    basis = (
                        f"selector {selector!r} is declared by "
                        + ", ".join(sorted(declaring))
                        + "; the observable is reachable and the refusal "
                        "is not verified"
                    )
                else:
                    is_verified = True
                    basis = (
                        f"no program in the envelope ({', '.join(programs)}) "
                        f"declares selector {selector!r}"
                        + (f" for jobtype {jobtype!r}" if jobtype else "")
                    )
            elif blocked_node_id:
                found = False
                for plan in self.scientific_toolchain_plans.values():
                    for node in plan.analysis_nodes:
                        if (
                            node.node_id == blocked_node_id
                            and node.support_state == "blocked_unsupported"
                            and any(
                                output.output_id == observable_id
                                for output in node.outputs
                            )
                        ):
                            found = True
                if found:
                    is_verified = True
                    basis = (
                        f"analysis node {blocked_node_id!r} is declared "
                        "blocked_unsupported in this session's plan and "
                        f"names {observable_id!r} as its output"
                    )
                else:
                    basis = (
                        f"no blocked_unsupported node {blocked_node_id!r} "
                        f"with output {observable_id!r} exists in this "
                        "session's plans; the refusal is not verified"
                    )
            else:
                # A precision no method in this envelope can reach is the
                # third way an observable is unreachable, and the only
                # one where the producer exists and the number was
                # computed. The verifier knew a missing selector and a
                # blocked node and nothing else, so a session facing this
                # would have had to invent an absence to say a true
                # thing -- and the menu offered refusal as a route it
                # could not walk. What the host checks is that the
                # requirement exists and is open on this goal's own
                # record; whether no conceivable calculation could reach
                # it is not claimed, by the host or by the word.
                assessment = self.requirement_assessments.get(observable_id)
                declaration = self.requested_observable_declarations.get(
                    observable_id, {}
                )
                state = str((assessment or {}).get("state") or "")
                if (
                    declaration.get("required_tolerance") is not None
                    and state in _OPEN_SUFFICIENCY_STATES
                ):
                    is_verified = True
                    basis = (
                        "the requirement "
                        f"{declaration['required_tolerance']} "
                        f"{declaration.get('unit', '')} stands open on this "
                        f"goal's own record ({state}) and the session "
                        "states the available evidence does not establish "
                        "it here; the host verifies the open requirement "
                        "and never that no calculation could reach it"
                    ).replace("  ", " ")
                elif declaration.get("required_tolerance") is not None:
                    basis = (
                        f"observable {observable_id!r} carries a required "
                        "tolerance and no open assessment stands against "
                        "it; claim it with the uncertainty you attribute "
                        "to it first, so what is refused is on the record"
                    )
                else:
                    basis = (
                        "no host-checkable producer was named (give "
                        "selector and jobtype, or blocked_node_id); the "
                        "refusal is stated, not verified"
                    )
            verified.append(
                {
                    "observable_id": observable_id,
                    "statement": statement,
                    "selector": selector,
                    "jobtype": jobtype,
                    "blocked_node_id": blocked_node_id,
                    "receipt_sha256s": receipts,
                    "verified": is_verified,
                    "basis": basis,
                }
            )
        return tuple(verified)

    def _verify_menu_route_dispositions(
        self, entries: Sequence[Mapping[str, Any]]
    ) -> tuple[dict[str, Any], ...]:
        """Check each menu disposition names an offered route and minted
        receipts; the choice itself is never graded.

        REACH-1 po3 cycle 3 (2026-09-06) rejected all four routes the
        repair menu offered, each with a mechanism, in prose the next
        cycle never saw: the wake carries the menu and not what was
        done with it. The owner ruled (R3) the disposition is a typed
        field of the decision. The host verifies only what it can
        know -- the route was offered this cycle, the receipts are its
        own -- and records the rest verbatim.
        """

        offered = tuple(self.offered_repair_routes)
        verified: list[dict[str, Any]] = []
        for entry in entries:
            route = str(entry.get("route") or "")
            disposition = str(entry.get("disposition") or "")
            reason = str(entry.get("reason") or "").strip()
            receipts = tuple(
                str(item) for item in (entry.get("receipt_sha256s") or ())
            )
            if route not in offered:
                raise RoutedContractError(
                    gate="decision.route_is_one_the_menu_offered",
                    invariant=(
                        "a disposition names a route the wake's repair "
                        "menu offered this cycle, by the menu's own key."
                    ),
                    diagnosis=(
                        f"{route!r} is not among the offered routes "
                        f"{list(offered)}."
                        if offered
                        else f"{route!r}: this wake carried no repair menu."
                    ),
                    route=(
                        "name one of the offered routes, or carry the "
                        "mechanism as an uncertainty of the decision."
                    ),
                )
            for receipt_sha256 in receipts:
                require_sha256(receipt_sha256, "disposition receipt_sha256")
                if not self._receipt_known(receipt_sha256):
                    raise RoutedContractError(
                        gate="decision.receipt_is_one_the_host_minted",
                        invariant=(
                            "a disposition cites only receipts this "
                            "session or a recorded run minted."
                        ),
                        diagnosis=(
                            f"{receipt_sha256[:8]} is "
                            + (
                                self._digest_names(receipt_sha256)
                                or "no digest this host minted"
                            )
                            + "."
                        ),
                        route=(
                            "cite the receipt_sha256 a tool returned, or "
                            "one inspect_run lists."
                        ),
                    )
            verified.append(
                {
                    "route": route,
                    "disposition": disposition,
                    "reason": reason,
                    "receipt_sha256s": receipts,
                }
            )
        return tuple(verified)

    #: Every registry the host keys by a receipt digest it minted. A
    #: registry opts in here **by name**: reflection over ``__dict__``
    #: would make a private cache keyed by a digest into citable
    #: evidence, which is not the same question. Authenticity -- did
    #: this host mint this digest -- is one predicate and lives here.
    #: Whether a receipt is *admissible* for a particular role is a
    #: second, narrower question each caller asks for itself.
    #:
    #: Five hand-written copies of a shorter version of this tuple had
    #: drifted apart (the decision gate and ``_receipt_known`` carried
    #: five names, the uncertainty resolver four, the post-hoc detector
    #: a different four), and the decision gate refused five digests
    #: this host had minted and printed -- a completion receipt, a
    #: stationary-point characterisation, an anomaly observation --
    #: while telling the session they were "no receipt of this
    #: session". The model wrote that falsehood into its permanent
    #: record. One authority, derived readers (SUFFICIENCY-5,
    #: 2026-09-10).
    _RECEIPT_REGISTRY_NAMES = (
        "analysis_claim_records",
        "analysis_completion_receipts",
        "anomaly_observations",
        "capabilities",
        "command_inspections",
        "environments",
        "functional_resolutions",
        "preflights",
        "project_renders",
        "project_validations",
        "quantity_expression_receipts",
        "quantity_extraction_bindings",
        "quantity_extraction_selectors",
        "quantity_extractions",
        "reached_geometries",
        "result_inspections",
        "safe_previews",
        "scientific_decisions",
        "scientific_validation_receipts",
        "stationary_point_characterisations",
        "substitutions",
        "thermochemistry_receipts",
        "validators",
    )

    #: Registries keyed by something that is *not* a minted receipt, and
    #: what they are keyed by. A digest found here is answered with what
    #: it actually names rather than with "no receipt": a reader who
    #: cites an artifact digest has made a different mistake from one who
    #: invented a digest, and the old message could not tell them apart.
    _NON_RECEIPT_REGISTRY_KEYS = {
        "atom_appends": "artifact",
        "consulted_skill_records": "document",
        "consulted_skills": "document",
        "database_extractions": "artifact",
        "engine_bindings": "binding",
        "geometry_edits": "artifact",
        "invocations": "invocation",
        "materialized_workflows": "plan",
        "mode_displacements": "artifact",
        "molecular_compositions": "artifact",
        "molecular_derivations": "artifact",
        "program_bindings": "binding",
        "project_documents": "document",
        "pubchem_geometries": "artifact",
        "quantity_expression_requests": "receipt",
        "scientific_identities": "binding",
        "scientific_toolchain_plans": "plan",
        "scientific_workflow_plans": "plan",
        "symmetry_breaks": "artifact",
        "workflow_drafts": "plan",
    }

    #: Every host-owned way a geometry can descend from another, and how
    #: to reach its predecessor. Two walks hand-listed subsets of this
    #: and disagreed: the review chain followed edits, appends,
    #: displacements and symmetry breaks and anchored on derivations and
    #: database extractions, while the original-bound-geometry walk
    #: followed three of the four and anchored on neither -- so a
    #: composed or fetched molecule reached the human page with no
    #: origin hop at all, and the panel's own docstring says the hop
    #: that decides what the molecule IS can sit at the root.
    #:
    #: A hop names the field carrying its predecessor's digest, because
    #: they genuinely differ: an edit, an append and a symmetry break
    #: come from a parent *geometry*, while a mode displacement comes
    #: from a *result*. Getting that wrong is how the two walks drifted.
    _GEOMETRY_HOPS: tuple[tuple[str, str, str], ...] = (
        ("geometry_edit", "geometry_edits", "parent_sha256"),
        ("atom_append", "atom_appends", "parent_sha256"),
        ("mode_displacement", "mode_displacements", "result_sha256"),
        ("symmetry_break", "symmetry_breaks", "parent_sha256"),
    )

    #: Origins anchor a chain: they have no predecessor geometry in this
    #: workspace, so the walk ends on them and the review names them.
    #: A composition is deliberately here rather than among the hops --
    #: it merges two parents and a linear walk cannot follow it, so it
    #: terminates the chain and reports both lineages instead of
    #: silently picking one.
    _GEOMETRY_ORIGINS: tuple[tuple[str, str], ...] = (
        ("derivation", "molecular_derivations"),
        ("database_extraction", "database_extractions"),
        ("pubchem_geometry", "pubchem_geometries"),
        ("composition", "molecular_compositions"),
    )

    def _geometry_provenance(
        self, sha256: str
    ) -> tuple[tuple[dict[str, Any], ...], tuple[str, Any] | None]:
        """The hops a geometry descends through, and the origin anchoring it.

        One reader for both walks. Returns the hops nearest-first and the
        origin, if this workspace holds one.
        """

        hops: list[dict[str, Any]] = []
        cursor = str(sha256 or "")
        seen: set[str] = set()
        while cursor and cursor not in seen and len(seen) < 64:
            seen.add(cursor)
            for kind, registry_name, predecessor in self._GEOMETRY_HOPS:
                registry = getattr(self, registry_name, None)
                receipt = (
                    registry.get(cursor)
                    if isinstance(registry, dict)
                    else None
                )
                if receipt is None:
                    continue
                hops.append({"kind": kind, **canonical_data(receipt)})
                cursor = str(getattr(receipt, predecessor, "") or "")
                break
            else:
                break
        origin: tuple[str, Any] | None = None
        for kind, registry_name in self._GEOMETRY_ORIGINS:
            registry = getattr(self, registry_name, None)
            if isinstance(registry, dict) and cursor in registry:
                origin = (kind, registry[cursor])
                break
        return tuple(hops), origin

    def _receipt_registries(self) -> dict[str, Any]:
        """Kind to registry, for every registry keyed by a receipt."""

        found: dict[str, Any] = {}
        for name in self._RECEIPT_REGISTRY_NAMES:
            registry = getattr(self, name, None)
            if isinstance(registry, dict):
                found[name] = registry
        return found

    def _resolve_receipt(self, receipt_sha256: str) -> tuple[str, Any] | None:
        """The registry a minted receipt lives in, and the receipt.

        A registry keyed by something other than the receipt still holds
        receipts: ``pubchem_geometries`` is keyed by the artifact digest,
        so the receipt digest the host minted, emitted on its own event
        and handed to the session was indexed nowhere -- and the gate
        that refused it said, truthfully of the registry and falsely of
        the ledger, that no receipt of that name existed. Whether a
        digest was minted here cannot depend on which key its family
        happens to use, so the value's own ``receipt_sha256`` answers
        too (SUFFICIENCY-5, 2026-09-10).
        """

        for name, registry in self._receipt_registries().items():
            if receipt_sha256 in registry:
                return name, registry[receipt_sha256]
        for name in self._NON_RECEIPT_REGISTRY_KEYS:
            registry = getattr(self, name, None)
            if not isinstance(registry, dict):
                continue
            for held in registry.values():
                if (
                    str(getattr(held, "receipt_sha256", "")) == receipt_sha256
                    and receipt_sha256
                ):
                    return name, held
        return None

    def _digest_names(self, receipt_sha256: str) -> str:
        """What a digest names here, for a refusal that tells the truth.

        A gate whose message says a digest was never minted, about a
        digest this host minted and returned, teaches the model a false
        fact about the host's own ledger -- and a live session wrote
        exactly that into its recorded decision. So the refusal names
        what it found instead.
        """

        resolved = self._resolve_receipt(receipt_sha256)
        if resolved is not None:
            return f"a {resolved[0]} receipt"
        for name, keyed_by in self._NON_RECEIPT_REGISTRY_KEYS.items():
            registry = getattr(self, name, None)
            if isinstance(registry, dict) and receipt_sha256 in registry:
                # The key, not a receipt: an artifact digest is a real
                # handle for a different question, and a reader who
                # cited one has made a different mistake from a reader
                # who invented a digest.
                return f"a {keyed_by} digest ({name}), not a receipt"
        if self._recorded_run_receipt(receipt_sha256):
            return "a receipt a recorded run of this workspace minted"
        return ""

    def _receipt_known(self, receipt_sha256: str) -> bool:
        """A receipt this session or a recorded run of the workspace minted."""

        if self._resolve_receipt(receipt_sha256) is not None:
            return True
        return self._recorded_run_receipt(receipt_sha256)

    def _refuse_plan_beyond_engine_budget(self, plan: Any) -> None:
        """Refuse, at plan time, more executable nodes than calls remain.

        The review builder refuses this at session end, which is the
        right contract at the wrong moment: a woken session planned
        twelve engine nodes against the five calls its wake context said
        remained, the frontier called the plan approvable, and the goal
        returned to the human (live, 2026-09-02). The plan alone proves
        the count, so the plan is where it is refused, against the
        smaller of the envelope and what the goal has left.
        """

        self._refuse_unfounded_excursions(plan)
        limits = []
        excursion_limits = []
        if self.bounded_execution_envelope is not None:
            limits.append(
                int(self.bounded_execution_envelope.max_engine_calls)
            )
            excursion_limits.append(
                int(self.bounded_execution_envelope.max_excursion_calls)
            )
        if self.engine_calls_remaining is not None:
            limits.append(int(self.engine_calls_remaining))
        if self.excursion_calls_remaining is not None:
            excursion_limits.append(int(self.excursion_calls_remaining))
        if not limits:
            return
        limit = min(limits)
        non_executable = self._release_non_executable_node_ids(plan)
        executable = tuple(
            node.node_id
            for node in plan.nodes
            if node.node_id not in non_executable and not node.excursion
        )
        excursions = tuple(
            node.node_id
            for node in plan.nodes
            if node.node_id not in non_executable and node.excursion
        )
        if len(executable) > limit:
            raise ContractError(
                "scientific workflow exceeds the engine-call budget: "
                f"{len(executable)} executable nodes for {limit} remaining "
                "calls; plan within the budget, reading registered results "
                "where they already stand instead of re-running them"
            )
        excursion_limit = min(excursion_limits) if excursion_limits else 0
        if len(excursions) > excursion_limit:
            raise ContractError(
                "scientific workflow exceeds the excursion grant: "
                f"{len(excursions)} excursion nodes for {excursion_limit} "
                "remaining excursion calls; an excursion is charged to its "
                "own displayed line and never to the engine-call budget"
            )

    def _refuse_unfounded_excursions(self, plan: Any) -> None:
        """An excursion cites a recorded anomaly and feeds no deliverable.

        The grant pays for investigation, not delivery: a tagged node
        whose output an untagged node consumes would buy the asked
        observable with the free line, and a tag citing no receipt the
        host minted is a model-authored anomaly.
        """

        tagged = {node.node_id: node.excursion for node in plan.nodes}
        tagged = {k: v for k, v in tagged.items() if v}
        if not tagged:
            return
        known = set(self.anomaly_observations) | {
            str(item.get("receipt_sha256") or "")
            for item in self.prior_anomaly_observations
        }
        for node_id, digest in sorted(tagged.items()):
            if digest not in known:
                raise ContractError(
                    f"node {node_id!r} is tagged as an excursion citing "
                    f"anomaly {digest[:8]}, which the host never recorded; "
                    "an excursion investigates an anomaly observation the "
                    "wake context or a validation receipt named"
                )
        for edge in getattr(plan, "edges", ()):
            source = str(getattr(edge, "source_node_id", "") or "")
            target = str(getattr(edge, "target_node_id", "") or "")
            if source in tagged and target and target not in tagged:
                raise ContractError(
                    f"excursion node {source!r} feeds {target!r}, which is "
                    "not an excursion; the grant may investigate an anomaly "
                    "and may never produce the asked observable"
                )

    def _refuse_excursion_producing_required_output(
        self, draft: Any, toolchain: Any
    ) -> None:
        """The grant investigates; it never delivers the asked observable.

        The edge rule catches a tagged node feeding another calculation
        node, and a probe showed what it misses: a tagged LEAF node,
        whose result the analysis chain reads and claims, buys the asked
        observable with the free line (2026-09-03). Reading an excursion
        is its whole point -- a replication receipt needs its numbers --
        so the refusal is narrower than "no analysis may read it": no
        required output may descend from a tagged node.
        """

        tagged = {
            node.node_id
            for node in getattr(draft, "nodes", ())
            if getattr(node, "excursion", "")
        }
        if not tagged:
            return
        analysis_nodes = {
            node.node_id: node
            for node in getattr(toolchain, "analysis_nodes", ())
        }
        produced_by = {
            output.output_id: node.node_id
            for node in analysis_nodes.values()
            for output in getattr(node, "outputs", ())
        }

        def _sources(node_id: str, seen: frozenset[str]) -> set[str]:
            """Calculation nodes an analysis node's outputs stand on."""

            node = analysis_nodes.get(node_id)
            if node is None or node_id in seen:
                return set()
            seen = seen | {node_id}
            found: set[str] = set()
            for item in getattr(node, "inputs", ()):
                producer = str(getattr(item, "producer_node_id", "") or "")
                if not producer:
                    continue
                if getattr(item, "source_kind", "") == "analysis_output":
                    found |= _sources(producer, seen)
                else:
                    found.add(producer)
            return found

        for output_id in getattr(toolchain, "required_output_ids", ()):
            owner = produced_by.get(str(output_id))
            if owner is None:
                continue
            offenders = sorted(_sources(owner, frozenset()) & tagged)
            if offenders:
                raise ContractError(
                    f"required output {output_id!r} descends from excursion "
                    f"node(s) {', '.join(offenders)}; the grant investigates "
                    "an anomaly and never produces the asked observable -- "
                    "produce it with an ordinary node and keep the tagged "
                    "node's own reading separate"
                )

    def _refuse_occupied_node_ids(self, node_ids: tuple[str, ...]) -> None:
        """Refuse, at plan time, a node id whose workspace holds outputs.

        Every node runs in ``<workspace>/nodes/<node_id>``, and the launch
        guard refuses a directory that already contains outputs so no
        evidence is ever overwritten. A revision that reused a failed
        node's id passed admission and met that guard only at launch,
        after the one-shot bundle was spent (live, 2026-09-02). The plan
        is where the id is chosen, so the plan is where the refusal
        names the route.
        """

        root = self.run_evidence_root
        if root is None:
            return
        occupied = tuple(
            node_id
            for node_id in node_ids
            if existing_node_branches(root, node_id)
        )
        if occupied:
            raise ContractError(
                "node id(s) "
                + ", ".join(repr(item) for item in occupied)
                + " name a node workspace that already holds outputs from "
                "an earlier run; a re-run takes a fresh node id, and the "
                "earlier directory stays as evidence"
            )

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
        self._refuse_occupied_node_ids(
            tuple(str(raw_node["node_id"]) for raw_node in values["nodes"])
        )
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
                internal_coordinates=(
                    canonical_data(raw_node.get("internal_coordinates"))
                    if raw_node.get("internal_coordinates")
                    else None
                ),
                excursion=str(raw_node.get("excursion") or ""),
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
            self._refuse_plan_beyond_engine_budget(scientific_plan)
            self.scientific_workflow_plans[scientific_plan.plan_sha256] = (
                scientific_plan
            )
            # The plan the session is now standing on, named rather than
            # inferred from insertion order: an amendment that restores
            # an earlier plan re-inserts its key in place, so ordering
            # points at the abandoned one.
            self.current_scientific_plan_sha256 = scientific_plan.plan_sha256
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
                                uncertainty_producer_node_id=str(
                                    item.get(
                                        "uncertainty_producer_node_id", ""
                                    )
                                ),
                                uncertainty_producer_output_id=str(
                                    item.get(
                                        "uncertainty_producer_output_id", ""
                                    )
                                ),
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
                    use_weighted_mass=raw_node.get("use_weighted_mass", False),
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
        self._refuse_excursion_producing_required_output(draft, plan)
        self.scientific_toolchain_plans[plan.plan_sha256] = plan
        self._scientific_toolchain_command_results[plan.plan_sha256] = (
            command_result
        )
        # When the calculation side could not project a scientific plan, the
        # resolver has already written down exactly why -- unbound identity,
        # unresolved inputs -- into the findings. Failing later on the
        # toolchain-binding invariant discarded that diagnosis and handed the
        # model a sentence about the host's own bookkeeping ("scientific
        # toolchain lacks its task-bound scientific plan"), which three
        # identical retries in one live session could not act on. Refuse here
        # instead, with the findings and the actions they call for.
        draft_for_findings = command_result.get("workflow_draft")
        if (
            getattr(draft_for_findings, "nodes", ())
            and command_result.get("scientific_workflow_plan") is None
        ):
            findings = tuple(command_result.get("findings") or ())
            named = (
                "; ".join(
                    f"{item.get('node_id')}: {item.get('rule_id')}"
                    for item in findings
                )
                or "no finding was recorded"
            )
            # A live session that had blocked every stage (as instructed) and
            # bound two identities (as instructed) was still refused twice,
            # because the general sentence below was already satisfied
            # vacuously: nothing in its draft consumed any artifact, so
            # nothing could anchor a molecule. Say that case in its own words.
            anchors_nothing = not any(
                getattr(intent, "artifact_id", "")
                for node in draft_for_findings.nodes
                for intent in getattr(node, "inputs", ())
            )
            if anchors_nothing:
                raise ContractError(
                    "the calculation nodes could not be bound into a "
                    f"scientific plan ({named}). No input anywhere in this "
                    "workflow names a workspace artifact, so the workflow "
                    "anchors no molecule: blocked stages document intent but "
                    "cannot anchor identity, and producer-fed inputs only "
                    "chain back to initial ones. At least one plannable "
                    "initial node must consume an identity-bound workspace "
                    "geometry artifact (bind_scientific_identity, then name "
                    "its artifact_id in that node's input)."
                )
            raise ContractError(
                "the calculation nodes could not be bound into a scientific "
                f"plan ({named}). Every initial node input must name a "
                "workspace geometry artifact whose scientific identity is "
                "bound (bind_scientific_identity), or the node must be "
                "declared blocked_unsupported with its reason; an input with "
                "neither an artifact_id nor a producer cannot anchor "
                "molecular identity."
            )
        self._bind_program_toolchain(plan, command_result)
        host_executed = None
        if (
            self.execute_analysis_only_plans
            and not draft.nodes
            and plan.analysis_nodes
        ):
            host_executed = self._execute_analysis_only_plan(plan)
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
            completed_analysis_node_ids=tuple(
                item["node_id"]
                for item in (host_executed or {}).get("executed_nodes", ())
                if item.get("state") == "executed"
            ),
            non_executable_calculation_node_ids=self._non_executable_reasons(
                command_result.get("scientific_workflow_plan")
            ),
        )
        result = {
            "scientific_toolchain_plan": plan,
            "calculation_plan": command_result,
            "workflow_frontier": frontier,
        }
        if host_executed is not None:
            result["host_executed_analysis"] = host_executed
        return result

    def _execute_analysis_only_plan(
        self, plan: ScientificToolchainPlanV1
    ) -> dict[str, Any]:
        """Walk an analysis-only plan now, on this host, and say so.

        Under a goal wake the plan's only inputs are results the goal's
        decision already covers; nothing is previewable and nothing is
        launched, so there is no review to wait for. A woken session that
        planned exactly this and waited for one lost six computed numbers
        (NOVEL-2 po2, 2026-09-04). The receipts land in this session's
        stream, the completion gate runs, and the result names what was
        delivered and what was not.
        """

        from chemsmart.agent.executor import execute_analysis_only_toolchain

        run_directory = (
            self.analysis_only_run_directory
            or Path(self.event_store.path).parent
        )
        workspace = self.analysis_only_workspace or (
            Path(self.run_evidence_root)
            if getattr(self, "run_evidence_root", None)
            else run_directory
        )
        record = execute_analysis_only_toolchain(
            host=self,
            toolchain=plan,
            run_directory=run_directory,
            task_spec_sha256=self._resolve_task_spec_reference(
                {}, "task_spec_sha256"
            ),
            workspace=workspace,
        )
        record["meaning"] = (
            "the host executed this analysis-only plan now, under the "
            "goal's standing decision, with no engine call; its claims "
            "and completion are in this session's receipts, so record "
            "the scientific decision over them rather than waiting for "
            "a review"
        )
        return record

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

    def _amend_scientific_workflow(self, turn_id: str, values: dict) -> Any:
        """Repair how one part of a planned DAG is expressed.

        A DAG arrives in one payload of roughly thirteen kilobytes over nine
        nodes, and a single mistyped identifier used to cost a resubmission of
        the whole graph.  Most rejections are about *expression* -- an
        identifier's case, a missing unit, a declared kind the operation does
        not derive, a selector the result does not resolve -- and repairing one
        node should not mean re-authoring the other eight.

        This is deliberately the same operation as rebinding a project, because
        it is the same act: clone the latest plan, change how one part is
        stated, append the result under its own digest.  Nothing is mutated in
        place, so the previous revision stays addressable, and every downstream
        digest still changes, so an amended plan goes through the normal
        materialise, review and approve cycle exactly as a fresh one would.

        What it will not do is change the science.  Molecular identity, state,
        program, job type, the producing node an input reads from, an analysis
        kind, the thermochemical conditions, and validation thresholds are all
        refused: those redefine the question rather than fix how it was
        written, and they belong in a new plan that a human reviews.
        """

        project_replacements = tuple(values.get("project_replacements") or ())
        analysis_repairs = tuple(values.get("analysis_repairs") or ())
        support_repairs = tuple(values.get("support_repairs") or ())
        if (
            not project_replacements
            and not analysis_repairs
            and not support_repairs
        ):
            raise ContractError(
                "an amendment must supply project_replacements, "
                "analysis_repairs, or support_repairs; it changes something "
                "or it is not an amendment"
            )
        return self._rebind_scientific_workflow_projects(
            turn_id,
            {
                "workflow_id": values["workflow_id"],
                "replacements": project_replacements,
                "analysis_repairs": analysis_repairs,
                "support_repairs": support_repairs,
            },
        )

    def _repaired_analysis_nodes(
        self, nodes: tuple, repairs: tuple
    ) -> tuple[tuple, dict[str, str]]:
        """Apply expression-level repairs to named analysis nodes.

        Each repair addresses an element that already exists -- an output, a
        selector, or an input -- and replaces only the field named.  Addressing
        something absent is refused with the inventory that node does carry,
        because a caller guessing at a name is the failure this exists to end.

        Returns the revised nodes and any output renames, so the caller can
        carry them into the required-output set.  Renaming an output otherwise
        orphans a requirement that names the old id, and asking the caller to
        restate the requirement is the kind of bookkeeping this whole change
        exists to stop demanding.
        """

        renames: dict[str, str] = {}
        if not repairs:
            return nodes, renames
        nodes_by_id = {node.node_id: node for node in nodes}
        unknown = sorted(
            {str(item["node_id"]) for item in repairs}.difference(nodes_by_id)
        )
        if unknown:
            raise ContractError(
                f"analysis repairs reference unknown nodes {unknown}; this "
                f"workflow has {sorted(nodes_by_id)}"
            )

        revised = dict(nodes_by_id)
        for repair in repairs:
            node = revised[str(repair["node_id"])]
            outputs = list(node.outputs)
            selectors = list(node.selectors)
            inputs = list(node.inputs)
            local_renames: dict[str, str] = {}

            for item in repair.get("outputs") or ():
                target = str(item["output_id"])
                index = next(
                    (
                        position
                        for position, output in enumerate(outputs)
                        if output.output_id == target
                    ),
                    None,
                )
                if index is None:
                    raise ContractError(
                        f"analysis node {node.node_id!r} has no output "
                        f"{target!r}; it declares "
                        f"{[output.output_id for output in outputs]}"
                    )
                current = outputs[index]
                renamed = str(item.get("new_output_id") or target)
                if renamed != target:
                    # Keyed by producing node as well as name: two nodes may
                    # each declare an output called the same thing, and a
                    # rename of one must not follow the other's consumers.
                    renames[(node.node_id, target)] = renamed
                    local_renames[target] = renamed
                outputs[index] = AnalysisOutputIntentV1(
                    output_id=renamed,
                    quantity_kind=str(
                        item.get("quantity_kind") or current.quantity_kind
                    ),
                    unit=str(item.get("unit") or current.unit),
                )

            for item in repair.get("selectors") or ():
                target = str(item["quantity_id"])
                index = next(
                    (
                        position
                        for position, selector in enumerate(selectors)
                        if selector.quantity_id == target
                    ),
                    None,
                )
                if index is None:
                    raise ContractError(
                        f"analysis node {node.node_id!r} has no selector for "
                        f"{target!r}; it declares "
                        f"{[item.quantity_id for item in selectors]}"
                    )
                selectors[index] = AnalysisSelectorIntentV1(
                    quantity_id=target,
                    selector=str(item["selector"]),
                )

            for item in repair.get("inputs") or ():
                target = str(item["input_id"])
                index = next(
                    (
                        position
                        for position, value in enumerate(inputs)
                        if getattr(value, "input_id", "") == target
                    ),
                    None,
                )
                if index is None:
                    raise ContractError(
                        f"analysis node {node.node_id!r} has no input "
                        f"{target!r}; it declares "
                        f"{[getattr(value, 'input_id', '') for value in inputs]}"
                    )
                current = inputs[index]
                if not isinstance(current, AnalysisInputIntentV1):
                    raise ContractError(
                        f"input {target!r} reads a registered result, whose "
                        "artifact is scientific intent rather than an "
                        "expression; plan that as a new workflow"
                    )
                # Only which named output of the *same* producer is read.
                # Repointing to another producer substitutes a different
                # calculation's result, which is a different question.
                inputs[index] = replace(
                    current,
                    producer_output_id=str(item["producer_output_id"]),
                )

            expression_nodes = node.expression_nodes
            expression_output_node_ids = node.expression_output_node_ids
            if node.analysis_kind == "quantity_expression" and local_renames:
                # The execution receipt keys produced quantities by
                # expression node id and the node contract requires the
                # declared outputs to name the exported expression nodes,
                # so an output rename follows into the expression DAG
                # rather than orphaning it.
                expression_output_node_ids = tuple(
                    local_renames.get(item, item)
                    for item in node.expression_output_node_ids
                )
                followed_nodes = []
                for item in node.expression_nodes:
                    entry = dict(item)
                    identifier = str(entry.get("node_id", ""))
                    if identifier in local_renames:
                        entry["node_id"] = local_renames[identifier]
                    for key in ("operand_ids", "input_ids"):
                        if key in entry:
                            entry[key] = tuple(
                                local_renames.get(str(value), str(value))
                                for value in entry[key]
                            )
                    followed_nodes.append(entry)
                expression_nodes = tuple(followed_nodes)
            revised[node.node_id] = replace(
                node,
                outputs=tuple(
                    sorted(outputs, key=lambda value: value.output_id)
                ),
                selectors=tuple(
                    sorted(selectors, key=lambda value: value.quantity_id)
                ),
                inputs=tuple(
                    sorted(
                        inputs,
                        key=lambda value: getattr(value, "input_id", ""),
                    )
                ),
                expression_nodes=expression_nodes,
                expression_output_node_ids=expression_output_node_ids,
            )

        # A renamed output is still read by whoever consumed it under the old
        # name, and leaving those references behind orphans them: the plan
        # builder refuses an input naming an output its producer no longer
        # declares.  Following the rename through every consumer is the same
        # bookkeeping this whole operation exists to stop demanding, so the
        # host does it rather than asking for a second repair per consumer.
        if renames:
            for node_id, node in list(revised.items()):
                followed = []
                changed = False
                for value in node.inputs:
                    key = (
                        getattr(value, "producer_node_id", ""),
                        getattr(value, "producer_output_id", ""),
                    )
                    if isinstance(value, AnalysisInputIntentV1) and (
                        key in renames
                    ):
                        followed.append(
                            replace(value, producer_output_id=renames[key])
                        )
                        changed = True
                    else:
                        followed.append(value)
                if changed:
                    revised[node_id] = replace(node, inputs=tuple(followed))
        return tuple(revised[node.node_id] for node in nodes), renames

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
            # A live session amended straight after a rejected plan call and
            # was told only that the ID was unknown -- true, but useless: a
            # rejected plan_scientific_workflow records nothing, so there was
            # never anything to amend. Say what is known and what the repair
            # actually is.
            known = sorted(
                {
                    plan.workflow_id
                    for plan in self.scientific_toolchain_plans.values()
                }
            )
            raise RoutedContractError(
                gate="plan.amend_needs_a_recorded_plan",
                invariant=(
                    "a rejected plan_scientific_workflow call records "
                    "nothing, so only a recorded plan can be amended."
                ),
                diagnosis=(
                    f"unknown scientific workflow ID {workflow_id!r}; "
                    "recorded workflow IDs: "
                    f"{', '.join(known) if known else 'none'}."
                ),
                route=(
                    "submit the corrected plan as a fresh "
                    "plan_scientific_workflow call."
                ),
            )
        current_plan = candidates[-1]
        current_result = self._scientific_toolchain_command_results[
            current_plan.plan_sha256
        ]
        current_draft = current_result["workflow_draft"]
        # A frozen approval carries the *v2* workflow digest, not this v1
        # toolchain digest.  Comparing the two compares hashes of disjoint
        # bodies whose schema_version alone differs, so the guard could never
        # fire; resolve the v2 plan for this workflow the way every other
        # frozen-approval check in this file does.
        approved_plan = self.scientific_plans.get(current_plan.workflow_id)
        if (
            self.frozen_workflow_approval is not None
            and approved_plan is not None
            and self.frozen_workflow_approval.plan_sha256
            == approved_plan.plan_sha256
        ):
            raise ContractError(
                "this workflow already carries a frozen execution approval; "
                "an amended plan is a new workflow and needs its own review"
            )
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
        # A support repair converts evidence learned mid-session into
        # declared intent. A real session planned the paper's third
        # functional, watched the program validator refuse it -- the
        # functional is not implemented in that program -- and correctly
        # said the stage should remain as non-executable scientific intent;
        # nothing could then make it so, and five nodes blocked an approval
        # they were never going to be part of.
        #
        # The repair is one-directional by design. Declaring a stage
        # blocked_unsupported only narrows what runs: it is displayed with
        # the workflow, excluded from the approval, and never launched.
        # The reverse direction would widen the executable partition and
        # therefore belongs to a fresh plan a human reviews, not an amend.
        support_repairs = tuple(values.get("support_repairs") or ())
        for repair in support_repairs:
            node_id = str(repair["node_id"])
            if node_id not in annotations:
                raise ContractError(
                    f"support repair references unknown node {node_id!r}; "
                    f"this workflow has {sorted(annotations)}"
                )
            reason = str(repair.get("blocked_reason") or "").strip()
            if not reason:
                raise ContractError(
                    "declaring a stage non-executable requires its reason"
                )
            if annotations[node_id]["support_state"] == "blocked_unsupported":
                continue
            annotations[node_id] = {
                **annotations[node_id],
                "support_state": "blocked_unsupported",
                "blocked_reason": reason,
            }
        # Recompiling the calculation side is only needed when a project
        # binding changed.  An analysis-only repair leaves every command,
        # binding and receipt on that side exactly as reviewed, so it is
        # carried by digest rather than rebuilt -- the mirror image of the way
        # this function has always carried the analysis side verbatim.
        command_result = (
            current_result
            if not replacements and not support_repairs
            else self._plan_command_workflow(
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
                            # Every per-node scientific fact has to be listed
                            # here by hand, and anything forgotten is dropped
                            # silently. A scan whose project role was amended
                            # came back out of this rebuild with no driven
                            # coordinate at all -- still typed `scan`, compiled
                            # without --coordinates, and therefore a plain
                            # optimisation wearing a scan's name. Changing a
                            # project must not change the chemistry.
                            **(
                                {
                                    "internal_coordinates": (
                                        node.internal_coordinates
                                    )
                                }
                                if node.internal_coordinates is not None
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
        )
        draft = command_result["workflow_draft"]
        observables = dict(current_plan.calculation_observables)
        revised_analysis, renames = self._repaired_analysis_nodes(
            current_plan.analysis_nodes,
            tuple(values.get("analysis_repairs") or ()),
        )
        revised_plan = build_scientific_toolchain_plan(
            plan_id=current_plan.plan_id,
            workflow_id=current_plan.workflow_id,
            command_workflow_draft_sha256=draft.draft_sha256,
            calculation_nodes=draft.nodes,
            calculation_observables=observables,
            analysis_nodes=revised_analysis,
            required_output_ids=tuple(
                _renamed_output_id(renames, output_id)
                for output_id in current_plan.required_output_ids
            ),
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
            # Name what exists: a woken session guessed the previous
            # cycle's workflow id twice in one goal, and the bare refusal
            # left it guessing (live, 2026-09-02).
            known = sorted(self._latest_program_workflows)
            if known:
                raise ContractError(
                    f"unknown scientific workflow ID {workflow_id!r}; this "
                    "session holds workflow(s) "
                    + ", ".join(repr(item) for item in known)
                )
            # Under a goal wake the only imperative here used to be
            # "plan one first", and a woken session that held every
            # number it needed obeyed it, planned a chain no review
            # could carry, and lost six computed numbers (NOVEL-2 po2,
            # 2026-09-04). The route depends on where the session stands.
            route = (
                "an earlier cycle's results are read through inspect_run "
                "and its workflow through inspect_run's run outcome; to "
                "deliver from receipts already in hand, call "
                "extract_result_quantities, derive_thermochemistry, "
                "evaluate_quantity_expression and record_analysis_claims "
                "yourself, or plan_scientific_workflow with no "
                "calculation_nodes -- the host executes that plan when "
                "planned under the goal's standing decision"
                if self.execute_analysis_only_plans
                else "plan_scientific_workflow plans one; an earlier "
                "cycle's workflow is read through inspect_run"
            )
            raise RoutedContractError(
                gate="workflow.id_names_a_plan_this_session_holds",
                invariant=(
                    "a workflow id resolves only against a plan this "
                    "session recorded."
                ),
                diagnosis=(
                    f"{workflow_id!r} names no workflow in this session, "
                    "which holds none yet."
                ),
                route=route,
            )
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

    #: Selectors that constitute reading a geometry's structure rather
    #: than its labels: positions is the gateway every distance, angle,
    #: and dihedral expression consumes, and connectivity is the
    #: perceived bond graph. The frontier's unread listing keys on these
    #: and nothing else.
    _STRUCTURAL_READ_SELECTORS = frozenset({"positions", "connectivity"})

    def _artifacts_without_structural_read(self) -> tuple[str, ...]:
        """Name structurally readable artifacts nothing has yet read.

        Host bookkeeping over state already held: every registered
        artifact whose kind a structural read could serve -- geometry
        files and each registered reader's result kind, so an identity
        audit over archived results is covered exactly like one over
        supplied geometries -- minus those with an extraction receipt
        whose request included a structural selector. The listing names
        artifacts, never what a reading would say -- measuring stays the
        session's act, and an id here is a question, not a verdict.
        Scope is this session's own reads: a rehydrated host starts with
        empty selector records because receipts do not carry selectors.
        """

        from chemsmart.agent.postprocessing import (
            typed_result_artifact_kind,
        )
        from chemsmart.analysis.result_readers import (
            registered_reader_programs,
        )

        readable_kinds = {"geometry_xyz"}
        for program in registered_reader_programs():
            try:
                readable_kinds.add(typed_result_artifact_kind(program))
            except Exception:
                continue
        read: set[str] = set()
        for sha, selectors in self.quantity_extraction_selectors.items():
            receipt = self.quantity_extractions.get(sha)
            if receipt is None:
                continue
            if any(
                selector in self._STRUCTURAL_READ_SELECTORS
                for selector in selectors
            ):
                read.add(receipt.artifact_id)
        return tuple(
            sorted(
                artifact_id
                for artifact_id, artifact in self.artifacts.items()
                if getattr(artifact, "kind", "") in readable_kinds
                and artifact_id not in read
            )
        )

    def _select_execution_wave(self, turn_id: str, values: dict) -> Any:
        """Judge the wave the Agent named against the host's frontier.

        The host owns readiness and owns concurrency; the Agent owns
        which ready calculations belong in one wave. So this asks
        ``_workflow_context`` -- the same projection
        ``inspect_workflow_frontier`` already serves -- rather than
        deriving a second frontier, and it never chooses, reorders or
        refuses. The order the Agent gave is the order the array elements
        take, because which calculation is element 0 is a scientific
        decision the host has no basis to overrule.

        A wave the host cannot dispatch comes back as a per-member
        verdict, not an exception: an exception is what teaches a session
        to carry workarounds for something the host should simply have
        reported (owner ruling, 2026-09-16).
        """

        del turn_id
        from chemsmart.agent.cohort import validate_wave

        resolved = self._resolve_program_workflow(values["workflow_id"])
        draft = resolved.draft
        scientific = getattr(resolved, "scientific_plan", None)
        context = self._workflow_context(
            draft,
            scientific_plan_sha256=(
                getattr(scientific, "plan_sha256", "") if scientific else ""
            ),
        )
        edges = tuple(
            (
                str(getattr(item, "producer_node_id", "") or ""),
                str(node.node_id),
            )
            for node in getattr(draft, "nodes", ()) or ()
            for item in getattr(node, "inputs", ()) or ()
            if getattr(item, "producer_node_id", "")
        )
        proposed = tuple(str(item) for item in values.get("node_ids") or ())
        verdict = validate_wave(
            proposed=proposed,
            ready=tuple(context.ready_node_ids),
            edges=edges,
            planned=tuple(
                str(node.node_id)
                for node in getattr(draft, "nodes", ()) or ()
            ),
        )
        dispatchable = bool(verdict.rows) and all(
            row.status == "ready" for row in verdict.rows
        )
        # Where the dispatcher reads it. A wave that lives only in a
        # tool reply is a wave the array never hears about -- and a wave
        # left standing after the Agent has moved on is worse, because
        # the driver would submit calculations it has already seen. The
        # last word is the selection, whether or not it is dispatchable.
        self.selected_execution_wave = (
            tuple(verdict.members) if dispatchable else ()
        )
        record = verdict.public_record()
        return {
            "status": "ready" if dispatchable else "not_dispatchable",
            "workflow_id": str(draft.workflow_id),
            "node_ids": list(verdict.members) if dispatchable else [],
            "members": record.get("rows", []),
            "summary": verdict.summary,
            "next_action": (
                "this wave is what will be submitted; every member runs "
                "and you are woken once, when all of them have ended"
                if dispatchable
                else "select again from the members the host reports "
                "ready, and choose the rest after reading this wave"
            ),
        }

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
        # The system prompt sends a session here for host-derived next
        # actions, and this was the one projection that never said whether the
        # workflow could actually be approved. `workflow_context` reports
        # dependency state -- a node stays `actionable` forever whether or not
        # it holds a green preview -- so a session could read the frontier,
        # see nothing outstanding, and stop while a node still blocked
        # approval. The readiness is already computed; it just was not offered
        # where the model was told to look.
        readiness = (
            self._approval_readiness(scientific_v2)
            if scientific_v2 is not None
            else None
        )
        plan = resolved.scientific_toolchain_plan
        if plan is None:
            result = {
                "scientific_workflow_plan": scientific_v2,
                "workflow_context": context,
                "artifacts_without_structural_read": (
                    self._artifacts_without_structural_read()
                ),
            }
            if readiness is not None:
                result["approval_readiness"] = readiness
            return result
        analysis_receipts = self._scientific_toolchain_analysis_receipts(
            plan,
            task_spec_sha256=draft.task_spec_id,
        )
        result = {
            "scientific_toolchain_plan": plan,
            "workflow_context": context,
            "artifacts_without_structural_read": (
                self._artifacts_without_structural_read()
            ),
            "workflow_frontier": project_scientific_toolchain_frontier(
                plan,
                actionable_calculation_node_ids=actionable,
                unresolved_calculation_node_ids=unresolved,
                completed_calculation_node_ids=context.completed_node_ids,
                completed_analysis_node_ids=analysis_receipts,
                non_executable_calculation_node_ids=(
                    self._non_executable_reasons(scientific_v2)
                ),
            ),
        }
        if readiness is not None:
            result["approval_readiness"] = readiness
        return result

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
        calculation_ids = set(plan.calculation_node_ids)

        def _calculation_dependency_satisfied(dependency: str) -> bool:
            receipt = self.execution_receipts.get(dependency)
            return receipt is not None and bool(
                getattr(receipt, "validated", False)
            )

        def _dependency_receipts(
            node: AnalysisNodeIntentV1,
        ) -> dict[str, set[str]]:
            # A dependency on a calculation node used to route through
            # ``matched`` like an analysis dependency, and calculation nodes
            # are never in ``matched`` -- so no producer-fed analysis node
            # could ever be reported complete, even after its producer had
            # executed and validated. Calculation dependencies are satisfied
            # by validated execution receipts; only analysis dependencies
            # carry receipt digests forward.
            receipts: dict[str, set[str]] = {}
            for dependency in node.dependencies:
                if dependency in calculation_ids:
                    if _calculation_dependency_satisfied(dependency):
                        receipts[dependency] = {"calculation-validated"}
                    else:
                        receipts[dependency] = set()
                else:
                    receipts[dependency] = set(matched.get(dependency, ()))
            return receipts

        def _producer_result_artifact_ids(
            node: AnalysisNodeIntentV1,
        ) -> set[str]:
            """Registered engine-result artifacts feeding this node."""

            artifact_ids: set[str] = set()
            for item in node.inputs:
                if isinstance(item, RegisteredResultInputIntentV1):
                    continue
                producer = item.producer_node_id
                if producer not in calculation_ids:
                    continue
                if not _calculation_dependency_satisfied(producer):
                    continue
                prefix = f"result.{producer}."
                artifact_ids.update(
                    artifact_id
                    for artifact_id in self.artifacts
                    if artifact_id.startswith(prefix)
                )
            return artifact_ids

        def _decision_evidence(
            decision: ScientificDecisionRecordV1,
        ) -> set[str]:
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
                if len(registered) == 1:
                    extraction_artifact_ids = {registered[0].artifact_id}
                elif not registered:
                    extraction_artifact_ids = _producer_result_artifact_ids(
                        node
                    )
                    if not extraction_artifact_ids:
                        continue
                else:
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
                    and receipt.artifact_id in extraction_artifact_ids
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
                        and receipt.artifact_id in extraction_artifact_ids
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
                        if (receipt := self.quantity_extractions.get(digest))
                        is not None
                    }
                    source_artifact_ids.update(
                        _producer_result_artifact_ids(node)
                    )
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
                        if (
                            tuple(quantities[quantity_kind].dimension)
                            != dimension
                        ):
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
                    if str(expression.get("node_id") or "") in required_outputs
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
                        and (producer := nodes.get(item.producer_node_id))
                        is not None
                        and producer.analysis_kind == "scientific_validation"
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
                                ).intersection(dependency_receipts)
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
                        sorted(
                            {digest for pair in candidates for digest in pair}
                        )
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
                    producer_stage = next(
                        (
                            candidate.stage
                            for candidate in scientific_v2.nodes
                            if candidate.node_id == item.producer_node_id
                        ),
                        "",
                    )
                    waiting_producers.append(
                        {
                            "binding_id": item.binding_id,
                            "producer_node_id": item.producer_node_id,
                            "producer_output_id": item.producer_output_id,
                            "producer_stage": producer_stage,
                            "deferrable_within_one_approval": (
                                producer_stage
                                in DEFERRABLE_GEOMETRY_PRODUCER_STAGES
                            ),
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
                **_undeferrable_producer_finding(waiting_producers),
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
            # A node names its project by artifact id, so a role that differs
            # from the promoted id only in spelling reads here as no project at
            # all -- and "render, promote, and validate this project role" is
            # then advice to redo work that is already done. Name what has been
            # promoted, the way needs_engine_selection above names the engines
            # it found, so the mismatch is visible rather than inferred.
            promoted = tuple(sorted(self.project_promotions))
            return {
                "status": "needs_project",
                "workflow_id": workflow_id,
                "node_id": node_id,
                "project_role": node.project_role,
                "promoted_project_roles": promoted,
                "next_action": (
                    "promote this project role, or amend the node to name one "
                    "of the promoted roles"
                    if promoted
                    else "render, promote, and validate this project role"
                ),
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
            job_option_values=_node_coordinates(node),
            live_schema=self.live_schema,
            server=(
                self.execution_server
                if self.surface.profile
                == "command_compiled_approved_execution"
                else self.preview_server
            ),
        )
        context = _CommandContext(
            internal_coordinates=getattr(node, "internal_coordinates", None),
            excursion=str(getattr(node, "excursion", "") or ""),
            proposal=proposal,
            capability=capability,
            program_binding=program_binding,
            engine_binding=engine_binding,
            project_artifact=project,
            project_validation=validation,
            input_artifact=input_artifact,
            root_artifact=self._root_artifact_for(
                proposal.node_id, input_artifact
            ),
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
            "observations": compile_time_observations(
                program=node.program,
                jobtype=node.jobtype,
                settings=validation.settings,
                atom_count=int(
                    _review_molecule_identity(input_artifact).get("atom_count")
                    or 0
                ),
                geometry=_geometry_for_observation(input_artifact),
            )
            # The program's own word on the bytes just written, in the
            # reply the model reads. It reached only the review page --
            # which under a standing approval nobody reads -- so a
            # 0.107 s uncharged probe that printed `RIJCOSX` by name
            # informed no decision, eight nodes launched over aborted
            # checks across two windows, and both cycles re-learned the
            # identical `engine_lines` from the dead runs. Verified:
            # po3-r18 cycle 2's model-visible transcript contains zero
            # occurrences of `RIJCOSX` and zero of `aborting the run`.
            + self._probe_observations_for(node.node_id),
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
                else "unresolved_future" if unresolved else "resolvable"
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
                    excursion=node.excursion,
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
        # An excursion was never a stage the task asked for; dropping it
        # next cycle is not a regression of the deliverable.
        dropped = sorted(
            node.node_id
            for node in previous.nodes
            if node.node_id not in current_ids
            and not getattr(node, "excursion", "")
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
        # Only the geometry edges are counted: admission keys each
        # producer edge by its consumer role, so an ORCA IRC or TS node
        # carrying a Hessian edge beside its geometry edge is one
        # candidate, not two. This predicate counted every data edge and
        # called po3's two IRC nodes blocking while the review resolved
        # and ran them (REACH-1, 2026-09-06) -- the frontier disagreeing
        # with the review in the other direction from ino3's. Whether
        # the auxiliary edge has a legal shape is the resolver's word,
        # which the frontier now asks before it calls a node deferred.
        geometry_counts: dict[str, int] = {}
        for edge in data_edges:
            if edge.artifact_class == "geometry_xyz":
                geometry_counts[edge.target_node_id] = (
                    geometry_counts.get(edge.target_node_id, 0) + 1
                )
        deferred = set()
        for edge in data_edges:
            producer = nodes.get(edge.source_node_id)
            target = nodes.get(edge.target_node_id)
            if (
                producer is None
                or target is None
                or geometry_counts.get(edge.target_node_id) != 1
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
        planned_ids = {node.node_id for node in getattr(plan, "nodes", ())}
        executable_ids = planned_ids - non_executable_ids
        # The review resolves every deferred node to one project, one
        # capability and one environment before it can be displayed. The
        # frontier used to admit a deferred node by its edge shape alone,
        # so it said approvable while the review then refused (REACH-1
        # ino3, 2026-09-06: two organs, opposite answers, and the refusal
        # never reached the session). The frontier now asks the review's
        # own resolver, without binding, and a node it would refuse is a
        # blocking node carrying that reason.
        review_target_ids = {
            edge.target_node_id
            for edge in getattr(plan, "edges", ())
            if getattr(edge, "edge_kind", "") == "data"
            and edge.target_node_id not in non_executable_ids
        }
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
            deferral_refusal = ""
            if deferred:
                try:
                    self._bounded_node_context(
                        plan=plan,
                        planned_node=node,
                        data_target_ids=review_target_ids,
                        bind=False,
                    )
                except ContractError as exc:
                    deferral_refusal = str(exc)
                    deferred = False
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
                    **(
                        {"deferral_refusal": deferral_refusal}
                        if deferral_refusal
                        else {}
                    ),
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
                node["node_id"] for node in nodes if node["non_executable"]
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
        # The executor rebuilds every approved node through this tool, so a
        # per-node coordinate that is not accepted here cannot be rebuilt at
        # all: two scans of one molecule over different ranges synthesised to
        # the same coordinate-free argv and were correctly refused as
        # differing from the reviewed operation.
        coordinates = values.get("internal_coordinates") or None
        invocation = compile_command(
            proposal,
            capability=capability,
            binding=binding,
            project=project,
            project_validation=validation,
            input_artifact=input_artifact,
            scientific_identity=identity,
            job_option_values=native_coordinate_options(
                values["program"], coordinates
            ),
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
            internal_coordinates=coordinates,
            excursion=str(values.get("excursion") or ""),
            proposal=proposal,
            capability=capability,
            program_binding=program_binding,
            engine_binding=binding,
            project_artifact=project,
            project_validation=validation,
            input_artifact=input_artifact,
            root_artifact=self._root_artifact_for(
                proposal.node_id, input_artifact
            ),
            scientific_identity=identity,
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
            retain_root=self.preview_retention_root,
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
                current_invocation, _context = self._plan_invocation_for_node(
                    plan=candidate_plan,
                    node_id=values["node_id"],
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
            plan = self.scientific_workflow_plans.get(invocation_plan_sha256)
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
                if (
                    observed.state != "pending"
                    or values["node_id"] not in ready
                ):
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
        self._probe_input_check(
            turn_id,
            node_id=values["node_id"],
            program=str(
                getattr(getattr(capability, "query", None), "program", "")
                or getattr(capability, "program", "")
                or ""
            ),
            safe_preview=safe_preview,
        )
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

    def _probe_observations_for(self, node_id: str) -> tuple[str, ...]:
        """The input-check probe's own lines for one node, if it ran.

        One accessor, so the compile reply, the frontier and the review
        all carry the same words instead of the review carrying them
        alone.
        """

        probe = getattr(self, "_input_check_by_node", {}).get(node_id)
        if probe is None:
            return ()
        from chemsmart.agent.input_check import probe_observation_lines

        return tuple(probe_observation_lines(probe))

    def _probe_input_check(
        self,
        turn_id: str,
        *,
        node_id: str,
        program: str,
        safe_preview: Any,
    ) -> Any:
        """Run the program's own input check on the previewed input.

        Two REACH-1 cycles died at ORCA's input check under green
        previews (2026-09-06). A green preview is ChemSmart's compile;
        the probe is ORCA's check, bounded and never charged (owner
        ruling R2). It runs only where the host-owned server profile
        names an ORCA executable, never inside a scheduler allocation
        -- the wake at a job's tail plans from one -- and only on a
        previewed input the preview retained by digest. Its word is an
        observation on the review beside the node, never a refusal.
        """

        from chemsmart.agent.input_check import (
            not_run_receipt,
            probe_orca_input_check,
        )
        from chemsmart.agent.preview import retained_preview_artifact

        executable = getattr(self, "input_check_executable", None)
        if program != "orca" or executable is None:
            return None
        if safe_preview is None or safe_preview.status != "previewed":
            return None
        inputs = tuple(
            artifact
            for artifact in safe_preview.artifacts
            if str(artifact.relative_path).endswith(".inp")
        )
        input_sha256 = inputs[0].sha256 if len(inputs) == 1 else ""
        cap = self.input_check_cap_seconds
        if any(key in os.environ for key in ("SLURM_JOB_ID", "PBS_JOBID")):
            receipt = not_run_receipt(
                node_id=node_id,
                program=program,
                input_sha256=input_sha256,
                reason="inside a scheduler allocation; the probe runs on "
                "the controller only",
                cap_seconds=cap,
            )
        elif len(inputs) != 1:
            receipt = not_run_receipt(
                node_id=node_id,
                program=program,
                input_sha256="",
                reason=f"the preview emitted {len(inputs)} .inp files",
                cap_seconds=cap,
            )
        else:
            retained = retained_preview_artifact(
                self.preview_retention_root, input_sha256
            )
            if retained is None:
                receipt = not_run_receipt(
                    node_id=node_id,
                    program=program,
                    input_sha256=input_sha256,
                    reason="the previewed input was not retained by digest",
                    cap_seconds=cap,
                )
            elif not executable.is_file():
                receipt = not_run_receipt(
                    node_id=node_id,
                    program=program,
                    input_sha256=input_sha256,
                    reason=f"no executable at {executable}",
                    cap_seconds=cap,
                )
            else:
                receipt = probe_orca_input_check(
                    node_id=node_id,
                    input_path=retained,
                    executable=executable,
                    env=self.input_check_env,
                    cap_seconds=cap,
                )
        self._input_check_by_node[node_id] = receipt
        self._emit(
            turn_id,
            EventKind.INPUT_CHECK_PROBED,
            receipt.receipt_sha256,
            node_id=node_id,
            program=program,
            status=receipt.status,
            reason=receipt.reason,
            input_sha256=receipt.input_sha256,
            engine_lines=receipt.engine_lines,
            wall_seconds=receipt.wall_seconds,
            cap_seconds=receipt.cap_seconds,
            charged=False,
        )
        return receipt

    def _current_scientific_plan(
        self, predicate=None
    ) -> ScientificWorkflowPlanV2 | None:
        """The plan the session is standing on, for every reader.

        One function, because three readers asked this question by
        taking the last value out of a digest-keyed dict and a dict
        re-insertion keeps its key in place: an amendment restoring an
        earlier plan left the frontier and the execution review
        answering about different plans, and the review refused what
        the frontier had just called approvable (po3-r17, 2026-09-11).
        Insertion order remains the fallback for a host restored from a
        record that predates the pointer.
        """

        # getattr, like the sibling registry a few readers below: a host
        # restored for preflight is built without running __init__, so
        # the pointer can be absent and insertion order is then all
        # there is.
        candidates = tuple(
            plan
            for plan in self.scientific_workflow_plans.values()
            if predicate is None or predicate(plan)
        )
        current = self.scientific_workflow_plans.get(
            getattr(self, "current_scientific_plan_sha256", "")
        )
        if current is not None and current in candidates:
            return current
        # The single place in this host where insertion order is ever
        # consulted, and only when the pointer names no admissible plan
        # -- a host restored for preflight is built without running
        # __init__. Every other reader calls this function, so a
        # digest-keyed container is reduced to one entry in exactly one
        # function and the lint can be absolute.
        return candidates[-1] if candidates else None

    def _materialize_scientific_workflow(
        self, *, turn_id: str, node_id: str
    ) -> MaterializedWorkflowV1 | None:
        """Ground the plan containing ``node_id`` from host receipts.

        The session's current plan first; only then the other
        registered plans, newest insertion first.
        """

        current = self._current_scientific_plan()
        candidates = ((current,) if current is not None else ()) + tuple(
            reversed(tuple(self.scientific_workflow_plans.values()))
        )
        plan = next(
            (
                candidate
                for candidate in candidates
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
        analysis_nodes = {node.node_id: node for node in plan.analysis_nodes}
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
        # The receipt states which required outputs the approved plan
        # itself declared blocked: a delivery with stated limitations
        # and a full delivery must not share one word.
        limitation_output_ids = tuple(
            sorted(
                {
                    output_id
                    for output_id in plan.required_output_ids
                    if output_producers[output_id] in missing_required
                }
            )
        )
        doubted_quantity_ids = self._claims_under_recorded_doubt(
            task_spec_sha256=draft.task_spec_id
        )
        failed_criterion_ids = self._claims_on_a_failed_criterion(
            plan, task_spec_sha256=draft.task_spec_id
        )
        status = "passed"
        findings: tuple[str, ...] = ()
        if failed_criterion_ids:
            # The number stays delivered and the word says what it
            # stands under: the plan's own acceptance criterion for the
            # result this claim descends from did not hold.
            status = "partial"
            findings = findings + tuple(
                f"analysis.claim_on_failed_criterion.{output_id}"
                for output_id in failed_criterion_ids
            )
        if doubted_quantity_ids:
            # A session that doubts a receipt and claims from it has said
            # both things in typed form; the completion word carries that
            # truth instead of certifying past it.
            status = "partial"
            findings = tuple(
                f"analysis.claim_under_recorded_doubt.{quantity_id}"
                for quantity_id in doubted_quantity_ids
            )
        return self._record_toolchain_completion(
            plan.plan_sha256,
            task_spec_sha256=draft.task_spec_id,
            source_receipt_sha256s=source_receipts,
            status=status,
            findings=findings,
            limitation_output_ids=limitation_output_ids,
        )

    def _claims_on_a_failed_criterion(
        self, plan: Any, *, task_spec_sha256: str
    ) -> tuple[str, ...]:
        """Claims that descend from a result whose own criterion failed.

        A plan may carry scientific_validation nodes -- an imaginary-mode
        count, a spin-manifold check, a scan ridge above the barrier it
        brackets -- and their verdicts reached nothing: REACH-1 po3
        planned exactly those rules and its delivery would have stood
        whatever they said (2026-09-06). The join needs no new field:
        both a claim and a criterion declare, in the plan the human
        approved, which producers they read, and a claim whose producers
        include one a failed criterion judged is named here.

        Never a refusal and never a silent drop. A failed criterion is a
        stated finding, the number stays delivered, and the reader is
        told which criterion it stands under.
        """

        nodes = {node.node_id: node for node in plan.analysis_nodes}
        if not nodes:
            return ()

        def producers(node_id: str, seen: frozenset[str] = frozenset()) -> set:
            node = nodes.get(node_id)
            if node is None or node_id in seen:
                return set()
            reached: set[str] = set()
            for item in node.inputs:
                producer = str(getattr(item, "producer_node_id", "") or "")
                artifact = str(getattr(item, "artifact_id", "") or "")
                if artifact:
                    reached.add(artifact)
                if not producer:
                    continue
                if producer in nodes:
                    reached |= producers(producer, seen | {node_id})
                else:
                    reached.add(producer)
            return reached

        failed_producers: set[str] = set()
        failed_by_node: dict[str, str] = {}
        for receipt in self.scientific_validation_receipts.values():
            if getattr(receipt, "all_rules_passed", True):
                continue
            node_id = str(getattr(receipt, "node_id", "") or "")
            if node_id not in nodes:
                continue
            reached = producers(node_id)
            failed_producers |= reached
            for producer in reached:
                failed_by_node[producer] = node_id
        if not failed_producers:
            return ()
        # Per claim, not per node. A claim is rendered under its own
        # input_id and the plan gate requires that id to be the node's
        # output, so each claimed value has its own lineage -- and
        # naming the whole node punished the ordinary, efficient shape
        # of batching eleven claims into one rendering node by
        # condemning all eleven for one.
        named: set[str] = set()
        for node in plan.analysis_nodes:
            if node.analysis_kind != "claim_rendering":
                continue
            by_id = {str(item.input_id): item for item in node.inputs}
            node_lineage = producers(node.node_id)
            for output in node.outputs:
                output_id = str(output.output_id)
                # A claim is rendered under its own input_id, and the
                # plan gate requires that id to be the node's output, so
                # each claimed value has its own lineage. Where a plan
                # does not pair them, the node's lineage still answers.
                item = by_id.get(output_id)
                if item is not None:
                    producer = str(getattr(item, "producer_node_id", "") or "")
                    lineage = (
                        {producer} | producers(producer) if producer else set()
                    )
                else:
                    lineage = node_lineage
                if lineage & failed_producers:
                    named.add(output_id)
        return tuple(sorted(named))

    def _claims_under_recorded_doubt(
        self, *, task_spec_sha256: str
    ) -> tuple[str, ...]:
        """Name claim quantities whose supporting receipt a doubt cites.

        Set intersection over digests the session itself minted: every
        ``doubt:{receipt_sha256}`` evidence reference on a task-bound
        scientific decision, against each rendered claim's
        source_receipt_sha256. No prose is read and no science judged;
        a prose-only doubt keeps the prior behavior, because the host
        cannot know which claim free text points at.
        """

        doubted = {
            str(reference)[len("doubt:") :]
            for record in self.scientific_decisions.values()
            if record.task_spec_sha256 == task_spec_sha256
            for reference in record.evidence_refs
            if str(reference).startswith("doubt:")
        }
        if not doubted:
            return ()
        return tuple(
            sorted(
                {
                    claim.quantity_id
                    for record in self.analysis_claim_records.values()
                    if getattr(record, "task_spec_sha256", "")
                    == task_spec_sha256
                    for claim in getattr(record, "claims", ())
                    if claim.source_receipt_sha256 in doubted
                }
            )
        )

    def _replication_receipt(
        self,
        *,
        node_id: str,
        context: Any,
        anomalies: Sequence[Mapping[str, Any]],
        source_receipt_sha256: str,
    ) -> AnomalyObservationV1 | None:
        """Replication before belief, as a superseding receipt.

        An excursion node cites the anomaly it investigates. When it
        ends with a verdict, the host compares the same sensor on the
        new result: the signal tripping again is ``replicated``, its
        silence is ``refuted``, and either is a second immutable receipt
        that supersedes the first. The host names the standing and
        never the meaning; whether a replicated anomaly is a discovery
        stays the scientist's claim.
        """

        cited = str(getattr(context, "excursion", "") or "")
        if not cited:
            return None
        known: dict[str, Mapping[str, Any]] = {
            str(item.get("receipt_sha256") or ""): item
            for item in self.prior_anomaly_observations
        }
        known.update(
            (digest, canonical_data(receipt))
            for digest, receipt in self.anomaly_observations.items()
        )
        source = known.get(cited)
        if source is None:
            return None
        signal_id = str(source.get("signal_id") or "")
        again = next(
            (
                item
                for item in anomalies
                if str(item.get("signal_id") or "") == signal_id
            ),
            None,
        )
        values = (
            {k: v for k, v in dict(again).items() if k != "signal_id"}
            if again is not None
            else {"signal_tripped": False}
        )
        return build_anomaly_observation(
            node_id=node_id,
            program=context.proposal.program,
            jobtype=context.proposal.jobtype,
            signal_id=signal_id,
            values=values,
            source_receipt_sha256=source_receipt_sha256,
            status="replicated" if again is not None else "refuted",
            supersedes_sha256=cited,
        )

    @staticmethod
    def _opened_result_output(receipt: Any) -> Any:
        """The parsed output behind a validation receipt, or None."""

        from chemsmart.analysis.result_readers import reader_for

        for artifact in getattr(receipt, "output_artifacts", ()) or ():
            try:
                return reader_for(
                    str(getattr(receipt, "program", ""))
                ).open_output(str(artifact.path))
            except Exception:  # noqa: BLE001 - not a readable result
                continue
        return None

    def _anomaly_output_ids(self) -> tuple[str, ...]:
        """Every anomaly the host recorded on this task, as output ids.

        The status rides in the id, so an unreplicated observation is a
        true word rather than a hidden one; the digest prefix names the
        receipt a decision may cite. A superseded receipt is history:
        the chain's latest receipt carries the standing.
        """

        records: list[Mapping[str, Any]] = [
            dict(item) for item in self.prior_anomaly_observations
        ]
        records.extend(
            canonical_data(receipt)
            for receipt in self.anomaly_observations.values()
        )
        records = list(anomaly_standing(records))
        return tuple(
            sorted(
                {
                    "anomaly:"
                    f"{record.get('signal_id')}:{record.get('status')}:"
                    f"{str(record.get('receipt_sha256') or '')[:8]}"
                    for record in records
                    if record.get("signal_id") and record.get("receipt_sha256")
                }
            )
        )

    def _record_toolchain_completion(
        self,
        policy_sha256: str,
        *,
        task_spec_sha256: str,
        source_receipt_sha256s: tuple[str, ...],
        status: str = "passed",
        findings: tuple[str, ...] = (),
        limitation_output_ids: tuple[str, ...] = (),
    ) -> tuple[str, ...]:
        """Record one toolchain-as-policy completion receipt and its event."""

        declared_misses, declared_limitations = (
            self._declared_observable_completion(
                task_spec_sha256=task_spec_sha256
            )
        )
        # Predictions ride the event beside the misses and never enter
        # the receipt body: a falsified expectation is a result, not a
        # defect, so it moves no status and no limitation.
        declared_predictions = self._declared_observable_predictions(
            task_spec_sha256=task_spec_sha256
        )
        limitation_output_ids = tuple(
            sorted(set(tuple(limitation_output_ids) + declared_limitations))
        )
        # A falsified pre-registration is a result, never a defect: it
        # rides the observation list under its own prefix, so the word
        # carries it and no status or limitation moves.
        falsified = tuple(
            f"falsified_expectation:{row.get('observable_id')}"
            for row in declared_predictions
            if row.get("agreement") == "diverged" and row.get("observable_id")
        )
        anomaly_output_ids = tuple(
            sorted(set(self._anomaly_output_ids()) | set(falsified))
        )
        body = {
            "schema_version": "chemsmart.analysis-completion-receipt.v1",
            # A scientific toolchain is already a visible, typed output
            # contract.  Its digest fills the existing aggregate gate's
            # policy identity without inventing a parallel policy file.
            "policy_sha256": policy_sha256,
            "task_spec_sha256": task_spec_sha256,
            "source_receipt_sha256s": source_receipt_sha256s,
            "status": status,
            "findings": tuple(findings),
        }
        if limitation_output_ids:
            body["limitation_output_ids"] = tuple(limitation_output_ids)
        if anomaly_output_ids:
            body["anomaly_output_ids"] = anomaly_output_ids
        completion = AnalysisCompletionReceiptV1(
            **body, receipt_sha256=canonical_sha256(body)
        )
        self.analysis_completion_receipts[completion.receipt_sha256] = (
            completion
        )
        completion_record = canonical_data(completion)
        completion_record.pop("receipt_sha256")
        if not completion.limitation_output_ids:
            # Mirror the digest body: the field is present only when
            # non-empty, so pre-field records and full deliveries share
            # one arithmetic.
            completion_record.pop("limitation_output_ids", None)
        if not completion.anomaly_output_ids:
            completion_record.pop("anomaly_output_ids", None)
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
                "critical_finding_count": len(completion.findings),
                "limitation_output_ids": completion.limitation_output_ids,
                "anomaly_output_ids": completion.anomaly_output_ids,
                "declared_observable_misses": declared_misses,
                "declared_observable_predictions": declared_predictions,
                "declared_observable_join_fields": dict(
                    sorted(self._declared_observable_join_fields.items())
                ),
                "completion_kind": "scientific_toolchain",
                "record": completion_record,
            },
            idempotency_key=(
                "scientific-toolchain-completion:" + completion.receipt_sha256
            ),
        )
        return (completion.receipt_sha256,)

    def evaluate_approved_toolchain_completion(
        self,
        plan: ScientificToolchainPlanV1,
        *,
        source_receipt_sha256s: tuple[str, ...],
    ) -> tuple[str, ...]:
        """Executor-side completion over an approved, fully executed chain.

        The session gate infers completion by matching receipts backward onto
        the plan; the executor needs no inference -- it dispatched every
        analysis call itself and holds the ledger of resulting receipt
        digests.  The gates that remain host-owned here: every approved
        calculation node must hold a validated execution receipt, and the
        ledger must be non-empty.  A scientific decision is deliberately not
        required -- interpretation stays a session act.
        """

        for node_id in plan.calculation_node_ids:
            receipt = self.execution_receipts.get(node_id)
            if receipt is None or not bool(
                getattr(receipt, "validated", False)
            ):
                raise ContractError(
                    "approved toolchain completion requires every "
                    f"calculation node validated; {node_id!r} is not"
                )
        receipts = tuple(sorted(set(source_receipt_sha256s)))
        if not receipts:
            raise ContractError(
                "approved toolchain completion requires at least one "
                "analysis receipt"
            )
        task_spec_sha256 = self._resolve_task_spec_reference(
            {}, "task_spec_sha256"
        )
        return self._record_toolchain_completion(
            plan.plan_sha256,
            task_spec_sha256=task_spec_sha256,
            source_receipt_sha256s=receipts,
        )

    def evaluate_partial_toolchain_completion(
        self,
        plan: ScientificToolchainPlanV1,
        *,
        source_receipt_sha256s: tuple[str, ...],
        findings: tuple[str, ...],
    ) -> tuple[str, ...]:
        """Bind what an interrupted approved chain actually produced.

        One benchmark lost 22 executed runs at this seam: engines validated,
        receipts landed, a late analysis node refused, and the reader got
        nothing.  A partial completion states the true outcome -- the
        receipts that exist plus findings naming every node that did not
        execute -- so the surviving evidence can be rendered without ever
        promoting an unvalidated number to a claim.  Findings are
        mandatory, and a calculation node that never validated is admitted
        only when the findings name it in the executor's fixed shape: in a
        batch, one record's failed engine run settles that record while
        the other records' receipts still deliver, and the non-delivery is
        rendered disclosure rather than a silent hole.
        """

        for node_id in plan.calculation_node_ids:
            receipt = self.execution_receipts.get(node_id)
            if receipt is None or not bool(
                getattr(receipt, "validated", False)
            ):
                expected = f"{node_id} (calculation):"
                if not any(item.startswith(expected) for item in findings):
                    raise ContractError(
                        "a partial toolchain completion must name every "
                        f"non-validated calculation node; {node_id!r} is "
                        "not named in its findings"
                    )
        if not findings:
            raise ContractError(
                "a partial toolchain completion must name its findings"
            )
        return self._record_toolchain_completion(
            plan.plan_sha256,
            task_spec_sha256=self._resolve_task_spec_reference(
                {}, "task_spec_sha256"
            ),
            source_receipt_sha256s=tuple(sorted(set(source_receipt_sha256s))),
            status="partial",
            findings=tuple(findings),
        )

    def completion_receipts_for_delivered_claims(self) -> tuple[str, ...]:
        """Certify a delivery made directly from registered results.

        A session may answer a task from results already in the
        workspace: extract, derive, evaluate, claim, decide. That route
        is ordinary and the charter names it, but finalisation knew only
        two ways to mint a certificate -- a task-owned analysis policy,
        or a preflighted workflow -- so a session with neither ended
        `blocked` however much it had delivered. SUFFICIENCY-1 recorded
        seventy claims and its scientific decision and returned to the
        human for want of a ceremony, which is the proximate reason that
        window produced no settlement anyone could read.

        The certificate is minted from what such a session actually has:
        the requirements it declared and the receipts its claims stand
        on. No plan is invented and no further model act is asked for --
        the same declared requirements and the same receipt graph as
        every other route, which is what makes the word comparable
        across them.
        """

        records = [
            record
            for record in self.analysis_claim_records.values()
            if getattr(record, "claims", ())
        ]
        if not records:
            raise ContractError(
                "no analysis claim has been recorded, so there is no "
                "delivery to certify"
            )
        if not self.scientific_decisions:
            raise ContractError(
                "claims are recorded and no scientific decision stands "
                "beside them; the delivery is not finished"
            )
        task_spec_sha256 = str(
            getattr(records[-1], "task_spec_sha256", "") or ""
        )
        sources = tuple(
            sorted(
                {
                    str(claim.source_receipt_sha256)
                    for record in records
                    for claim in record.claims
                }
            )
        )
        # The policy identity is the contract this delivery answers: the
        # goal's own declarations, in the order they were first made.
        policy_sha256 = canonical_sha256(
            {
                "schema_version": "chemsmart.delivered-claims-policy.v1",
                "declarations": tuple(
                    {
                        "observable_id": str(
                            record.get("observable_id") or ""
                        ),
                        "unit": str(record.get("unit") or ""),
                        "required_tolerance": (
                            float(record["required_tolerance"])
                            if record.get("required_tolerance") is not None
                            else -1.0
                        ),
                    }
                    for record in (
                        self.requested_observable_declarations.values()
                    )
                ),
            }
        )
        return self._record_toolchain_completion(
            policy_sha256,
            task_spec_sha256=task_spec_sha256,
            source_receipt_sha256s=sources,
        )

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
                stage_receipts["quantity_extraction"][
                    artifact_sha256
                ] = receipt.receipt_sha256
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
                stage_receipts["thermochemistry"][
                    artifact_sha256
                ] = receipt.receipt_sha256
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
        toolchain = self.scientific_toolchain_plans.get(
            completion.policy_sha256
        )
        if toolchain is not None:
            # The toolchain itself was the completion policy.  Claims are
            # rendered when a claim stage recorded them; a scientific decision
            # is a session act and may legitimately not exist yet.
            #
            # Several claim records are ordinary, not evidence of a break.
            # A live chain delivered two reduction potentials and their
            # difference from three claim stages, every node executed and
            # every verdict passed, and refusing the render here downgraded
            # that whole delivery to "partial" and manufactured a finding
            # whose text was about a binding rule rather than about any
            # chemistry -- while the report it then wrote rendered all
            # three records anyway.  The refusal prevented nothing and
            # mislabelled everything, so it is gone: each record renders
            # under its own label whatever the status.  One decision still
            # binds, because the renderer interprets through a single
            # decision and silently picking the first of several is the
            # ambiguity worth refusing.
            if completion.status == "passed" and len(decisions) > 1:
                raise ContractError(
                    "a toolchain completion binds at most one decision"
                )
            return self._render_toolchain_analysis_report(
                completion=completion,
                toolchain=toolchain,
                claim_records=tuple(
                    sorted(claim_records, key=lambda item: item.receipt_sha256)
                ),
                decision=decisions[0] if decisions else None,
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
            HOST_REPORT_TITLE,
            "",
            f"{COMPLETION_RECEIPT_LABEL}: `{completion.receipt_sha256}`",
            f"{CLAIM_RECORD_LABEL}: `{claims.receipt_sha256}`",
        ]
        if policy.required_conditions:
            lines.extend(
                (
                    "",
                    CONDITIONS_HEADING,
                    "",
                    f"| Condition | Value | Unit | Origin | {EVIDENCE_COLUMN} |",
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
                CLAIMS_HEADING,
                "",
                f"| Claim | Value | Unit | {SOURCE_RECEIPT_COLUMN} |",
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
        sections = tuple(
            zip(
                DECISION_SECTIONS,
                (
                    (decision.method_rationale,),
                    decision.assumptions,
                    decision.diagnostics,
                    decision.uncertainties,
                    decision.alternatives,
                ),
            )
        )
        for title, values in sections:
            entries = tuple(value for value in values if value)
            if not entries:
                continue
            lines.extend(("", f"## {title}", ""))
            lines.extend(f"- {value}" for value in entries)
        return "\n".join(lines)

    @staticmethod
    def _divergent_source_note(item: Any, phrase: str) -> str:
        """The receipt's display form, when it differs from canonical."""
        source_unit = str(getattr(item, "source_unit", "") or "")
        if not source_unit or source_unit == str(item.unit):
            return ""
        display = json.dumps(
            canonical_data(item.source_value),
            ensure_ascii=False,
            separators=(",", ":"),
        )
        return f" (= `{display}` {source_unit} {phrase})"

    def _render_toolchain_analysis_report(
        self,
        *,
        completion: AnalysisCompletionReceiptV1,
        toolchain: ScientificToolchainPlanV1,
        claim_records: tuple[Any, ...],
        decision: Any,
    ) -> str:
        """Render an approved toolchain's executed analysis, decision or not.

        A partial completion renders too: the headline states how much of
        the chain executed, the findings name what did not and why, and the
        receipts that survived are shown as evidence at their rung -- never
        as claims.  The host renders no number without a typed receipt, and
        an unvalidated value never gains claim standing by appearing here.
        """

        lines = [
            HOST_REPORT_TITLE,
            "",
            f"{COMPLETION_RECEIPT_LABEL}: `{completion.receipt_sha256}`",
            f"{TOOLCHAIN_PLAN_LABEL}: `{toolchain.plan_sha256}`",
        ]
        if completion.status == "partial":
            # The renderer cannot know how many nodes executed -- a finding
            # may describe an unexecuted node or a post-execution refusal --
            # so the header states only what the receipt actually carries
            # and the findings lines name each state.
            lines.extend(
                (
                    "",
                    f"{PARTIAL_STATUS_PREFIX}: "
                    f"{len(completion.findings)} finding(s) over "
                    f"{len(toolchain.analysis_nodes)} analysis nodes; "
                    f"{len(claim_records)} claim record(s) rendered.",
                    "",
                    FINDINGS_HEADING,
                    "",
                )
            )
            lines.extend(f"- {finding}" for finding in completion.findings)
        # Rendered whether the run reads completed or partial.  A chain can
        # deliver every output it declared it was for and still have been
        # refused a quantity along the way; saying so is what keeps "absence
        # is meaning" true once a refusal stops ending the whole extraction.
        # The reasons are read straight off the extraction receipts the
        # completion already binds -- the host states what it could not give
        # and never quietly gives something else.
        absences = tuple(
            (quantity_id, selector, reason)
            for digest in sorted(completion.source_receipt_sha256s)
            for quantity_id, selector, reason in getattr(
                self.quantity_extractions.get(digest), "absent", ()
            )
        )
        if absences:
            lines.extend(
                (
                    "",
                    ABSENCES_HEADING,
                    "",
                    "| Quantity | Selector | Why the result does not carry it |",
                    "|---|---|---|",
                )
            )
            lines.extend(
                f"| {quantity_id} | {selector} | {reason} |"
                for quantity_id, selector, reason in absences
            )
        conditions = tuple(
            node
            for node in toolchain.analysis_nodes
            if node.analysis_kind == "thermochemistry"
            and node.temperature_k is not None
        )
        if conditions:
            # Every convention the receipt binds is displayed: a 1 mol/L
            # solution free energy used to render as "1.0 atm" because the
            # concentration standard state, entropy model, and scale factor
            # were receipt-bound but invisible -- a 1.89 kcal/mol-per-species
            # convention hidden from the reader.
            lines.extend(
                (
                    "",
                    THERMO_CONDITIONS_HEADING,
                    "",
                    "| Stage | Temperature (K) | Standard state | "
                    "Entropy model | Frequency scale |",
                    "|---|---:|---|---|---:|",
                )
            )
            for node in conditions:
                if node.concentration_mol_l is not None:
                    standard_state = f"{node.concentration_mol_l:g} mol/L"
                else:
                    standard_state = f"{node.pressure_atm:g} atm"
                entropy_model = node.entropy_method
                if node.entropy_cutoff_cm1 is not None:
                    entropy_model += (
                        f" (cutoff {node.entropy_cutoff_cm1:g} cm^-1)"
                    )
                lines.append(
                    f"| {node.node_id} | `{node.temperature_k}` | "
                    f"{standard_state} | {entropy_model} | "
                    f"`{node.frequency_scale_factor:g}` |"
                )
        constant_names: list[str] = []
        for node in toolchain.analysis_nodes:
            if node.analysis_kind != "quantity_expression":
                continue
            for item in node.expression_nodes:
                if str(item.get("operation", "")) != "constant":
                    continue
                name = str(item.get("constant_name", ""))
                if name and name not in constant_names:
                    constant_names.append(name)
        if constant_names:
            # Every host-owned value the chain selected, with the
            # standard-state convention that gives it meaning.  The number
            # a reader must not have to take on faith is exactly the one a
            # model was never allowed to write.
            lines.extend(
                (
                    "",
                    LITERATURE_CONSTANTS_HEADING,
                    "",
                    "| Constant | Value | Unit | Family | Convention |",
                    "|---|---:|---|---|---|",
                )
            )
            families: list[str] = []
            for name in constant_names:
                try:
                    entry = literature_constant(name)
                except UnknownLiteratureConstantError:
                    # Reachable only if a registry entry was removed after
                    # this chain was approved; the report must still say
                    # why rather than fail to render.
                    lines.append(
                        f"| {name} | — | — | — | no longer registered |"
                    )
                    continue
                if (
                    entry.convention_family != "independent"
                    and entry.convention_family not in families
                ):
                    families.append(entry.convention_family)
                lines.append(
                    f"| {entry.name} | `{entry.value:g}` | {entry.unit} | "
                    f"{entry.convention_family} | {entry.convention} |"
                )
            if len(families) > 1:
                # An absolute electrode potential means one thing beside the
                # proton solvation free energy determined on the same scale
                # and another beside a different one, and the literature
                # circulates the halves separately. A crossed pair is silent:
                # both values are correct, the units are right, nothing
                # diverges, and the answer moves by more than method error.
                lines.append("")
                lines.append(
                    "This analysis combined "
                    f"{len(families)} convention families "
                    f"({', '.join(sorted(families))}). Constants determined "
                    "on different scales are not interchangeable even when "
                    "each is correct on its own."
                )
        # What the host made of each precision requirement, on the page a
        # human reads. `grep -rn sufficiency chemsmart/agent/tui/`
        # returned nothing and this report had no row either, so the
        # 2026-09-10 boundary's claim that `met` is "auditable where a
        # human meets it" held on the settlement's achieved branch and
        # on no rendered surface at all (SUFFICIENCY-5).
        if self.requirement_assessments:
            lines.extend(
                (
                    "",
                    "## Precision requirements",
                    "",
                    "| Observable | Required | Stated | Basis | Whose | "
                    "Combined by | State | The host observed |",
                    "|---|---:|---:|---|---|---|---|---|",
                )
            )
            for observable_id, row in sorted(
                self.requirement_assessments.items()
            ):
                observed = ", ".join(
                    str(item)
                    for item in row.get("uncertainty_observations") or ()
                )
                unquantified = tuple(row.get("unquantified_components") or ())
                if unquantified:
                    observed = (observed + "; " if observed else "") + (
                        f"{len(unquantified)} term(s) named unquantified"
                    )
                unit = str(row.get("unit") or "")
                stated = row.get("uncertainty")
                lines.append(
                    f"| {observable_id} "
                    f"| {row.get('required_tolerance')} {unit} "
                    f"| {'—' if stated is None else stated} {unit} "
                    f"| {row.get('uncertainty_basis') or '—'} "
                    f"| {row.get('tolerance_origin') or 'unstated'} "
                    f"| {(row.get('uncertainty_combination') or {}).get('rule') or 'not stated'} "
                    f"| {row.get('state')} "
                    f"| {observed or '—'} |"
                )
            lines.extend(
                (
                    "",
                    "Only `met` discharges a requirement, and it states "
                    "that the host resolved a stated uncertainty inside a "
                    "declared tolerance -- never that the estimate is "
                    "adequate, which is the session's claim to defend.",
                )
            )

        for claims in claim_records:
            lines.extend(
                (
                    "",
                    f"{CLAIM_RECORD_LABEL}: `{claims.receipt_sha256}`",
                    "",
                    CLAIMS_HEADING,
                    "",
                    f"| Claim | Value | Unit | {SOURCE_RECEIPT_COLUMN} |",
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
                    f"| {claim.claim_id} | `{display}` | {claim.display_unit} "
                    f"| `{claim.source_receipt_sha256}` |"
                )
        sources = set(completion.source_receipt_sha256s)
        if completion.status == "partial":
            # Each row carries the canonical value and, when the receipt
            # holds a different display form, that form beside it. A live
            # report rendered a quantity named e-diff-kjmol as its
            # canonical 0.1545 hartree while the 405.7 kJ/mol the plan's
            # convert had produced sat unprinted in the receipt -- a
            # reader trusting the row name over the unit column quotes
            # the wrong number. An expression's source form is the one
            # the plan requested; a parsed quantity's is the one the
            # program printed; thermochemistry's source unit is the
            # engine's internal representation and is not annotated.
            survivors: list[tuple[str, Any, str, str, str, str]] = []
            for digest in sorted(sources):
                extraction = self.quantity_extractions.get(digest)
                if extraction is not None:
                    survivors.extend(
                        (
                            item.quantity_id,
                            item.value,
                            item.unit,
                            "parsed",
                            digest,
                            self._divergent_source_note(item, "as printed"),
                        )
                        for item in extraction.quantities
                    )
                thermochemistry = self.thermochemistry_receipts.get(digest)
                if thermochemistry is not None:
                    survivors.extend(
                        (
                            item.quantity_id,
                            item.value,
                            item.unit,
                            "derived",
                            digest,
                            "",
                        )
                        for item in thermochemistry.quantities
                    )
                expression = self.quantity_expression_receipts.get(digest)
                if expression is not None:
                    survivors.extend(
                        (
                            item.quantity_id,
                            item.value,
                            item.unit,
                            "derived",
                            digest,
                            self._divergent_source_note(item, "as requested"),
                        )
                        for item in expression.outputs
                    )
            if survivors:
                lines.extend(
                    (
                        "",
                        SURVIVING_HEADING,
                        "",
                        "| Quantity | Value | Unit | Evidence rung | "
                        f"{SOURCE_RECEIPT_COLUMN} |",
                        "|---|---:|---|---|---|",
                    )
                )
                for quantity_id, value, unit, rung, digest, note in survivors:
                    display = json.dumps(
                        canonical_data(value),
                        ensure_ascii=False,
                        separators=(",", ":"),
                    )
                    lines.append(
                        f"| {quantity_id} | `{display}`{note} | {unit} "
                        f"| {rung} | `{digest}` |"
                    )
        verdict_rows = []
        for digest in sorted(sources):
            validation = self.scientific_validation_receipts.get(digest)
            if validation is None:
                continue
            for result in validation.rule_results:
                verdict_rows.append(
                    (
                        validation.node_id,
                        result.rule_id,
                        result.predicate,
                        result.observed_value,
                        (
                            result.threshold
                            if result.threshold is not None
                            else result.expected_count
                        ),
                        result.unit,
                        "passed" if result.passed else "failed",
                    )
                )
        if verdict_rows:
            lines.extend(
                (
                    "",
                    VERDICTS_HEADING,
                    "",
                    "| Node | Rule | Predicate | Observed | Bound | Unit | "
                    "Verdict |",
                    "|---|---|---|---:|---:|---|---|",
                )
            )
            for row in verdict_rows:
                node_id, rule_id, predicate, observed, bound, unit, word = row
                lines.append(
                    f"| {node_id} | {rule_id} | {predicate} | `{observed}` "
                    f"| `{bound if bound is not None else ''}` | {unit} | "
                    f"{word} |"
                )
        prediction_rows = self._declared_observable_predictions(
            task_spec_sha256=completion.task_spec_sha256
        )
        if prediction_rows:
            # Beside the verdicts and never among them: a verdict is a
            # criterion the delivery had to meet, an expectation is what
            # a scientist thought beforehand, and only the first can
            # fail.
            lines.extend(
                (
                    "",
                    PREDICTIONS_HEADING,
                    "",
                    "| Observable | Expected | Delivered | Unit | "
                    "Agreement | Basis |",
                    "|---|---|---:|---|---|---|",
                )
            )
            for row in prediction_rows:
                expected = row["expected_sign"]
                if row["expected_low"] != "":
                    band = f"{row['expected_low']}..{row['expected_high']}"
                    expected = f"{expected} {band}".strip()
                lines.append(
                    f"| {row['observable_id']} | {expected} | "
                    f"`{row['delivered_value']}` | {row['delivered_unit']} "
                    f"| {row['agreement']} | {row['expectation_basis']} |"
                )
            lines.extend(
                (
                    "",
                    "An expectation is displayed, never scored: a "
                    "diverging row settles nothing and means the "
                    "chemistry disagreed with the reasoning, which is a "
                    "result the reader owns.",
                )
            )
        if decision is not None:
            sections = tuple(
                zip(
                    DECISION_SECTIONS,
                    (
                        (decision.method_rationale,),
                        decision.assumptions,
                        decision.diagnostics,
                        decision.uncertainties,
                        decision.alternatives,
                    ),
                )
            )
            for title, values in sections:
                entries = tuple(value for value in values if value)
                if not entries:
                    continue
                lines.extend(("", f"## {title}", ""))
                lines.extend(f"- {value}" for value in entries)
        else:
            lines.extend(
                (
                    "",
                    f"{NO_DECISION_PREFIX} -- interpretation is "
                    "a session act; this run executed extraction, "
                    "thermochemistry, expressions, validation verdicts, and "
                    "claim rendering only.",
                )
            )
        if completion.status == "partial":
            lines.extend(
                (
                    "",
                    f"{RECOVERY_PREFIX}: the engine outputs and every "
                    "receipt above are unchanged and remain readable. A "
                    "later explicit analysis request -- a new session over "
                    "this workspace planning an analysis-only toolchain on "
                    "the registered results -- may re-extract, validate, "
                    "and claim them without re-running any engine.",
                )
            )
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
                {
                    role: self.hessian_handoffs[node_id]
                    for role in ("hess_filename", "inhess_filename")
                    if role in dict(context.job_artifact_options)
                    and node_id in self.hessian_handoffs
                }
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
                approved = execution_invocation.environment_identity_sha256 in (
                    set(
                        self.frozen_workflow_approval.environment_identity_sha256s
                    )
                )
            if not approved:
                raise ContractError(
                    "execution environment differs from the exact frozen "
                    "node approval"
                )
        node_workspace = node_branch_directory(
            self.approved_workspace, node_id, cycle_label=self.cycle_label
        )
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
        # The approved node bindings live on the one-shot approval, not
        # on the frozen approval; reading them off the wrong object
        # crashed every launch of a whole sealed window (void E4 window
        # 1, 2026-09-03) while the suite stayed green, because no test
        # drives a node through this method past the approval checks.
        excursion_node_ids = frozenset(
            binding.node_id
            for binding in approval.node_bindings
            if getattr(binding, "excursion", "")
        )
        effective_timeout_seconds = self._require_bounded_launch_budget(
            excursion=node_id in excursion_node_ids,
            excursion_node_ids=excursion_node_ids,
        )
        fence = self.event_store.reserve_workflow_node_launch(
            turn_id=turn_id,
            plan=scientific_plan,
            materialized_workflow=materialized,
            approval=frozen_approval,
            invocation=execution_invocation,
            run_id=v2_run_id,
            timestamp=started,
            # The host's own bound on how long this engine may live, so a
            # concurrent sibling holding a node can be told from a
            # process that died holding one. Plus the postprocessing
            # reserve, because the host's own evaluation happens inside
            # the same reservation.
            lease_seconds=_launch_lease_seconds(
                self.execution_resources, self.bounded_execution_envelope
            ),
            reserver=_launch_reserver(),
            # The approved grant, checked inside the fence's own lock.
            # The process-local count below stays as a cheap early
            # refusal; this is the one that is exact when several array
            # elements reserve at the same moment.
            max_engine_calls=(
                None
                if self.bounded_execution_envelope is None
                else self.bounded_execution_envelope.max_engine_calls
            ),
            max_excursion_calls=(
                None
                if self.bounded_execution_envelope is None
                else self.bounded_execution_envelope.max_excursion_calls
            ),
            excursion_node_ids=frozenset(excursion_node_ids),
            excursion=node_id in excursion_node_ids,
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
        # The branch is written before the engine starts, so a node that
        # dies still says what was asked of it. A chemist opening this
        # folder sees the exact CHEMSMART command and the exact project
        # configuration it ran from, beside the engine's own files.
        _write_branch_request(node_workspace, command)
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
                    # The wrapper traps SIGTERM to stop the engine and copy
                    # its partial outputs back from scratch before exiting.
                    # The default one-second grace SIGKILLed it mid-copy, so a
                    # timed-out scan's converged points -- data the parsers
                    # deliberately read as a truncated surface -- never
                    # reached the job folder. Thirty seconds bounds the
                    # salvage without changing what the timeout means.
                    termination_grace_seconds=30.0,
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
        # The engine's window ends here; what follows is the host's own
        # reading, timed separately and never charged to the engine wall.
        finished = datetime.now(timezone.utc).isoformat()
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
            expected_root_artifact=context.root_artifact,
            expected_project_artifact=context.project_artifact,
            output_artifacts=outputs,
            exit_status=wrapper_exit_status,
            expected_environment_receipt_sha256=(
                execution_invocation.environment_receipt_sha256
            ),
            capability_environment_receipt=capability_environment,
            pyscf_engine_observation=pyscf_engine,
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
            # Retired mechanism, retained digest field (always empty).
            stationary_point_policy_sha256="",
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
        sibling_observations = ()
        if result_validation_receipt.state == "valid":
            sibling_observations = _same_structure_observations(
                self.result_validation_receipts,
                node_id,
                self._opened_result_output(result_validation_receipt),
                handoffs=self.handoffs,
                input_sha256=str(
                    getattr(
                        result_validation_receipt,
                        "input_artifact_sha256",
                        "",
                    )
                    or ""
                ),
                output_sha256s=tuple(
                    item.sha256
                    for item in result_validation_receipt.output_artifacts
                ),
            )
        for anomaly in (*evaluation.anomalies, *sibling_observations):
            observation_receipt = build_anomaly_observation(
                node_id=node_id,
                program=context.proposal.program,
                jobtype=context.proposal.jobtype,
                signal_id=str(anomaly["signal_id"]),
                values={
                    key: value
                    for key, value in dict(anomaly).items()
                    if key != "signal_id"
                },
                source_receipt_sha256=(
                    result_validation_receipt.receipt_sha256
                ),
                flagged_artifact_sha256s=tuple(
                    item.sha256
                    for item in result_validation_receipt.output_artifacts
                ),
            )
            self.anomaly_observations[observation_receipt.receipt_sha256] = (
                observation_receipt
            )
            self._emit(
                turn_id,
                EventKind.ANOMALY_OBSERVED,
                observation_receipt.receipt_sha256,
                status=observation_receipt.status,
                node_id=node_id,
                signal_id=observation_receipt.signal_id,
                record=canonical_data(observation_receipt),
            )
        replication = self._replication_receipt(
            node_id=node_id,
            context=context,
            anomalies=evaluation.anomalies,
            source_receipt_sha256=result_validation_receipt.receipt_sha256,
        )
        if replication is not None:
            self.anomaly_observations[replication.receipt_sha256] = replication
            self._emit(
                turn_id,
                EventKind.ANOMALY_OBSERVED,
                replication.receipt_sha256,
                status=replication.status,
                node_id=node_id,
                signal_id=replication.signal_id,
                record=canonical_data(replication),
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
                else "engine_complete" if engine_complete else "failed"
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
        evaluated = datetime.now(timezone.utc).isoformat()
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
            evaluated_at=evaluated,
        )
        self.event_store.record_program_execution_receipt(
            turn_id=turn_id,
            workflow_id=scientific_plan.workflow_id,
            run_id=v2_run_id,
            receipt=receipt,
            excursion=node_id in excursion_node_ids,
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
                        item.selection_rule != "validated_optimized_geometry",
                        item.consumer_node_id,
                    ),
                )
            )
            for edge in outgoing_edges:
                if edge.producer_node_id != node_id:
                    continue
                consumer_binding = approval.node(edge.consumer_node_id)
                hessian_role = hessian_role_for_rule(edge.selection_rule)
                if edge.selection_rule == "validated_final_orca_ts_hessian":
                    if (
                        context.proposal.program
                        != hessian_role.producer_program
                        or context.proposal.jobtype
                        not in hessian_role.producer_stages
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
                elif edge.selection_rule == "validated_producer_orca_hessian":
                    if (
                        context.proposal.program
                        != hessian_role.producer_program
                        or context.proposal.jobtype
                        not in hessian_role.producer_stages
                    ):
                        raise ContractError(
                            "a producer ORCA Hessian handoff requires a "
                            "frequency-bearing ORCA producer"
                        )
                    result_candidates = tuple(
                        item for item in outputs if item.kind == "orca_output"
                    )
                    hessian_candidates = tuple(
                        item for item in outputs if item.kind == "orca_hessian"
                    )
                    if len(result_candidates) != 1 or not hessian_candidates:
                        raise ContractError(
                            "a validated ORCA producer requires one output "
                            "and at least one native Hessian candidate"
                        )
                    artifact, observed = (
                        handoff_validated_orca_producer_hessian(
                            producer_receipt=receipt,
                            result_artifact=result_candidates[0],
                            hessian_candidates=hessian_candidates,
                            producer_edge=edge,
                            approved_workspace=self.approved_workspace,
                            hessian_artifact_id=(
                                f"hessian.{edge.producer_node_id}-to-"
                                f"{edge.consumer_node_id}"
                            ),
                            expected_charge=(
                                context.scientific_identity.charge
                            ),
                            expected_multiplicity=(
                                context.scientific_identity.multiplicity
                            ),
                        )
                    )
                    if observed.consumer_state != (
                        consumer_binding.charge,
                        consumer_binding.multiplicity,
                    ):
                        raise ContractError(
                            "an ORCA TS search must start from a Hessian on "
                            "its own charge and multiplicity surface"
                        )
                    self.artifacts[artifact.artifact_id] = artifact
                    self.hessian_handoffs[edge.consumer_node_id] = observed
                elif edge.selection_rule == "validated_scan_minimum_geometry":
                    if context.proposal.program != "orca":
                        raise ContractError(
                            "a scan-minimum geometry handoff requires an "
                            "ORCA scan producer"
                        )
                    candidates = tuple(
                        item for item in outputs if item.kind == "orca_output"
                    )
                    if len(candidates) != 1:
                        raise ContractError(
                            "a validated ORCA scan requires exactly one "
                            "orca_output result"
                        )
                    artifact, observed = handoff_scan_minimum_geometry(
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
                if edge.selection_rule in {
                    "validated_optimized_geometry",
                    "validated_scan_minimum_geometry",
                }:
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
                    self.scientific_identities[identity.binding_sha256] = (
                        identity
                    )
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
                if edge.selection_rule in {
                    "validated_optimized_geometry",
                    "validated_scan_minimum_geometry",
                }:
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

    def _require_bounded_launch_budget(
        self,
        *,
        excursion: bool = False,
        excursion_node_ids: frozenset[str] = frozenset(),
    ) -> float:
        """Return this launch's timeout while preserving analysis time.

        The envelope's node timeout is an upper bound, not a requirement that
        the full duration remain available before every node.  Planning and
        earlier nodes may consume part of the episode; a later launch receives
        the smaller live window left after reserving postprocessing time.
        An excursion launch is counted against the grant's own line, and
        the engine-call line never sees it.
        """

        envelope = self.bounded_execution_envelope
        if envelope is None:
            return float(self.execution_resources.node_timeout_seconds)
        excursions_used = sum(
            1
            for launched in self.execution_receipts
            if launched in excursion_node_ids
        )
        plain_used = len(self.execution_receipts) - excursions_used
        if excursion:
            if excursions_used >= envelope.max_excursion_calls:
                raise ContractError(
                    "bounded execution excursion grant exhausted"
                )
        elif plain_used >= envelope.max_engine_calls:
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
        """Plan stages that will not execute, and must not block the rest.

        Two ways a stage earns this, and both are narrowings. Release maturity
        is a host fact: a family outside the executable matrix cannot run
        whatever the plan asked for. And a *declared* ``blocked_unsupported``
        stage is the scientist retaining intent the program cannot realise --
        which can be true even when the job family itself is executable. A
        real session hit exactly that: it planned the paper's third functional,
        ORCA's validator refused it (the functional is not implemented there),
        and the session correctly declared those stages non-executable intent
        -- but this predicate honoured the declaration only for non-executable
        *families*, so five nodes riding on plain ``orca/opt`` and
        ``orca/scan`` kept blocking an approval they were never going to be
        part of.

        Honouring the declaration is safe by construction: a blocked stage is
        displayed with the workflow, excluded from the approval, and never
        launched, so the only thing the declaration can do is narrow what
        runs. The node contract already requires a stated reason.
        """

        non_executable: set[str] = set()
        for node in plan.nodes:
            # A family outside the executable matrix that is *not* declared
            # blocked stays a blocker on purpose: the planner is required to
            # state the blockage, keeping this state explicit, never inferred.
            if node.support_state == "blocked_unsupported":
                non_executable.add(node.node_id)
        # Non-executability cascades along data edges. A consumer of a
        # blocked producer can never receive its input, so it is
        # non-executable by implication -- a pure narrowing the host can
        # derive, sparing the planner from hand-marking every downstream
        # node of a stage it already declared. Without this, the readiness
        # projection called such a workflow approvable while the review
        # builder refused it ("an executed node cannot consume the output of
        # a stage this release cannot execute"), and a session ended
        # honest-looking with no packet. Observed live on a cluster-continuum
        # plan whose complex-optimisation stage was declared blocked and fed
        # the TD stage.
        while True:
            grew = False
            for edge in getattr(plan, "edges", ()):
                if (
                    getattr(edge, "edge_kind", "") == "data"
                    and edge.source_node_id in non_executable
                    and edge.target_node_id not in non_executable
                ):
                    non_executable.add(edge.target_node_id)
                    grew = True
            if not grew:
                break
        return frozenset(non_executable)

    def _non_executable_reasons(
        self, plan: ScientificWorkflowPlanV2 | None
    ) -> dict[str, str]:
        """Per-node reasons for the same set the readiness projection uses.

        The frontier used to compute actionability without consulting
        support state at all, so it told a session "compile_and_preview" for
        a stage the readiness was simultaneously excluding as
        non-executable, and the session had to flag the contradiction
        itself. One truth source: the declared-plus-cascaded set above,
        joined with the declared reason where one exists and the naming of
        the blocked producer where the state is inherited.
        """

        if plan is None:
            return {}
        non_executable = self._release_non_executable_node_ids(plan)
        declared = {
            node.node_id: node.blocked_reason
            for node in plan.nodes
            if node.support_state == "blocked_unsupported"
        }
        reasons: dict[str, str] = {}
        for node_id in non_executable:
            if node_id in declared:
                reasons[node_id] = (
                    declared[node_id] or "declared non-executable intent"
                )
            else:
                sources = sorted(
                    edge.source_node_id
                    for edge in getattr(plan, "edges", ())
                    if getattr(edge, "edge_kind", "") == "data"
                    and edge.target_node_id == node_id
                    and edge.source_node_id in non_executable
                )
                reasons[node_id] = (
                    "consumes the output of non-executable stage "
                    + (", ".join(sources) if sources else "upstream")
                )
        return reasons

    def execution_review_ineligibility_reason(
        self,
        *,
        plan: ScientificWorkflowPlanV2,
        planned_node: ScientificWorkflowNodeV2,
    ) -> str:
        """Explain why one valid plan node cannot enter human execution review.

        A stage this release declares preview-only is deferred rather than
        ineligible: see ``_release_non_executable_node_ids``.

        This method answers a question; it never raises to ask it.  A
        ``ContractError`` is ChemSmart's refusal mechanism and its message
        already names the boundary that was hit, so here it *is* the answer.
        Letting one escape turned a stated refusal into exit 1: the caller in
        ``live_session`` builds its ineligible-node list by asking this for
        every node, so a single unresolvable input discarded the report for
        the whole session -- plan, previews and all -- and the run looked like
        a crash rather than the honest "this cannot run because ..." it was.
        Other exception types still propagate, because a programming fault
        must not be able to disguise itself as a scientific reason.
        """

        try:
            return self._execution_review_ineligibility_reason(
                plan=plan, planned_node=planned_node
            )
        except ContractError as exc:
            return str(exc)

    def _execution_review_ineligibility_reason(
        self,
        *,
        plan: ScientificWorkflowPlanV2,
        planned_node: ScientificWorkflowNodeV2,
    ) -> str:
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
                return "ORCA transition-state review lacks validated project"
            settings = dict(context.project_validation.settings)
            if not bool(settings.get("freq") or settings.get("numfreq")):
                return (
                    "ORCA transition-state execution requires a requested "
                    "frequency analysis to establish a first-order saddle"
                )
        return ""

    def execution_review_wanted(self) -> bool:
        """Whether this session's ending should carry the host's own review.

        True when a bounded envelope requested an inert review, a
        calculation plan exists and is materialised, and the host knows
        the workspace the request is filed against.
        """

        if self.bounded_execution_envelope is None:
            return False
        if self.run_evidence_root is None:
            return False
        if not any(
            plan.nodes for plan in self.scientific_workflow_plans.values()
        ):
            return False
        return self.bounded_review_is_materialized()

    def prepare_execution_review(self) -> WorkflowExecutionReviewV1:
        """Build the review the readiness gates promise, while the runtime
        stream is still open.

        A planning session that reached "host readiness gates passed"
        used to have its review built after the loop had sealed the
        stream; when the builder refused, the refusal could not be
        recorded (the store is absorbing), the ContractError escaped,
        and the goal settled on a Python error with the reason lost --
        with its whole grant unspent (REACH-1 ino3, 2026-09-06; the same
        ending on 2026-09-03). The same checks the session runner made
        run here, and a refusal names what it refused.
        """

        # The plan the session is standing on, through the one function
        # every reader shares. This was
        # `[p for p in ...values() if p.nodes][-1]`, a *fourth* reader
        # resolving the question by insertion order -- and dd29df23
        # repaired only the three spelled `tuple(...)[-1]`, so the
        # eligibility loop below judged the plan a restoring amendment
        # had abandoned while `build_execution_review` built the
        # reviewed packet from the restored one. Proven by probe
        # (2026-09-11): after A -> B -> A the gate answered B and the
        # packet answered A, which puts an unchecked node into the
        # single human decision.
        plan = self._current_scientific_plan()
        if plan is None or not plan.nodes:
            raise ContractError(
                "execution review requires a scientific workflow with "
                "at least one node"
            )
        ineligible = []
        for node in plan.nodes:
            reason = self.execution_review_ineligibility_reason(
                plan=plan, planned_node=node
            )
            if reason:
                ineligible.append(
                    f"{node.node_id} ({node.program}/{node.engine}/"
                    f"{node.stage}: {reason})"
                )
        try:
            if ineligible:
                raise ContractError(
                    "execution review is not eligible: "
                    + "; ".join(ineligible)
                )
            review = self.build_execution_review(
                workspace=self.run_evidence_root
            )
        except ContractError as exc:
            self.execution_review_refusal = {
                "workflow_id": plan.workflow_id,
                "reason": str(exc),
            }
            raise
        self.prepared_execution_review = review
        return review

    def bounded_review_is_materialized(self) -> bool:
        """Whether the plan a review would use has been materialised.

        ``build_execution_review`` refuses an unmaterialised plan, which is
        the right contract at the wrong moment.  A session that inspected
        capability, bound an environment, rendered, promoted and validated
        project YAML, planned a workflow and then ran out of turns has
        produced real evidence a human should still see; raising instead of
        reporting discards the whole record, transcript included.  Asking
        first lets the caller report the session it actually got, using the
        preview-only status the result vocabulary already carries.

        This answers exactly the question the review builder would fail on,
        so a deliberate refusal -- an unsupported node, an exceeded engine
        budget -- still surfaces as before.
        """

        plan = self._current_scientific_plan()
        if plan is None:
            return False
        return any(
            workflow.plan_sha256 == plan.plan_sha256
            for workflow in self.materialized_workflows.values()
        )

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
            raise ContractError(
                "review resources differ from execution envelope"
            )
        plan = self._current_scientific_plan()
        if plan is None:
            raise ContractError(
                "execution review requires a scientific workflow"
            )
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
            analysis_only = (
                "; this envelope's budget is 0, so it authorises analysis "
                "over registered results only"
                if envelope.max_engine_calls == 0
                else ""
            )
            raise ContractError(
                "scientific workflow exceeds bounded engine-call budget: "
                f"{len(executable_nodes)} nodes for "
                f"{envelope.max_engine_calls} calls" + analysis_only
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
            elif is_validated_scan_minimum_geometry_edge(plan, edge):
                selection_rule = "validated_scan_minimum_geometry"
            elif is_validated_orca_ts_hessian_edge(plan, edge):
                selection_rule = "validated_final_orca_ts_hessian"
            elif is_validated_producer_orca_hessian_edge(plan, edge):
                selection_rule = "validated_producer_orca_hessian"
            else:
                raise ContractError(
                    "execution review has no exact selection rule for data "
                    f"edge {edge.edge_id!r}; expected optimized geometry, "
                    "an ORCA scan minimum-energy point geometry, an ORCA "
                    "final-TS Hessian for IRC, or an ORCA producer "
                    "Hessian for a TS inhess_filename role"
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
            if edge.selection_rule
            in {
                "validated_optimized_geometry",
                "validated_scan_minimum_geometry",
            }
        )
        edge_by_target = {
            edge.consumer_node_id: edge for edge in geometry_edges
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
        node_observations: list[dict[str, Any]] = []
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
            # A node may deliberately reuse a producer geometry on another
            # charge/multiplicity surface, which is how a redox or
            # hydrogen-transfer series is written.  That freedom is exactly
            # where an impossible pair enters, so the state is checked against
            # this molecule's electron count before it reaches a review.  No
            # producer rule changes the atom set, so the input artifact's
            # symbols are the consumer's symbols.
            try:
                _node_symbols = _review_molecule_identity(
                    context.input_artifact
                )["atom_order"]
            except Exception:
                _node_symbols = ()
            if _node_symbols:
                refuse_impossible_electronic_state(
                    _node_symbols,
                    target_charge,
                    target_multiplicity,
                    context=f"node {planned_node.node_id!r}",
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
                invocation, invocation_context = (
                    self._plan_invocation_for_node(
                        plan=plan, node_id=planned_node.node_id
                    )
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
                    # A node whose geometry comes from a producer compiles
                    # here rather than through the prepare path, so a driven
                    # coordinate omitted here is dropped for exactly the
                    # workflows that need it most: a scan of a structure some
                    # earlier stage optimised. Observed on a real paper task,
                    # where two correctly specified torsion scans reached the
                    # human review as a bare `scan` with no range at all.
                    job_option_values=native_coordinate_options(
                        context.proposal.program, context.internal_coordinates
                    ),
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
            composition = self.molecular_compositions.get(
                context.input_artifact.sha256
            )
            if composition is not None:
                # The human review displays the composed arrangement's full
                # lineage: which approved parents, which contact, at what
                # distance -- the single /approve covers exactly this
                # displayed construction.
                approved_parent_sha256s = {
                    identity.geometry_sha256
                    for identity in (
                        self.approved_molecular_identities.values()
                    )
                }
                parents_approved = {
                    composition.fragment_a_sha256,
                    composition.fragment_b_sha256,
                }.issubset(approved_parent_sha256s)
                composition_record = canonical_data(composition)
                molecular_identity = {
                    "identity_evidence_status": (
                        "composed-from-approved-parents"
                        if parents_approved
                        else "composed-task-bound"
                    ),
                    **_review_molecule_identity(context.input_artifact),
                    "composition": composition_record,
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
            elif approved_identity is None:
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
            # Host-owned input lineage travels into the displayed review
            # regardless of which identity branch applied.  A derived
            # species' kept/removed atoms are exactly the per-record facts
            # a reviewer needs on the page (a wrong-atom removal was once
            # caught only from the artifact bytes), and a database-record
            # extraction carries its stored fields as labelled
            # observations beside the explicitly bound state.
            # A built geometry may be a CHAIN: the cis rotamer with both
            # methyls staggered is three edits deep, and the hop that
            # decides which rotamer the molecule is sits at the root.  The
            # first chained build displayed only the final hop, hiding
            # exactly the edit a reviewer would need to catch, so the walk
            # follows parent digests through every edit and append and
            # attaches the whole chain root-first, ending on a derivation
            # or extraction when one anchors it.
            # Derived from `_GEOMETRY_HOPS` and `_GEOMETRY_ORIGINS`
            # rather than hand-listed here, because the sibling walk
            # hand-listed a different subset and they drifted: a
            # composed or PubChem-fetched molecule reached this page
            # with no origin hop, while a derived one carried its
            # panel. `pubchem_geometries` had no reader anywhere in the
            # tree (SUFFICIENCY-5, 2026-09-10).
            hops, origin = self._geometry_provenance(
                context.input_artifact.sha256
            )
            geometry_lineage = list(hops)
            if origin is not None:
                origin_kind, origin_receipt = origin
                molecular_identity[origin_kind] = canonical_data(
                    origin_receipt
                )
            if geometry_lineage:
                molecular_identity["geometry_lineage"] = tuple(
                    reversed(geometry_lineage)
                )
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
            review_atom_count = int(molecular_identity.get("atom_count") or 0)
            if not review_atom_count:
                try:
                    review_atom_count = int(
                        _review_molecule_identity(context.input_artifact).get(
                            "atom_count"
                        )
                        or 0
                    )
                except Exception:
                    review_atom_count = 0
            stated: tuple[str, ...] = ()
            if review_atom_count:
                stated = compile_time_observations(
                    program=planned_node.program,
                    jobtype=planned_node.stage,
                    settings=validation.settings,
                    atom_count=review_atom_count,
                    geometry=_geometry_for_observation(context.input_artifact),
                )
            probe = getattr(self, "_input_check_by_node", {}).get(
                planned_node.node_id
            )
            if probe is not None:
                from chemsmart.agent.input_check import (
                    probe_observation_lines,
                )

                stated = tuple(stated) + probe_observation_lines(probe)
            if stated:
                node_observations.append(
                    {
                        "node_id": planned_node.node_id,
                        "observations": stated,
                    }
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
                    internal_coordinates=context.internal_coordinates,
                    excursion=context.excursion,
                    root_artifact_id=(
                        context.root_artifact.artifact_id
                        if context.root_artifact is not None
                        else ""
                    ),
                    root_artifact_sha256=(
                        context.root_artifact.sha256
                        if context.root_artifact is not None
                        else ""
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
                f"{node_id}: {reason}"
                for node_id, reason in sorted(unbindable)
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
            str(request_id).strip() or "review-" + plan.plan_sha256[:16]
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
            non_executable_node_ids=tuple(sorted(non_executable_ids)),
            scientific_toolchain_plan=self._toolchain_plan_for_review(plan),
            requested_observable_declarations=tuple(
                record
                for _observable_id, record in sorted(
                    self.requested_observable_declarations.items()
                )
            ),
            node_observations=tuple(
                sorted(node_observations, key=lambda item: item["node_id"])
            ),
            consulted_domain_knowledge=tuple(
                sorted(
                    self.consulted_skill_records.values(),
                    key=lambda item: (
                        str(item.get("skill_id", "")),
                        str(item.get("document_sha256", "")),
                    ),
                )
            ),
        )

    def _toolchain_plan_for_review(
        self, plan: ScientificWorkflowPlanV2
    ) -> ScientificToolchainPlanV1 | None:
        """The analysis chain bound to exactly this plan, if one exists.

        The toolchain plan used to live only in session RAM, so the human
        approved a packet that never displayed the analysis nodes and the
        executor had nothing to walk. Attach it to the review only when it
        belongs to this exact workflow resolution -- the same exactness rule
        _resolve_program_workflow applies -- and stay silent otherwise so a
        calculation-only review keeps its historical bytes.
        """

        resolved = self._latest_program_workflows.get(plan.workflow_id)
        if resolved is None:
            return None
        if resolved.scientific_plan is not plan:
            return None
        toolchain = resolved.scientific_toolchain_plan
        if toolchain is None:
            return None
        if (
            toolchain.command_workflow_draft_sha256
            != resolved.draft.draft_sha256
        ):
            return None
        return toolchain

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

        def _carries_node(candidate) -> bool:
            return any(
                item.node_id == node_id for item in candidate.nodes
            ) and (not plan_sha256 or candidate.plan_sha256 == plan_sha256)

        plan = self._current_scientific_plan(predicate=_carries_node)
        if plan is None:
            raise ContractError(
                "bounded execution requires a current scientific workflow "
                "containing the requested node"
            )
        # Node IDs are workflow-local.  Choose the current exact plan first,
        # then require its own materialization.  Falling back to an older plan
        # with the same node name is what turned a missing identity/preparation
        # step in an edge-free diagnostic into an unrelated future-edge error.
        #
        # "The current exact plan" was implemented as `plans[-1]` --
        # insertion order over a digest-keyed dict -- so the comment
        # stated the requirement and the code approximated it. After an
        # amendment restores an earlier plan, the restored key keeps its
        # original position and this, the bounded *execution admission*
        # path, would admit the node against the plan the session
        # abandoned. Resolve through the one owner, and keep the
        # filtered fallback for a host restored without the pointer.

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
            analysis_only = (
                "; this envelope's budget is 0, so it authorises analysis "
                "over registered results only"
                if envelope.max_engine_calls == 0
                else ""
            )
            raise ContractError(
                "scientific workflow exceeds bounded engine-call budget: "
                f"{len(plan.nodes)} nodes for {envelope.max_engine_calls} "
                "calls" + analysis_only
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
            elif is_validated_scan_minimum_geometry_edge(plan, edge):
                selection_rule = "validated_scan_minimum_geometry"
            elif is_validated_orca_ts_hessian_edge(plan, edge):
                selection_rule = "validated_final_orca_ts_hessian"
            elif is_validated_producer_orca_hessian_edge(plan, edge):
                selection_rule = "validated_producer_orca_hessian"
            else:
                raise ContractError(
                    "bounded execution has no exact selection rule for data "
                    f"edge {edge.edge_id!r}; expected optimized geometry, "
                    "an ORCA scan minimum-energy point geometry, an ORCA "
                    "final-TS Hessian for IRC, or an ORCA producer "
                    "Hessian for a TS inhess_filename role"
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
            if edge.selection_rule
            in {
                "validated_optimized_geometry",
                "validated_scan_minimum_geometry",
            }
        )
        edge_by_target = {
            edge.consumer_node_id: edge for edge in geometry_edges
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
            # A node may deliberately reuse a producer geometry on another
            # charge/multiplicity surface, which is how a redox or
            # hydrogen-transfer series is written.  That freedom is exactly
            # where an impossible pair enters, so the state is checked against
            # this molecule's electron count before it reaches a review.  No
            # producer rule changes the atom set, so the input artifact's
            # symbols are the consumer's symbols.
            try:
                _node_symbols = _review_molecule_identity(
                    context.input_artifact
                )["atom_order"]
            except Exception:
                _node_symbols = ()
            if _node_symbols:
                refuse_impossible_electronic_state(
                    _node_symbols,
                    target_charge,
                    target_multiplicity,
                    context=f"node {planned_node.node_id!r}",
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
                    internal_coordinates=context.internal_coordinates,
                    excursion=context.excursion,
                    root_artifact_id=(
                        context.root_artifact.artifact_id
                        if context.root_artifact is not None
                        else ""
                    ),
                    root_artifact_sha256=(
                        context.root_artifact.sha256
                        if context.root_artifact is not None
                        else ""
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
        bind: bool = True,
    ) -> _CommandContext:
        """Resolve current or future context without model-carried receipts."""

        if planned_node.node_id not in data_target_ids:
            _invocation, context = self._plan_invocation_for_node(
                plan=plan, node_id=planned_node.node_id, bind=bind
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
            item
            for item in (node.inputs if node is not None else ())
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
        # One declared Hessian role per auxiliary input, read from the
        # role table rather than spelled out per program and stage.
        valid_hessian_auxiliary = bool(
            len(auxiliary_inputs) == 1
            and len(geometry_inputs) == 1
            and geometry_inputs[0].binding_id == "filename"
            and any(
                planned_node.program == role.consumer_program
                and planned_node.stage in role.consumer_stages
                and auxiliary_inputs[0].binding_id == role.consumer_input_id
                and auxiliary_inputs[0].artifact_class == role.artifact_class
                for role in HESSIAN_CONSUMER_ROLES.values()
            )
        )
        if (
            node is None
            or len(geometry_inputs) != 1
            or (auxiliary_inputs and not valid_hessian_auxiliary)
        ):
            raise ContractError(
                "future bounded node requires one filename/geometry_xyz input; "
                "ORCA IRC may additionally consume one "
                "hess_filename/orca_hessian input, and an ORCA TS one "
                "inhess_filename/orca_hessian starting-Hessian input"
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
        # The producer edge is matched on the draft's own geometry input
        # binding, never on a literal name: the executor hands the
        # produced geometry to whatever binding the draft declared, and
        # a resolver that answered only to one name was program-specific
        # in a program-neutral layer.
        geometry_binding_id = str(geometry_inputs[0].binding_id)
        producer = next(
            (
                edge.source_node_id
                for edge in plan.edges
                if edge.edge_kind == "data"
                and edge.target_node_id == planned_node.node_id
                and edge.artifact_class == "geometry_xyz"
                and edge.consumer_input_id == geometry_binding_id
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
                f"into its {geometry_binding_id!r} input"
            )
        if producer in data_target_ids:
            # A producer that is itself a deferred data target has no
            # compiled invocation yet. The chain neutral opt -> cation
            # opt -> single point is ordinary science, and REACH-1 ino3's
            # review refused exactly its two second-hop nodes with "node
            # has no compiled command invocation" while the frontier had
            # called the same plan approvable (replayed 2026-09-06). The
            # producer's context is resolved as its own would be; the
            # plan is acyclic, so the recursion ends at a compiled root.
            producer_node = next(
                (item for item in plan.nodes if item.node_id == producer),
                None,
            )
            if producer_node is None:
                raise ContractError(
                    f"future node {planned_node.node_id!r} names a producer "
                    f"{producer!r} the plan does not carry"
                )
            producer_context = self._bounded_node_context(
                plan=plan,
                planned_node=producer_node,
                data_target_ids=data_target_ids,
                bind=bind,
            )
        else:
            _producer_invocation, producer_context = (
                self._plan_invocation_for_node(
                    plan=plan, node_id=producer, bind=bind
                )
            )
        return _CommandContext(
            internal_coordinates=getattr(node, "internal_coordinates", None),
            excursion=str(getattr(node, "excursion", "") or ""),
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
        plan_bindings = getattr(self, "_invocation_workflow_plan_sha256s", {})
        for invocation in reversed(tuple(invocations.values())):
            context = contexts[invocation.invocation_sha256]
            if context.proposal.node_id == node_id and (
                not plan_sha256
                or plan_bindings.get(invocation.invocation_sha256)
                == plan_sha256
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
        bind: bool = True,
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
        if bind:
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
        if (
            self.execution_resources is None
            or canonical_data(self.execution_resources)
            != review.execution_resources
        ):
            raise ContractError("execution resources differ from human review")
        if (
            context.project_artifact is None
            or context.project_validation is None
        ):
            raise ContractError(
                "reviewed execution requires a validated project"
            )
        reviewed_settings = json.loads(review.project_settings_text)
        effective_settings = dict(context.project_validation.settings)
        if canonical_data(effective_settings) != canonical_data(
            reviewed_settings
        ):
            raise ContractError(
                "effective project settings differ from human review"
            )
        current_environment = self.environments.get(
            context.engine_binding.environment_receipt_sha256
        )
        if (
            current_environment is None
            or environment_review_summary(current_environment)
            != review.environment_summary
        ):
            raise ContractError(
                "execution environment facts differ from human review"
            )
        molecular = review.molecular_identity
        if context.scientific_identity.charge != molecular.get(
            "charge"
        ) or context.scientific_identity.multiplicity != molecular.get(
            "multiplicity"
        ):
            raise ContractError(
                "molecular electronic state differs from human review"
            )
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
            ),
        }
        if coordinate.get("kind") == "exact-input-artifact":
            expected_input = str(
                coordinate.get("geometry_artifact_sha256", "")
            )
            if context.input_artifact.sha256 != expected_input:
                raise ContractError(
                    "molecular input bytes differ from human review"
                )
            if (
                molecular.get("scientific_identity_sha256")
                != context.scientific_identity.binding_sha256
            ):
                raise ContractError(
                    "molecular identity differs from human review"
                )
            input_role = "molecular-input"
            input_digest = expected_input
        elif coordinate.get("kind") == "validated-producer-output":
            binding = self.workflow_execution_approval.node(node_id)
            input_digest = str(coordinate.get("producer_edge_sha256", ""))
            if (
                binding.input_mode != "producer"
                or binding.producer_edge_sha256 != input_digest
            ):
                raise ContractError(
                    "producer geometry edge differs from human review"
                )
            if self.handoffs.get(node_id) is None:
                raise ContractError(
                    "producer geometry lacks validated handoff"
                )
            input_role = "producer-geometry"
        else:
            raise ContractError("unknown reviewed coordinate identity")
        path_bindings[context.input_artifact.cli_value] = (
            input_role,
            input_digest,
        )
        if self.execution_server:
            server_path = Path(self.execution_server)
            if not server_path.is_file() or server_path.is_symlink():
                raise ContractError("execution server profile is unavailable")
            path_bindings[self.execution_server] = (
                "server-profile",
                review.server_profile_sha256,
            )
        auxiliary_by_name = dict(context.job_artifact_options)
        for binding in invocation.auxiliary_input_bindings:
            artifact = auxiliary_by_name.get(binding.parameter_name)
            if artifact is None or artifact.sha256 != binding.artifact_sha256:
                raise ContractError(
                    "auxiliary input differs from human review"
                )
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
            kind = _output_artifact_kind(program, path)
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
        expected_root_artifact: TrustedArtifactRefV1 | None = None,
        expected_project_artifact: TrustedArtifactRefV1 | None = None,
        output_artifacts: tuple[TrustedArtifactRefV1, ...],
        exit_status: int | None,
        expected_environment_receipt_sha256: str = "",
        capability_environment_receipt: (
            EnvironmentCapabilityReceiptV1 | None
        ) = None,
        pyscf_engine_observation: _PySCFEngineObservation | None = None,
        process_observation: ProcessObservationV1 | None = None,
    ) -> _ExecutionValidationEvaluation:
        findings: list[str] = []
        sensor_inputs: dict[str, Any] = {}
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
                # The runner certifies the artifact's invariants and says
                # validated; the order of a stationary point is the host's
                # program-neutral verdict below, for every program alike.
                # A deferral to a "downstream classification" that no organ
                # performed used to admit an unclassified Hessian here.
                if run_receipt.get("scientifically_validated") is not True:
                    findings.append(
                        "pyscf.run_receipt.scientific_validation_failed"
                    )
                elif runner_findings:
                    findings.append(
                        "pyscf.run_receipt.validation_state_inconsistent"
                    )
                if run_receipt.get("state") != "validated":
                    findings.append("pyscf.run_receipt.state_not_validated")
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
                    from chemsmart.io.native_failure import (
                        summarize_pyscf_native_failure,
                    )
                    from chemsmart.io.pyscf.output import read_pyscf_h5
                    from chemsmart.jobs.pyscf.validation import (
                        FIXED_GEOMETRY_JOBTYPES,
                        validate_pyscf_result,
                    )

                    result_spec, _provenance, result_status, _results = (
                        read_pyscf_h5(
                            _current_artifact_path(
                                result, field_name="PySCF HDF5 result"
                            )
                        )
                    )
                    # The driver's own typed account -- the stage that
                    # raised, or a stage that quietly returned
                    # unconverged -- is summarised before validation so
                    # a validator refusal cannot mask it.
                    native = summarize_pyscf_native_failure(result_status)
                    if native is not None:
                        pyscf_observation = observation.setdefault("pyscf", {})
                        pyscf_observation["native_failure"] = native.as_dict()
                        findings.append(
                            f"pyscf.native_failure.{native.error_class}"
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
                        # The validator holds a fixed-geometry job type
                        # to the geometry it was handed; which job types
                        # those are is the validator's own set, never a
                        # second list here: this list read {"sp", "hess"}
                        # while the validator's had grown to td, and E1 of
                        # PySCF round 2 was typed failed on roots that were
                        # right to the digit (2026-09-13).
                        expected_positions=(
                            expected_positions
                            if jobtype in FIXED_GEOMETRY_JOBTYPES
                            else None
                        ),
                        expected_receipt=expected_receipt,
                    )
                    observation["result_validation"] = canonical_data(
                        result_validation
                    )
                    findings.extend(
                        str(item.rule_id)
                        for item in result_validation["findings"]
                    )
                    result_state = result_validation.get("state")
                    if (
                        result_state != "validated"
                        and not result_validation.get("findings")
                    ):
                        findings.append(
                            "pyscf.result.validation_state_inconsistent"
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
        elif exit_status != 0 and program not in {
            "orca",
            "gaussian",
            "xtb",
        }:
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
                    sensor_inputs.update(
                        _imaginary_mode_sensor_inputs(output, frequencies)
                    )
                    sensor_inputs["basin"] = _basin_sensor_inputs(
                        expected_input_artifact, output, jobtype
                    )
                    if expected_root_artifact is not None and (
                        expected_input_artifact is None
                        or expected_root_artifact.sha256
                        != expected_input_artifact.sha256
                    ):
                        # The walk from the goal's own start, beside the
                        # walk from this node's input: a repair that
                        # restarts from a distorted geometry must not
                        # extinguish the anomaly it answers.
                        sensor_inputs["basin_root"] = {
                            **_basin_sensor_inputs(
                                expected_root_artifact, output, jobtype
                            ),
                            "reference": "goal_root",
                            "reference_artifact_id": (
                                expected_root_artifact.artifact_id
                            ),
                        }
                    if jobtype == "scan":
                        sensor_inputs["scan_profile"] = tuple(
                            {
                                "coordinate": float(point["coordinate"]),
                                "energy": float(point["energy"]),
                            }
                            for point in (output.scan_profile or ())
                        )
                    finite_frequencies = bool(frequencies) and all(
                        math.isfinite(float(value)) for value in frequencies
                    )
                    consequential_imaginary_modes = tuple(
                        float(value)
                        for value in frequencies
                        if math.isfinite(float(value)) and float(value) < -20.0
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
                            "imaginary_frequencies_cm1": tuple(
                                float(value)
                                for value in frequencies
                                if float(value) < 0.0
                            ),
                            "consequential_imaginary_mode_count": (
                                len(consequential_imaginary_modes)
                                if frequencies
                                else None
                            ),
                            **_spin_square_observation(
                                output, output.multiplicity
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
                        jobtype in GEOMETRY_SEARCH_JOBTYPES
                        and output.converged is not True
                    ):
                        findings.append(
                            "orca.result.optimization_not_converged"
                        )
                    # A relaxed scan prints the convergence marker per
                    # step, so an exhausted step turns the tri-state
                    # False even when the run terminates normally --
                    # and such a scan used to validate silently, its
                    # surface partial with nothing saying so. Absence
                    # of any marker stays an absent fact, not a
                    # verdict.
                    if jobtype == "scan" and output.converged is False:
                        findings.append("orca.result.scan_step_not_converged")
                    if output.final_energy is None:
                        findings.append("orca.result.energy_missing")
                    requested_frequency_analysis = any(
                        bool((expected_settings or {}).get(field))
                        for field in ("freq", "numfreq", "vpt2")
                    )
                    if requested_frequency_analysis and not frequencies:
                        findings.append("orca.result.frequencies_missing")
                    elif (
                        requested_frequency_analysis and not finite_frequencies
                    ):
                        findings.append("orca.result.frequencies_invalid")
                    if jobtype == "ts":
                        if not finite_frequencies:
                            if (
                                "orca.result.frequencies_invalid"
                                not in findings
                            ):
                                findings.append(
                                    "orca.result.frequencies_invalid"
                                )
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
            # A failed xTB run used to bail out above with the single
            # generic process finding, so the 68-code receipt audit and
            # the native-failure account below were reachable only when
            # xTB exited 0 -- the whole vocabulary vanished exactly when
            # it was needed. xTB now takes the same path as ORCA: the
            # process finding is recorded and its own audit still runs.
            if exit_status != 0:
                findings.append("execution.process.nonzero_or_unknown")
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
            # The receipt says whether the run satisfied its contract; it
            # does not say what xTB itself complained about.  xTB is one of
            # the three programs this release executes, so a failed run
            # deserves the same account as ORCA gets -- and most of all
            # when the crash left no receipt, which is why this block
            # sits outside the receipt branch.
            # The xTB log itself is registered as kind "xtb_output";
            # "program_output" covers the charges/wbo/xtbrestart/topology
            # sidecars, which never contain the "finished run" sentinel.
            # Summarizing over those called every healthy live run an
            # incomplete_output failure -- observed on the first real
            # xtb execution through agent run, where all three
            # normally terminated single points failed validation.
            xtb_logs = tuple(
                artifact
                for artifact in output_artifacts
                if artifact.kind == "xtb_output"
                and Path(artifact.path).suffix.lower() != ".err"
            )
            # The receipt carries no frequencies; the log does, and the
            # typed layer's own parser reads them. Observed live (W1): a
            # hess node's count came back None from the receipt and the
            # stationary-point rule never reached xTB.
            observation.setdefault("xtb", {})[
                "consequential_imaginary_mode_count"
            ] = consequential_imaginary_mode_count(
                _xtb_log_frequencies(xtb_logs)
            )
            if xtb_logs:
                from chemsmart.io.native_failure import (
                    summarize_xtb_native_failure,
                )

                failure_summary = summarize_xtb_native_failure(
                    _iter_native_diagnostic_lines(xtb_logs),
                    diagnostic_lines=_iter_native_diagnostic_lines(
                        _native_diagnostic_artifacts(output_artifacts)
                    ),
                )
                if failure_summary is not None:
                    xtb_observation = observation.setdefault("xtb", {})
                    xtb_observation["native_failure"] = (
                        _bound_native_failure_summary(
                            failure_summary,
                            artifacts=xtb_logs[:1],
                        )
                    )
                    findings.append(
                        "xtb.native_failure." + failure_summary.error_class
                    )
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
                                "imaginary_frequencies_cm1": tuple(
                                    float(value)
                                    for value in frequencies
                                    if float(value) < 0.0
                                ),
                                "consequential_imaginary_mode_count": (
                                    consequential_imaginary_mode_count(
                                        frequencies
                                    )
                                ),
                                "transition_count": len(transitions),
                                "wavefunction_stability_history": (
                                    stability_history
                                ),
                                **_spin_square_observation(
                                    output, output.multiplicity
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
                            expected_result_jobtype in GEOMETRY_SEARCH_JOBTYPES
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
        # The facts the sensors below consume, read once through the
        # program's own reader; a branch above that already wrote a key
        # keeps its value, and a program whose branch wrote nothing (PySCF,
        # every one of them until 2026-09-12) is fed the same way.
        neutral_block, neutral_inputs = _neutral_sensor_facts(
            program=program,
            jobtype=jobtype,
            multiplicity=multiplicity,
            output_artifacts=output_artifacts,
            expected_input_artifact=expected_input_artifact,
            expected_root_artifact=expected_root_artifact,
        )
        if neutral_block:
            program_block = observation.setdefault(program, {})
            if isinstance(program_block, dict):
                for key, value in neutral_block.items():
                    program_block.setdefault(key, value)
        for key, value in neutral_inputs.items():
            sensor_inputs.setdefault(key, value)
        # One program-neutral verdict on the order of the stationary point,
        # from the frequencies the program itself printed and the jobtype
        # the human approved. ORCA's transition-state check above stays;
        # this one also covers a minimum that landed on a saddle, and the
        # other programs, so the same physics gets the same word.
        observed_order = _observed_imaginary_mode_count(observation, program)
        order_finding = stationary_point_order_finding(jobtype, observed_order)
        if order_finding:
            findings.append(order_finding)
        anomalies: list[dict[str, Any]] = []
        basin = dict(sensor_inputs.pop("basin", {}) or {})
        rmsd = basin.get("heavy_atom_rmsd_angstrom")
        if rmsd is not None and float(rmsd) >= 0.3:
            anomalies.append(
                {
                    "signal_id": "geometry.heavy_atom_rmsd_ge_0.3",
                    **basin,
                }
            )
        if basin.get("bonds_made") or basin.get("bonds_broken"):
            anomalies.append(
                {"signal_id": "geometry.connectivity_changed", **basin}
            )
        root = dict(sensor_inputs.pop("basin_root", {}) or {})
        fired = {item["signal_id"] for item in anomalies}
        root_rmsd = root.get("heavy_atom_rmsd_angstrom")
        if (
            root_rmsd is not None
            and float(root_rmsd) >= 0.3
            and "geometry.heavy_atom_rmsd_ge_0.3" not in fired
        ):
            anomalies.append(
                {"signal_id": "geometry.heavy_atom_rmsd_ge_0.3", **root}
            )
        if (
            root.get("bonds_made") or root.get("bonds_broken")
        ) and "geometry.connectivity_changed" not in fired:
            anomalies.append(
                {"signal_id": "geometry.connectivity_changed", **root}
            )
        scan_boundary = _scan_boundary_sensor(
            sensor_inputs.pop("scan_profile", ()) or ()
        )
        if scan_boundary is not None:
            anomalies.append(scan_boundary)
        deviation = _observed_spin_deviation(observation, program)
        if deviation is not None and abs(deviation) >= 0.2:
            # ⟨S²⟩ is an observation, never a gate; a deviation this size
            # is the surprise a scientist weighs (a wrong state, a
            # multireference character), recorded with its numbers.
            anomalies.append(
                {
                    "signal_id": "spin.s2_deviation_ge_0.2",
                    "spin_square_deviation": float(f"{deviation:.4f}"),
                    "bound_multiplicity": multiplicity,
                }
            )
        if order_finding:
            # The verdict says the run failed its promise; the observation
            # says what the structure is. Both are true, and the second
            # used to have nowhere to live.
            anomalies.append(
                {
                    "signal_id": "stationary_point.unexpected_order",
                    "expected_imaginary_modes": expected_imaginary_mode_count(
                        jobtype
                    ),
                    "observed_imaginary_modes": observed_order,
                    **sensor_inputs,
                }
            )
        stationarity_gradient = sensor_inputs.pop(
            "stationarity_gradient", None
        )
        if (
            stationarity_gradient is not None
            and float(stationarity_gradient)
            > HESS_STATIONARITY_GRADIENT_EH_PER_BOHR
        ):
            anomalies.append(_gradient_anomaly(float(stationarity_gradient)))
        soft_mode = _observed_soft_imaginary_mode(
            observation, program, jobtype=jobtype
        )
        if soft_mode is not None:
            anomalies.append(
                {
                    "signal_id": "stationary_point.imaginary_mode_lt_50",
                    "imaginary_mode_cm1": float(f"{soft_mode:.2f}"),
                    "noise_convention_cm1": 20.0,
                    "soft_band_cm1": SOFT_IMAGINARY_MODE_BAND_CM1,
                    "expected_imaginary_modes": 1,
                }
            )
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
            anomalies=tuple(anomalies),
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

    def unapproved_workflow_summary(self) -> dict[str, Any] | None:
        """Say whether the newest plan could be approved, and what stops it.

        A session ends on its first message without tool calls, and at that
        moment nothing re-checked whether the workflow it built could actually
        be approved.  The readiness was computed once, on the receipt of the
        plan itself, before any node had been prepared -- so a session could
        do the work, never see the result of it, and stop believing it had
        delivered an approvable plan while the host knew otherwise.

        This reports the state, it does not judge the chemistry: the node ids
        are the host's own bookkeeping, already computed by
        ``_approval_readiness``.  Returns ``None`` when there is no plan to
        speak about, which keeps an analysis-only session unaffected.
        """

        if not self.scientific_plans:
            return None
        plan = next(reversed(self.scientific_plans.values()))
        readiness = self._approval_readiness(plan)
        if readiness.get("approvable"):
            return None
        return {
            "workflow_id": getattr(plan, "workflow_id", ""),
            "blocking_node_ids": tuple(readiness.get("blocking_node_ids", ())),
            "deferred_node_ids": tuple(readiness.get("deferred_node_ids", ())),
            "workflow_blocked_reason": readiness.get(
                "workflow_blocked_reason", ""
            ),
        }

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

    def _inspect_run_outcome(self, turn_id: str, values: dict) -> Any:
        """Return how a recorded run ended, as typed terminal states.

        The durable stream has always carried the facts; a session could
        read none of them -- a failed run's account lived in terminal
        text a human had to retype into the next task. This tool serves
        the same derivation the goal loop's wake context uses, so what a
        session can ask for and what the loop acts on cannot drift.

        Without a ``run`` reference it lists the runs the workspace
        records; with one it returns that run's full outcome. The host
        resolves run references inside the user workspace's own
        ``.chemsmart-agent`` records only -- never a caller-supplied
        path, and never the session's private preview root, where no
        run is ever recorded.
        """

        if self.run_evidence_root is None:
            raise ContractError(
                "run inspection requires the workspace whose "
                ".chemsmart-agent directory records runs; this host "
                "was constructed without one"
            )
        from chemsmart.agent.terminal_states import (
            derive_run_outcome,
            read_run_events,
        )

        records_root = self.run_evidence_root / ".chemsmart-agent"
        streams: dict[str, Path] = {}
        for pattern in (
            "replays/*/run/events.jsonl",
            "executions/*/events.jsonl",
            "goals/*/runs/*/events.jsonl",
            "runs/*/events.jsonl",
        ):
            for path in sorted(records_root.glob(pattern)):
                if path.is_symlink() or not path.is_file():
                    continue
                reference = str(path.parent.relative_to(records_root))
                streams[reference] = path

        requested = str(values.get("run", "") or "").strip()
        if not requested:
            listing = []
            for reference, path in streams.items():
                try:
                    outcome = derive_run_outcome(read_run_events(path))
                except (ContractError, ValueError, TypeError):
                    listing.append({"run": reference, "readable": False})
                    continue
                listing.append(
                    {
                        "run": reference,
                        "readable": True,
                        "workflow_id": outcome.workflow_id,
                        "workflow_state": outcome.workflow_state,
                        "engine_calls_consumed": (
                            outcome.engine_calls_consumed
                        ),
                        "node_states": {
                            node.node_id: node.state for node in outcome.nodes
                        },
                    }
                )
            return {
                "workspace_runs": tuple(listing),
                "note": (
                    "name one run to read its full typed outcome, "
                    "including per-node terminal facts and the evidence "
                    "digests a revision cites"
                ),
            }
        path = streams.get(requested)
        resolved_from = ""
        if path is None:
            # The wake's previous_run_outcome names its run_id; the
            # reference the stream is filed under ends in it. A session
            # asked by run_id and was refused twice (NOVEL-3 ino3).
            by_last_segment = [
                reference
                for reference in streams
                if reference.rsplit("/", 1)[-1] == requested
            ]
            if len(by_last_segment) == 1:
                resolved_from = requested
                requested = by_last_segment[0]
                path = streams[requested]
            elif by_last_segment:
                raise RoutedContractError(
                    gate="inspect.run_reference_names_one_recorded_run",
                    invariant=(
                        "a run outcome is read from exactly one recorded "
                        "stream."
                    ),
                    diagnosis=(
                        f"run_id {requested!r} names "
                        f"{len(by_last_segment)} recorded runs: "
                        f"{sorted(by_last_segment)}."
                    ),
                    route="name one by its full reference.",
                )
            else:
                raise RoutedContractError(
                    gate="inspect.run_reference_names_one_recorded_run",
                    invariant=(
                        "a run outcome is read from exactly one recorded "
                        "stream."
                    ),
                    diagnosis=(
                        f"this workspace records no run {requested!r}; "
                        f"recorded runs: {sorted(streams)}."
                    ),
                    route=(
                        "name one of the recorded references, or a run_id "
                        "that names exactly one of them; inspect_run with no "
                        "arguments lists them."
                    ),
                )
        stream_bytes = path.read_bytes()
        outcome = derive_run_outcome(read_run_events(path))
        # A named read is a typed act: the durable event binds the run
        # reference to the exact stream bytes served, so the revision
        # gate can verify which run's outcome entered this session --
        # a bare listing proves nothing and records nothing.
        stream_sha256 = hashlib.sha256(stream_bytes).hexdigest()
        self.event_store.append(
            turn_id=turn_id,
            kind=EventKind.RUN_OUTCOME_INSPECTED.value,
            payload={
                "run": requested,
                "workflow_id": outcome.workflow_id,
                "workflow_state": outcome.workflow_state,
                "engine_calls_consumed": outcome.engine_calls_consumed,
                "stream_sha256": stream_sha256,
            },
            idempotency_key=(
                f"run-outcome-inspected:{requested}:{stream_sha256}"
            ),
        )
        record = outcome.public_record()
        record["run"] = requested
        if resolved_from:
            record["resolved_from"] = resolved_from
        record["nodes"] = tuple(
            {
                **node.public_record(),
                "evidence_event_hashes": node.evidence_event_hashes,
            }
            for node in outcome.nodes
        )
        return record

    def _inspect_result_selectors(self, turn_id: str, values: dict) -> Any:
        """Return what one completed result resolves, by probing its parser.

        The reader has always been able to answer this, and until now the
        answer only ever reached a session appended to a refusal -- one
        selector at a time, after an extraction had already been planned and
        run.  Probing costs one attempted read per accessor on an already
        opened result and no re-parse.

        It states what the artifact holds and chooses nothing: which of these
        the scientific question needs is not the host's judgement, and a
        selector missing here is a property of the calculation that produced
        this result rather than a gap in the parser.
        """

        artifact = self._artifact(values["artifact_id"])
        program = str(values["program"]).strip().lower()
        reader = reader_for(program)
        if reader is None:
            raise ContractError(
                "no typed result reader is registered for "
                f"{program!r}; registered: "
                f"{registered_reader_programs()}"
            )
        if artifact.kind != reader.artifact_kind:
            raise ContractError(
                f"{program} selector inspection requires a bound "
                f"{reader.artifact_kind} artifact, not {artifact.kind!r}"
            )
        output = reader.open_output(Path(artifact.path))
        available = reader.available_selectors(output)
        jobtype = str(getattr(output, "jobtype", "") or "").casefold()
        declared = reader.selectors_for_jobtype(jobtype)
        requestable = tuple(
            selector
            for selector in available
            if declared is None or selector in declared
        )
        return {
            "artifact_id": artifact.artifact_id,
            "program": program,
            "parser_id": reader.parser_id,
            "jobtype": jobtype,
            "available_selectors": available,
            # A selector this artifact resolves is still refused unless the
            # job type declares it, because a declaration is a claim about
            # what the value means for that job type.
            "requestable_selectors": requestable,
            # Which molecular state each value belongs to -- supplied,
            # reached, the thermochemistry reference -- so a session asks
            # for the role it needs instead of learning it from a receipt.
            "structural_states": {
                selector: reader.structural_state(selector)
                for selector in requestable
            },
            # Whose density or method each value belongs to -- the
            # reference, an excited root, the correlated method -- resolved
            # against this artifact, because a geometry identity says
            # nothing about whose density a dipole is: an excited-root
            # optimisation carries the reference's dipole beside the root's
            # energy at one structure.
            "electronic_provenance": {
                selector: reader.electronic_provenance_for_output(
                    output, selector
                )
                for selector in requestable
            },
            # The level this artifact computed at, from its own record --
            # method, basis, frozen core, the response and the followed
            # root -- so a session names it beside the number it delivers
            # instead of inferring it from a project it may not hold.
            "level": reader.level_for_output(output),
        }

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
        record = canonical_data(
            canonical_extraction_receipt_body(
                schema_version=receipt.schema_version,
                artifact_id=receipt.artifact_id,
                artifact_sha256=receipt.artifact_sha256,
                program=receipt.program,
                parser_id=receipt.parser_id,
                quantities=receipt.quantities,
                status=receipt.status,
                absent=receipt.absent,
                derived_adjacency=receipt.derived_adjacency,
                electronic_provenance=receipt.electronic_provenance,
            )
        )
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
        try:
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
                reaction_coordinate_mode=self._reaction_coordinate_mode(
                    values
                ),
                frequency_scale_factor=float(
                    values.get("frequency_scale_factor", 1.0)
                ),
            )
        except ValueError as exc:
            # The thermochemistry kernel refuses scientifically meaningless
            # requests with a bare ValueError -- observed live: an acetate
            # optimization that converged onto its methyl-torsion saddle
            # (-64.7 cm^-1), where a pure-RRHO Gibbs correction does not
            # exist. Typing the refusal here keeps the delivery envelope's
            # rule honest: typed errors settle the node and reach the
            # partial report as findings. The other half of this comment
            # used to read "while a bare exception stays what it should
            # be -- a defect that crashes", and that half is withdrawn.
            # An escaping TypeError killed a goal unsettled and lost
            # every sibling receipt with it, while its traceback went to
            # a log nothing reads. The executor now settles an
            # unexpected exception as a host defect, which is both
            # louder and recoverable (executor._run_analysis_phase).
            raise ContractError(str(exc)) from exc
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
        low_modes = self._low_frequency_mode_observation(
            artifact=artifact,
            temperature_k=float(values["temperature_k"]),
        )
        self._emit(
            turn_id,
            EventKind.THERMOCHEMISTRY_DERIVED,
            receipt.receipt_sha256,
            status=receipt.status,
            artifact_sha256=receipt.artifact_sha256,
            temperature_k=receipt.temperature_k,
            pressure_atm=receipt.pressure_atm,
            record=record,
            **({"observations": (low_modes,)} if low_modes else {}),
        )
        if low_modes:
            self._reply_observations = (low_modes,)
        return receipt

    def _low_frequency_mode_observation(
        self, *, artifact: TrustedArtifactRefV1, temperature_k: float
    ) -> dict[str, Any] | None:
        """Name the modes below 50 cm-1 and the RRHO entropy they carry.

        Read from the result's own printed frequencies with the same
        engine the derivation used; an observation beside the receipt,
        never inside it and never a verdict.
        """

        from chemsmart.analysis.result_quantities import (
            low_frequency_mode_entropy,
        )
        from chemsmart.analysis.thermochemistry import Thermochemistry

        try:
            engine = Thermochemistry(
                filename=str(artifact.path), temperature=temperature_k
            )
            frequencies = tuple(engine.real_frequencies or ())
        except Exception:  # noqa: BLE001 -- an observation never fails a call
            return None
        summary = low_frequency_mode_entropy(
            frequencies, temperature_k=temperature_k
        )
        if not summary["low_modes_cm1"]:
            return None
        return {
            "kind": "low_frequency_modes_under_rrho",
            **summary,
            "meaning": (
                f"{len(summary['low_modes_cm1'])} real mode(s) below "
                f"{summary['threshold_cm1']:.0f} cm-1 contribute "
                f"{summary['entropy_term_kj_per_mol']:.2f} kJ/mol of T*S "
                "under the harmonic oscillator, which describes such modes "
                "worst; two enantiomeric rotamers have differed by 0.4 "
                "kJ/mol in Gibbs energy on this term alone. Compare it "
                "with the difference you deliver; entropy_method and "
                "entropy_cutoff_cm1 select a quasi-harmonic treatment"
            ),
        }

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
            expression_node_from_plan(item) for item in values["nodes"]
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
        geometry_observations = self._geometry_operation_observations(
            values, nodes
        )
        self._emit(
            turn_id,
            EventKind.QUANTITY_EXPRESSION_EVALUATED,
            receipt.receipt_sha256,
            status=receipt.status,
            output_ids=tuple(item.quantity_id for item in receipt.outputs),
            semantic_signature_sha256=receipt.semantic_signature_sha256,
            record=record,
            **(
                {"geometry_observations": geometry_observations}
                if geometry_observations
                else {}
            ),
        )
        if geometry_observations:
            self._reply_observations = geometry_observations
        return receipt

    def _geometry_operation_observations(
        self, values: Mapping[str, Any], nodes: Sequence[Any]
    ) -> tuple[dict[str, Any], ...]:
        """Say when an internal coordinate was measured off a non-bond.

        distance, angle and dihedral promise "bonded order a-b-c-d" and
        nothing checked it. A woken session measured every F-C-C-S
        torsion to refute its own rotamer labels -- exactly the right
        move -- with the sulfone's indices one atom off, so it reported
        the well at 78-82 deg (F-C-C-O) where F-C-C-S sits at 62-71
        deg, and the number was relayed as a finding (NOVEL-2 po2,
        2026-09-04). The adjacency the same receipt carried would have
        caught it. An observation, never a refusal: a torsion over a
        non-bonded chain is sometimes exactly what a scientist wants.
        """

        by_id = {node.node_id: node for node in nodes}
        source_by_input: dict[str, str] = {}
        for item in values.get("inputs") or ():
            source_by_input[str(item.get("input_id") or "")] = str(
                item.get("receipt_sha256") or ""
            )

        def _atom_of(node_id: str) -> tuple[str, int] | None:
            node = by_id.get(node_id)
            if node is None or node.operation != "ref":
                return None
            if len(node.indices) != 1:
                return None
            receipt = source_by_input.get(node.reference, "")
            return (receipt, int(node.indices[0])) if receipt else None

        observations: list[dict[str, Any]] = []
        for node in nodes:
            if node.operation not in {"distance", "angle", "dihedral"}:
                continue
            atoms = [_atom_of(input_id) for input_id in node.input_ids]
            if any(atom is None for atom in atoms):
                continue
            receipts = {receipt for receipt, _index in atoms}
            if len(receipts) != 1:
                continue
            extraction = self.quantity_extractions.get(next(iter(receipts)))
            adjacency = getattr(extraction, "derived_adjacency", None)
            if not isinstance(adjacency, Mapping):
                continue
            bonds = {
                tuple(sorted((int(a), int(b))))
                for a, b in adjacency.get("bond_atom_pairs") or ()
            }
            if not bonds:
                continue
            indices = [index for _receipt, index in atoms]
            unbonded = [
                [first, second]
                for first, second in zip(indices, indices[1:])
                if tuple(sorted((first, second))) not in bonds
            ]
            if not unbonded:
                continue
            symbols = ()
            for quantity in getattr(extraction, "quantities", ()):
                if getattr(quantity, "quantity_id", "") == "symbols" or (
                    "symbols" in str(getattr(quantity, "evidence_ref", ""))
                ):
                    symbols = tuple(
                        str(item) for item in (quantity.value or ())
                    )
                    break
            observations.append(
                {
                    "kind": "geometry_operation_over_non_bond",
                    "node_id": node.node_id,
                    "operation": node.operation,
                    "indices": list(indices),
                    "elements": (
                        [
                            symbols[index] if index < len(symbols) else "?"
                            for index in indices
                        ]
                        if symbols
                        else []
                    ),
                    "unbonded_pairs": unbonded,
                    "meaning": (
                        f"{node.operation} promises bonded order and the "
                        "receipt's own perceived bonds do not join "
                        + ", ".join(f"{a}-{b}" for a, b in unbonded)
                        + "; check the indices against the atom order "
                        "before reading this as an internal coordinate"
                    ),
                }
            )
        return tuple(observations)

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
                # Absence and ambiguity are different faults with different
                # repairs, and one message for both names neither.  Say which
                # it is, and for an absent quantity name what the receipt does
                # carry, so the reason does not have to be guessed from the
                # shape of the failure.
                carried = ", ".join(
                    sorted(
                        {
                            quantity.quantity_id
                            for quantity in getattr(receipt, collection)
                        }
                    )
                )
                if not matches:
                    raise ContractError(
                        f"{operation} input {quantity_id!r} is not in this "
                        f"{source_kind} receipt, which carries: {carried}"
                    )
                raise ContractError(
                    f"{operation} input {quantity_id!r} appears "
                    f"{len(matches)} times in this {source_kind} receipt, so "
                    "the receipt alone cannot say which one is meant"
                )
            return source_kind, receipt, matches[0]
        raise ContractError(f"{operation} references an unknown receipt")

    def _evaluate_scientific_validation(
        self, turn_id: str, values: dict
    ) -> ScientificValidationReceiptV1:
        """Evaluate one exact validation node without model-authored rules."""

        workflow_id = str(values["workflow_id"])
        node_id = str(values["node_id"])
        approved = self.approved_scientific_toolchain_plan
        if (
            approved is not None
            and workflow_id == approved.workflow_id
            and workflow_id not in self._latest_program_workflows
        ):
            # Executor host: no planning session ran in this process, so the
            # sealed rules come from the approved bundle's digest-bound plan
            # and the task from the single approved task spec.
            plan = approved
            task_spec_sha256 = self._resolve_task_spec_reference(
                {}, "task_spec_sha256"
            )
        else:
            resolved = self._resolve_program_workflow(workflow_id)
            plan = resolved.scientific_toolchain_plan
            task_spec_sha256 = resolved.draft.task_spec_id
        if plan is None:
            raise ContractError(
                "workflow has no scientific analysis toolchain"
            )
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
            task_spec_sha256=task_spec_sha256,
        )
        input_intents = {
            item.input_id: item
            for item in node.inputs
            if isinstance(item, AnalysisInputIntentV1)
        }
        supplied = tuple(values["inputs"])
        supplied_ids = tuple(
            sorted(str(item["input_id"]) for item in supplied)
        )
        if len(supplied_ids) != len(
            set(supplied_ids)
        ) or supplied_ids != tuple(sorted(input_intents)):
            raise ContractError(
                "scientific validation requires exactly its planned inputs"
            )
        producer_nodes = {item.node_id: item for item in plan.analysis_nodes}
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
                    (
                        "outputs"
                        if source_kind
                        in {"quantity_expression", "scientific_validation"}
                        else "quantities"
                    ),
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
                "source_receipt_sha256s": (receipt.source_receipt_sha256s),
                "all_rules_passed": receipt.all_rules_passed,
                "status": receipt.status,
                "record": record,
            },
            idempotency_key=(
                "scientific-validation:" + receipt.receipt_sha256
            ),
        )
        return receipt

    def _expression_is_model_authored(
        self,
        receipt_sha256: str,
        quantity_id: str = "",
        seen: set[str] | None = None,
    ) -> bool:
        """Whether an expression's value rests on a literal the model wrote.

        An expression receipt is host-minted arithmetic, but its leaves
        need not be host-owned: the vocabulary separates a ``literal``
        the session supplies from a ``constant`` the registry owns, and
        a second expression over a literal-rooted first inherits none of
        that provenance. A value the model supplied stays the model's
        however many receipts sit above it.

        The answer is the expression layer's own, read off the receipt.
        ``model_authored_constants`` already names every number a node
        contributed that no measurement produced -- a ``literal``, a
        ``power`` exponent, a ``scale_factor``, an extrapolation
        exponent -- accumulated to each output through its own sources.
        A second scan here read only ``literal`` operations, so a
        registered constant multiplied by a model-authored
        ``scale_factor`` of 1e-5 was called host-owned evidence and
        discharged a tolerance, while the receipt beside it named that
        factor as the model's. One part of ChemSmart must not call a
        contribution model-authored while another treats its result as
        wholly the host's. Reading the receipt also survives a restart,
        which the request dictionary did not.
        """

        seen = set() if seen is None else seen
        key = f"{receipt_sha256}:{quantity_id}"
        if key in seen:
            return False
        seen.add(key)
        receipt = self.quantity_expression_receipts.get(receipt_sha256)
        if receipt is None:
            return False
        rows = tuple(getattr(receipt, "output_dependencies", ()) or ())
        if quantity_id:
            # The citation names an output; the receipt's other outputs
            # are other numbers. Reading every row made one expression's
            # provenance the whole receipt's: a session that put the
            # task's own three oxidant potentials in the same expression
            # as its uncertainty budget had a clean, fully derived
            # functional term refused because those literals shared its
            # receipt, and the remedy it learned was to partition
            # receipts rather than to derive anything. Forty lines below,
            # the two-source guard already scopes by output_id; these
            # are two readers of one citation and only one was right.
            rows = tuple(
                row
                for row in rows
                if str(getattr(row, "output_id", "")) == quantity_id
            )
        for dependency in rows:
            if getattr(dependency, "model_authored_constants", ()):
                return True
            for source in (
                getattr(dependency, "source_receipt_sha256s", ()) or ()
            ):
                # Across receipts the constants are not accumulated, so
                # the upstream walk stays receipt-scoped: conservative
                # where the host cannot see which output was read.
                if source in self.quantity_expression_receipts:
                    if self._expression_is_model_authored(source, "", seen):
                        return True
        return False

    def _cited_uncertainty_quantity(
        self, receipt_sha256: str, quantity_id: str
    ) -> Any | None:
        """The typed quantity a citation names inside a receipt."""

        registries = (
            (self.quantity_extractions, "quantities"),
            (self.thermochemistry_receipts, "quantities"),
            (self.quantity_expression_receipts, "outputs"),
            (self.scientific_validation_receipts, "outputs"),
        )
        for registry, attribute in registries:
            receipt = registry.get(receipt_sha256)
            if receipt is None:
                continue
            for quantity in getattr(receipt, attribute, ()) or ():
                if str(getattr(quantity, "quantity_id", "")) == quantity_id:
                    return quantity
        return None

    def _uncertainty_evidence(
        self, reference: str, quantity_id: str = ""
    ) -> tuple[bool, str]:
        """Whether a stated uncertainty rests on evidence the host owns.

        The reference names a receipt this host minted or a constant the
        registry owns. Nothing here judges whether the number is a good
        estimate of anything -- that is chemistry and it stays the
        session's -- only whether the host can say where it came from.
        """

        ref = str(reference or "").strip()
        if not ref:
            return False, "no reference was given"
        try:
            from chemsmart.analysis.literature_constants import (
                literature_constant,
            )

            literature_constant(ref)
            return True, f"registered constant {ref!r}"
        except Exception:
            pass
        registries = (
            ("quantity_extraction", self.quantity_extractions),
            ("thermochemistry", self.thermochemistry_receipts),
            ("quantity_expression", self.quantity_expression_receipts),
            ("scientific_validation", self.scientific_validation_receipts),
        )
        for kind, registry in registries:
            if ref in registry:
                # A literal-rooted chain used to be refused here. The
                # refusal read spelling rather than value: the operation
                # vocabulary is rational-complete over any non-zero
                # quantity, so `divide(x, x)` summed to any integer
                # reaches the same number a `literal` would, carrying no
                # model-authored constant at all. Thirteen firings over
                # six windows produced operand rewriting and receipt
                # partitioning and not one derivation, while the same
                # refusal rejected coefficients a definition fixes --
                # the electron count of a one-electron couple, the half
                # that makes a half-range. Provenance is what the host
                # owns, so the authorship is recorded and reported by
                # `_uncertainty_observations`; whether an equivalently
                # spelled coefficient is worth less is a scientific
                # judgement and it is not the host's to make (owner
                # ruling, 2026-09-10).
                return True, f"{kind} receipt"
        return False, "no receipt or registered constant of that name"

    def _expression_source_closure(
        self,
        receipt_sha256: str,
        quantity_id: str = "",
        seen: set[str] | None = None,
    ) -> set[str]:
        """Every receipt a cited output actually descends from.

        Counting the immediate edge called a spread over eight receipts
        a spread over one: the claim cited an output of a receipt whose
        only job was restating hartree per electron as volts, so its
        single dependency was the receipt that had done the comparing.
        An algebraically identity transformation must not change what
        the host says the evidence is (SUFFICIENCY-5, 2026-09-10).
        """

        seen = set() if seen is None else seen
        key = f"{receipt_sha256}:{quantity_id}"
        if key in seen:
            return set()
        seen.add(key)
        receipt = self.quantity_expression_receipts.get(receipt_sha256)
        if receipt is None:
            return set()
        rows = tuple(getattr(receipt, "output_dependencies", ()) or ())
        if quantity_id:
            rows = tuple(
                row
                for row in rows
                if str(getattr(row, "output_id", "")) == quantity_id
            )
        closure: set[str] = set()
        for dependency in rows:
            for source in (
                getattr(dependency, "source_receipt_sha256s", ()) or ()
            ):
                if source in self.quantity_expression_receipts:
                    upstream = self._expression_source_closure(
                        source, "", seen
                    )
                    # A restatement contributes its own parents; a real
                    # comparison contributes itself as well.
                    closure |= upstream or {source}
                else:
                    closure.add(source)
        return closure

    def _authored_constants_in_chain(
        self,
        receipt_sha256: str,
        quantity_id: str = "",
        seen: set[str] | None = None,
    ) -> tuple[tuple[str, str, str], ...]:
        """The coefficients the receipt itself names, with role and value.

        The host resolved these and then reduced them to one bare id, so
        a reader could not tell the one half that *defines* a half-range
        from a 0.1 V solvation budget somebody typed. Both are the
        session's; they are not the same statement, and the receipt has
        the difference in hand.
        """

        seen = set() if seen is None else seen
        key = f"{receipt_sha256}:{quantity_id}"
        if key in seen:
            return ()
        seen.add(key)
        receipt = self.quantity_expression_receipts.get(receipt_sha256)
        if receipt is None:
            return ()
        rows = tuple(getattr(receipt, "output_dependencies", ()) or ())
        if quantity_id:
            rows = tuple(
                row
                for row in rows
                if str(getattr(row, "output_id", "")) == quantity_id
            )
        found: list[tuple[str, str, str]] = []
        for dependency in rows:
            for constant in (
                getattr(dependency, "model_authored_constants", ()) or ()
            ):
                node = str(getattr(constant, "node_id", "") or "")
                role = str(getattr(constant, "role", "") or "")
                value = str(getattr(constant, "value", "") or "")
                if not (node or role or value):
                    # A receipt that names a coefficient without the
                    # three fields still names one; say so rather than
                    # printing an empty pair.
                    role = str(constant)
                found.append((node, role, value))
            for source in (
                getattr(dependency, "source_receipt_sha256s", ()) or ()
            ):
                if source in self.quantity_expression_receipts:
                    # Receipt-scoped upstream, because the edge records
                    # which receipts were read and never which of their
                    # outputs. That imprecision is reported rather than
                    # hidden: see `uncertainty.upstream_scope`.
                    found.extend(
                        self._authored_constants_in_chain(source, "", seen)
                    )
        ordered = sorted(set(found))
        return tuple(ordered)

    def _reaction_coordinate_mode(self, values: Mapping[str, Any]) -> int:
        """The mode a session names as the reaction coordinate.

        Admitted only over a stationary-point characterisation this host
        minted for the same bytes, so the order the treatment assumes is
        one the host checked against the program's own printed
        frequencies rather than one the session asserted. That receipt
        became citable in the same round; this is the reader that makes
        the charter's "a claim standing on the characterised result says
        so" mean something.
        """

        mode = int(values.get("reaction_coordinate_mode", 0) or 0)
        if mode <= 0:
            return 0
        artifact = self._artifact(values["artifact_id"])
        characterised = tuple(
            receipt
            for receipt in self.stationary_point_characterisations.values()
            if str(getattr(receipt, "result_artifact_sha256", ""))
            == str(getattr(artifact, "sha256", ""))
        )
        if not characterised:
            raise ContractError(
                "naming a reaction coordinate says this structure is a "
                "saddle, and the host will not take that on your word: "
                "characterise_stationary_point on this result first, "
                "which checks the order you state against the "
                "program's own printed frequencies. Then name the mode"
            )
        orders = {
            int(getattr(receipt, "order_claimed", 0))
            for receipt in characterised
        }
        if orders and mode > max(orders):
            raise ContractError(
                f"mode {mode} is past the order the characterisation "
                f"established for this result ({sorted(orders)}): name a "
                "mode the structure actually has"
            )
        return mode

    def _component_observations(
        self, component: Mapping[str, Any], display_unit: str = ""
    ) -> dict[str, Any]:
        """What the host saw about one component's cited magnitude.

        Including what the cited quantity actually reads. The claim's
        own ``measured`` magnitude is checked against its citation, and
        a component's is deliberately not -- a term restated at two
        sigma from a one sigma receipt is ordinary and grading it would
        be the host judging a statistical convention. But the reader
        could not see the difference either: ino3-r17 cycle 3
        (2026-09-10) stated a 0.214 V component whose cited quantity
        holds -2.103 V, the SVP potential itself rather than any
        spread, and nothing on the row said so. So the host reports
        what it read and rules on none of it, which is the trade this
        whole surface was rebuilt on (owner ruling, 2026-09-10).
        """

        reference = str(component.get("reference") or "")
        if not reference:
            return {}
        receipt, _, quantity_id = reference.partition(":")
        observed = list(self._uncertainty_observations(receipt, quantity_id))
        if quantity_id and display_unit:
            cited = self._cited_uncertainty_quantity(receipt, quantity_id)
            if cited is not None:
                try:
                    reads = convert_normalized_value(
                        cited.value, cited.dimension, display_unit
                    )
                except Exception:
                    reads = None
                if isinstance(reads, (int, float)):
                    observed.append(
                        f"uncertainty.cited_quantity_reads[{quantity_id}]"
                        f"={reads:.6g} {display_unit}"
                    )
        return {"observations": observed} if observed else {}

    def _uncertainty_observations(
        self, reference: str, quantity_id: str = ""
    ) -> tuple[str, ...]:
        """What the host saw about a cited magnitude, and did not judge.

        Three of these were refusals until 2026-09-10, and each was
        removed for the same reason: it was computable and
        scientifically arbitrary. A spread of exactly zero is a real
        observation -- three treatments agreeing to printed precision,
        an equality symmetry enforces -- and refusing it is the host
        asserting that a measurement cannot be zero. A variance over
        many samples inside one receipt compares plenty. And an
        equivalent coefficient does not become worth less because
        operators spelled it. None of the three stopped what it was
        aimed at, because a substitution walks past all of them.

        What replaced them was too coarse to pay for them, and the two
        audits of that trade agreed on why: the ids named an
        implementation feature rather than what the host saw. So each
        now carries the numbers behind it -- which coefficients, what
        role, what value, how many receipts the value really descends
        from, how much of it the toolkit's own vocabulary produced --
        and the report says where its own reading is imprecise.

        The ids name a measurement and never a verdict, as an anomaly id
        does. Nothing here is refused.
        """

        ref = str(reference or "").strip()
        if not ref:
            return ()
        observations: list[str] = []
        authored = self._authored_constants_in_chain(ref, quantity_id)
        if authored:
            named = "; ".join(
                (
                    f"{role or 'constant'}={value}@{node}"
                    if node
                    else f"{role or 'constant'}={value}"
                )
                for node, role, value in authored
            )
            observations.append(
                f"uncertainty.model_authored_constants[{named}]"
            )
        composed = self.quantity_expression_receipts.get(ref)
        if composed is not None and quantity_id:
            closure = self._expression_source_closure(ref, quantity_id)
            observations.append(f"uncertainty.source_receipts={len(closure)}")
            for dependency in (
                getattr(composed, "output_dependencies", ()) or ()
            ):
                if str(getattr(dependency, "output_id", "")) != quantity_id:
                    continue
                nodes = int(
                    getattr(dependency, "arithmetic_node_count", 0) or 0
                )
                conventions = tuple(
                    getattr(dependency, "convention_operations", ()) or ()
                )
                # The analysis layer already answers "did the toolkit's
                # vocabulary produce this number, or did the model
                # assemble it" and nothing outside that layer read the
                # answer. It is a better-shaped question than any of the
                # three ids above and it was already being computed.
                observations.append(
                    "uncertainty.arithmetic_nodes="
                    f"{nodes} conventions="
                    + (",".join(conventions) if conventions else "none")
                )
                if any(
                    source in self.quantity_expression_receipts
                    for source in (
                        getattr(dependency, "source_receipt_sha256s", ()) or ()
                    )
                ):
                    observations.append("uncertainty.upstream_scope=receipt")
                break
        if quantity_id:
            cited = self._cited_uncertainty_quantity(ref, quantity_id)
            value = getattr(cited, "value", None)
            if isinstance(value, (int, float)) and not isinstance(value, bool):
                if float(value) == 0.0:
                    observations.append("uncertainty.zero_magnitude")
        return tuple(observations)

    def _declarations_for_claim(
        self, claim_id: str, quantity_id: str
    ) -> tuple[Mapping[str, Any], ...]:
        """Every declaration this claim delivers, joined as the gate joins.

        The gate credits a claim under both its ``claim_id`` and its
        host-minted ``quantity_id``, and this join returned only the
        first match. So a claim on a tolerance-bearing quantity, made
        under a ``claim_id`` that names a second, tolerance-free
        declaration, delivered *both* ids and was assessed against
        neither: the completion certified, the record carried no
        sufficiency row, and the requirement disappeared from the
        question rather than standing open. Every obligation a claim
        discharges is assessed, so a delivered tolerance-bearing
        declaration with nothing said about it is ``unstated``.
        """

        found: list[Mapping[str, Any]] = []
        seen: set[str] = set()
        for key in (claim_id, quantity_id):
            if not key or key in seen:
                continue
            declaration = self.requested_observable_declarations.get(key)
            if declaration is None:
                continue
            seen.add(key)
            found.append(declaration)
        return tuple(found)

    def _restated_sufficiency_row(
        self,
        declaration: Mapping[str, Any] | None,
        *,
        display_unit: str,
        display_value: float,
        uncertainty: Any,
        uncertainty_basis: str,
        evidence_backed: bool = False,
        unquantified_components: tuple[str, ...] = (),
        uncertainty_observations: tuple[str, ...] = (),
        uncertainty_combination: Mapping[str, Any] | None = None,
    ) -> dict[str, Any]:
        """The claim's numbers in the unit its declaration asked for.

        A conversion the host cannot make is a refusal, not a silent
        comparison: an uncertainty in a unit of another dimension says
        nothing about a tolerance, and stating it would be a host word
        that is false.
        """

        row = {
            "display_value": display_value,
            "uncertainty": uncertainty,
            "uncertainty_basis": uncertainty_basis,
            "uncertainty_evidence_backed": evidence_backed,
            "unquantified_components": list(unquantified_components),
            # What the host observed and did not rule on. It rides the
            # assessment rather than only the claim because the
            # assessment is what the session, the workspace record and
            # the settlement all read: a fact computed where nothing
            # consumes it is the pattern this laboratory has paid for
            # five times, and three refusals became this report on the
            # condition that a reader of `met` can see what it rests on
            # (owner ruling, 2026-09-10).
            "uncertainty_observations": list(uncertainty_observations),
            # Carried, never graded: which rule produced the total is
            # what decided the word in three windows running.
            "uncertainty_combination": (
                dict(uncertainty_combination)
                if uncertainty_combination
                else None
            ),
        }
        if (
            declaration is None
            or declaration.get("required_tolerance") is None
        ):
            return row
        declared_unit = str(declaration.get("unit") or "")
        if not declared_unit or declared_unit == display_unit:
            return row
        # Only the uncertainty. The value was restated too, for a
        # `decision_boundary` comparison that is retired, and
        # judge_sufficiency never read it -- so a claim in the wrong
        # dimension was refused outright when its declaration carried a
        # tolerance and recorded when it did not, decided by a field
        # nothing in the comparison uses. Under a tolerance the number
        # never reached the record at all, which is where the charter's
        # own repair route starts: declare the corrected observable and
        # name the one it retires, with both on the record.
        value = row["uncertainty"]
        if value is None:
            return row
        restated = _restate_display_value(
            float(value), display_unit, declared_unit
        )
        if restated is None:
            raise ContractError(
                f"claim {declaration.get('observable_id')!r} states its "
                f"uncertainty in {display_unit!r} while its declaration "
                f"asks for {declared_unit!r}, and the host owns no "
                "conversion between them; state it in the declared unit"
            )
        row["uncertainty"] = restated
        return row

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
        sufficiency_rows: list[Any] = []
        for item in values["claims"]:
            # Shape before lookup: a pairing error in the claim itself is
            # the model's to fix and costs nothing to find, so it is
            # named before the receipt registry is searched.
            if (
                item.get("uncertainty") is not None
                or str(item.get("uncertainty_reference", "")).strip()
            ) and not str(item.get("uncertainty_basis", "")).strip():
                raise ContractError(
                    "an uncertainty needs uncertainty_basis: measured, "
                    "inferred, or asserted. All three count in full; the "
                    "word says which it is, and a reader cannot tell a "
                    "computed spread from a recalled one without it"
                )
            # A word that names evidence must name the evidence. The
            # basis was a free string, so "measured" cost nothing to
            # write and the host could not tell it from a judgement.
            shape_basis = str(item.get("uncertainty_basis", "")).strip()
            shape_reference = str(
                item.get("uncertainty_reference", "")
            ).strip()
            if shape_basis in {"measured", "inferred"} and not shape_reference:
                raise ContractError(
                    f"uncertainty_basis {shape_basis!r} names evidence, so "
                    "it needs uncertainty_reference. For 'measured' that is "
                    "'<receipt_sha256>:<quantity_id>' -- the host reads that "
                    "quantity and checks your number against it, exactly as "
                    "it copies the value you claim. For 'inferred' it is the "
                    "receipt or registered constant your judgement rests on. "
                    "State 'asserted' if it is your judgement alone -- an "
                    "assertion is never penalised, it simply does not "
                    "discharge a tolerance on its own"
                )
            if shape_basis == "asserted" and shape_reference:
                raise ContractError(
                    "uncertainty_basis 'asserted' is your own judgement and "
                    "takes no uncertainty_reference; name the basis the "
                    "reference supports instead"
                )
            for component in item.get("uncertainty_components") or ():
                component_basis = str(component.get("basis", "")).strip()
                component_reference = str(
                    component.get("reference", "")
                ).strip()
                if component_basis in {"measured", "inferred"}:
                    # The schema promises a reader that a component's
                    # reference is "the receipt or registered constant",
                    # and the host checked only that the string was not
                    # empty -- so a component could name evidence that
                    # does not exist under a word that claims it does.
                    # A sentence at a point of use naming a check the
                    # host does not make is a capability claim.
                    if not component_reference:
                        raise ContractError(
                            "an uncertainty component with basis "
                            f"{component_basis!r} needs its reference"
                        )
                    # A component may cite the quantity too, in the
                    # same shape the claim's own reference uses; the
                    # host resolves the receipt half. It does not check
                    # a component's magnitude, because how the terms add
                    # is the science and the total is what a tolerance
                    # is compared against.
                    component_receipt, _, component_quantity = (
                        component_reference.partition(":")
                    )
                    resolved, detail = self._uncertainty_evidence(
                        component_receipt, component_quantity
                    )
                    if resolved and component_quantity:
                        # The schema promises a reader that a component's
                        # reference is "resolved by the host as the
                        # claim's own reference is", and the host
                        # resolved only the receipt half -- so a
                        # component could name a quantity that does not
                        # exist under a word that says it does. The
                        # magnitude is deliberately not graded: how the
                        # terms add, and at what coverage each is
                        # stated, is the science. Existence is
                        # provenance and the host owns that.
                        if (
                            self._cited_uncertainty_quantity(
                                component_receipt, component_quantity
                            )
                            is None
                        ):
                            resolved = False
                            detail = (
                                f"receipt {component_receipt[:8]} carries no "
                                f"quantity {component_quantity!r}"
                            )
                    if not resolved:
                        raise ContractError(
                            "uncertainty component reference "
                            f"{component_reference!r} does not resolve: "
                            f"{detail}"
                        )
                elif component_basis == "asserted" and component_reference:
                    raise ContractError(
                        "an uncertainty component with basis 'asserted' "
                        "is your own judgement and takes no reference"
                    )
                if component_basis == "unquantified":
                    if component.get("magnitude") is not None:
                        raise ContractError(
                            "an unquantified component carries no "
                            "magnitude; give it a basis that matches the "
                            "number, or drop the number and keep the term"
                        )
                elif component.get("magnitude") is None:
                    raise ContractError(
                        "an uncertainty component needs a magnitude, or "
                        "basis 'unquantified' if you cannot give one"
                    )
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
                # Naming neither the quantity nor the alternatives made this
                # refusal unactionable: a campaign lost four analysis chains
                # to it and the logs could not say why without opening the
                # receipts by hand. Every other ChemSmart refusal names the
                # boundary it hit; this one now does too.
                available = sorted(
                    str(quantity.quantity_id) for quantity in quantities
                )
                trouble = "is ambiguous" if matches else "is absent"
                raise ContractError(
                    f"analysis claim quantity {quantity_id!r} {trouble} in "
                    f"the cited {source_kind} receipt; it carries "
                    f"{available}"
                )
            quantity = matches[0]
            if quantity.data_kind in {"text", "text_vector"}:
                raise ContractError("analysis claims must be numerical")
            display_unit = str(item["display_unit"])
            display_value = convert_normalized_value(
                quantity.value, quantity.dimension, display_unit
            )
            uncertainty = item.get("uncertainty")
            uncertainty_basis = str(item.get("uncertainty_basis", "")).strip()
            uncertainty_reference = str(
                item.get("uncertainty_reference", "")
            ).strip()
            # Only a magnitude the host can check discharges an
            # evidence obligation (owner ruling, 2026-09-09). `met` had
            # meant "a sufficiently small number stated beside a
            # resolvable citation": an uncertainty of 0.01 kJ/mol citing
            # a receipt whose quantity was an unrelated standard-state
            # correction reached it. Binding a number to its numerical
            # source is provenance and the host owns that; judging
            # whether it estimates the relevant scientific error is
            # chemistry and stays the session's. ChemSmart already draws
            # that line for the claimed value, which the model never
            # types -- the uncertainty is held to the same rule.
            evidence_backed = False
            uncertainty_observations: tuple[str, ...] = ()
            approximates = item.get("approximates")
            if approximates is not None and not isinstance(
                approximates, Mapping
            ):
                raise ContractError(
                    "approximates is a mapping naming the declared "
                    "observable this number approximates, the "
                    "relationship in your own words, and its basis"
                )
            approximates_record = (
                {
                    "observable_id": str(
                        approximates.get("observable_id") or ""
                    ),
                    "relationship": str(
                        approximates.get("relationship") or ""
                    ),
                    "basis": str(approximates.get("basis") or ""),
                }
                if isinstance(approximates, Mapping)
                else None
            )
            if approximates_record is not None and not (
                approximates_record["observable_id"]
                and approximates_record["relationship"]
            ):
                raise ContractError(
                    "approximates needs the observable_id it stands in "
                    "for and the relationship in your own words: the "
                    "host ships no vocabulary of approximation kinds and "
                    "will not infer one from an identifier"
                )
            combination = item.get("uncertainty_combination")
            if combination is not None and not isinstance(
                combination, Mapping
            ):
                raise ContractError(
                    "uncertainty_combination is a mapping stating the rule "
                    "you combined your components with, what its inputs "
                    "mean, the coverage the total claims, and what you "
                    "assume about dependence between the terms"
                )
            combination_record = (
                {
                    "rule": str(combination.get("rule") or ""),
                    "input_meaning": str(
                        combination.get("input_meaning") or ""
                    ),
                    "coverage": str(combination.get("coverage") or ""),
                    "dependence": str(combination.get("dependence") or ""),
                }
                if isinstance(combination, Mapping)
                else None
            )
            if uncertainty_reference:
                cited_receipt, _, cited_quantity = (
                    uncertainty_reference.partition(":")
                )
                resolved, evidence_detail = self._uncertainty_evidence(
                    cited_receipt, cited_quantity
                )
                if not resolved:
                    raise ContractError(
                        f"uncertainty_reference {uncertainty_reference!r} "
                        f"does not resolve: {evidence_detail}. If a "
                        "number of yours entered the chain, derive it "
                        "instead of typing it -- an electron count is the "
                        "difference of the two states' own charge "
                        "selectors, a threshold the host owns is a "
                        "convention rather than your value -- and the whole "
                        "chain stays the host's. Otherwise name a receipt "
                        "this host minted or a registered constant, or "
                        "state the basis as 'asserted'"
                    )
                if uncertainty_basis == "measured":
                    if not cited_quantity:
                        raise ContractError(
                            "uncertainty_basis 'measured' means the host can "
                            "read your number where you got it, so the "
                            "reference names the quantity too: "
                            f"'{cited_receipt}:<quantity_id>'. Cite the "
                            "quantity, or state 'inferred' if the number is "
                            "your judgement standing on that receipt"
                        )
                    # Its own name: this loop's `quantity` is the one
                    # being claimed, and rebinding it here made every
                    # claim on the `measured` path record the
                    # *uncertainty's* quantity_id, canonical value and
                    # value digest -- so the claim said 0.219 eV and
                    # 0.350 eV at once and its digest no longer hashed
                    # the delivered number. The one path that can
                    # discharge a tolerance was the one that corrupted
                    # the claim.
                    cited = self._cited_uncertainty_quantity(
                        cited_receipt, cited_quantity
                    )
                    if cited is None:
                        raise ContractError(
                            f"receipt {cited_receipt} carries no quantity "
                            f"{cited_quantity!r}; name one it does"
                        )
                    stated = convert_normalized_value(
                        cited.value, cited.dimension, display_unit
                    )
                    if not isinstance(stated, (int, float)):
                        raise ContractError(
                            f"quantity {cited_quantity!r} is not a number "
                            "this host can read as an uncertainty"
                        )
                    # A zero magnitude and a single-receipt spread
                    # were refused here until the boundary was drawn.
                    # Both checks were formal and neither was
                    # scientifically defensible: a spread of exactly
                    # zero is a real observation, a variance over many
                    # samples inside one receipt compares plenty, and
                    # neither refusal stopped the substitution it was
                    # aimed at -- 1e-9 walks past the first and a split
                    # receipt past the second. What replaces them is a
                    # report: the host says what it saw about the number
                    # and the session owns the claim that the number is
                    # adequate (owner ruling, 2026-09-10). `met` is
                    # therefore cheaper than it was, and it is now
                    # auditable at the point a human reads it, which the
                    # refusals never made it.
                    if uncertainty is None:
                        # Cite the quantity and the host supplies the
                        # magnitude, exactly as it copies the value you
                        # claim: the model never writes the number. This
                        # is the route a planned estimator takes -- the
                        # executor names the output its own analysis
                        # chain produced and the host reads it.
                        uncertainty = abs(float(stated))
                    elif not (
                        _reads_as(float(uncertainty), abs(stated))
                        or _reads_as(abs(stated), float(uncertainty))
                    ):
                        raise ContractError(
                            f"quantity {cited_quantity!r} reads {stated} "
                            f"{display_unit} and your uncertainty states "
                            f"{uncertainty}. 'measured' means the host reads "
                            "your number where you got it; state that "
                            "number, omit it and let the host copy it, or "
                            "use 'inferred'"
                        )
                    evidence_backed = True
                uncertainty_observations = self._uncertainty_observations(
                    cited_receipt, cited_quantity
                )
            unquantified_components = tuple(
                str(component.get("meaning") or "")
                for component in item.get("uncertainty_components") or ()
                if str(component.get("basis", "")).strip() == "unquantified"
            )
            # A tolerance is declared in the observable's own unit and a
            # claim states its uncertainty in the claim's display unit,
            # and nothing converted between them: a 2 kJ/mol tolerance
            # certified a 1 kcal/mol uncertainty as met, and printed it
            # as "1.0 kJ/mol". The comparison is the host's to own, so
            # the numbers are restated into the declared unit here --
            # through the same function the expectation row uses -- and
            # judge_sufficiency stays a comparison over normalised
            # inputs.
            # Either field answers, claim_id first, exactly as the
            # completion gate joins: a plan carries the host-minted
            # quantity_id and its own short input label in claim_id, and
            # a gate reading claim_id alone lost six observables three
            # fields away from the id it wanted. Reading one field here
            # while the gate reads two meant a legitimate alternate
            # label delivered the observable and produced no assessment
            # of it at all.
            # Every declaration this claim delivers, not the first:
            # the gate credits both ids, and assessing only one let a
            # claim deliver a tolerance-bearing observable under a
            # tolerance-free label and be assessed against neither.
            for declaration in self._declarations_for_claim(
                str(item["claim_id"]), str(quantity.quantity_id)
            ):
                sufficiency_rows.append(
                    judge_sufficiency(
                        declaration,
                        self._restated_sufficiency_row(
                            declaration,
                            display_unit=display_unit,
                            display_value=display_value,
                            uncertainty=uncertainty,
                            uncertainty_basis=uncertainty_basis,
                            evidence_backed=evidence_backed,
                            unquantified_components=unquantified_components,
                            uncertainty_observations=(
                                uncertainty_observations
                            ),
                            uncertainty_combination=combination_record,
                        ),
                    )
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
                    uncertainty=(
                        None if uncertainty is None else float(uncertainty)
                    ),
                    uncertainty_basis=uncertainty_basis,
                    uncertainty_reference=uncertainty_reference,
                    uncertainty_observations=uncertainty_observations,
                    uncertainty_combination=combination_record,
                    approximates=approximates_record,
                    uncertainty_components=tuple(
                        {
                            "meaning": str(component.get("meaning") or ""),
                            "magnitude": (
                                None
                                if component.get("magnitude") is None
                                else float(component["magnitude"])
                            ),
                            "basis": str(component.get("basis") or ""),
                            "reference": str(component.get("reference") or ""),
                            # The report reached the claim's own citation
                            # and stopped there, so moving a magnitude
                            # into a component escaped it entirely --
                            # and a component is the shape the charter
                            # actively encourages. A live claim carried
                            # a 0.10 V solvation term as an `asserted`
                            # component and a reference term composed
                            # from three literature constants with no
                            # engine output at all, and neither was
                            # observed anywhere (SUFFICIENCY-5).
                            **self._component_observations(
                                component,
                                str(item.get("display_unit") or ""),
                            ),
                        }
                        for component in item.get("uncertainty_components")
                        or ()
                    ),
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
                # How each delivered number stands against the precision
                # its declaration asked for. The driver reads this from
                # the stream: the same judgement, from one function, for
                # the gate and the settlement alike.
                "sufficiency": tuple(
                    row for row in sufficiency_rows if row is not None
                ),
            },
            idempotency_key="analysis-claims:" + record.receipt_sha256,
        )
        # And tell the session what the host made of its own numbers, in
        # the turn it made them. The verdict was written to the stream
        # and never to the reply, so a session could not learn in-turn
        # that its number missed the tolerance it had itself declared --
        # it had to spend a whole cycle to be told through the wake,
        # while two of the three routes that answer cost no engine call
        # at all. The refusal and the reply are where this model is
        # actually taught.
        self._reply_observations = tuple(
            row for row in sufficiency_rows if row is not None
        )
        # The current assessment of each requirement, latest wins, so
        # the refusal verifier can check that a precision a session
        # says it cannot establish is one this goal actually has open.
        for row in sufficiency_rows:
            if row is not None and row.get("observable_id"):
                self.requirement_assessments[str(row["observable_id"])] = row
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
            if label == "trusted artifact":
                # The refusal a session meets when it names evidence the
                # host printed: a digest off a run outcome, a node id off
                # the record. Ten such refusals in 22 seconds named neither
                # the id rule nor the nearest bound id (REACH-1 ino3,
                # 2026-09-06).
                import difflib

                nearest = difflib.get_close_matches(
                    str(key), [str(item) for item in known], n=3, cutoff=0.5
                )
                results = sorted(
                    item for item in known if "-result-" in str(item)
                )
                if re.fullmatch(r"[0-9a-f]{64}", str(key) or ""):
                    shape = (
                        f"{key!r} is a content digest, not an id; a result "
                        "registered from it carries the id "
                        f"<program>-result-{str(key)[:16]}."
                    )
                else:
                    shape = f"{key!r} is not a registered artifact id."
                raise RoutedContractError(
                    gate="artifact.id_is_registered",
                    invariant=(
                        "a reading tool opens only an artifact the host "
                        "registered under its id; ids are host-minted, never "
                        "typed from a digest or a node name."
                    ),
                    diagnosis=shape
                    + (f" Nearest bound ids: {nearest}." if nearest else "")
                    + (
                        f" Registered results: {results[:8]}"
                        + (" ..." if len(results) > 8 else "")
                        + "."
                        if results
                        else " No result is registered in this session."
                    ),
                    route=(
                        "open the result by the artifact_id the workspace "
                        "record or inspect_run shows "
                        "(<program>-result-<16 hex of its digest>); a digest "
                        "is citable in a decision, never an argument."
                    ),
                ) from exc
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


def _renamed_output_id(
    renames: dict[tuple[str, str], str], output_id: str
) -> str:
    """Follow an output rename into the required-output set.

    A required output names an id, not the node that produces it, so the
    rename is matched on the name alone.  Renaming the same name on two
    different producers is therefore ambiguous here and is left alone rather
    than guessed at; the plan builder will say so if the requirement no longer
    resolves.
    """

    matches = {
        renamed
        for (_node_id, old_id), renamed in renames.items()
        if old_id == output_id
    }
    return matches.pop() if len(matches) == 1 else output_id


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
    findings: list[str] = []
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
        findings.append(f"{tool_name} requires {described}")
    unknown = sorted(set(arguments).difference(properties))
    if unknown:
        findings.append(
            f"{tool_name} does not accept {unknown}; it accepts "
            f"{sorted(properties)}"
        )
    for name, value in arguments.items():
        # An unknown name has already been reported and has no schema to
        # check it against.
        if name in properties:
            _collect_json_violations(name, value, properties[name], findings)
    if findings:
        raise ContractError(_combined_violation_text(findings))


#: How many violations one rejection may enumerate.  A systematic mistake in a
#: large DAG can produce a violation per node; the caller needs enough to fix a
#: batch in one revision, not a message longer than the payload it describes.
_REPORTED_VIOLATIONS = 12


def _combined_violation_text(findings: list[str]) -> str:
    """Render one rejection covering every violation that was found.

    A single violation keeps its exact wording, so nothing that reads these
    messages has to learn a second shape for the common case.
    """

    if len(findings) == 1:
        return findings[0]
    shown = findings[:_REPORTED_VIOLATIONS]
    body = "; ".join(
        f"({index}) {message}" for index, message in enumerate(shown, 1)
    )
    text = f"{len(findings)} arguments are invalid: {body}"
    remaining = len(findings) - len(shown)
    if remaining:
        text += f"; and {remaining} more not listed"
    return text


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

    This raises on the first violation it finds.  Whole-payload validation goes
    through :func:`_collect_json_violations` instead, so one submission can be
    told about every independent slip at once.
    """

    findings: list[str] = []
    _collect_json_violations(name, value, schema, findings)
    if findings:
        raise ContractError(findings[0])


def _collect_json_violations(
    name: str,
    value: Any,
    schema: Mapping[str, Any],
    findings: list[str],
) -> None:
    """Append every violation under ``name`` rather than raising the first.

    A workflow DAG is authored in one payload, so a single bad identifier in
    node ten used to cost a full resubmission of all ten nodes -- and then the
    next independent slip cost another.  Across the recorded sessions, most
    multi-failure streaks carried violations on entirely different fields, so
    reporting them one per turn was pure round-trip tax rather than a cascade
    that had to be unwound in order.

    Descent stops at a value whose *type* is already wrong: checking a pattern
    or a minimum against it would only add noise about a value the caller is
    going to replace anyway.
    """

    alternatives = schema.get("oneOf")
    if isinstance(alternatives, list):
        matches = 0
        for alternative in alternatives:
            if not isinstance(alternative, Mapping):
                continue
            probe: list[str] = []
            _collect_json_violations(name, value, alternative, probe)
            if not probe:
                matches += 1
        if matches != 1:
            shapes = ", ".join(
                str(item.get("type", "?"))
                for item in alternatives
                if isinstance(item, Mapping)
            )
            findings.append(
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
            findings.append(
                f"tool argument {name} must be one of {allowed}, but got "
                f"{type(value).__name__} {_offending_text(value)}"
            )
            return
        if value is None:
            # No further keyword applies to an explicit null.
            return
    elif not _matches(expected):
        findings.append(
            f"tool argument {name} must be {expected}, but got "
            f"{type(value).__name__} {_offending_text(value)}"
        )
        return
    if "enum" in schema and value not in schema["enum"]:
        allowed = list(schema["enum"])
        findings.append(
            f"tool argument {name} is {_offending_text(value)}, which is not "
            f"one of {allowed}"
        )
    if isinstance(value, str) and schema.get("pattern"):
        pattern = str(schema["pattern"])
        if re.fullmatch(pattern, value) is None:
            finding = (
                f"tool argument {name} is {_offending_text(value)}, which "
                f"does not match the required pattern {pattern}"
            )
            if re.fullmatch(r"[0-9a-f]{64}", value):
                finding += (
                    "; that is a content digest, and a result is opened by "
                    f"its registered id <program>-result-{value[:16]}, which "
                    "the workspace record and inspect_run show"
                )
            findings.append(finding)
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        if "minimum" in schema and value < float(schema["minimum"]):
            findings.append(
                f"tool argument {name} is {value}, below its minimum "
                f"{schema['minimum']}"
            )
        if "maximum" in schema and value > float(schema["maximum"]):
            findings.append(
                f"tool argument {name} is {value}, above its maximum "
                f"{schema['maximum']}"
            )
        if "exclusiveMinimum" in schema and value <= float(
            schema["exclusiveMinimum"]
        ):
            findings.append(
                f"tool argument {name} is {value}, which must be greater than "
                f"{schema['exclusiveMinimum']}"
            )
        if "exclusiveMaximum" in schema and value >= float(
            schema["exclusiveMaximum"]
        ):
            findings.append(
                f"tool argument {name} is {value}, which must be less than "
                f"{schema['exclusiveMaximum']}"
            )
    if isinstance(value, list):
        if "minItems" in schema and len(value) < int(schema["minItems"]):
            findings.append(
                f"tool argument {name} has {len(value)} items, fewer than the "
                f"required {schema['minItems']}"
            )
        if "maxItems" in schema and len(value) > int(schema["maxItems"]):
            findings.append(
                f"tool argument {name} has {len(value)} items, more than the "
                f"allowed {schema['maxItems']}"
            )
    if isinstance(value, list) and isinstance(schema.get("items"), Mapping):
        for index, item in enumerate(value):
            # Index the path: "outputs[]" tells a caller with eight outputs
            # nothing about which one to change.
            _collect_json_violations(
                f"{name}[{index}]", item, schema["items"], findings
            )
    if isinstance(value, dict):
        properties = schema.get("properties", {})
        additional = schema.get("additionalProperties")
        required = set(schema.get("required", ()))
        missing = sorted(required.difference(value))
        if missing:
            findings.append(
                f"tool argument {name} is missing {missing}; it supplied "
                f"{sorted(value)}"
            )
        if schema.get("additionalProperties") is False:
            unknown = sorted(set(value).difference(properties))
            if unknown:
                findings.append(
                    f"tool argument {name} supplied {unknown}, which this "
                    f"object does not accept; it accepts {sorted(properties)}"
                )
        for child_name, child_value in value.items():
            child_schema = properties.get(child_name)
            if child_schema is None and isinstance(additional, Mapping):
                child_schema = additional
            if child_schema is not None:
                _collect_json_violations(
                    f"{name}.{child_name}",
                    child_value,
                    child_schema,
                    findings,
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


def _launch_lease_seconds(resources: Any, envelope: Any) -> int:
    """How long a launch reservation may still be a live engine.

    The approved node timeout is the host's own bound on an engine's
    life, plus the postprocessing reserve the human granted, because the
    host's own reading of the result happens while the reservation is
    still held. Nothing is invented: both numbers are already approved.

    Zero when no bound is known, which reads as "nothing can be
    concluded" rather than "alive" -- the same reading every reservation
    written before leases existed gets.
    """

    seconds = 0.0
    node_timeout = getattr(resources, "node_timeout_seconds", None)
    try:
        seconds += float(node_timeout or 0.0)
    except (TypeError, ValueError):
        return 0
    reserve = getattr(envelope, "postprocess_reserve_seconds", None)
    try:
        seconds += float(reserve or 0.0)
    except (TypeError, ValueError):
        pass
    return int(seconds) if seconds > 0 else 0


def _launch_reserver() -> str:
    """Who took this reservation, for a human reading the record.

    Diagnosis only: liveness is decided by the lease, never by this
    string, because a host name and a pid mean nothing to a reader on
    another node.
    """

    parts = [socket.gethostname(), str(os.getpid())]
    for variable in ("SLURM_JOB_ID", "SLURM_ARRAY_TASK_ID"):
        value = os.environ.get(variable)
        if value:
            parts.append(f"{variable}={value}")
    return " ".join(parts)


def _write_branch_request(node_workspace: Path, command: list[str]) -> None:
    """Record what was asked of this calculation, before it is asked.

    Two files, both host-owned and both written before launch so a node
    killed mid-engine still carries them:

    ``command.txt`` -- the exact CHEMSMART command, shell-quoted, as it
    was handed to the operating system. Not a reconstruction: the same
    list the process was started with.

    ``project.yaml`` -- a copy of the project configuration that command
    names, so the branch does not depend on a file elsewhere in the
    workspace still existing, or still saying what it said.

    A missing or unreadable project file is recorded as such rather than
    raising: the engine's own run is the evidence, and the host does not
    fail a calculation over its own bookkeeping.
    """

    _write_host_execution_artifact(
        node_workspace / "command.txt",
        shlex.join(str(item) for item in command) + "\n",
    )
    project = ""
    for flag, value in zip(command, command[1:]):
        if str(flag) in {"-p", "--project"}:
            project = str(value)
            break
    if not project:
        return
    source = Path(project)
    try:
        payload = source.read_text(encoding="utf-8")
    except OSError as exc:
        _write_host_execution_artifact(
            node_workspace / "project.yaml.missing",
            f"{source}: {exc}\n",
        )
        return
    _write_host_execution_artifact(
        node_workspace / "project.yaml",
        f"# copied from {source}\n{payload}",
    )


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
    if observation.state.startswith("external_signal"):
        findings.add("execution.process.external_signal")
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
        producer_output_id = str(getattr(item, "producer_output_id", "") or "")
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
            excursion=node.excursion,
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
