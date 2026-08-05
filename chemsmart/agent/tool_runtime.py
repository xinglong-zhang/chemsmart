"""Host-owned dispatcher for the command-compiled model tool surface."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys
from typing import Any, Mapping

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
    ProgramComponentConformanceReceiptV1,
    ProgramSupportOverlayV1,
    ProgramCapabilityRegistryV1,
    ResolvedEngineBindingV1,
    ResolvedProgramBindingV1,
    TrustedComputeEnvironmentReceiptV1,
    build_approved_execution_overlay,
    build_command_compiled_preview_overlay,
    consume_pyscf_compute_environment_receipt,
    load_program_capabilities,
    query_capability,
    query_environment,
    resolve_engine_binding,
    resolve_program_binding,
)
from chemsmart.agent.cli_schema import LiveClickSchemaV1, build_live_click_schema
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
    ExecutionResourceSpecV1,
    FrozenWorkflowApprovalV1,
    OptimizedGeometryHandoffV1,
    ProgramExecutionReceiptV1,
    ProgramResultValidationReceiptV1,
    ProjectArtifactPromotionV1,
    ScientificDecisionRecordV1,
    WorkflowExecutionApprovalV1,
    WorkflowNodeRunStateV1,
    WorkflowRunStateV1,
    bind_project_promotion_validation,
    build_program_execution_invocation,
    build_program_execution_receipt,
    build_program_result_validation_receipt,
    build_scientific_decision_record,
    build_validated_data_edge_binding,
    handoff_optimized_pyscf_geometry,
    promote_project_candidate,
)
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
    project_scientific_materializations,
    read_project_yaml,
    render_project_yaml,
    validate_project_yaml,
)
from chemsmart.agent.runtime.event_store import RuntimeEventStore
from chemsmart.agent.runtime.events import EventKind
from chemsmart.agent.tool_specs import (
    AgentToolSurfaceV1,
    build_approved_execution_tool_surface,
    build_command_compiled_tool_surface,
    build_single_agent_baseline_tool_surface,
)
from chemsmart.analysis.quantity_expressions import (
    QuantityExpressionNodeV1,
    QuantityExpressionRequestV1,
    convert_normalized_value,
    quantity_expression_receipt_from_record,
)
from chemsmart.analysis.result_quantities import (
    QuantityExtractionReceiptV1,
    QuantitySelectorV1,
    QuantityValueV1,
    ThermochemistryReceiptV1,
    make_quantity_value,
    quantity_extraction_receipt_from_record,
    thermochemistry_receipt_from_record,
)
from chemsmart.agent.workflow_context import (
    project_workflow_context,
    workflow_context_enabled,
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
    launch_failure_observation,
    observe_process,
)


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
            item.evidence_ref
            for item in materializations
            if item.evidence_ref
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
    receipt: EnvironmentCapabilityReceiptV1 | TrustedComputeEnvironmentReceiptV1,
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


def _pyscf_environment_evidence(
    *,
    output_artifacts: tuple[TrustedArtifactRefV1, ...],
    run_receipt: Mapping[str, Any] | None,
    capability_environment: EnvironmentCapabilityReceiptV1 | None,
) -> tuple[Mapping[str, Any], tuple[str, ...]]:
    """Compare different environment receipt types by stable semantics."""

    findings: list[str] = []
    candidates: list[tuple[TrustedArtifactRefV1, Mapping[str, Any] | None]] = []
    for artifact in output_artifacts:
        if artifact.kind != "json":
            continue
        try:
            path = _current_artifact_path(
                artifact, field_name="PySCF environment receipt"
            )
            raw = json.loads(path.read_text(encoding="utf-8"))
        except (ContractError, OSError, TypeError, ValueError, json.JSONDecodeError):
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
    if run_receipt is None or run_receipt.get(
        "environment_receipt_sha256"
    ) != run_environment_sha256:
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
            observed_facts = _environment_semantic_facts(adapted)
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
        if result_artifact is not None and run_receipt.get(
            "result_sha256"
        ) != result_artifact.sha256:
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
        raise ContractError("stationary-point binding requires a real approval")
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
        raise ContractError(
            "frozen stationary-point policy is unavailable"
        )
    if plan is None or not isinstance(plan, ScientificWorkflowPlanV2):
        raise ContractError(
            "stationary-point policy requires its exact scientific plan"
        )
    if plan.plan_sha256 != frozen_approval.plan_sha256:
        raise ContractError("stationary-point policy plan differs from approval")
    if policy.task_spec_sha256 != frozen_approval.task_spec_sha256 or (
        policy.task_spec_sha256 != plan.task_spec_sha256
    ):
        raise ContractError("stationary-point policy task differs from approval")
    if policy.hessian_node_id not in frozen_approval.approved_node_ids:
        raise ContractError("stationary-point Hessian node is not approved")
    matching_nodes = tuple(
        node for node in plan.nodes if node.node_id == policy.hessian_node_id
    )
    if len(matching_nodes) != 1 or matching_nodes[0].stage != "hess":
        raise ContractError("stationary-point policy must bind a Hessian node")
    if hessian_node_id and policy.hessian_node_id != hessian_node_id:
        raise ContractError("stationary-point policy targets another Hessian node")
    if policy.require_finite_modes is not True:
        raise ContractError("stationary-point policy must require finite modes")
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
    """Admit only an exact policy-bound HESS classification handoff.

    The CLI runner owns engine/artifact invariants but does not own an
    imaginary-mode expectation.  An approved agent execution may therefore
    accept its ``engine_complete``/``unclassified`` intermediate receipt only
    when the exact frozen policy is present.  The parent result validator must
    still apply that policy and return ``validated`` before the node succeeds.
    """

    if jobtype != "hess" or not isinstance(
        stationary_point_policy, StationaryPointValidationPolicyV1
    ):
        return False
    if (
        not hessian_node_id
        or stationary_point_policy.hessian_node_id != hessian_node_id
    ):
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
        and stationary_point_policy.require_finite_modes is True
        and stationary_point_policy.require_symmetric_hessian is True
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
        program_binding_receipts: Mapping[
            str, ResolvedProgramBindingV1
        ] = {},
        engine_binding_receipts: Mapping[
            str, ResolvedEngineBindingV1
        ] = {},
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
        stationary_point_policy: (
            StationaryPointValidationPolicyV1 | None
        ) = None,
        scientific_workflow_plan: ScientificWorkflowPlanV2 | None = None,
        execution_server: str = "",
        execution_environment: Mapping[str, str] = {},
    ) -> None:
        self.event_store = event_store
        self.registry = registry or load_program_capabilities()
        self.live_schema = live_schema or build_live_click_schema()
        preview_overlay = build_command_compiled_preview_overlay(
            self.registry,
            conformance_receipts=component_conformance_receipts,
            live_schema=self.live_schema,
        )
        h1_surface = build_command_compiled_tool_surface(self.registry)
        h0_surface = build_single_agent_baseline_tool_surface(self.registry)
        execution_surface = build_approved_execution_tool_surface(self.registry)
        if tool_surface is not None and tool_surface.tool_schema_sha256 not in {
            h0_surface.tool_schema_sha256,
            h1_surface.tool_schema_sha256,
            execution_surface.tool_schema_sha256,
        }:
            raise ContractError("injected tool surface is not a canonical profile")
        self.surface = tool_surface or h1_surface
        self.artifacts = dict(artifacts)
        self.task_spec_sha256s = frozenset(task_spec_sha256s)
        self.approved_workspace = (
            Path(approved_workspace).resolve()
            if approved_workspace is not None
            else None
        )
        self.execution_resources = execution_resources
        self.workflow_execution_approval = workflow_execution_approval
        self.frozen_workflow_approval = frozen_workflow_approval
        self.stationary_point_policy = stationary_point_policy
        _validate_stationary_point_policy_binding(
            self.frozen_workflow_approval,
            self.stationary_point_policy,
            plan=scientific_workflow_plan,
        )
        self.execution_server = str(execution_server)
        self.execution_environment = {
            str(key): str(value) for key, value in execution_environment.items()
        }
        if self.surface.profile == "command_compiled_approved_execution":
            if self.approved_workspace is None:
                raise ContractError("execution profile requires an approved workspace")
            if self.execution_resources is None:
                raise ContractError("execution profile requires host-owned resources")
            if self.workflow_execution_approval is None:
                raise ContractError("execution profile requires workflow approval")
            if self.frozen_workflow_approval is not None:
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
            if (
                Path(self.workflow_execution_approval.workspace).resolve()
                != self.approved_workspace
            ):
                raise ContractError("workflow approval targets another workspace")
            execution_evidence_sha256 = canonical_sha256(
                {
                    "approval_sha256": (
                        self.workflow_execution_approval.approval_sha256
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
            evidence_overlay = build_approved_execution_overlay(
                registry=self.registry,
                preview_overlay=preview_overlay,
                approved_nodes=(
                    (item.program, item.jobtype, item.engine)
                    for item in self.workflow_execution_approval.node_bindings
                ),
                execution_evidence_sha256=execution_evidence_sha256,
            )
        else:
            evidence_overlay = preview_overlay
        if support_overlay is not None:
            if support_overlay.base_registry_sha256 != self.registry.registry_sha256:
                raise ContractError("injected support overlay uses another registry")
            if support_overlay.overlay_sha256 != evidence_overlay.overlay_sha256:
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
            raise ContractError("approved molecular identity registry key mismatch")
        self.environment_targets = tuple(environment_targets)
        self.compute_environment_receipts = tuple(
            compute_environment_receipts
        )
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
                self.functional_resolutions[
                    resolution.receipt_sha256
                ] = resolution
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
        self.quantity_extractions: dict[
            str, QuantityExtractionReceiptV1
        ] = {}
        self.quantity_extraction_selectors: dict[str, tuple[str, ...]] = {}
        self.quantity_extraction_bindings: dict[str, dict[str, str]] = {}
        self.thermochemistry_receipts: dict[
            str, ThermochemistryReceiptV1
        ] = {}
        self.quantity_expression_receipts: dict[str, Any] = {}
        self.quantity_expression_requests: dict[
            str, QuantityExpressionRequestV1
        ] = {}
        self.analysis_claim_records: dict[str, Any] = {}
        self.analysis_completion_receipts: dict[str, Any] = {}
        self.counterexamples: dict[str, CommandCounterexampleV1] = {}
        self.workflow_drafts: dict[str, CommandWorkflowDraftV1] = {}
        #: Last accepted scientific plan per workflow, so a later repair can be
        #: checked against the question and not only against its own findings.
        self.scientific_plans: dict[str, Any] = {}
        self.scientific_workflow_plans: dict[
            str, ScientificWorkflowPlanV2
        ] = {}
        self.materialized_workflows: dict[
            str, MaterializedWorkflowV1
        ] = {}
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
        self._command_contexts: dict[str, _CommandContext] = {}
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
                or receipt.live_cli_schema_sha256 != self.live_schema.schema_sha256
                or receipt.overlay_sha256 != self.overlay.overlay_sha256
            ):
                raise ContractError("seeded capability receipt is stale")
        for environment in self.environments.values():
            self._latest_environment_by_capability[
                environment.capability_receipt_sha256
            ] = environment
        self._rehydrate_analysis_event_records()

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
        """Persist host-prebound H0 evidence before any model action."""

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
            "synthesize_command": self._synthesize_command,
            "repair_command": self._repair_command,
            "preview_command": self._preview_command,
            "preflight_program_node": self._preflight_program_node,
            "inspect_calculation_artifact": self._inspect_calculation_artifact,
            "extract_result_quantities": self._extract_result_quantities,
            "derive_thermochemistry": self._derive_thermochemistry,
            "evaluate_quantity_expression": self._evaluate_quantity_expression,
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
        task_spec_sha256 = values["task_spec_sha256"]
        if task_spec_sha256 not in self.task_spec_sha256s:
            raise ContractError("scientific identity targets an unknown task spec")
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
        self._latest_environment_by_capability[
            capability.receipt_sha256
        ] = environment
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
                raise ContractError("functional equivalence targets another request")
            observed_claim_receipts = tuple(
                sorted(item.receipt_sha256 for item in verified_claims)
            )
            if (
                equivalence.status == "verified"
                and equivalence.claim_evidence_receipt_sha256s
                != observed_claim_receipts
            ):
                raise ContractError("functional equivalence uses other claim evidence")
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
            self.program_bindings[program_binding.binding_sha256] = program_binding
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
                self.engine_bindings[
                    engine_binding.binding_sha256
                ] = engine_binding
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
        self.project_validations[receipt.receipt_sha256] = receipt
        materializations = project_scientific_materializations(receipt)
        for materialization in materializations:
            self.functional_resolutions[
                materialization.receipt_sha256
            ] = materialization
        promotion = self.project_promotions.get(artifact.artifact_id)
        if promotion is not None and promotion.validation_status == "pending":
            self.project_promotions[artifact.artifact_id] = (
                bind_project_promotion_validation(
                    promotion, artifact, receipt
                )
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
        return {
            **canonical_data(receipt),
            "scientific_materializations": tuple(
                item.public_record() for item in materializations
            ),
            "decision_binding": _scientific_decision_binding_requirement(
                materializations
            ),
        }

    def _record_scientific_decision(self, turn_id: str, values: dict) -> Any:
        task_spec_sha256 = values["task_spec_sha256"]
        if task_spec_sha256 not in self.task_spec_sha256s:
            raise ContractError("scientific decision targets an unknown task spec")
        postprocessing_registries = (
            self.quantity_extractions,
            self.thermochemistry_receipts,
            self.quantity_expression_receipts,
            self.analysis_claim_records,
        )
        postprocessing_receipt_sha256s = tuple(
            str(item)
            for item in values.get("postprocessing_receipt_sha256s", ())
        )
        for receipt_sha256 in postprocessing_receipt_sha256s:
            require_sha256(
                receipt_sha256, "postprocessing_receipt_sha256"
            )
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
        if re.search(
            r"(?i)(?<![a-z0-9])(?:vwn\s*[35]|b3lypg|b3lyp5)(?![a-z0-9])",
            convention_narrative,
        ) and not functional_resolution_refs:
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

    def _plan_command_workflow(self, turn_id: str, values: dict) -> Any:
        """Record a broad DAG before execution-grade evidence is available."""

        nodes = []
        findings: list[dict[str, str]] = []
        declared_programs = {item.program: item for item in self.registry.programs}
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
                    {"node_id": node.node_id, "rule_id": "workflow.program.not_declared"}
                )
            elif node.jobtype not in capability.jobtypes:
                findings.append(
                    {"node_id": node.node_id, "rule_id": "workflow.job.not_declared"}
                )
            for item in node.inputs:
                if (
                    not item.producer_node_id
                    and (
                        not item.artifact_id
                        or item.artifact_id not in self.artifacts
                    )
                ):
                    findings.append(
                        {
                            "node_id": node.node_id,
                            "rule_id": "workflow.input.unresolved",
                        }
                    )
            nodes.append(node)
        draft = build_command_workflow_draft(
            workflow_id=values["workflow_id"],
            task_spec_id=values["task_spec_id"],
            nodes=tuple(nodes),
        )
        self.workflow_drafts[draft.draft_sha256] = draft
        scientific_plan = self._scientific_plan_from_draft(
            draft, findings=findings
        )
        if scientific_plan is not None:
            self.scientific_workflow_plans[
                scientific_plan.plan_sha256
            ] = scientific_plan
        unresolved_nodes = {
            item["node_id"] for item in findings
        } | {
            node.node_id
            for node in draft.nodes
            if node.dependencies or node.unresolved_fields
        }
        actionable = tuple(
            node.node_id
            for node in draft.nodes
            if node.node_id not in unresolved_nodes
        )
        unresolved = tuple(
            node.node_id for node in draft.nodes if node.node_id in unresolved_nodes
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
        result = {
            "workflow_draft": draft,
            "scientific_workflow_plan": scientific_plan,
            "actionable_node_ids": actionable,
            "unresolved_node_ids": unresolved,
            "findings": tuple(findings),
        }
        context = self._workflow_context(draft)
        if context is not None:
            result["workflow_context"] = context
        return result

    def _workflow_context(self, draft: CommandWorkflowDraftV1) -> Any:
        """Derive the dependency context the model would otherwise reconstruct.

        Host-derived and read-only: the model is told which nodes are runnable
        and what each waiting node is waiting for, and can never assert it.
        """

        if not workflow_context_enabled():
            return None
        return project_workflow_context(
            workflow_id=draft.workflow_id,
            nodes=draft.nodes,
            materialized_artifact_ids=self.artifacts,
        )

    def _scientific_plan_from_draft(
        self,
        draft: CommandWorkflowDraftV1,
        *,
        findings: list[dict[str, str]],
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
            else canonical_sha256(
                {"scientific_identity_sha256s": identities}
            )
        )
        scientific_nodes = []
        for node in draft.nodes:
            matching_capabilities = tuple(
                receipt
                for receipt in self.capabilities.values()
                if receipt.query.program == node.program
                and receipt.query.jobtype == node.jobtype
            )
            engines = tuple(
                sorted({receipt.query.engine for receipt in matching_capabilities})
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
                    receipt.receipt_sha256
                    == binding.capability_receipt_sha256
                    for receipt in matching_capabilities
                )
            }
            requested_program = (
                next(iter(requested_programs))
                if len(requested_programs) == 1
                else node.program
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
                        edge_id=(
                            "control." + dependency + "." + node.node_id
                        ),
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
        )
        if self.frozen_workflow_approval is not None and (
            self.frozen_workflow_approval.plan_sha256 != plan.plan_sha256
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
        validation_digest = values.get(
            "project_validation_receipt_sha256", ""
        )
        validation = (
            self._get(
                self.project_validations,
                validation_digest,
                "project validation receipt",
            )
            if validation_digest
            else None
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
            scientific_identity_sha256=values[
                "scientific_identity_sha256"
            ],
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
            live_schema=self.live_schema,
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
            critical_finding_count=len(
                validator.critical_finding_sha256s
            ),
            source_receipt_sha256=validator.source_receipt_sha256,
        )
        return {"safe_preview": receipt, "validator": validator}

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
        project_digest = values.get(
            "project_validation_receipt_sha256", ""
        )
        project = (
            self._get(
                self.project_validations,
                project_digest,
                "project validation receipt",
            )
            if project_digest
            else None
        )
        preview_digest = values.get("safe_preview_receipt_sha256", "")
        safe_preview = (
            self._get(
                self.safe_previews, preview_digest, "safe preview receipt"
            )
            if preview_digest
            else None
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
            if {
                item.receipt_sha256 for item in supplied_validators
            } != {
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
            safe_preview_receipt_sha256=(
                receipt.safe_preview_receipt_sha256
            ),
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
            previewed = (
                preflight is not None
                and preflight.plan_state == "previewed"
                and not preflight.critical_finding_sha256s
            )
            materialized_nodes.append(
                MaterializedNodeV1(
                    node_id=planned_node.node_id,
                    input_artifact_sha256=context.input_artifact.sha256,
                    project_artifact_sha256=project_artifact.sha256,
                    project_validation_receipt_sha256=project.receipt_sha256,
                    environment_receipt_sha256=environment_sha256,
                    invocation_sha256=invocation.invocation_sha256,
                    preflight_receipt_sha256=(
                        preflight.receipt_sha256 if previewed else ""
                    ),
                    state="previewed" if previewed else "compiled",
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
        self.materialized_workflows[
            workflow.materialized_sha256
        ] = workflow
        self.event_store.record_materialized_workflow(
            turn_id=turn_id, workflow=workflow
        )
        return workflow

    def completion_receipts_for_latest_preflight(self) -> tuple[str, ...]:
        """Return host-derived gates only for a green previewed preflight."""

        if not self.preflights:
            raise ContractError("no node preflight has been observed")
        latest = next(reversed(self.preflights.values()))
        if latest.plan_state != "previewed" or latest.critical_finding_sha256s:
            raise ContractError("latest node preflight is not green")
        completion = list(self._completion_sets[latest.receipt_sha256])
        if self.surface.profile == "command_compiled_approved_execution":
            if not self.execution_receipts:
                raise ContractError("approved workflow has not executed a node")
            approval = self.workflow_execution_approval
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
                preferred = tuple(
                    receipt
                    for receipt in matches
                    if receipt.receipt_sha256 in downstream_source_receipts
                ) or matches
                receipt = sorted(
                    preferred, key=lambda item: item.receipt_sha256
                )[-1]
                stage_receipts["quantity_extraction"][artifact_sha256] = (
                    receipt.receipt_sha256
                )
                selected.append(receipt.receipt_sha256)
                evidence_receipts.append(receipt.receipt_sha256)

        if "thermochemistry" in policy.required_stages:
            required_ids = set(
                policy.required_thermochemistry_quantity_ids
            )
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
                        or claims[
                            requirement.claim_id
                        ].source_receipt_sha256
                        in {
                            stage_receipts[requirement.source_kind][artifact]
                            for artifact in requirement.source_artifact_sha256s
                        }
                    )
                    and (
                        not requirement.source_expression_id
                        or claims[
                            requirement.claim_id
                        ].source_receipt_sha256
                        == expression_stage_receipts.get(
                            requirement.source_expression_id
                        )
                    )
                    and (
                        not requirement.source_selector
                        or self.quantity_extraction_bindings.get(
                            claims[
                                requirement.claim_id
                            ].source_receipt_sha256,
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
                sorted(matches, key=lambda item: item.record_sha256)[-1]
                .record_sha256
            )

        completion = build_analysis_completion_receipt(
            policy=policy,
            source_receipt_sha256s=tuple(sorted(set(selected))),
        )
        self.analysis_completion_receipts[
            completion.receipt_sha256
        ] = completion
        completion_record = canonical_data(completion)
        completion_record.pop("receipt_sha256")
        self.event_store.append(
            turn_id=turn_id,
            kind=EventKind.ANALYSIS_COMPLETION_EVALUATED.value,
            payload={
                "receipt_sha256": completion.receipt_sha256,
                "policy_sha256": completion.policy_sha256,
                "task_spec_sha256": completion.task_spec_sha256,
                "source_receipt_sha256s": (
                    completion.source_receipt_sha256s
                ),
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

    def _execute_approved_program_node(self, turn_id: str, values: dict) -> Any:
        """Resolve and execute one approved node without model-authored argv."""

        if self.surface.profile != "command_compiled_approved_execution":
            raise ContractError("execution tool is absent from this profile")
        node_id = values["node_id"]
        approval = self.workflow_execution_approval
        scientific_plan = self._execution_scientific_plan()
        frozen_approval = getattr(self, "frozen_workflow_approval", None)
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
            return {
                "execution": replayed,
                "idempotent_replay": True,
                "handoff": self.handoffs.get(node_id),
            }
        existing = self.execution_receipts.get(node_id)
        if existing is not None:
            raise ContractError(
                "process-local execution receipt lacks replay evidence"
            )
        approved_node = approval.node(node_id)
        if approved_node.jobtype == "hess":
            _validate_stationary_point_policy_binding(
                self.frozen_workflow_approval,
                self.stationary_point_policy,
                plan=scientific_plan,
                hessian_node_id=node_id,
                require_for_hessian=True,
            )
        invocation, context = self._latest_invocation_for_node(node_id)
        if context.project_artifact is None or context.project_validation is None:
            raise ContractError("execution requires a validated project artifact")
        if context.project_validation.status != "valid":
            raise ContractError("execution project loader gate is red")
        if context.project_validation.settings_sha256 != approved_node.settings_sha256:
            raise ContractError("effective project settings differ from approval")
        if context.project_artifact.sha256 != approved_node.project_artifact_sha256:
            raise ContractError("project bytes differ from workflow approval")
        if not self._invocation_has_green_preflight(invocation.invocation_sha256):
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
        )
        if self.frozen_workflow_approval is not None:
            future_rules = self.frozen_workflow_approval.producer_rules_for(
                node_id
            )
            expected_environments = (
                {item.environment_receipt_sha256 for item in future_rules}
                if future_rules
                else set(
                    self.frozen_workflow_approval.environment_receipt_sha256s
                )
            )
            if execution_invocation.environment_receipt_sha256 not in (
                expected_environments
            ):
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
                raise ContractError("terminal replay lacks an execution receipt")
            self.execution_receipts[node_id] = replayed
            return {
                "execution": replayed,
                "idempotent_replay": True,
                "handoff": self.handoffs.get(node_id),
            }
        command = [sys.executable, "-m", "chemsmart", *real_argv[1:]]
        environment = os.environ.copy()
        environment.update(self.execution_environment)
        source_root = str(Path(__file__).resolve().parents[2])
        current_pythonpath = environment.get("PYTHONPATH", "")
        environment["PYTHONPATH"] = (
            source_root
            if not current_pythonpath
            else source_root + os.pathsep + current_pythonpath
        )
        launch_ambiguous = False
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
                timeout_seconds=(
                    self.execution_resources.node_timeout_seconds
                ),
                memory_limit_mb=self.execution_resources.memory_gb * 1024.0,
                error_type=type(exc).__name__,
            )
            stdout_text = ""
            stderr_text = type(exc).__name__
        else:
            process_result = observe_process(
                process,
                timeout_seconds=(
                    self.execution_resources.node_timeout_seconds
                ),
                memory_limit_mb=self.execution_resources.memory_gb * 1024.0,
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
        outputs = self._execution_output_artifacts(node_id, node_workspace)
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
                if approved_node.jobtype == "hess"
                else ""
            ),
            approved_hessian_node_id=(
                node_id if approved_node.jobtype == "hess" else ""
            ),
            process_observation=process_observation,
        )
        result_validation_receipt = build_program_result_validation_receipt(
            invocation=execution_invocation,
            validator_id=evaluation.validator_id,
            validator_schema_version=(
                evaluation.validator_schema_version
            ),
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
        if receipt.validated and context.proposal.program == "pyscf":
            for edge in approval.producer_edges:
                if edge.producer_node_id != node_id:
                    continue
                result = next(
                    (item for item in outputs if item.kind == "pyscf_hdf5"),
                    None,
                )
                if result is None:
                    raise ContractError("validated PySCF OPT lacks an HDF5 result")
                geometry, observed = handoff_optimized_pyscf_geometry(
                    producer_receipt=receipt,
                    result_artifact=result,
                    producer_edge=edge,
                    approved_workspace=self.approved_workspace,
                    geometry_artifact_id=(
                        f"geometry.{edge.producer_node_id}-to-{edge.consumer_node_id}"
                    ),
                    expected_charge=context.scientific_identity.charge,
                    expected_multiplicity=context.scientific_identity.multiplicity,
                )
                identity = build_scientific_identity_binding(
                    task_spec_sha256=approval.task_spec_sha256,
                    geometry_artifact=geometry,
                    charge=observed.charge,
                    multiplicity=observed.multiplicity,
                )
                self.artifacts[geometry.artifact_id] = geometry
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
                        "optimized handoff lacks an exact scientific data edge"
                    )
                produced_handoffs.append(
                    {
                        "handoff": observed,
                        "artifact": geometry,
                        "scientific_identity": identity,
                    }
                )
                pending_data_edges.append(
                    (scientific_edge, edge, observed, identity)
                )
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
                        result_validation_receipt=(
                            result_validation_receipt
                        ),
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
        }

    def _execution_scientific_plan(self) -> ScientificWorkflowPlanV2:
        """Resolve the exact execution DAG without inferring tuple order."""

        approval = self.workflow_execution_approval
        frozen_approval = getattr(self, "frozen_workflow_approval", None)
        plans = getattr(self, "scientific_workflow_plans", None)
        if plans is None:
            plans = {}
            self.scientific_workflow_plans = plans
        if frozen_approval is not None:
            plan = plans.get(
                frozen_approval.plan_sha256
            )
            if plan is None:
                raise ContractError(
                    "frozen workflow approval has no registered scientific plan"
                )
            return plan
        plan = _scientific_plan_from_v1_approval(approval)
        plans.setdefault(plan.plan_sha256, plan)
        return plan

    def _latest_invocation_for_node(
        self, node_id: str
    ) -> tuple[CanonicalCommandInvocationV1, _CommandContext]:
        invocations = getattr(self, "invocations", {})
        contexts = getattr(self, "_command_contexts", {})
        for invocation in reversed(tuple(invocations.values())):
            context = contexts[invocation.invocation_sha256]
            if context.proposal.node_id == node_id:
                return invocation, context
        raise ContractError("node has no compiled command invocation")

    def _invocation_has_green_preflight(self, invocation_sha256: str) -> bool:
        return any(
            invocation_sha256 in self._completion_sets.get(
                receipt.receipt_sha256, ()
            )
            and receipt.plan_state == "previewed"
            and receipt.execution_ready
            and not receipt.critical_finding_sha256s
            for receipt in self.preflights.values()
        )

    def _real_execution_argv(
        self, invocation: CanonicalCommandInvocationV1
    ) -> tuple[str, ...]:
        if invocation.command_path[0] != "run":
            raise ContractError("local execution requires a run command")
        program = invocation.command_path[1]
        try:
            program_index = invocation.argv.index(program, 2)
        except ValueError as exc:
            raise ContractError("compiled argv does not contain its program") from exc
        resources = self.execution_resources
        root = ["chemsmart", "run", "--no-fake"]
        root.append("--no-scratch" if resources.scratch_policy == "none" else "--scratch")
        root.extend(("--num-cores", str(resources.cores)))
        root.extend(("--num-gpus", str(resources.gpu_count)))
        memory = (
            str(int(resources.memory_gb))
            if resources.memory_gb.is_integer()
            else str(resources.memory_gb)
        )
        root.extend(("--mem-gb", memory))
        if self.execution_server:
            root.extend(("--server", self.execution_server))
        return tuple(root + list(invocation.argv[program_index:]))

    def _execution_output_artifacts(
        self, node_id: str, workspace: Path
    ) -> tuple[TrustedArtifactRefV1, ...]:
        artifacts = []
        ordinal = 0
        for path in sorted(workspace.rglob("*")):
            if path.is_symlink():
                raise ContractError("execution emitted a symbolic link")
            if not path.is_file() or path.name.startswith("controller."):
                continue
            ordinal += 1
            kind = (
                "pyscf_hdf5"
                if path.suffix.lower() == ".h5"
                else "geometry_xyz"
                if path.suffix.lower() == ".xyz"
                else "json"
                if path.suffix.lower() == ".json"
                else "program_output"
            )
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
            findings.extend(
                _process_observation_findings(process_observation)
            )
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
                observation["runner_validation_delegation"] = (
                    "approved_stationary_point_policy"
                    if deferred_hessian_classification
                    else "none"
                )
                if not deferred_hessian_classification:
                    if (
                        run_receipt.get("scientifically_validated")
                        is not True
                    ):
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
                if (
                    expected_input_artifact is not None
                    and run_receipt.get("input_artifact_sha256")
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

                    expected_symbols, expected_positions = _pyscf_input_geometry(
                        expected_input_artifact
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
                    if result_validation.get("state") != "validated" and not (
                        result_validation.get("findings")
                    ):
                        findings.append(
                            "pyscf.result.validation_state_inconsistent"
                        )
                    from chemsmart.io.pyscf.output import read_pyscf_h5

                    result_spec, _provenance, _status, _results = read_pyscf_h5(
                        _current_artifact_path(
                            result, field_name="PySCF HDF5 program binding"
                        )
                    )
                    expected_engine = (
                        capability_environment_receipt.engine
                        if capability_environment_receipt is not None
                        else str((expected_settings or {}).get("engine") or "")
                    )
                    if result_spec.get("program") != "pyscf":
                        findings.append("pyscf.result.program_binding_mismatch")
                    if expected_engine and result_spec.get("engine") != (
                        expected_engine
                    ):
                        findings.append("pyscf.result.engine_binding_mismatch")
                except Exception as exc:
                    observation["pyscf_result_error_type"] = type(exc).__name__
                    findings.append("pyscf.result.unreadable")
            elif run_receipt is not None and observation.get(
                "runner_validation_delegation"
            ) == "approved_stationary_point_policy":
                findings.append("pyscf.result.policy_validation_unavailable")
        elif exit_status != 0:
            findings.append("execution.process.nonzero_or_unknown")
        elif program == "xtb":
            receipts: list[Path] = []
            for artifact in output_artifacts:
                if (
                    artifact.kind != "json"
                    or "xtb-result-receipt"
                    not in Path(artifact.path).name
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
        elif program not in {"pyscf", "xtb"}:
            findings.append("execution.program.validator_unavailable")
        normalized_findings = tuple(sorted(set(findings)))
        validator_schema_version = str(
            (
                observation.get("result_validation") or {}
            ).get("schema_version")
            or "chemsmart.generic-result-validation.v1"
        )
        environment_validation = observation.get("environment_validation") or {}
        return _ExecutionValidationEvaluation(
            validator_id=(
                "pyscf-result-validator"
                if program == "pyscf"
                else "xtb-result-validator"
                if program == "xtb"
                else "program-result-validator"
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
    def _validate_execution_outputs(**values: Any) -> tuple[
        bool, tuple[str, ...], tuple[str, ...]
    ]:
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
            quantity_ids=tuple(item.quantity_id for item in receipt.quantities),
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
        )
        self.thermochemistry_receipts[receipt.receipt_sha256] = receipt
        record = canonical_data(receipt)
        record.pop("receipt_sha256")
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
            receipt = self.quantity_extractions.get(receipt_sha256)
            if receipt is None:
                receipt = self.thermochemistry_receipts.get(receipt_sha256)
            if receipt is None:
                raise ContractError("quantity expression references an unknown receipt")
            quantity_id = str(item["quantity_id"])
            matches = tuple(
                quantity
                for quantity in receipt.quantities
                if quantity.quantity_id == quantity_id
            )
            if len(matches) != 1:
                raise ContractError(
                    "quantity expression input is absent or ambiguous in its receipt"
                )
            source = matches[0]
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
                input_ids=tuple(str(value) for value in item.get("input_ids", ())),
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
            )
            for item in values["nodes"]
        )
        request = QuantityExpressionRequestV1(
            schema_version="chemsmart.quantity-expression-request.v1",
            expression_id=str(values["expression_id"]),
            inputs=tuple(inputs),
            nodes=nodes,
            output_node_ids=tuple(str(item) for item in values["output_node_ids"]),
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

    def _record_analysis_claims(self, turn_id: str, values: dict) -> Any:
        """Render reportable values from exact typed receipt outputs."""

        task_spec_sha256 = str(values["task_spec_sha256"])
        if task_spec_sha256 not in self.task_spec_sha256s:
            raise ContractError("analysis claims target an unknown task spec")
        registries = (
            ("quantity_extraction", self.quantity_extractions, "quantities"),
            ("thermochemistry", self.thermochemistry_receipts, "quantities"),
            ("quantity_expression", self.quantity_expression_receipts, "outputs"),
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
        inspection = inspect_command(
            invocation, live_schema=self.live_schema
        )
        self.invocations[invocation.invocation_sha256] = invocation
        self.command_inspections[inspection.receipt_sha256] = inspection
        self._command_contexts[invocation.invocation_sha256] = context
        self._emit(
            turn_id,
            EventKind.COMMAND_COMPILED,
            invocation.invocation_sha256,
            status=invocation.status,
            input_sha256=invocation.input_sha256,
            project_sha256=invocation.project_sha256,
            scientific_identity_sha256=(
                invocation.scientific_identity_sha256
            ),
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
                detail = f"no {label} is bound yet"
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
    missing = sorted(required.difference(arguments))
    if missing:
        raise ContractError("tool arguments missing: " + ", ".join(missing))
    properties = schema.get("properties", {})
    unknown = sorted(set(arguments).difference(properties))
    if unknown:
        raise ContractError("tool arguments not allowed: " + ", ".join(unknown))
    for name, value in arguments.items():
        _validate_json_value(name, value, properties[name])


def _validate_json_value(name: str, value: Any, schema: Mapping[str, Any]) -> None:
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
            raise ContractError(
                f"tool argument {name} must match exactly one allowed shape"
            )
        return
    expected = schema.get("type")
    type_ok = {
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
    }.get(expected, True)
    if not type_ok:
        raise ContractError(f"tool argument {name} has the wrong type")
    if "enum" in schema and value not in schema["enum"]:
        raise ContractError(f"tool argument {name} is outside its enum")
    if isinstance(value, str) and schema.get("pattern"):
        if re.fullmatch(str(schema["pattern"]), value) is None:
            raise ContractError(f"tool argument {name} does not match its pattern")
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        if "minimum" in schema and value < float(schema["minimum"]):
            raise ContractError(f"tool argument {name} is below its minimum")
        if "maximum" in schema and value > float(schema["maximum"]):
            raise ContractError(f"tool argument {name} is above its maximum")
        if "exclusiveMinimum" in schema and value <= float(
            schema["exclusiveMinimum"]
        ):
            raise ContractError(
                f"tool argument {name} is not above its exclusive minimum"
            )
        if "exclusiveMaximum" in schema and value >= float(
            schema["exclusiveMaximum"]
        ):
            raise ContractError(
                f"tool argument {name} is not below its exclusive maximum"
            )
    if isinstance(value, list):
        if "minItems" in schema and len(value) < int(schema["minItems"]):
            raise ContractError(f"tool argument {name} has too few items")
        if "maxItems" in schema and len(value) > int(schema["maxItems"]):
            raise ContractError(f"tool argument {name} has too many items")
    if isinstance(value, list) and isinstance(schema.get("items"), Mapping):
        for item in value:
            _validate_json_value(name + "[]", item, schema["items"])
    if isinstance(value, dict):
        properties = schema.get("properties", {})
        additional = schema.get("additionalProperties")
        required = set(schema.get("required", ()))
        missing = sorted(required.difference(value))
        if missing:
            raise ContractError(
                f"tool argument {name} missing: " + ", ".join(missing)
            )
        if schema.get("additionalProperties") is False:
            unknown = sorted(set(value).difference(properties))
            if unknown:
                raise ContractError(
                    f"tool argument {name} not allowed: " + ", ".join(unknown)
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
                validator_receipt_sha256s=(
                    receipt.validator_receipt_sha256s
                ),
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
        receipt.started_at for receipt in receipts.values() if receipt.started_at
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
            if workflow_state in {"validated", "failed"}
            and finished_values
            else ""
        ),
    }
    return WorkflowRunStateV1(
        **body, run_state_sha256=canonical_sha256(body)
    )


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
