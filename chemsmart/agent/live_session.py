"""Plan-first live ChemSmart Agent sessions.

This module is the small composition root behind ``chemsmart agent plan``
and the terminal interface.  It binds a provider session to exact,
pre-existing XYZ artifacts, the live Click schema, observed program
conformance, and a private Runtime V2 event stream.  It never generates
coordinates, engine inputs, shell commands, or scientific readiness
decisions.

A planning session ends at ``waiting_for_approval``.  Real execution is the
provider-free ``chemsmart agent run`` over a one-shot approval bundle, or the
single ``/approve`` action inside the terminal interface.
"""

from __future__ import annotations

import hashlib
import inspect
import json
import logging
import math
import os
import re
import shutil
import sys
import uuid
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import TYPE_CHECKING, Any, Callable, Iterable, Mapping

if TYPE_CHECKING:
    from chemsmart.agent.scientific_toolchain import (
        ScientificToolchainPlanV1,
    )

from chemsmart.agent._contracts import (
    AuxiliaryArtifactBindingV1,
    ContractError,
    TrustedArtifactRefV1,
    canonical_data,
    canonical_json,
    canonical_sha256,
    file_sha256,
    require_sha256,
)
from chemsmart.agent.analysis_completion import load_analysis_completion_policy
from chemsmart.agent.api_access import (
    DEFAULT_KEY_LABELS,
    PROVIDER_KEY_LABEL_TOKENS,
    load_secret_lease,
    normalize_key_label,
)
from chemsmart.agent.bootstrap import (
    bootstrap_program_conformance,
    probe_python_compute_environment,
)
from chemsmart.agent.capabilities import (
    EnvironmentTargetV1,
    ProgramCapabilityRegistryV1,
    ProgramComponentConformanceReceiptV1,
    TrustedComputeEnvironmentReceiptV1,
    build_program_component_conformance_receipt,
    load_program_capabilities,
)
from chemsmart.agent.cli_schema import build_live_click_schema
from chemsmart.agent.commands import build_scientific_identity_binding
from chemsmart.agent.dependency_context import (
    ContextSelectionReceiptV1,
    TaskDependencyContextV2,
    build_dependency_context_public_projection,
)
from chemsmart.agent.execution import (
    ApprovedNodeBindingV1,
    ExecutionResourceSpecV1,
    FrozenMaterializedNodePreviewV1,
    FrozenProducerEdgeRuleV1,
    FrozenWorkflowApprovalV1,
    ProducerEdgeRuleV1,
    WorkflowApprovalRequestV1,
    WorkflowEnvironmentBindingV1,
    WorkflowExecutionApprovalBundleV1,
    WorkflowExecutionApprovalV1,
    WorkflowExecutionNodeReviewV1,
    WorkflowExecutionReviewV1,
    WorkflowReviewResolutionV1,
    workflow_execution_approval_bundle_json,
    workflow_execution_review_json,
)
from chemsmart.agent.execution_envelope import (
    BoundedExecutionEnvelopeV1,
    load_bounded_execution_envelope,
)
from chemsmart.agent.identity import (
    ApprovedMolecularIdentityV1,
    ApprovedMolecularInputV1,
    validate_identity_for_geometry,
)
from chemsmart.agent.knowledge_packs import (
    BUILTIN_PROGRAM_PACKS,
    activate_program_knowledge,
    skills_for_activation,
)
from chemsmart.agent.projects import project_document, render_project_yaml
from chemsmart.agent.provider_config import (
    AgentProviderProfileV1,
    load_agent_provider_selection,
)
from chemsmart.agent.request_context import (
    ProviderNetworkBudgetV1,
    RequestContextProvenanceV1,
    build_request_context_provenance,
)
from chemsmart.agent.runtime.contracts import (
    ResourceBudgetV1,
    TaskEnvelopeV1,
    TaskPhase,
)
from chemsmart.agent.runtime.event_store import RuntimeEventStore
from chemsmart.agent.runtime.transport import ProviderTurnDeadlinesV1
from chemsmart.agent.services.unified_session import UnifiedSessionRunner
from chemsmart.agent.skills import (
    SkillDocumentV1,
    resolve_skills,
    skills_enabled,
)
from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1
from chemsmart.agent.tool_specs import (
    build_approved_execution_tool_surface,
    build_command_compiled_tool_surface,
)
from chemsmart.agent.workflows import (
    MaterializedNodeV1,
    MaterializedWorkflowV1,
    ScientificWorkflowEdgeV2,
    ScientificWorkflowNodeV2,
    ScientificWorkflowPlanV2,
    StationaryPointValidationPolicyV1,
)
from chemsmart.analysis.result_quantities import (
    QuantityExtractionError,
    validate_pyscf_analysis_artifact,
)

_SESSION_WALL_TIME_SECONDS = 90 * 60
_MAX_TOOL_CALLS = 256


def _declared_pyscf_interpreter() -> Path | None:
    """Return the interpreter the server YAML declares for PySCF, if any.

    The server YAML is this file's stated canonical source for which programs
    an installation controls, but PySCF was resolved from ``sys.executable``
    instead.  On any host that installs PySCF in its own environment -- the
    ordinary arrangement, since it pins numpy and h5py -- that answers about
    the controller interpreter, the probe reports the required packages
    missing, and the engine is declared unavailable while it is installed and
    working.  Reading the operator's declaration is honouring it; substituting
    the controller environment for it is the implicit replacement the project
    forbids.
    """

    try:
        from chemsmart.io.yaml import YAMLFile
        from chemsmart.settings.user import CHEMSMARTUserSettings

        settings = CHEMSMARTUserSettings()
        available = list(settings.all_available_servers or ())
        preferred = os.environ.get("CHEMSMART_AGENT_SERVER") or "local"
        name = (
            preferred
            if preferred in available
            else (available[0] if available else "")
        )
        if not name:
            return None
        path = Path(settings.user_server_dir) / f"{name}.yaml"
        if not path.is_file():
            return None
        block = YAMLFile(filename=str(path)).yaml_contents_dict.get("PYSCF")
        folder = str((block or {}).get("EXEFOLDER") or "").strip()
        if not folder:
            return None
        candidate = Path(folder).expanduser() / "python"
        if candidate.is_file():
            return candidate.resolve()
    except Exception:  # noqa: BLE001 - fall back rather than fail to start
        return None
    return None


def _resolve_pyscf_interpreter() -> Path:
    """Explicit override first, then the declaration, then this process."""

    override = os.environ.get("CHEMSMART_PYSCF_INTERPRETER")
    if override:
        return Path(override).expanduser().resolve()
    declared = _declared_pyscf_interpreter()
    if declared is not None:
        return declared
    return Path(sys.executable).expanduser().resolve()


_PYSCF_INTERPRETER = _resolve_pyscf_interpreter()
logger = logging.getLogger(__name__)
_PRIVATE_ROOT_NAME = ".chemsmart-agent"


@dataclass(frozen=True)
class _XyzObservation:
    artifact: TrustedArtifactRefV1
    atom_count: int
    symbols: tuple[str, ...]

    def public_record(self) -> dict[str, Any]:
        return {
            "artifact_id": self.artifact.artifact_id,
            "artifact_class": self.artifact.kind,
            "sha256": self.artifact.sha256,
            "size_bytes": self.artifact.size_bytes,
            "atom_count": self.atom_count,
            "symbols": self.symbols,
            "provenance_status": "workspace_exact_user_approved",
        }


@dataclass(frozen=True)
class _PySCFResultObservation:
    artifact: TrustedArtifactRefV1
    jobtype: str
    method: str
    applied_method: str
    basis: str
    engine: str
    charge: int
    multiplicity: int
    project_yaml_sha256: str
    input_artifact_sha256: str
    validation_receipt_sha256: str
    scientific_validation_state: str

    @property
    def program(self) -> str:
        return "pyscf"

    def public_record(self) -> dict[str, Any]:
        return {
            "artifact_id": self.artifact.artifact_id,
            "artifact_class": self.artifact.kind,
            "sha256": self.artifact.sha256,
            "size_bytes": self.artifact.size_bytes,
            "program": "pyscf",
            "jobtype": self.jobtype,
            "requested_method": self.method,
            "applied_method": self.applied_method,
            "basis": self.basis,
            "engine": self.engine,
            "charge": self.charge,
            "multiplicity": self.multiplicity,
            "project_yaml_sha256": self.project_yaml_sha256,
            "input_artifact_sha256": self.input_artifact_sha256,
            "validation_receipt_sha256": self.validation_receipt_sha256,
            "scientific_validation_state": self.scientific_validation_state,
            "functional_evidence_ref": (
                "result_functional_resolution:"
                f"{self.validation_receipt_sha256}"
            ),
            "provenance_status": "workspace_exact_structured_result",
        }


@dataclass(frozen=True)
class _LoggedResultObservation:
    """One validated native-output artifact exposed through a shared reader."""

    artifact: TrustedArtifactRefV1
    program: str
    jobtype: str
    method: str
    basis: str | None
    engine: str
    charge: int
    multiplicity: int
    project_yaml_sha256: str
    input_artifact_sha256: str
    validation_receipt_sha256: str
    scientific_validation_state: str
    provenance_status: str = "workspace_exact_validated_native_result"
    provenance_limitations: tuple[str, ...] = ()
    #: The dispersion correction the run actually applied, when the native
    #: output names one.  A functional and its dispersion correction are not
    #: the same method, and for some observables the correction is the larger
    #: part of the answer, so reporting only the functional understates what
    #: was run rather than merely abbreviating it.
    dispersion: str | None = None

    def public_record(self) -> dict[str, Any]:
        return {
            "artifact_id": self.artifact.artifact_id,
            "artifact_class": self.artifact.kind,
            "sha256": self.artifact.sha256,
            "size_bytes": self.artifact.size_bytes,
            "program": self.program,
            "jobtype": self.jobtype,
            "requested_method": self.method,
            "applied_method": self.method,
            "dispersion": self.dispersion,
            "basis": self.basis,
            "engine": self.engine,
            "charge": self.charge,
            "multiplicity": self.multiplicity,
            "project_yaml_sha256": self.project_yaml_sha256,
            "input_artifact_sha256": self.input_artifact_sha256,
            "validation_receipt_sha256": self.validation_receipt_sha256,
            "scientific_validation_state": self.scientific_validation_state,
            "provenance_status": self.provenance_status,
            "provenance_limitations": self.provenance_limitations,
        }


_ResultObservation = _PySCFResultObservation | _LoggedResultObservation


@dataclass(frozen=True)
class LiveAgentSessionResultV1:
    """Path-free public projection of one local live session."""

    schema_version: str
    session_id: str
    task_spec_sha256: str
    terminal_state: str
    execution_requested: bool
    execution_profile_status: str
    final_text: str
    artifact_records: tuple[dict[str, Any], ...]
    conformance_records: tuple[dict[str, Any], ...]
    public_transcript: tuple[dict[str, Any], ...]
    successful_tool_calls: int
    failed_tool_calls: int
    execution_review: dict[str, Any]
    event_stream_head_sha256: str
    result_sha256: str
    prepared_execution: WorkflowExecutionReviewV1 | None = None

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.live-agent-session-result.v1":
            raise ContractError("unsupported live session result schema")
        if self.terminal_state not in {
            "complete",
            "planned",
            "failed",
            "blocked",
            "cancelled",
            "waiting_for_approval",
        }:
            raise ContractError("invalid live session terminal state")
        body = self._body()
        if self.result_sha256 != canonical_sha256(body):
            raise ContractError("live session result digest mismatch")

    def _body(self) -> dict[str, Any]:
        return {
            key: value
            for key, value in asdict(self).items()
            if key not in {"result_sha256", "prepared_execution"}
        }

    def public_summary_json(self) -> str:
        """Return the exact visible result; no local path is included."""

        return canonical_json(
            {**self._body(), "result_sha256": self.result_sha256}
        )


def run_live_agent_session(
    *,
    task: str,
    provider: str | None,
    provider_config_file: str | Path | None = None,
    secret_file: str | Path | None = None,
    workspace: str | Path,
    execution_enabled: bool,
    approval_file: str | Path | None,
    execution_envelope_file: str | Path | None = None,
    analysis_completion_file: str | Path | None = None,
    approved_molecular_identity: ApprovedMolecularIdentityV1 | None = None,
    approved_molecular_identities: Iterable[ApprovedMolecularIdentityV1] = (),
    approved_molecular_inputs: Iterable[ApprovedMolecularInputV1] = (),
    dependency_context: TaskDependencyContextV2 | None = None,
    dependency_context_selection_receipt: (
        ContextSelectionReceiptV1 | None
    ) = None,
    dependency_public_records: Mapping[str, Mapping[str, Any]] | None = None,
    review_file: str | Path | None = None,
    on_run_directory: Callable[[Path], None] | None = None,
    should_stop: Callable[[], bool] | None = None,
) -> LiveAgentSessionResultV1:
    """Run one agent.yaml-selected session over exact workspace artifacts.

    A real provider request is made only after the workspace, coordinate
    artifacts, live schema, program previews, and provider budget have been
    bound.  The private run directory contains the append-only event stream
    and sanitized public transcript; it is intentionally outside the source
    tree when the CLI is used as documented.
    """

    selection = load_agent_provider_selection(
        provider_config_file,
        requested_profile=str(provider).strip() if provider else None,
    )
    profile = selection.active_profile
    normalized_provider = profile.provider
    if (dependency_context is None) != (
        dependency_context_selection_receipt is None
    ):
        raise ContractError(
            "dependency context and selection receipt must be supplied together"
        )
    if dependency_context is None:
        if dependency_public_records:
            raise ContractError(
                "dependency public records require a selected dependency context"
            )
        dependency_context_public_projection: dict[str, Any] = {}
    else:
        if dependency_public_records is None:
            raise ContractError(
                "selected dependency context requires exact public record payloads"
            )
        dependency_context_public_projection = (
            build_dependency_context_public_projection(
                context=dependency_context,
                selection_receipt=dependency_context_selection_receipt,
                records=dependency_public_records,
            )
        )
    task = str(task).strip()
    if not task:
        raise ContractError("live agent task must not be empty")
    if approval_file is not None:
        raise ContractError(
            "live provider sessions cannot consume approvals; use the "
            "model-free 'chemsmart agent run' command"
        )
    bounded_envelope = (
        load_bounded_execution_envelope(execution_envelope_file)
        if execution_envelope_file is not None
        else None
    )
    bounded_review_requested = bounded_envelope is not None
    if execution_enabled:
        raise ContractError(
            "live provider sessions are planning and safe-preview only"
        )
    if review_file is not None and bounded_envelope is None:
        raise ContractError(
            "an execution review requires --execution-envelope"
        )
    workspace_path = _validated_workspace(workspace)
    observations = _scan_xyz_artifacts(workspace_path)
    result_observations = _scan_result_artifacts(workspace_path)
    molecular_inputs = _coerce_approved_molecular_inputs(
        approved_molecular_inputs
    )
    direct_identities = tuple(approved_molecular_identities)
    if molecular_inputs and (
        approved_molecular_identity is not None or direct_identities
    ):
        raise ContractError(
            "use approved molecular inputs or direct molecular identities, "
            "not both"
        )
    identities = (
        tuple(item.molecular_identity for item in molecular_inputs)
        if molecular_inputs
        else _coerce_approved_identities(
            approved_molecular_identity, direct_identities
        )
    )
    identity_records = _validated_identity_records(observations, identities)
    task_spec_sha256 = _task_spec_sha256(
        task,
        observations,
        identities,
        result_observations=result_observations,
        approved_molecular_inputs=molecular_inputs,
    )
    prebound_scientific_identities = _approved_input_state_bindings(
        molecular_inputs,
        observations=observations,
        task_spec_sha256=task_spec_sha256,
    )
    analysis_completion_policy = (
        load_analysis_completion_policy(
            analysis_completion_file,
            task_spec_sha256=task_spec_sha256,
        )
        if analysis_completion_file is not None
        else None
    )
    analysis_only_session = bool(
        not observations and result_observations and not execution_enabled
    )
    session_id = _session_id(task_spec_sha256)
    run_directory = _private_run_directory(workspace_path, session_id)

    if not observations and not analysis_only_session:
        return _local_result(
            session_id=session_id,
            task_spec_sha256=task_spec_sha256,
            terminal_state="blocked",
            execution_requested=execution_enabled,
            execution_profile_status="not_started",
            final_text=(
                "No exact XYZ or supported completed-result artifact is "
                "present in the approved workspace. Add user-approved input "
                "or result data; coordinates and results were not generated."
            ),
        )

    if analysis_only_session:
        registry = load_program_capabilities()
        live_schema = build_live_click_schema()
        conformance = ()
        environment_targets, compute_receipts, environment_records = (
            _observe_environments()
        )
        conformance_records = tuple(
            sorted(environment_records, key=_record_sort_key)
        )
    else:
        registry = load_program_capabilities()
        live_schema = build_live_click_schema()
        conformance, conformance_records = _bootstrap_conformance(
            run_directory=run_directory,
            input_artifact=observations[0].artifact,
            registry_sha256=registry.registry_sha256,
            live_schema=live_schema,
            resources=(
                bounded_envelope.resources
                if bounded_envelope is not None
                else None
            ),
        )
        environment_targets, compute_receipts, environment_records = (
            _observe_environments()
        )
        conformance_records = tuple(
            sorted(
                (*conformance_records, *environment_records),
                key=_record_sort_key,
            )
        )

    use_execution_surface = False
    approved_project_records: tuple[dict[str, Any], ...] = ()
    approved_workflow_record: dict[str, Any] = {}
    surface = (
        build_approved_execution_tool_surface(registry)
        if use_execution_surface
        else build_command_compiled_tool_surface(registry)
    )

    event_store = RuntimeEventStore(
        run_directory / "events.jsonl", session_id=session_id
    )
    if on_run_directory is not None:
        # A view-only announcement: an embedding interface may tail the
        # append-only event stream from the moment it exists.  It grants no
        # authority and changes no session behavior.
        on_run_directory(run_directory)
    preview_server_path = _ensure_preview_server(
        run_directory,
        resources=(
            bounded_envelope.resources
            if bounded_envelope is not None
            else None
        ),
    )
    host_kwargs: dict[str, Any] = {
        "event_store": event_store,
        "artifacts": {
            item.artifact.artifact_id: item.artifact
            for item in (*observations, *result_observations)
        },
        "environment_targets": environment_targets,
        "compute_environment_receipts": compute_receipts,
        "component_conformance_receipts": conformance,
        "tool_surface": surface,
        "registry": registry,
        "live_schema": live_schema,
        "task_spec_sha256s": (task_spec_sha256,),
        "approved_molecular_identities": (
            {identity.identity_sha256: identity for identity in identities}
        ),
        "scientific_identities": prebound_scientific_identities,
        # Preview candidates stay inside this session's private workspace.
        # The execution composition below replaces this with the exact
        # user-approved workflow workspace.
        "approved_workspace": run_directory,
        "preview_server": str(preview_server_path),
        "result_functional_evidence": {
            item.validation_receipt_sha256: item.public_record()
            for item in result_observations
            if item.program == "pyscf"
        },
        "analysis_completion_policy": analysis_completion_policy,
        "dependency_context": dependency_context,
        "dependency_context_selection_receipt": (
            dependency_context_selection_receipt
        ),
    }
    if bounded_envelope is not None:
        _validate_bounded_envelope_against_registry(
            bounded_envelope, registry=registry
        )
        host_kwargs["execution_resources"] = bounded_envelope.resources
        host_kwargs["bounded_execution_envelope"] = bounded_envelope
    if use_execution_surface and approval_file is not None:
        execution_inputs = _execution_composition_inputs(
            host_type=CommandCompiledToolHostV1,
            workspace=workspace_path,
            run_directory=run_directory,
            approval_file=Path(approval_file),
            task_spec_sha256=task_spec_sha256,
        )
        approved_plan = execution_inputs.pop("approved_scientific_plan")
        approved_workflow_record = _public_workflow_approval(
            execution_inputs["workflow_execution_approval"],
            approved_plan=approved_plan,
        )
        if approved_plan is not None:
            # Execution resolves the exact reviewed DAG by its frozen digest.
            host_kwargs["scientific_workflow_plan"] = approved_plan
        approved_materialization = execution_inputs.pop(
            "approved_materialized_workflow"
        )
        if approved_materialization is not None:
            # Execution resolves the reviewed materialization by its digest.
            host_kwargs["materialized_workflow"] = approved_materialization
        approved_projects = execution_inputs.pop("approved_project_artifacts")
        host_kwargs["artifacts"].update(
            {item.artifact_id: item for item in approved_projects}
        )
        approved_project_records = tuple(
            {
                "artifact_id": item.artifact_id,
                "artifact_class": item.kind,
                "sha256": item.sha256,
                "size_bytes": item.size_bytes,
                "approval_status": "workflow_bound",
            }
            for item in approved_projects
        )
        host_kwargs.update(execution_inputs)
    elif use_execution_surface and bounded_envelope is not None:
        _validate_bounded_envelope_against_registry(
            bounded_envelope, registry=registry
        )
        host_kwargs.update(
            _bounded_execution_composition_inputs(
                workspace=workspace_path,
                run_directory=run_directory,
                envelope=bounded_envelope,
            )
        )
    if use_execution_surface:
        provider_secret_labels = {
            profile.api_key_env,
            *(
                label
                for labels in DEFAULT_KEY_LABELS.values()
                for label in labels
            ),
        }
        host_kwargs["execution_environment_remove"] = tuple(
            sorted(provider_secret_labels)
        )
    host = CommandCompiledToolHostV1(**host_kwargs)

    context = _public_context(
        task=task,
        task_spec_sha256=task_spec_sha256,
        observations=observations,
        result_observations=result_observations,
        conformance_records=conformance_records,
        registry_sha256=registry.registry_sha256,
        live_schema_sha256=live_schema.schema_sha256,
        execution_requested=execution_enabled,
        execution_available=use_execution_surface,
        execution_review_requested=bounded_review_requested,
        bounded_execution_record=(
            bounded_envelope.public_record()
            if bounded_envelope is not None
            else None
        ),
        approved_project_records=approved_project_records,
        approved_workflow_record=approved_workflow_record,
        provider_record=_provider_public_record(
            profile=profile,
            fallback_profiles=selection.fallback_profiles,
        ),
        approved_identity_records=identity_records,
        approved_input_records=tuple(
            item.public_record() for item in molecular_inputs
        ),
        analysis_completion_record=(
            analysis_completion_policy.public_record()
            if analysis_completion_policy is not None
            else None
        ),
    )
    if dependency_context is not None:
        context = {
            **context,
            **dependency_context_public_projection,
        }
    base_messages = _coordinator_base_messages(
        context=context,
        approved_workflow=approved_workflow_record,
        bounded_review_requested=bounded_review_requested,
        task=task,
    )
    messages = [
        base_messages[0],
        {"role": "user", "content": canonical_json(context)},
    ]
    provider_budget = _provider_budget(
        profile,
        max_concurrency=1,
        wall_time_seconds=(
            bounded_envelope.episode_wall_time_seconds
            if bounded_envelope is not None
            else _SESSION_WALL_TIME_SECONDS
        ),
    )
    request_context = _request_context(
        messages=messages,
        task_spec_sha256=task_spec_sha256,
        tool_schema_sha256=surface.tool_schema_sha256,
        provider_budget=provider_budget,
        provider_profile_sha256=profile.profile_sha256,
        execution_requested=execution_enabled,
        dependency_context=dependency_context,
    )
    envelope = _task_envelope(
        session_id=session_id,
        messages=messages,
        task_spec_sha256=task_spec_sha256,
        workspace_records=tuple(
            item.public_record()
            for item in (*observations, *result_observations)
        ),
        tool_schema_sha256=surface.tool_schema_sha256,
        execution_enabled=use_execution_surface,
        provider_profile=profile,
        wall_time_seconds=(
            bounded_envelope.episode_wall_time_seconds
            if bounded_envelope is not None
            else _SESSION_WALL_TIME_SECONDS
        ),
        chemistry_engine_calls=0,
    )
    lease = load_secret_lease(
        provider=normalized_provider,
        path=secret_file,
        # The profile already says which key it bills; selecting a different
        # one here is how a second key for the same provider gets charged by
        # accident.
        label=profile.api_key_env,
        ttl_seconds=(
            bounded_envelope.episode_wall_time_seconds
            if bounded_envelope is not None
            else _SESSION_WALL_TIME_SECONDS
        )
        + 60,
    )
    loop_result = UnifiedSessionRunner(
        host=host,
        event_store=event_store,
        credential_lease=lease,
        provider_config=profile.runtime_config(),
    ).run(
        messages=messages,
        envelope=envelope,
        request_context=request_context,
        provider_budget=provider_budget,
        should_stop=should_stop,
    )

    execution_review_record: dict[str, Any] = {}
    prepared_execution: WorkflowExecutionReviewV1 | None = None
    # A resource envelope constrains real engine work; it must not turn an
    # analysis-only session into an execution-approval request.  The normal
    # ChemSmart result readers and typed analysis DAG are sufficient when the
    # scientific plan contains no calculation nodes.
    calculation_plans = tuple(
        plan for plan in host.scientific_workflow_plans.values() if plan.nodes
    )
    latest_calculation_plan = (
        calculation_plans[-1] if calculation_plans else None
    )
    execution_ineligible_nodes: tuple[str, ...] = ()
    if bounded_envelope is not None and latest_calculation_plan is not None:
        ineligible = []
        for node in latest_calculation_plan.nodes:
            reason = host.execution_review_ineligibility_reason(
                plan=latest_calculation_plan,
                planned_node=node,
            )
            if reason:
                ineligible.append(
                    f"{node.node_id} ({node.program}/{node.engine}/"
                    f"{node.stage}: {reason})"
                )
        execution_ineligible_nodes = tuple(ineligible)
    if (
        bounded_envelope is not None
        and latest_calculation_plan is not None
        and not execution_ineligible_nodes
        # A cancelled session must never mint pending authority: the human
        # withdrew the request, so no review packet is built and the result
        # keeps its cancelled state.
        and loop_result.terminal_state != "cancelled"
        # A planned but unmaterialised workflow is a session that stopped
        # early, not a broken contract.  Report it as the preview-only run it
        # is rather than losing the transcript to an unhandled refusal.
        and host.bounded_review_is_materialized()
    ):
        try:
            review = host.build_execution_review(workspace=workspace_path)
        except ContractError as exc:
            # The refusal itself is correct -- a plan that exceeds the engine
            # budget, or reaches a stage this release cannot execute, must not
            # become an approvable packet.  Losing the whole session to it is
            # not: the model has planned, compiled and previewed real work, and
            # an operator who sees only a traceback cannot tell whether to
            # widen the envelope or narrow the plan.  Record the reason and let
            # the existing preview-only status carry it.
            logger.warning("execution review refused: %s", exc)
            execution_ineligible_nodes = execution_ineligible_nodes + (
                f"workflow ({latest_calculation_plan.workflow_id}: {exc})",
            )
            review = None
        if review is not None:
            prepared_execution = review
            if review_file is not None:
                write_workflow_execution_review(review, review_file)
            execution_review_record = {
                "schema_version": "chemsmart.public-workflow-execution-review.v1",
                "status": review.status,
                "workflow_id": review.scientific_plan.workflow_id,
                "plan_sha256": review.scientific_plan.plan_sha256,
                "materialized_workflow_sha256": (
                    review.materialized_workflow.materialized_sha256
                ),
                "resource_sha256": review.execution_resources.resource_sha256,
                "review_sha256": review.review_sha256,
                "non_executable_node_ids": review.non_executable_node_ids,
                "next_action": "review the displayed ChemSmart workflow and approve it explicitly",
            }

    execution_status = (
        "review_packet_ready"
        if execution_review_record
        else (
            "preview_only_not_execution_eligible"
            if execution_ineligible_nodes
            else (
                "preview_only_resource_bounds_loaded"
                if bounded_envelope is not None
                else "not_requested"
            )
        )
    )
    terminal_state = loop_result.terminal_state
    final_text = loop_result.final_text
    if execution_review_record:
        terminal_state = "waiting_for_approval"
        suffix = (
            " Planning, command compilation, and safe preview are complete. "
            "An inert exact review packet was written for human approval. "
            "No chemistry engine was launched."
        )
        non_executable = tuple(
            execution_review_record.get("non_executable_node_ids") or ()
        )
        if non_executable:
            suffix += (
                " The following planned stages remain in the workflow as "
                "declared scientific intent but are not executable in this "
                "release, so they are displayed with the review and will not "
                "run: " + ", ".join(non_executable) + "."
            )
        final_text = (final_text.rstrip() + suffix).strip()
    elif execution_ineligible_nodes:
        suffix = (
            " Planning, project-YAML validation, command compilation, and "
            "safe preview remain valid. No chemistry engine was launched. "
            "The following nodes are not eligible for Agent execution on "
            "this release or fall outside the execution envelope: "
            + ", ".join(execution_ineligible_nodes)
            + "."
        )
        final_text = (final_text.rstrip() + suffix).strip()
    body = {
        "schema_version": "chemsmart.live-agent-session-result.v1",
        "session_id": session_id,
        "task_spec_sha256": task_spec_sha256,
        "terminal_state": terminal_state,
        "execution_requested": bool(execution_enabled),
        "execution_profile_status": execution_status,
        "final_text": final_text,
        "artifact_records": (
            *(item.public_record() for item in observations),
            *(item.public_record() for item in result_observations),
            *identity_records,
        ),
        "conformance_records": conformance_records,
        "public_transcript": loop_result.public_transcript,
        "successful_tool_calls": loop_result.successful_tool_calls,
        "failed_tool_calls": loop_result.failed_tool_calls,
        "execution_review": execution_review_record,
        "event_stream_head_sha256": loop_result.event_stream_head_sha256,
    }
    return LiveAgentSessionResultV1(
        **body,
        result_sha256=canonical_sha256(body),
        prepared_execution=prepared_execution,
    )


def _provider_public_record(
    *,
    profile: AgentProviderProfileV1,
    fallback_profiles: Iterable[AgentProviderProfileV1],
) -> dict[str, Any]:
    record = {
        "profile_name": profile.profile_name,
        "provider": profile.provider,
        "model": profile.model,
        "endpoint_origin": profile.endpoint,
        "reasoning_effort": profile.reasoning_effort,
        "preserve_thinking": profile.preserve_thinking,
        "profile_sha256": profile.profile_sha256,
        "fallback_profiles": tuple(
            item.profile_name for item in fallback_profiles
        ),
        "fallback_policy": "explicit_attributed_provider_unavailability_only",
    }
    deadlines = getattr(profile, "transport_deadlines", None)
    record["transport_deadlines"] = (
        deadlines.configuration_record()
        if deadlines is not None
        else ProviderTurnDeadlinesV1().configuration_record()
    )
    record["transport_deadline_source"] = (
        "agent_provider_profile_v2"
        if deadlines is not None
        else "runtime_immutable_defaults"
    )
    return record


def cross_program_skill_ids() -> tuple[str, ...]:
    """Skills every program pack carries, and which therefore belong to none.

    Derived from the packs rather than listed separately: a skill that all of
    them advertise is program-neutral by construction, so adding or removing
    one cannot leave a second list behind to drift.
    """

    advertised = [set(pack.skill_ids) for pack in BUILTIN_PROGRAM_PACKS]
    if not advertised:
        return ()
    return tuple(sorted(set.intersection(*advertised)))


def activated_skill_documents(
    task: str,
) -> tuple[tuple[str, ...], tuple[SkillDocumentV1, ...]]:
    """Resolve the advisory skills available for a task.

    Returns the activated pack digests and the resolved documents.

    Program-specific advice stays text-gated by each pack's own
    ``activation_terms``.  Program-neutral skills do not, and that distinction
    is the whole point: this function's own docstring used to promise that "a
    cross-program conventions skill is reachable from any request" while
    delivering the opposite.  Every pack's terms are program brand names, so a
    task saying "a DFT setup and a cheaper semi-empirical one" -- correct
    chemistry, naming no product -- matched nothing, the index came back empty,
    and the prompt then listed no skills at all.  The tool's own description
    tells the model to consult a skill "listed in the system prompt", so an
    empty list is an instruction not to consult anything.

    Three recorded observations of a frozen task were read as evidence that the
    model would not consult; they were measuring this gate.  General chemistry
    knowledge is not gated behind mentioning a vendor.
    """

    if not skills_enabled():
        return (), ()
    targets = sorted(
        {
            (pack.target_program, pack.target_engine)
            for pack in BUILTIN_PROGRAM_PACKS
        }
    )
    pack_sha256s: set[str] = set()
    skill_ids: set[str] = set(cross_program_skill_ids())
    for program, engine in targets:
        receipt = activate_program_knowledge(
            request=task, program=program, engine=engine
        )
        # Only genuinely matched packs are recorded as activated; the neutral
        # skills above are available without claiming any pack fired.
        pack_sha256s.update(receipt.activated_pack_sha256s)
        skill_ids.update(skills_for_activation(receipt))
    return tuple(sorted(pack_sha256s)), resolve_skills(
        tuple(sorted(skill_ids))
    )


def _coordinator_base_messages(
    *,
    context: Mapping[str, Any],
    approved_workflow: Mapping[str, Any] | None,
    bounded_review_requested: bool = False,
    task: str = "",
) -> list[dict[str, str]]:
    _, documents = activated_skill_documents(task)
    return [
        {
            "role": "system",
            "content": _system_prompt(
                approved_workflow,
                bounded_review_requested=bounded_review_requested,
                skill_index=tuple(item.index_entry() for item in documents),
            ),
        },
        {"role": "user", "content": canonical_json(context)},
    ]


def _validated_workspace(value: str | Path) -> Path:
    path = Path(value).expanduser()
    if not path.is_absolute():
        raise ContractError("agent workspace must be an absolute path")
    if not path.is_dir() or path.is_symlink():
        raise ContractError("agent workspace must be a regular directory")
    return path.resolve()


def _private_run_directory(workspace: Path, session_id: str) -> Path:
    private_root = workspace / _PRIVATE_ROOT_NAME
    if private_root.is_symlink():
        raise ContractError(
            "private agent directory cannot be a symbolic link"
        )
    private_root.mkdir(exist_ok=True, mode=0o700)
    if not private_root.is_dir():
        raise ContractError("private agent path is not a directory")
    private_root.chmod(0o700)
    root = private_root / "runs"
    if root.is_symlink():
        raise ContractError("private run root cannot be a symbolic link")
    root.mkdir(exist_ok=True, mode=0o700)
    if not root.is_dir():
        raise ContractError("private run root is not a directory")
    root.chmod(0o700)
    target = root / session_id
    target.mkdir(mode=0o700)
    target.chmod(0o700)
    return target


def _scan_xyz_artifacts(workspace: Path) -> tuple[_XyzObservation, ...]:
    observations: dict[str, _XyzObservation] = {}
    private_root = workspace / _PRIVATE_ROOT_NAME
    host_artifact_root = workspace / "artifacts"
    host_node_root = workspace / "nodes"
    for candidate in sorted(workspace.rglob("*.xyz")):
        if any(
            root in candidate.parents
            for root in (private_root, host_artifact_root, host_node_root)
        ):
            continue
        if candidate.is_symlink() or not candidate.is_file():
            continue
        resolved = candidate.resolve()
        try:
            resolved.relative_to(workspace)
        except ValueError as exc:
            raise ContractError(
                "XYZ artifact escapes the approved workspace"
            ) from exc
        digest = file_sha256(resolved)
        atom_count, symbols = _inspect_xyz(resolved)
        artifact_id = f"geometry-{digest[:16]}"
        artifact = TrustedArtifactRefV1(
            artifact_id=artifact_id,
            kind="geometry_xyz",
            sha256=digest,
            size_bytes=resolved.stat().st_size,
            path=str(resolved),
            cli_value=str(resolved),
        )
        observations.setdefault(
            digest,
            _XyzObservation(
                artifact=artifact,
                atom_count=atom_count,
                symbols=symbols,
            ),
        )
    return tuple(
        sorted(
            observations.values(), key=lambda item: item.artifact.artifact_id
        )
    )


def _scan_pyscf_result_artifacts(
    workspace: Path,
) -> tuple[_PySCFResultObservation, ...]:
    """Bind user-placed structured results without exposing model paths.

    Only normally terminated PySCF HDF5 artifacts with finite, typed metadata
    enter the tool host. Invalid or legacy result files are ignored rather than
    being misrepresented as analysis-ready evidence.
    """

    observations: dict[str, _PySCFResultObservation] = {}
    private_root = workspace / _PRIVATE_ROOT_NAME
    host_artifact_root = workspace / "artifacts"
    for candidate in sorted(workspace.rglob("*.h5")):
        if any(
            root in candidate.parents
            for root in (private_root, host_artifact_root)
        ):
            continue
        if candidate.is_symlink() or not candidate.is_file():
            continue
        resolved = candidate.resolve()
        try:
            resolved.relative_to(workspace)
        except ValueError as exc:
            raise ContractError(
                "result artifact escapes the approved workspace"
            ) from exc
        digest = file_sha256(resolved)
        try:
            output, receipt = validate_pyscf_analysis_artifact(
                resolved, expected_sha256=digest
            )
            receipt_sha256 = str(receipt["receipt_sha256"])
            jobtype = str(output.jobtype or "").strip().lower()
            method = str(output.method or "").strip()
            applied_method = str(output.spec.get("xc") or method).strip()
            basis = str(output.basis or "").strip()
            engine = str(output.engine or "").strip().lower()
            charge = output.charge
            multiplicity = output.multiplicity
            project_yaml_sha256 = str(output.project_yaml_digest or "").strip()
            input_artifact_sha256 = str(
                output.spec.get("input_artifact_sha256") or ""
            ).strip()
            if (
                jobtype not in {"sp", "opt", "hess"}
                or not method
                or not applied_method
                or not basis
                or engine not in {"cpu", "gpu"}
                or not isinstance(charge, int)
                or not isinstance(multiplicity, int)
                or multiplicity < 1
            ):
                continue
        except (
            OSError,
            UnicodeDecodeError,
            json.JSONDecodeError,
            ValueError,
            KeyError,
            TypeError,
            QuantityExtractionError,
        ):
            continue
        artifact = TrustedArtifactRefV1(
            artifact_id=f"pyscf-result-{digest[:16]}",
            kind="pyscf_hdf5",
            sha256=digest,
            size_bytes=resolved.stat().st_size,
            path=str(resolved),
            cli_value=str(resolved),
        )
        observations.setdefault(
            digest,
            _PySCFResultObservation(
                artifact=artifact,
                jobtype=jobtype,
                method=method,
                applied_method=applied_method,
                basis=basis,
                engine=engine,
                charge=charge,
                multiplicity=multiplicity,
                project_yaml_sha256=project_yaml_sha256,
                input_artifact_sha256=input_artifact_sha256,
                validation_receipt_sha256=receipt_sha256,
                scientific_validation_state=str(
                    receipt["scientific_validation_state"]
                ),
            ),
        )
    return tuple(
        sorted(
            observations.values(), key=lambda item: item.artifact.artifact_id
        )
    )


def _scan_xtb_result_artifacts(
    workspace: Path,
) -> tuple[_LoggedResultObservation, ...]:
    """Admit validated xTB native outputs to the shared quantity tool chain.

    The xTB runner already emits a complete result receipt next to its native
    output.  This adapter reuses that program contract and exposes the output
    under the registered ``xtb_output`` reader kind; it does not add a second
    parser or a paper-specific result path.
    """

    from chemsmart.jobs.xtb.validation import audit_xtb_result_receipt

    observations: dict[str, _LoggedResultObservation] = {}
    private_root = workspace / _PRIVATE_ROOT_NAME
    host_artifact_root = workspace / "artifacts"
    for receipt_path in sorted(workspace.rglob("*.xtb-result-receipt.json")):
        if any(
            root in receipt_path.parents
            for root in (private_root, host_artifact_root)
        ):
            continue
        if receipt_path.is_symlink() or not receipt_path.is_file():
            continue
        try:
            receipt_path.resolve().relative_to(workspace)
            receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
            requested = receipt["requested_settings"]
            jobtype = str(requested["jobtype"]).strip().lower()
            charge = requested["charge"]
            multiplicity = requested["multiplicity"]
            observation, findings = audit_xtb_result_receipt(
                receipt_path,
                expected_jobtype=jobtype,
                expected_charge=charge,
                expected_multiplicity=multiplicity,
                expected_settings={
                    key: requested.get(key)
                    for key in (
                        "gfn_version",
                        "solvent_model",
                        "solvent_id",
                        "grad",
                    )
                },
                audit_mode="archive",
            )
            if findings:
                continue
            artifacts = receipt["artifacts"]
            receipt_suffix = ".xtb-result-receipt.json"
            main_output_name = (
                receipt_path.name[: -len(receipt_suffix)] + ".out"
            )
            if main_output_name not in artifacts:
                continue
            output_path = (receipt_path.parent / main_output_name).resolve()
            output_path.relative_to(workspace)
            output_record = artifacts[main_output_name]
            digest = file_sha256(output_path)
            if (
                output_record.get("sha256") != digest
                or output_record.get("size") != output_path.stat().st_size
                or observation.get("final_energy_hartree") is None
                or not isinstance(charge, int)
                or not isinstance(multiplicity, int)
                or multiplicity < 1
            ):
                continue
            project_record = receipt.get("project_artifact") or {}
            source_record = receipt.get("source_artifact") or {}
            validation_receipt_sha256 = str(receipt["receipt_sha256"])
            artifact = TrustedArtifactRefV1(
                artifact_id=f"xtb-result-{digest[:16]}",
                kind="xtb_output",
                sha256=digest,
                size_bytes=output_path.stat().st_size,
                path=str(output_path),
                cli_value=str(output_path),
            )
            observations.setdefault(
                digest,
                _LoggedResultObservation(
                    artifact=artifact,
                    program="xtb",
                    jobtype=jobtype,
                    method=str(requested["gfn_version"]).strip().lower(),
                    basis=None,
                    engine="cpu",
                    charge=charge,
                    multiplicity=multiplicity,
                    project_yaml_sha256=str(
                        project_record.get("sha256") or ""
                    ),
                    input_artifact_sha256=str(
                        source_record.get("sha256") or ""
                    ),
                    validation_receipt_sha256=validation_receipt_sha256,
                    scientific_validation_state=str(
                        receipt["validation_state"]
                    ),
                    provenance_status=str(observation["provenance_status"]),
                    provenance_limitations=tuple(
                        str(item)
                        for item in observation["provenance_limitations"]
                    ),
                ),
            )
        except (
            ContractError,
            KeyError,
            OSError,
            RuntimeError,
            TypeError,
            ValueError,
            json.JSONDecodeError,
        ):
            continue
    return tuple(
        sorted(
            observations.values(), key=lambda item: item.artifact.artifact_id
        )
    )


def _scan_orca_result_artifacts(
    workspace: Path,
) -> tuple[_LoggedResultObservation, ...]:
    """Expose normally terminated ORCA outputs to the typed result reader."""

    from chemsmart.io.orca.output import ORCAOutput

    observations: dict[str, _LoggedResultObservation] = {}
    private_root = workspace / _PRIVATE_ROOT_NAME
    host_artifact_root = workspace / "artifacts"
    for candidate in sorted(workspace.rglob("*.out")):
        if any(
            root in candidate.parents
            for root in (private_root, host_artifact_root)
        ):
            continue
        if candidate.is_symlink() or not candidate.is_file():
            continue
        resolved = candidate.resolve()
        try:
            resolved.relative_to(workspace)
            output = ORCAOutput(filename=resolved)
            energy = output.final_energy
            molecule = output.molecule
            positions = molecule.positions
            method = str(output.method or "").strip().lower()
            basis = str(output.basis or "").strip().lower()
            dispersion = str(output.dispersion or "").strip().lower()
            semiempirical = str(output.semiempirical or "").strip().lower()
            jobtype = str(output.jobtype or "").strip().lower()
            charge = output.charge
            multiplicity = output.multiplicity
            if (
                not output.normal_termination
                or energy is None
                or not math.isfinite(float(energy))
                or not method
                or (not basis and not semiempirical)
                # A hand-listed set like this is how a new family gets
                # silently dropped: when `scan` joined the surface, completed
                # scans stopped being registered here, so a session opening a
                # workspace that held two validated 72-point surfaces was
                # shown neither and rationally replanned both from scratch.
                # (Before the reader learned to classify a %geom Scan block,
                # scan outputs slipped through this filter mislabelled as
                # `opt` -- the exclusion only started biting when the label
                # became correct.)
                or jobtype
                not in {
                    "sp",
                    "opt",
                    "freq",
                    "td",
                    "ts",
                    "irc",
                    "scan",
                    "modred",
                }
                or not isinstance(charge, int)
                or not isinstance(multiplicity, int)
                or multiplicity < 1
                or not molecule.chemical_symbols
                or not all(
                    math.isfinite(float(value))
                    for row in positions
                    for value in row
                )
            ):
                continue
        except (
            AttributeError,
            IndexError,
            OSError,
            TypeError,
            ValueError,
        ):
            continue
        digest = file_sha256(resolved)
        artifact = TrustedArtifactRefV1(
            artifact_id=f"orca-result-{digest[:16]}",
            kind="orca_output",
            sha256=digest,
            size_bytes=resolved.stat().st_size,
            path=str(resolved),
            cli_value=str(resolved),
        )
        observations.setdefault(
            digest,
            _LoggedResultObservation(
                artifact=artifact,
                program="orca",
                jobtype=jobtype,
                method=method,
                basis=basis or None,
                dispersion=dispersion or None,
                engine="cpu",
                charge=charge,
                multiplicity=multiplicity,
                project_yaml_sha256="",
                input_artifact_sha256="",
                validation_receipt_sha256="",
                scientific_validation_state="parsed_normal_termination",
            ),
        )
    return tuple(
        sorted(
            observations.values(), key=lambda item: item.artifact.artifact_id
        )
    )


def _scan_gaussian_result_artifacts(
    workspace: Path,
) -> tuple[_LoggedResultObservation, ...]:
    """Expose normally terminated Gaussian results to the typed reader.

    Gaussian remains an operator-installed, separately licensed program.  A
    result that a user deliberately places in the workspace is nevertheless a
    first-class scientific artifact: the existing Gaussian parser can recover
    geometries, energies, frequencies, thermochemistry, excited states, spin
    diagnostics, trajectories, and stability evidence without launching the
    executable.  Admission here is based on native normal termination and
    finite molecular/energy data, just as for the ORCA native-output path.
    """

    from chemsmart.io.gaussian.output import Gaussian16Output

    observations: dict[str, _LoggedResultObservation] = {}
    private_root = workspace / _PRIVATE_ROOT_NAME
    host_artifact_root = workspace / "artifacts"
    candidates = {
        *workspace.rglob("*.log"),
        *workspace.rglob("*.out"),
    }
    for candidate in sorted(candidates):
        if any(
            root in candidate.parents
            for root in (private_root, host_artifact_root)
        ):
            continue
        if candidate.is_symlink() or not candidate.is_file():
            continue
        resolved = candidate.resolve()
        try:
            resolved.relative_to(workspace)
            output = Gaussian16Output(filename=str(resolved))
            energies = tuple(output.energies or ())
            molecule = output.molecule
            positions = molecule.positions
            method = str(output.method or "").strip().lower()
            basis = str(output.basis or "").strip().lower() or None
            jobtype = str(output.jobtype or "").strip().lower()
            charge = output.charge
            multiplicity = output.multiplicity
            if (
                not output.normal_termination
                or not energies
                or not math.isfinite(float(energies[-1]))
                or not method
                or jobtype
                not in {"sp", "opt", "freq", "td", "ts", "irc", "link"}
                or not isinstance(charge, int)
                or not isinstance(multiplicity, int)
                or multiplicity < 1
                or not molecule.chemical_symbols
                or not all(
                    math.isfinite(float(value))
                    for row in positions
                    for value in row
                )
            ):
                continue
        except (
            AttributeError,
            IndexError,
            OSError,
            TypeError,
            ValueError,
        ):
            continue
        digest = file_sha256(resolved)
        artifact = TrustedArtifactRefV1(
            artifact_id=f"gaussian-result-{digest[:16]}",
            kind="gaussian_output",
            sha256=digest,
            size_bytes=resolved.stat().st_size,
            path=str(resolved),
            cli_value=str(resolved),
        )
        observations.setdefault(
            digest,
            _LoggedResultObservation(
                artifact=artifact,
                program="gaussian",
                jobtype=jobtype,
                method=method,
                basis=basis,
                engine="cpu",
                charge=charge,
                multiplicity=multiplicity,
                project_yaml_sha256="",
                input_artifact_sha256="",
                validation_receipt_sha256="",
                scientific_validation_state="parsed_normal_termination",
            ),
        )
    return tuple(
        sorted(
            observations.values(), key=lambda item: item.artifact.artifact_id
        )
    )


def _scan_result_artifacts(workspace: Path) -> tuple[_ResultObservation, ...]:
    """Collect every registered, analysis-ready program result."""

    return tuple(
        sorted(
            (
                *_scan_pyscf_result_artifacts(workspace),
                *_scan_xtb_result_artifacts(workspace),
                *_scan_orca_result_artifacts(workspace),
                *_scan_gaussian_result_artifacts(workspace),
            ),
            key=lambda item: item.artifact.artifact_id,
        )
    )


def _inspect_xyz(path: Path) -> tuple[int, tuple[str, ...]]:
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
        count = int(lines[0].strip())
    except (IndexError, UnicodeDecodeError, ValueError) as exc:
        raise ContractError(
            "workspace XYZ has a malformed atom count"
        ) from exc
    if count < 1 or len(lines) < count + 2:
        raise ContractError("workspace XYZ is truncated")
    symbols = []
    for row in lines[2 : count + 2]:
        fields = row.split()
        if len(fields) < 4 or not fields[0].isalpha():
            raise ContractError("workspace XYZ has a malformed atom row")
        try:
            coordinates = tuple(float(item) for item in fields[1:4])
        except ValueError as exc:
            raise ContractError(
                "workspace XYZ coordinates are not numeric"
            ) from exc
        if not all(math.isfinite(item) for item in coordinates):
            raise ContractError("workspace XYZ coordinates must be finite")
        symbols.append(fields[0])
    return count, tuple(symbols)


def _task_spec_sha256(
    task: str,
    observations: Iterable[_XyzObservation],
    approved_molecular_identity: (
        ApprovedMolecularIdentityV1
        | Iterable[ApprovedMolecularIdentityV1]
        | None
    ) = None,
    *,
    result_observations: Iterable[_ResultObservation] = (),
    approved_molecular_inputs: Iterable[ApprovedMolecularInputV1] = (),
) -> str:
    identities = _coerce_approved_identities(
        (
            approved_molecular_identity
            if isinstance(
                approved_molecular_identity, ApprovedMolecularIdentityV1
            )
            else None
        ),
        (
            approved_molecular_identity
            if approved_molecular_identity is not None
            and not isinstance(
                approved_molecular_identity, ApprovedMolecularIdentityV1
            )
            else ()
        ),
    )
    body = {
        "schema_version": "chemsmart.live-scientific-task.v1",
        "task": task,
        "coordinate_artifacts": tuple(
            {
                "artifact_id": item.artifact.artifact_id,
                "sha256": item.artifact.sha256,
                "atom_count": item.atom_count,
                "symbols": item.symbols,
            }
            for item in observations
        ),
        "result_artifacts": tuple(
            {
                "artifact_id": item.artifact.artifact_id,
                "sha256": item.artifact.sha256,
                "program": item.program,
                "jobtype": item.jobtype,
            }
            for item in result_observations
        ),
    }
    if len(identities) <= 1:
        # Preserve the existing single-identity task digest contract.
        body["approved_molecular_identity_sha256"] = (
            identities[0].identity_sha256 if identities else ""
        )
    else:
        body["approved_molecular_identity_sha256s"] = tuple(
            identity.identity_sha256 for identity in identities
        )
    molecular_inputs = _coerce_approved_molecular_inputs(
        approved_molecular_inputs
    )
    if molecular_inputs:
        body["approved_molecular_input_sha256s"] = tuple(
            item.assignment_sha256 for item in molecular_inputs
        )
    return canonical_sha256(body)


def _coerce_approved_molecular_inputs(
    values: Iterable[ApprovedMolecularInputV1],
) -> tuple[ApprovedMolecularInputV1, ...]:
    approved = tuple(values)
    if any(
        not isinstance(item, ApprovedMolecularInputV1) for item in approved
    ):
        raise ContractError("approved molecular inputs must be typed records")
    input_ids = tuple(item.input_id for item in approved)
    assignment_sha256s = tuple(item.assignment_sha256 for item in approved)
    if len(input_ids) != len(set(input_ids)) or len(assignment_sha256s) != len(
        set(assignment_sha256s)
    ):
        raise ContractError("approved molecular inputs must be unique")
    return tuple(sorted(approved, key=lambda item: item.input_id))


def _approved_input_state_bindings(
    molecular_inputs: tuple[ApprovedMolecularInputV1, ...],
    *,
    observations: tuple[_XyzObservation, ...],
    task_spec_sha256: str,
) -> dict[str, Any]:
    """Prebind only states explicitly approved for exact input bytes."""

    observations_by_sha256 = {
        item.artifact.sha256: item for item in observations
    }
    bindings = {}
    for molecular_input in molecular_inputs:
        observation = observations_by_sha256.get(
            molecular_input.molecular_identity.geometry_sha256
        )
        if observation is None:
            raise ContractError(
                f"approved molecular input {molecular_input.input_id!r} has no "
                "matching coordinate artifact"
            )
        binding = build_scientific_identity_binding(
            task_spec_sha256=task_spec_sha256,
            geometry_artifact=observation.artifact,
            charge=molecular_input.charge,
            multiplicity=molecular_input.multiplicity,
        )
        bindings[binding.binding_sha256] = binding
    return bindings


def _coerce_approved_identities(
    single: ApprovedMolecularIdentityV1 | None,
    multiple: Iterable[ApprovedMolecularIdentityV1],
) -> tuple[ApprovedMolecularIdentityV1, ...]:
    values = tuple(multiple)
    if single is not None:
        if values:
            raise ContractError(
                "use approved_molecular_identity or approved_molecular_identities, "
                "not both"
            )
        values = (single,)
    if any(
        not isinstance(item, ApprovedMolecularIdentityV1) for item in values
    ):
        raise ContractError(
            "approved molecular identities must be typed records"
        )
    identity_ids = tuple(item.identity_id for item in values)
    identity_sha256s = tuple(item.identity_sha256 for item in values)
    if len(identity_ids) != len(set(identity_ids)) or len(
        identity_sha256s
    ) != len(set(identity_sha256s)):
        raise ContractError("approved molecular identities must be unique")
    return tuple(sorted(values, key=lambda item: item.identity_id))


def _validated_identity_records(
    observations: tuple[_XyzObservation, ...],
    identity: (
        ApprovedMolecularIdentityV1
        | Iterable[ApprovedMolecularIdentityV1]
        | None
    ),
) -> tuple[dict[str, Any], ...]:
    identities = _coerce_approved_identities(
        (
            identity
            if isinstance(identity, ApprovedMolecularIdentityV1)
            else None
        ),
        (
            identity
            if identity is not None
            and not isinstance(identity, ApprovedMolecularIdentityV1)
            else ()
        ),
    )
    if not identities:
        return ()
    observations_by_sha256: dict[str, _XyzObservation] = {}
    for observation in observations:
        previous = observations_by_sha256.setdefault(
            observation.artifact.sha256, observation
        )
        if previous.symbols != observation.symbols:
            raise ContractError(
                "identical coordinate bytes produced inconsistent atom-order records"
            )
    records = []
    for approved in identities:
        observation = observations_by_sha256.get(approved.geometry_sha256)
        if observation is None:
            raise ContractError(
                f"approved molecular identity {approved.identity_id!r} has no "
                "matching coordinate artifact"
            )
        validate_identity_for_geometry(
            approved,
            geometry_sha256=observation.artifact.sha256,
            atom_order=observation.symbols,
        )
        records.append(approved.public_record())
    return tuple(records)


def _session_id(task_spec_sha256: str) -> str:
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    return f"live-{timestamp}-{task_spec_sha256[:10]}-{uuid.uuid4().hex[:8]}"


#: The core scientific stages a hierarchical workflow is built from.
#: Conformance covers the intersection of this set with what each program
#: declares, so bringing a new program to ``preview_only`` needs a ChemSmart
#: declaration and a project fixture -- not a new branch in this function.
_CONFORMANCE_CORE_STAGES = frozenset({"hess", "opt", "sp", "td"})

#: Section shape of each program's project YAML.  ChemSmart does not declare
#: this anywhere, because it is a property of each program's settings loader:
#: PySCF keys sections by jobtype, while the route-building programs split
#: ``gas`` from ``solv``.  It is data here rather than an if-ladder so a new
#: program is one entry.
_CONFORMANCE_PROJECT_SHAPES = {
    "gaussian": "route_gas_solv",
    "orca": "route_gas_solv",
    "pyscf": "pyscf_stage_keyed",
}


def _conformance_jobtypes(program: str, engine: str) -> tuple[str, ...]:
    """Return the declared stages ``program`` previews on ``engine``.

    An explicit engine/job matrix is the program owner's bounded preview
    surface.  The legacy core-stage intersection remains only for older
    declarations whose independent engine and job lists would otherwise form
    an unverified Cartesian product.
    """

    from chemsmart.settings.capabilities import ENGINE_CAPABILITIES

    capability = ENGINE_CAPABILITIES.get(program)
    if capability is None:
        return ()
    pairs = capability.engine_job_capabilities
    if pairs:
        # A program that declares per-engine capability is authoritative about
        # it; GPU PySCF, for example, previews no excited-state stage.
        declared = {
            item.jobtype
            for item in pairs
            if item.engine == engine and item.preview_supported
        }
        return tuple(sorted(declared))
    declared = set(capability.jobtypes)
    return tuple(sorted(declared & _CONFORMANCE_CORE_STAGES))


#: Stages of a route-building program that read their own project section
#: instead of ``gas``/``solv``.  Their project sections include the same typed
#: scientific settings that a user supplies in a normal project YAML; the
#: conformance preview must exercise that canonical path rather than rely on
#: hidden command-specific defaults.
_ROUTE_PROGRAM_STAGE_SECTIONS = frozenset({"link", "neb", "td"})


def _conformance_project_sections(
    program: str,
) -> dict[str, dict[str, Any]] | None:
    """Return deterministic conformance sections, or ``None`` when unneeded."""

    shape = _CONFORMANCE_PROJECT_SHAPES.get(program)
    if shape == "pyscf_stage_keyed":
        return _pyscf_conformance_sections()
    if shape == "route_gas_solv":
        # A deliberately plain, cheap method: this fixture exists to exercise
        # the compile/preview path, and must never look like a scientific
        # recommendation for any particular system.
        common = {"basis": "def2-SVP", "functional": "B3LYP"}
        sections = {"gas": dict(common), "solv": {**common, "freq": False}}
        # ``gas``/``solv`` feed ordinary stages, but a few jobtypes read
        # their own section and the loader returns ``None`` when it is absent.
        # A fixture must therefore declare every stage it claims to cover, or
        # conformance fails inside the CLI rather than reporting a gap.
        for stage in _conformance_jobtypes(program, "cpu"):
            if stage in _ROUTE_PROGRAM_STAGE_SECTIONS:
                if stage == "neb":
                    sections[stage] = {
                        **common,
                        "joboption": "NEB-TS",
                        "nimages": 2,
                        "preopt_ends": False,
                    }
                elif stage == "link":
                    sections[stage] = {
                        **common,
                        "freq": False,
                        "guess": "mix",
                        "jobtype": "opt",
                        "stable": "opt",
                    }
                elif stage == "td":
                    sections[stage] = (
                        {
                            **common,
                            "nstates": 3,
                            "states": "singlets",
                        }
                        if program == "gaussian"
                        else {
                            **common,
                            "nstates": 3,
                            "response_method": "tda",
                            "state_manifold": "singlet",
                        }
                    )
        return sections
    return None


def _bootstrap_conformance(
    *,
    run_directory: Path,
    input_artifact: TrustedArtifactRefV1,
    registry_sha256: str,
    live_schema: Any,
    resources: ExecutionResourceSpecV1 | None = None,
) -> tuple[
    tuple[ProgramComponentConformanceReceiptV1, ...],
    tuple[dict[str, Any], ...],
]:
    """Bootstrap every declared program through its real fake-preview path.

    Conformance never invokes a program binary: the generated input file is the
    conformance artifact.  A receipt therefore enables planning and safe
    preview only, and a program with no installed executable can still reach
    ``preview_only``.
    """

    bootstrap_directory = run_directory / "bootstrap"
    bootstrap_directory.mkdir(mode=0o700)
    server_path = bootstrap_directory / "preview-server.yaml"
    _write_private_exact(
        server_path,
        _preview_server_profile(resources=resources).encode("utf-8"),
    )

    receipts: list[ProgramComponentConformanceReceiptV1] = []
    records: list[dict[str, Any]] = []
    for program in _conformance_programs():
        project_path: Path | None = None
        sections = _conformance_project_sections(program)
        if sections is not None:
            project_path = bootstrap_directory / f"{program}-fixture.yaml"
            rendered = render_project_yaml(
                project_document(program=program, sections=sections)
            )
            _write_private_exact(
                project_path, rendered.rendered_yaml.encode("utf-8")
            )
        parts: list[ProgramComponentConformanceReceiptV1] = []
        for engine in _conformance_engines(program):
            jobtypes = _conformance_jobtypes(program, engine)
            if not jobtypes:
                continue
            try:
                receipt = bootstrap_program_conformance(
                    program=program,
                    engine=engine,
                    jobtypes=jobtypes,
                    input_path=input_artifact.path,
                    project_path=project_path,
                    server_path=server_path,
                    charge=0,
                    multiplicity=1,
                    registry_sha256=registry_sha256,
                    live_schema=live_schema,
                )
            except Exception as exc:  # Stays observable, fail-closed.
                records.append(
                    {
                        "record_kind": "program_conformance",
                        "program": program,
                        "engine": engine,
                        "status": "failed",
                        "error_class": type(exc).__name__,
                        "covered_jobtypes": (),
                    }
                )
                continue
            parts.append(receipt)
            records.append(_conformance_record(receipt, engine=engine))
        if parts:
            receipts.append(
                parts[0]
                if len(parts) == 1
                else _combine_program_conformance(parts)
            )
    return tuple(receipts), tuple(records)


def _ensure_preview_server(
    run_directory: Path,
    *,
    resources: ExecutionResourceSpecV1 | None = None,
) -> Path:
    """Materialize the host-selected server used by every safe preview.

    Bootstrap conformance already passes this file explicitly.  Normal model
    commands must do the same; otherwise Click falls back to a machine-local
    ``local.yaml`` and the observed program environment silently differs from
    the one shown to the model.
    """

    bootstrap_directory = run_directory / "bootstrap"
    bootstrap_directory.mkdir(mode=0o700, exist_ok=True)
    server_path = bootstrap_directory / "preview-server.yaml"
    _write_private_exact(
        server_path,
        _preview_server_profile(resources=resources).encode("utf-8"),
    )
    return server_path


def _conformance_engines(program: str) -> tuple[str, ...]:
    """Return the v1 agent preview engines ChemSmart declares."""

    from chemsmart.settings.capabilities import (
        AGENT_PROGRAM_PREVIEW_ENGINES,
    )

    return tuple(AGENT_PROGRAM_PREVIEW_ENGINES.get(program, ()))


def _conformance_programs() -> tuple[str, ...]:
    """Return the programs ChemSmart declares as executable.

    Deriving this from the declaration rather than a hand-written list is what
    makes a newly added program visible to the agent without touching the
    session code.
    """

    from chemsmart.settings.capabilities import AGENT_PROGRAMS

    return tuple(
        program
        for program in sorted(AGENT_PROGRAMS)
        if _conformance_jobtypes(program, "cpu")
    )


def _combine_program_conformance(
    receipts: Iterable[ProgramComponentConformanceReceiptV1],
) -> ProgramComponentConformanceReceiptV1:
    rows = tuple(receipts)
    if not rows or len({item.program for item in rows}) != 1:
        raise ContractError("combined conformance requires one program")
    passed = tuple(
        item
        for item in rows
        if all(
            value == "passed"
            for value in (
                item.compiler_status,
                item.preview_status,
                item.preflight_status,
            )
        )
    )
    covered_engine_job_pairs = tuple(
        sorted(
            {
                pair
                for item in passed
                for pair in item.effective_engine_job_pairs
            }
        )
    )
    covered_engines = tuple(
        sorted({engine for engine, _jobtype in covered_engine_job_pairs})
    )
    covered_jobtypes = tuple(
        sorted({jobtype for _engine, jobtype in covered_engine_job_pairs})
    )
    status = "passed" if covered_engine_job_pairs else "failed"

    def aggregate(field: str) -> str:
        return canonical_sha256(tuple(getattr(item, field) for item in rows))

    verifier_status = (
        "failed"
        if any(item.verifier_status == "failed" for item in rows)
        else (
            "passed"
            if all(item.verifier_status == "passed" for item in rows)
            else "not_observed"
        )
    )

    return build_program_component_conformance_receipt(
        program=rows[0].program,
        registry_sha256=rows[0].registry_sha256,
        live_cli_schema_sha256=rows[0].live_cli_schema_sha256,
        fixture_bundle_sha256=canonical_sha256(
            tuple(item.fixture_bundle_sha256 for item in rows)
        ),
        covered_jobtypes=covered_jobtypes,
        covered_engines=covered_engines,
        covered_engine_job_pairs=covered_engine_job_pairs,
        compiler_receipt_sha256=aggregate("compiler_receipt_sha256"),
        preview_receipt_sha256=aggregate("preview_receipt_sha256"),
        preflight_receipt_sha256=aggregate("preflight_receipt_sha256"),
        verifier_receipt_sha256=(
            aggregate("verifier_receipt_sha256")
            if verifier_status != "not_observed"
            else ""
        ),
        compiler_status=status,
        preview_status=status,
        preflight_status=status,
        verifier_status=verifier_status,
    )


def _pyscf_conformance_sections() -> dict[str, dict[str, Any]]:
    common = {
        "basis": "def2-svp",
        "defgrid": "defgrid2",
        "density_fit": False,
        "functional": "b3lyp",
        "scf_maxiter": 100,
        "scf_tol": 1e-9,
    }
    return {
        "sp": dict(common),
        "opt": {**common, "opt_maxsteps": 100, "opt_solver": "geometric"},
        "hess": dict(common),
        "td": {
            **common,
            "nstates": 5,
            "response_method": "tda",
            "state_manifold": "singlet",
        },
    }


#: Keys in a server YAML that describe the scheduler rather than a program.
_NON_PROGRAM_SERVER_KEYS = frozenset({"SERVER"})


#: A discovery stub announces itself with this marker so the harness can tell a
#: real program from a placeholder that only satisfies program discovery.
_DISCOVERY_STUB_MARKER = b"ChemSmart agent-harness DISCOVERY STUB"


def _executable_is_discovery_stub(path: str) -> bool:
    """Return whether ``path`` is a discovery stub, not a real program."""

    if not path:
        return False
    try:
        with open(path, "rb") as handle:
            head = handle.read(2048)
    except OSError:
        return False
    return _DISCOVERY_STUB_MARKER in head


def _declared_executable_path(program: str, folder: str) -> str:
    """Return the executable ChemSmart itself would run from ``folder``.

    The binary is not always named after the program -- Gaussian's is ``g16``.
    Asking ChemSmart's own executable registry keeps the agent's view of the
    environment identical to the one a real job would use, instead of a second,
    guessed convention that silently disagrees with it.
    """

    from chemsmart.settings.executable import Executable

    for subclass in Executable.subclasses():
        if str(subclass.PROGRAM or "").lower() != program:
            continue
        try:
            resolved = subclass(executable_folder=folder).get_executable()
        except Exception:
            return ""
        return str(resolved or "")
    return ""


def _declared_server_programs() -> tuple[tuple[str, str], ...]:
    """Return ``(program, executable_folder)`` for every declared program.

    ChemSmart's server YAML is the canonical statement of which programs this
    installation controls.  Deriving the agent's view from it -- instead of
    hard-coding a list here -- is what keeps the agent's environment the same
    environment ChemSmart itself uses.  Adding a program to the server YAML is
    then sufficient to make the agent see it.

    A declared program whose folder is absent is still returned: the agent
    should be able to tell "ChemSmart knows about this program and it is not
    installed" apart from "this program does not exist".
    """

    from chemsmart.io.yaml import YAMLFile
    from chemsmart.settings.user import CHEMSMARTUserSettings

    settings = CHEMSMARTUserSettings()
    available = list(settings.all_available_servers or ())
    preferred = os.environ.get("CHEMSMART_AGENT_SERVER") or "local"
    name = (
        preferred
        if preferred in available
        else (available[0] if available else "")
    )
    if not name:
        return ()
    path = Path(settings.user_server_dir) / f"{name}.yaml"
    if not path.is_file():
        return ()
    try:
        content = YAMLFile(filename=str(path)).yaml_contents_dict
    except Exception:  # A malformed server file must not break planning.
        return ()
    rows = []
    for key, value in (content or {}).items():
        if key in _NON_PROGRAM_SERVER_KEYS or not isinstance(value, dict):
            continue
        folder = value.get("EXEFOLDER")
        if not folder:
            continue
        rows.append((str(key).lower(), str(Path(str(folder)).expanduser())))
    return tuple(sorted(set(rows)))


def _active_server_program_blocks() -> dict[str, dict[str, Any]]:
    """Read complete program blocks from the active user server YAML.

    Bounded real execution must retain ORCA MPI and Gaussian environment
    setup, not reduce an installed program to its executable directory.
    Scheduler keys remain host-owned and are deliberately excluded.
    """

    from chemsmart.io.yaml import YAMLFile
    from chemsmart.settings.user import CHEMSMARTUserSettings

    settings = CHEMSMARTUserSettings()
    available = list(settings.all_available_servers or ())
    preferred = os.environ.get("CHEMSMART_AGENT_SERVER") or "local"
    name = (
        preferred
        if preferred in available
        else (available[0] if available else "")
    )
    if not name:
        return {}
    path = Path(settings.user_server_dir) / f"{name}.yaml"
    if not path.is_file() or path.is_symlink():
        return {}
    try:
        content = YAMLFile(filename=str(path)).yaml_contents_dict
    except Exception as exc:
        raise ContractError("active ChemSmart server YAML is invalid") from exc
    return {
        str(key).upper(): dict(value)
        for key, value in (content or {}).items()
        if str(key).upper() not in _NON_PROGRAM_SERVER_KEYS
        and isinstance(value, dict)
    }


def _preview_server_profile(
    *, resources: ExecutionResourceSpecV1 | None = None
) -> str:
    """Scheduler-shaped profile for non-submitting run/sub conformance."""

    cores = resources.cores if resources is not None else 4
    memory_gb = resources.memory_gb if resources is not None else 4
    gpu_count = resources.gpu_count if resources is not None else 0

    return (
        "SERVER:\n"
        "  SCHEDULER: PBS\n"
        "  QUEUE_NAME: preview\n"
        "  NUM_HOURS: 1\n"
        f"  MEM_GB: {memory_gb:g}\n"
        f"  NUM_CORES: {cores}\n"
        f"  NUM_GPUS: {gpu_count}\n"
        f"  NUM_THREADS: {cores}\n"
        "  SUBMIT_COMMAND: true\n"
        "  SCRATCH_DIR: null\n"
        "  PROJECT: preview\n"
        "  USE_HOSTS: false\n" + _local_program_server_blocks()
    )


def _local_program_server_blocks(
    *, preserve_active: bool = False, scratch_root: Path | None = None
) -> str:
    """Return the installed-program blocks shared by preview and execution."""

    if preserve_active:
        active = _active_server_program_blocks()
        known_provider_labels = {
            normalize_key_label(label)
            for labels in DEFAULT_KEY_LABELS.values()
            for label in labels
        }

        def is_provider_secret_assignment(line: str) -> bool:
            match = re.match(
                r"^\s*(?:export\s+)?([A-Za-z_][A-Za-z0-9_-]*)\s*=",
                line,
            )
            if match is None:
                return False
            label = normalize_key_label(match.group(1))
            return label in known_provider_labels or any(
                token in label for token in PROVIDER_KEY_LABEL_TOKENS.values()
            )

        for block in active.values():
            for field in ("ENVARS", "SCRIPTS"):
                raw = block.get(field)
                if not isinstance(raw, str):
                    continue
                retained = tuple(
                    line
                    for line in raw.splitlines()
                    if not is_provider_secret_assignment(line)
                )
                if retained:
                    block[field] = "\n".join(retained) + "\n"
                else:
                    block.pop(field, None)
        pyscf = active.get("PYSCF")
        if pyscf is not None:
            pyscf.setdefault("LOCAL_RUN", True)
            pyscf.setdefault("SCRATCH", False)
        pyscf_override = os.environ.get("CHEMSMART_PYSCF_INTERPRETER")
        if pyscf_override:
            block = active.setdefault("PYSCF", {})
            block["EXEFOLDER"] = str(
                Path(pyscf_override).expanduser().resolve().parent
            )
            block.setdefault("LOCAL_RUN", True)
            block.setdefault("SCRATCH", False)
        xtb_override = os.environ.get("CHEMSMART_XTB_EXECUTABLE")
        if xtb_override:
            block = active.setdefault("XTB", {})
            block["EXEFOLDER"] = str(
                Path(xtb_override).expanduser().resolve().parent
            )
            block.setdefault("LOCAL_RUN", True)
        elif "XTB" in active:
            active["XTB"].setdefault("LOCAL_RUN", True)
        if scratch_root is not None:
            scratch_assignment = re.compile(
                r"^\s*export\s+(?:SCRATCH|GAUSS_SCRDIR)\s*="
            )
            for program, block in active.items():
                raw_envars = str(block.get("ENVARS") or "")
                lines = tuple(
                    line
                    for line in raw_envars.splitlines()
                    if not scratch_assignment.match(line)
                )
                if program == "GAUSSIAN":
                    lines = (
                        *lines,
                        f"export GAUSS_SCRDIR={scratch_root}",
                    )
                if lines:
                    block["ENVARS"] = "\n".join(lines) + "\n"
                else:
                    block.pop("ENVARS", None)
        # Safe dumping preserves multiline ENVARS/MODULES/SCRIPTS/CONDA_ENV
        # as data and never evaluates them in this controller.
        import yaml

        return yaml.safe_dump(active, sort_keys=True)

    declared = dict(_declared_server_programs())
    xtb = os.environ.get("CHEMSMART_XTB_EXECUTABLE") or shutil.which("xtb")
    if "xtb" not in declared and xtb:
        declared["xtb"] = str(Path(xtb).expanduser().parent)
    # PySCF is a library backend whose executable is a Python interpreter, so
    # it needs no server declaration to be real; the rest come from the server
    # YAML with preview-safe values layered on top.
    blocks = [
        "PYSCF:\n"
        f"  EXEFOLDER: {str(_PYSCF_INTERPRETER.parent)!r}\n"
        "  LOCAL_RUN: true\n"
        "  SCRATCH: false\n"
    ]
    for program, folder in sorted(declared.items()):
        if program == "pyscf":
            continue
        blocks.append(
            f"{program.upper()}:\n"
            f"  EXEFOLDER: {folder!r}\n"
            "  LOCAL_RUN: true\n"
            "  SCRATCH: false\n"
        )
    return "".join(blocks)


def _write_private_exact(path: Path, payload: bytes) -> None:
    if path.exists():
        if path.is_symlink() or path.read_bytes() != payload:
            raise ContractError(
                "private bootstrap artifact conflicts with existing bytes"
            )
        return
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    try:
        pending = memoryview(payload)
        while pending:
            written = os.write(descriptor, pending)
            pending = pending[written:]
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _observe_environments() -> tuple[
    tuple[EnvironmentTargetV1, ...],
    tuple[TrustedComputeEnvironmentReceiptV1, ...],
    tuple[dict[str, Any], ...],
]:
    targets = [
        EnvironmentTargetV1(
            program="pyscf",
            engine="cpu",
            target_kind="compute_receipt",
            locator="registered-pyscf-cpu-interpreter",
            required_dependencies=("h5py", "numpy", "pyscf"),
            required_dependency_versions=(("pyscf", "2.14.0"),),
        ),
        EnvironmentTargetV1(
            program="pyscf",
            engine="gpu",
            target_kind="compute_receipt",
            locator="registered-pyscf-gpu-interpreter",
            required_dependencies=(
                "cupy",
                "gpu4pyscf",
                "h5py",
                "numpy",
                "pyscf",
            ),
            required_dependency_versions=(
                ("gpu4pyscf", "1.8.0"),
                ("pyscf", "2.14.0"),
            ),
            required_gpu_facts=("device_available", "gpu4pyscf_distribution"),
        ),
    ]
    # Every program ChemSmart declares becomes a discovery target, so the agent
    # observes the same installation ChemSmart controls rather than a list
    # maintained here.
    for program, folder in _declared_server_programs():
        if program == "pyscf":
            continue
        executable = _declared_executable_path(program, folder)
        # A discovery stub is useful for loader/preview conformance, but it is
        # not an execution environment.  Do not hand it to the environment
        # resolver as an executable target: doing so makes ``which`` succeed
        # and invites the agent to choose a program that cannot produce a
        # scientific result.  With no target, the ordinary binding remains
        # preview-only while project rendering and safe preview still work.
        if executable and _executable_is_discovery_stub(executable):
            continue
        targets.append(
            EnvironmentTargetV1(
                program=program,
                engine="cpu",
                target_kind="executable",
                # Observe the executable ChemSmart itself resolves from the
                # server profile.  Re-running ``which(program)`` creates a
                # second environment model and is wrong for names such as
                # Gaussian, whose executable is ``g16``.
                locator=executable or program,
            )
        )
    if not any(item.program == "xtb" for item in targets):
        targets.append(
            EnvironmentTargetV1(
                program="xtb",
                engine="cpu",
                target_kind="executable",
                locator="xtb",
            )
        )
    receipts = []
    records = []
    for engine in ("cpu", "gpu"):
        try:
            receipt = probe_python_compute_environment(
                _PYSCF_INTERPRETER, engine=engine
            )
            receipts.append(receipt)
            gpu = dict(receipt.gpu_evidence)
            records.append(
                {
                    "record_kind": "program_environment",
                    "program": "pyscf",
                    "engine": engine,
                    "status": (
                        "available"
                        if engine == "cpu"
                        else (
                            "available"
                            if gpu.get("device_available") is True
                            and gpu.get("gpu4pyscf_distribution") is True
                            else "missing"
                        )
                    ),
                    "dependency_versions": receipt.dependency_versions,
                    "gpu_available": bool(gpu.get("device_available", False)),
                    "receipt_sha256": receipt.evidence_sha256,
                }
            )
        except Exception as exc:
            records.append(
                {
                    "record_kind": "program_environment",
                    "program": "pyscf",
                    "engine": engine,
                    "status": "missing",
                    "error_class": type(exc).__name__,
                }
            )
    for program, folder in _declared_server_programs():
        if program == "pyscf":
            continue
        candidate = _declared_executable_path(program, folder)
        located = (
            candidate
            if candidate and Path(candidate).exists()
            else (shutil.which(program) or "")
        )
        is_discovery_stub = _executable_is_discovery_stub(located)
        records.append(
            {
                "record_kind": "program_environment",
                "program": program,
                "engine": "cpu",
                "status": (
                    "preview_only"
                    if is_discovery_stub
                    else ("available" if located else "missing")
                ),
                "declared_folder": folder,
                "observation_method": "declared_server_exefolder",
                # A stub satisfies discovery and planning but cannot compute.
                # Recording that keeps "a binary was found" distinct from "a
                # scientific result is obtainable" in the evidence chain.
                "is_discovery_stub": is_discovery_stub,
            }
        )
    return tuple(targets), tuple(receipts), tuple(records)


def _conformance_record(
    receipt: ProgramComponentConformanceReceiptV1, *, engine: str
) -> dict[str, Any]:
    status = (
        "passed"
        if all(
            item == "passed"
            for item in (
                receipt.compiler_status,
                receipt.preview_status,
                receipt.preflight_status,
                receipt.verifier_status,
            )
        )
        else "failed"
    )
    return {
        "record_kind": "program_conformance",
        "program": receipt.program,
        "engine": engine,
        "status": status,
        "covered_jobtypes": receipt.covered_jobtypes,
        "covered_engine_job_pairs": receipt.effective_engine_job_pairs,
        "receipt_sha256": receipt.receipt_sha256,
    }


def _record_sort_key(value: Mapping[str, Any]) -> tuple[str, str, str]:
    return (
        str(value.get("record_kind", "")),
        str(value.get("program", "")),
        str(value.get("engine", "")),
    )


def _public_context(
    *,
    task: str,
    task_spec_sha256: str,
    observations: tuple[_XyzObservation, ...],
    result_observations: tuple[_ResultObservation, ...] = (),
    conformance_records: tuple[dict[str, Any], ...],
    registry_sha256: str,
    live_schema_sha256: str,
    execution_requested: bool,
    execution_available: bool,
    execution_review_requested: bool = False,
    bounded_execution_record: Mapping[str, Any] | None = None,
    approved_project_records: tuple[dict[str, Any], ...] = (),
    approved_workflow_record: Mapping[str, Any] | None = None,
    provider_record: Mapping[str, Any] | None = None,
    approved_identity_records: tuple[dict[str, Any], ...] = (),
    approved_input_records: tuple[dict[str, Any], ...] = (),
    analysis_completion_record: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    public_workflow = dict(approved_workflow_record or {})
    return {
        "schema_version": "chemsmart.live-agent-public-context.v1",
        "task": task,
        "task_spec_sha256": task_spec_sha256,
        "registry_sha256": registry_sha256,
        "live_cli_schema_sha256": live_schema_sha256,
        "artifacts": tuple(
            item.public_record()
            for item in (*observations, *result_observations)
        ),
        "approved_molecular_identities": approved_identity_records,
        "approved_molecular_inputs": approved_input_records,
        "approved_project_artifacts": approved_project_records,
        "approved_workflow": public_workflow,
        "provider": dict(provider_record or {}),
        "analysis_completion_policy": dict(analysis_completion_record or {}),
        "approved_execution_contract": _approved_execution_context(
            public_workflow
        ),
        "conformance_observations": conformance_records,
        "execution_requested": bool(execution_requested),
        "execution_tool_available": bool(execution_available),
        "execution_review_requested": bool(execution_review_requested),
        "bounded_execution_envelope": dict(bounded_execution_record or {}),
        "authority": (
            "Use artifact IDs and typed tool fields only. The host owns paths, "
            "CLI argv, project materialization, validation, approval, and execution."
        ),
    }


def _approved_execution_context(
    public_workflow: Mapping[str, Any],
) -> dict[str, Any]:
    """Describe execution mechanics without inventing scientific settings."""

    if public_workflow.get("authorization_mode") in {
        "bounded_local",
        "bounded_continuous",
    }:
        return {
            "authorization_mode": "bounded_local",
            "scientific_settings_authority": "task_and_project_yaml",
            "producer_input_policy": "validated_exact_data_edge_only",
            "operating_bounds": dict(
                public_workflow.get("operating_bounds") or {}
            ),
        }
    nodes = tuple(public_workflow.get("nodes", ()))
    if not nodes:
        return {}
    return {
        "programs": tuple(
            sorted({str(item.get("program", "")) for item in nodes})
        ),
        "node_order": tuple(public_workflow["node_order"]),
        "scientific_settings_authority": (
            "task_and_approved_project_artifacts"
        ),
        "producer_input_policy": "validated_approved_edge_only",
    }


def _workflow_context_sentence() -> str:
    """State dependency structure so the model need not recall it."""

    return (
        "The workflow planning tools return a host-derived workflow_context or "
        "workflow_frontier: per "
        "node it says whether that node is ready, waiting, or blocked, which "
        "upstream node and output each waiting input needs, and which nodes "
        "depend on it. Read it instead of reconstructing the DAG from memory, "
        "and re-read it after each artifact appears. It is an observation, "
        "not a claim you can make: never assert a node is ready. "
    )


def _system_prompt(
    approved_workflow: Mapping[str, Any] | None,
    *,
    bounded_review_requested: bool = False,
    skill_index: tuple[str, ...] = (),
) -> str:
    execution_available = bool(approved_workflow)
    bounded_execution = bool(
        approved_workflow
        and approved_workflow.get("authorization_mode")
        in {"bounded_local", "bounded_continuous"}
    )
    bounded_review = bool(bounded_review_requested or bounded_execution)
    execution_sentence = (
        "Execution is exposed as execute_approved_program_node(node_id) under "
        "science-free operating bounds. First express the complete scientific "
        "DAG through the normal workflow tools. Validate projects and obtain a "
        "green preview/preflight for every initial calculation node before the "
        "first execution request. The host then freezes only that current "
        "ChemSmart plan and its exact evidence; it refuses programs, engines, "
        "resources, dependencies, or calls outside operating_bounds. Execute "
        "producer nodes before their consumers and use only validated returned "
        "handoffs. Continue parsing and postprocessing in this same conversation."
        if bounded_execution
        else (
            "Execution is exposed only as execute_approved_program_node(node_id). "
            "Use the exact approved_workflow node IDs and listed order, and stop on "
            "any red result. When an approved producer edge exists, execute its "
            "producer first; after validation, consume the returned produced_handoffs "
            "artifact and scientific identity before compiling, previewing, "
            "preflighting, and executing that edge's consumer. "
            "Reuse the listed approved project artifact rather than promoting a new "
            "copy. Treat the execution receipt's deterministic validators as the live "
            "result gate; do not call inspect_calculation_artifact for newly executed "
            "nodes because this composition has no separate legacy settings/run IDs."
            if execution_available
            else (
                "Execution is not exposed. Finish with a useful planned or "
                "previewed state."
            )
        )
    )
    approval_readiness_sentence = (
        "The supplied operating bounds request an inert exact workflow review "
        "after this planning turn; they do not expose engine execution or "
        "human approval to the provider. Every initially runnable node needs "
        "a green preview before bounded "
        "execution. An exact producer-data target may instead appear as "
        "deferred_admissible until its validated upstream geometry exists; "
        "keep that causal stage in the workflow and execute its producer "
        "first. A release-unsupported stage may instead appear as "
        "non_executable: retain it as scientific intent, but it needs no "
        "green preview and will not be approved or launched. At least one "
        "release-executable stage is required for human execution review. "
        "Read approval_readiness for preview_required, deferred, and "
        "non-executable nodes. "
        if bounded_review
        else (
            "Every currently materialized node needs a green preview before "
            "an exact approval can execute. Read approval_readiness for the "
            "nodes that still require preview. If a named program remains "
            "unsupported after repairing its preview findings, use a "
            "scientifically defensible supported alternative while preserving "
            "any required unavailable observable as blocked_unsupported. "
        )
    )
    skill_sentence = ""
    if skill_index:
        listing = " ".join(f"({item})" for item in skill_index)
        skill_sentence = (
            " Domain-knowledge skills are available through "
            "consult_domain_skill(skill_id). Available now: "
            f"{listing}. Consult the relevant skill before you commit to a "
            "reporting convention, an electronic-state assignment, or a "
            "workflow shape it covers; before you judge whether a method, "
            "basis, solvation model or conformer sample can answer the "
            "question asked, or compare a computed value with experiment; and "
            "before you repair a rejected analysis node. Say when a stated "
            "fact came from a "
            "skill. A skill is advisory knowledge only: it never establishes "
            "readiness, approval, terminal state, or an accuracy claim, and it "
            "never replaces a typed host receipt."
        )
    return (
        "You are a professional computational-chemistry planning agent operating "
        "ChemSmart 3.1.4. Work plan-first through typed tools. Inspect program "
        "capability and environment, bind exact artifact identity, render and promote "
        "stage-specific project YAML, validate it, build a scientific tool-chain DAG, compile safe "
        "commands, and preview every currently resolvable node. Keep every future "
        "producer input unresolved until its validated upstream artifact exists. "
        "For every request that ends in a calculated or derived value, use "
        "plan_scientific_workflow to record any required calculations, result "
        "extraction, validation, mathematics, and claim rendering in one "
        "connected DAG. For analysis of registered results, leave "
        "calculation_nodes empty instead of inventing a placeholder program call. "
        "Its analysis inputs name future producer node/output pairs, so do not wait "
        "for artifact hashes before planning postprocessing. Preserve an unavailable "
        "parser or external analysis as blocked_unsupported instead of deleting the "
        "requested observable. Use plan_command_workflow only for a calculation-only "
        "compatibility task, and use inspect_workflow_frontier for host-derived "
        "next actions. "
        "When the task names a program, plan that program. If its preview is "
        "refused, repair it from the findings preview_command returns, which "
        "name the field, the expected value and the observed one. Only when a "
        "named program still cannot preview green should you use a "
        "scientifically defensible supported alternative. "
        + approval_readiness_sentence
        + _workflow_context_sentence()
        + "Never author native "
        "Gaussian, ORCA, xTB, or PySCF input/script text. Never invent coordinates, "
        "paths, shell syntax, evidence, readiness, or terminal state. "
        "The execution target is host policy: preview planning compiles local run, "
        "and an approved execution profile uses its frozen resource target. Never "
        "choose or infer run versus scheduler submission. "
        "Explain method rationale, alternatives, uncertainty, and diagnostics in "
        "concise public "
        "English. A molecular or state-specific geometry name is authorized only "
        "when public context contains its approved_molecular_identity record. "
        "Use only one of that record's approved_names, bind it to the record's "
        "exact geometry_sha256, and cite its evidence_ref in the scientific "
        "decision. File names, "
        "XYZ comments, element lists, project settings, and preview artifacts do "
        "not establish molecular identity. An approved molecular identity never "
        "establishes charge or multiplicity. A public approved_molecular_input "
        "record separately establishes its stated geometry_role, charge, and "
        "multiplicity only for the exact geometry_sha256 it names. The host has "
        "already bound that declared state; do not infer another initial state. "
        "Preserve explicit scientific dependencies, but do not convert a "
        "presentation sequence into a control edge. SP(initial geometry) and "
        "OPT(initial geometry) are siblings unless the request supplies a separate "
        "control dependency; only a producer output consumed downstream creates a "
        "data edge. Distinguish loader-supported, preview-conformant, "
        "environment-ready, and scientifically suitable. Never infer one status "
        "from another. Do not assert quantitative accuracy, cost, or "
        "density-fitting effects without typed evidence, and do not claim an RI/DF "
        "path unless the exact project explicitly enables density_fit. A project "
        "functional literal is the requested value, not proof of the applied XC "
        "interpretation. Use only the functional-resolution record returned by "
        "project validation for an alias or correlation-convention claim, cite "
        "its exact functional_resolution evidence_ref, and treat exact LibXC "
        "components as unknown until target-runtime materialization. That host "
        "resolution is not environment-readiness or scientific-suitability "
        "evidence. When project validation returns decision_binding, call "
        "record_scientific_decision after validation with its exact "
        "evidence_refs before rendering any applied XC alias or correlation "
        "convention; an earlier task-level decision is insufficient. "
        "Present an alternative as "
        "runnable only when the current project loader, command preview, and observed "
        "environment support it; otherwise label it as a scientifically relevant "
        "but unmaterialized alternative. PySCF project stage keys are exactly "
        "sp, opt, hess, and preview-only td; xTB project stage keys are exactly "
        "sp, opt, and hess. Gaussian and ORCA projects retain gas/solv phase "
        "sections: SP consumes solv when present, otherwise gas, and an "
        "explicit sp override takes precedence; physical solvation is enabled "
        "only by the solvent settings themselves. "
        "workflow node IDs may separately express initial or optimized geometry. "
        "For each job, pass the exact receipt_sha256 returned by that job's "
        "inspect_program_capability call into environment inspection and project "
        "validation, then use the engine binding returned by environment inspection. "
        "Do not substitute conformance, joined-capability, or environment receipt "
        "digests for those typed fields. Keep project artifact IDs distinct from "
        "geometry artifact IDs. Bind scientific identity only to a geometry_xyz "
        "artifact, never to a project, and do this before planning the workflow. "
        "Every workflow node must declare at least one expected output. If "
        "plan_command_workflow or plan_scientific_workflow returns findings or a null scientific_workflow_plan, "
        "repair the binding or DAG and call it again; a workflow_draft alone is "
        "not the typed scientific DAG. In workflow inputs, represent an initial artifact "
        "with empty producer_node_id and producer_output_id strings; represent a "
        "future optimized input with its producer IDs and no invented artifact ID. "
        "Omit absent optional settings instead of encoding them as the string none. "
        "When a planned workflow needs repairing, use amend_scientific_workflow "
        "rather than resubmitting the whole DAG: it repairs how a named part is "
        "expressed, including a corrected project promoted under a new artifact "
        "ID, an identifier, a unit, a declared quantity kind, or a selector, and "
        "preserves every node you do not name. Do not leave a repaired project "
        "detached from the final workflow. "
        "When an approved project artifact is supplied, read and validate that exact "
        "artifact instead of rerendering an equivalent project. "
        "If critical evidence is missing, identify it and block honestly. "
        "When public context contains a host-bound structured result, finish the "
        "scientific data path rather than stopping at a calculation plan: use "
        "extract_result_quantities for raw observables, derive_thermochemistry "
        "with explicit temperature and pressure for RRHO quantities, and "
        "evaluate_quantity_expression for requested arithmetic or geometric "
        "derivations. Local input and intermediate node IDs are presentation-only; "
        "the host grades an identifier-independent symbolic DAG. When two inputs "
        "carry the same source quantity ID, the host tells them apart by each "
        "input's own id, which is already unique within the expression, so "
        "semantic_role is optional; supply one only to label an occurrence more "
        "readably than its input id does, and never use a receipt hash as a "
        "role. When a numerical condition already exists "
        "as a quantity on a typed receipt, reference that receipt quantity instead "
        "of duplicating it as a literal; use a literal only when no typed source "
        "quantity exists. "
        "When the planned frontier exposes a scientific_validation node, use "
        "evaluate_scientific_validation with the exact upstream typed receipt "
        "quantities. The host evaluates the already-declared rules and returns "
        "a typed verdict; a prose decision does not execute validation. "
        "Use record_analysis_claims to bind each requested reported "
        "number and display unit to an exact receipt quantity; the host, not the "
        "model, supplies the value. The host renders the authoritative final numeric "
        "section from that claim record. Report only those host-rendered claim values. "
        "Keep receipt IDs, digests, and artifact hashes internal unless the user "
        "explicitly asks for an audit; the public answer should explain the chemistry, "
        "evidence stage, and limitations rather than reciting bookkeeping. "
        "Never copy a "
        "paper's hidden target value into a tool call, and never replace a required "
        "target-producing calculation or postprocessing step by deleting it from "
        "the plan. If a result artifact is absent, leave postprocessing planned and "
        "state exactly which producer artifact is required. A structured result's "
        "requested/applied functional distinction may be cited only through its "
        "exact result_functional_resolution evidence_ref from public context; do "
        "not require a new project-validation receipt merely to analyze an existing "
        "result. When public context contains analysis_completion_policy, complete "
        "every listed stage and cite each extraction, thermochemistry, expression, "
        "and analysis-claim "
        "receipt in the final scientific decision by passing its exact digest in "
        "postprocessing_receipt_sha256s rather than constructing a free-form receipt "
        "label; the host, not the model, decides "
        "whether that task-owned policy passed. "
        + execution_sentence
        + skill_sentence
    )


def _provider_budget(
    provider_profile: AgentProviderProfileV1,
    *,
    max_concurrency: int = 1,
    wall_time_seconds: float = _SESSION_WALL_TIME_SECONDS,
) -> ProviderNetworkBudgetV1:
    body = {
        "schema_version": "chemsmart.provider-network-budget.v1",
        "allowed_provider": provider_profile.provider,
        "endpoint_origin": provider_profile.endpoint,
        "max_concurrency": int(max_concurrency),
        "max_input_tokens_per_request": provider_profile.context_tokens,
        "max_output_tokens_per_request": provider_profile.max_output_tokens,
        "task_wall_time_seconds": float(wall_time_seconds),
    }
    return ProviderNetworkBudgetV1(
        **body, budget_sha256=canonical_sha256(body)
    )


def _request_context(
    *,
    messages: list[dict[str, str]],
    task_spec_sha256: str,
    tool_schema_sha256: str,
    provider_budget: ProviderNetworkBudgetV1,
    provider_profile_sha256: str,
    execution_requested: bool,
    dependency_context: TaskDependencyContextV2 | None = None,
) -> RequestContextProvenanceV1:
    configuration_sha256 = canonical_sha256(
        {
            "provider_profile_sha256": provider_profile_sha256,
            "execution_requested": bool(execution_requested),
            "dependency_context_sha256": (
                dependency_context.context_sha256
                if dependency_context is not None
                else ""
            ),
        }
    )
    return build_request_context_provenance(
        task_spec_sha256=task_spec_sha256,
        prompt_sha256=canonical_sha256(messages),
        tool_schema_sha256=tool_schema_sha256,
        configuration_sha256=configuration_sha256,
        provider_budget_sha256=provider_budget.budget_sha256,
    )


def _task_envelope(
    *,
    session_id: str,
    messages: list[dict[str, str]],
    task_spec_sha256: str,
    workspace_records: tuple[dict[str, Any], ...],
    tool_schema_sha256: str,
    execution_enabled: bool,
    provider_profile: AgentProviderProfileV1,
    wall_time_seconds: float = _SESSION_WALL_TIME_SECONDS,
    chemistry_engine_calls: int | None = None,
) -> TaskEnvelopeV1:
    engine_calls = (
        int(chemistry_engine_calls)
        if chemistry_engine_calls is not None
        else (4 if execution_enabled else 0)
    )
    resource = ResourceBudgetV1(
        max_input_tokens_per_request=provider_profile.context_tokens,
        max_output_tokens_per_request=provider_profile.max_output_tokens,
        max_tool_calls=_MAX_TOOL_CALLS,
        wall_time_seconds=float(wall_time_seconds),
        chemistry_engine_calls=engine_calls,
        hpc_calls=0,
    )
    body = {
        "schema_version": "chemsmart.task-envelope.v1",
        "task_id": f"live-task-{task_spec_sha256[:16]}",
        "session_id": session_id,
        "turn_id": session_id + ".turn-1",
        "request_sha256": canonical_sha256(messages),
        "cwd_sha256": canonical_sha256(workspace_records),
        "phase": TaskPhase.ROUTE,
        "budget": resource,
        "tool_schema_sha256": tool_schema_sha256,
    }
    return TaskEnvelopeV1(**body, envelope_sha256=canonical_sha256(body))


def _execution_composition_available() -> bool:
    parameters = inspect.signature(CommandCompiledToolHostV1).parameters
    required = {
        "approved_workspace",
        "workflow_execution_approval",
        "execution_resources",
        "execution_server",
        "execution_environment",
    }
    return required.issubset(parameters)


def _validate_bounded_envelope_against_registry(
    envelope: BoundedExecutionEnvelopeV1,
    *,
    registry: ProgramCapabilityRegistryV1,
) -> None:
    """Reject an allowlist that the installed ChemSmart cannot represent."""

    result_validators = {"gaussian", "orca", "pyscf", "xtb"}
    for program, engines in envelope.allowed_program_engines:
        capability = registry.get(program)
        if capability is None:
            raise ContractError(
                f"bounded execution allows unknown program {program!r}"
            )
        unknown = sorted(set(engines).difference(capability.engines))
        if unknown:
            raise ContractError(
                f"bounded execution allows unsupported {program} engine(s): "
                + ", ".join(unknown)
            )
        if program not in result_validators:
            raise ContractError(
                "bounded local execution lacks a result validator for "
                f"program {program!r}"
            )
        executable_engines = {
            engine
            for engine, _jobtype in capability.execution_engine_job_pairs
        }
        unavailable = sorted(set(engines).difference(executable_engines))
        if unavailable:
            raise ContractError(
                f"bounded execution has no executable jobs for {program} "
                "engine(s): " + ", ".join(unavailable)
            )


def _bounded_execution_composition_inputs(
    *,
    workspace: Path,
    run_directory: Path,
    envelope: BoundedExecutionEnvelopeV1,
) -> dict[str, Any]:
    """Compose a deferred execution host from user-owned operating bounds."""

    requested_scratch_root = Path(envelope.scratch_root)
    if requested_scratch_root.is_symlink():
        raise ContractError(
            "bounded execution scratch root cannot be a symlink"
        )
    scratch_root = requested_scratch_root.resolve()
    scratch_root.mkdir(parents=True, exist_ok=True)
    server_profile = _write_execution_server_profile(
        run_directory,
        envelope.resources,
        scratch_root=scratch_root,
    )
    path_value = os.environ.get("PATH", "")
    xtb_executable = os.environ.get(
        "CHEMSMART_XTB_EXECUTABLE"
    ) or shutil.which("xtb")
    executable_directory = (
        str(Path(xtb_executable).expanduser().parent) if xtb_executable else ""
    )
    environment = {
        "PATH": (
            path_value
            if not executable_directory
            else (
                executable_directory
                if not path_value
                else executable_directory + os.pathsep + path_value
            )
        ),
        "PYTHONNOUSERSITE": "1",
    }
    return {
        "approved_workspace": workspace,
        "execution_resources": envelope.resources,
        "execution_server": str(server_profile),
        "execution_environment": environment,
        "bounded_execution_envelope": envelope,
    }


def _execution_composition_inputs(
    *,
    host_type: type[CommandCompiledToolHostV1],
    workspace: Path,
    run_directory: Path,
    approval_file: Path,
    task_spec_sha256: str,
) -> dict[str, Any]:
    """Load one exact workflow approval and its explicit resource profile."""

    if not _execution_composition_available():
        raise ContractError("approval-bound execution host is unavailable")
    expected_parameters = {
        "approved_workspace",
        "execution_resources",
        "workflow_execution_approval",
        "execution_server",
        "execution_environment",
    }
    if not expected_parameters.issubset(
        inspect.signature(host_type).parameters
    ):
        raise ContractError("execution host composition API is incomplete")
    if not approval_file.is_file() or approval_file.is_symlink():
        raise ContractError("approval file must be a current regular file")
    try:
        payload = json.loads(approval_file.read_text(encoding="utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ContractError("approval file is not valid JSON") from exc
    if not isinstance(payload, Mapping):
        raise ContractError("approval file root must be an object")
    if "workflow_execution_approval_bundle" in payload:
        bundle = load_workflow_execution_approval_bundle(approval_file)
        approval = bundle.workflow_approval
        frozen = bundle.frozen_workflow_approval
        resources = bundle.execution_resources
        if Path(approval.workspace).resolve() != workspace:
            raise ContractError(
                "execution approval bundle targets another workspace"
            )
        if task_spec_sha256 and approval.task_spec_sha256 != task_spec_sha256:
            raise ContractError(
                "execution approval bundle targets another task spec"
            )
        envelope = _parse_bounded_execution_envelope_record(
            bundle.execution_envelope,
            resources=resources,
        )
        requested_scratch_root = Path(envelope.scratch_root)
        if requested_scratch_root.is_symlink():
            raise ContractError("approved scratch root cannot be a symlink")
        scratch_root = requested_scratch_root.resolve()
        scratch_root.mkdir(parents=True, exist_ok=True)
        server_profile = _write_execution_server_profile(
            run_directory,
            resources,
            scratch_root=scratch_root,
        )
        path_value = os.environ.get("PATH", "")
        xtb_executable = os.environ.get(
            "CHEMSMART_XTB_EXECUTABLE"
        ) or shutil.which("xtb")
        executable_directory = (
            str(Path(xtb_executable).expanduser().parent)
            if xtb_executable
            else ""
        )
        environment = {
            "PATH": (
                path_value
                if not executable_directory
                else (
                    executable_directory
                    if not path_value
                    else executable_directory + os.pathsep + path_value
                )
            ),
            "PYTHONNOUSERSITE": "1",
        }
        approved_projects = _approved_project_artifacts(workspace, approval)
        return {
            "approved_workspace": workspace,
            "execution_resources": resources,
            "workflow_execution_approval": approval,
            "frozen_workflow_approval": frozen,
            "execution_server": str(server_profile),
            "execution_environment": environment,
            "approved_project_artifacts": approved_projects,
            "approved_scientific_plan": bundle.approved_scientific_plan,
            "approved_materialized_workflow": (
                bundle.approved_materialized_workflow
            ),
            "approved_environment_identities": (
                bundle.approved_environment_identities
            ),
            "stationary_point_policy": bundle.stationary_point_policy,
        }
    raw_approval = payload.get("workflow_approval", payload)
    if not isinstance(raw_approval, Mapping):
        raise ContractError("workflow approval must be an object")
    approval = _parse_workflow_approval(raw_approval)
    raw_resources = payload.get("execution_resources")
    if not isinstance(raw_resources, Mapping):
        raise ContractError(
            "approval file requires explicit execution_resources"
        )
    resources = _parse_execution_resources(raw_resources)
    if Path(approval.workspace).resolve() != workspace:
        raise ContractError("workflow approval targets another workspace")
    if approval.task_spec_sha256 != task_spec_sha256:
        raise ContractError("workflow approval targets another task spec")
    if approval.resource_sha256 != resources.resource_sha256:
        raise ContractError(
            "workflow approval resources differ from supplied resources"
        )
    server_profile = _write_execution_server_profile(run_directory, resources)
    path_value = os.environ.get("PATH", "")
    xtb_executable = os.environ.get(
        "CHEMSMART_XTB_EXECUTABLE"
    ) or shutil.which("xtb")
    executable_directory = (
        str(Path(xtb_executable).expanduser().parent) if xtb_executable else ""
    )
    environment = {
        "PATH": (
            path_value
            if not executable_directory
            else (
                executable_directory
                if not path_value
                else executable_directory + os.pathsep + path_value
            )
        ),
        "PYTHONNOUSERSITE": "1",
    }
    approved_projects = _approved_project_artifacts(workspace, approval)
    composed: dict[str, Any] = {
        "approved_workspace": workspace,
        "execution_resources": resources,
        "workflow_execution_approval": approval,
        "execution_server": str(server_profile),
        "execution_environment": environment,
        "approved_project_artifacts": approved_projects,
    }
    # The execution tool requires the Runtime V2 body as well; without it the
    # V1 approval is preview-only and no node can run. Carrying it is what
    # makes ``agent run --approval-file`` a reachable path rather than an
    # advertised one. Its absence stays legal so an approval that deliberately
    # authorises preview alone keeps working.
    raw_frozen = payload.get("frozen_workflow_approval")
    composed["approved_scientific_plan"] = None
    composed["approved_materialized_workflow"] = None
    if raw_frozen is not None:
        if not isinstance(raw_frozen, Mapping):
            raise ContractError("frozen workflow approval must be an object")
        frozen = _parse_frozen_workflow_approval(raw_frozen)
        if frozen.task_spec_sha256 != task_spec_sha256:
            raise ContractError(
                "frozen workflow approval targets another task spec"
            )
        composed["frozen_workflow_approval"] = frozen
        # The plan body is optional, and self-verifying: it is disclosed only
        # when it hashes to the digest the frozen approval already pins, so a
        # reviewer cannot show the session one plan while approving another.
        raw_plan = payload.get("approved_scientific_plan")
        if raw_plan is not None:
            if not isinstance(raw_plan, Mapping):
                raise ContractError(
                    "approved scientific plan must be an object"
                )
            plan = _parse_scientific_workflow_plan(raw_plan)
            if plan.plan_sha256 != frozen.plan_sha256:
                raise ContractError(
                    "approved scientific plan is not the plan that was frozen"
                )
            composed["approved_scientific_plan"] = plan
        raw_material = payload.get("approved_materialized_workflow")
        if raw_material is not None:
            if not isinstance(raw_material, Mapping):
                raise ContractError(
                    "approved materialized workflow must be an object"
                )
            material = _parse_materialized_workflow(raw_material)
            if material.materialized_sha256 != (
                frozen.materialized_workflow_sha256
            ):
                raise ContractError(
                    "approved materialization is not the one that was frozen"
                )
            composed["approved_materialized_workflow"] = material
        # The environment the reviewer approved, identified by machine rather
        # than by receipt digest. A receipt digest folds in its capability
        # receipt, which changes with the active overlay, so the digest a plan
        # session records is never the one execution computes even on the same
        # interpreter. The reviewer holds the plan session's receipts and can
        # state their identity; without it, only an exact digest match passes.
        raw_identities = payload.get("approved_environment_identities")
        if raw_identities is not None:
            if not isinstance(raw_identities, (list, tuple)):
                raise ContractError(
                    "approved environment identities must be a list"
                )
            composed["approved_environment_identities"] = tuple(
                require_sha256(item, "approved environment identity")
                for item in raw_identities
            )
    return composed


def _parse_materialized_workflow(
    value: Mapping[str, Any],
) -> MaterializedWorkflowV1:
    """Rebuild the materialization the frozen approval pins."""

    raw = dict(value)
    try:
        nodes = []
        for item in raw.get("nodes", ()):
            node = dict(item)
            node["auxiliary_input_bindings"] = tuple(
                AuxiliaryArtifactBindingV1(**dict(binding))
                for binding in node.get("auxiliary_input_bindings", ())
            )
            nodes.append(MaterializedNodeV1(**node))
        raw["nodes"] = tuple(nodes)
        raw["unresolved_node_ids"] = tuple(raw.get("unresolved_node_ids", ()))
        return MaterializedWorkflowV1(**raw)
    except (KeyError, TypeError, ValueError) as exc:
        raise ContractError(
            "approved materialized workflow does not match the v1 schema"
        ) from exc


def _parse_scientific_workflow_plan(
    value: Mapping[str, Any],
) -> ScientificWorkflowPlanV2:
    """Rebuild the approved plan a session is required to reproduce."""

    def _tuples(mapping: Mapping[str, Any]) -> dict[str, Any]:
        return {
            key: tuple(item) if isinstance(item, list) else item
            for key, item in mapping.items()
        }

    raw = dict(value)
    try:
        raw["nodes"] = tuple(
            ScientificWorkflowNodeV2(**_tuples(item))
            for item in raw.get("nodes", ())
        )
        raw["edges"] = tuple(
            ScientificWorkflowEdgeV2(**_tuples(item))
            for item in raw.get("edges", ())
        )
        raw["complexity_factors"] = tuple(raw.get("complexity_factors", ()))
        raw["required_observables"] = tuple(
            raw.get("required_observables", ())
        )
        return ScientificWorkflowPlanV2(**raw)
    except (KeyError, TypeError, ValueError) as exc:
        raise ContractError(
            "approved scientific plan does not match the v2 schema"
        ) from exc


def _approved_project_artifacts(
    workspace: Path, approval: WorkflowExecutionApprovalV1
) -> tuple[TrustedArtifactRefV1, ...]:
    """Resolve workflow-bound project bytes without model rematerialization."""

    required = {
        item.project_artifact_sha256 for item in approval.node_bindings
    }
    matches: dict[str, Path] = {}
    project_roots = [workspace / "projects"]
    run_root = workspace / ".chemsmart-agent" / "runs"
    if run_root.is_dir() and not run_root.is_symlink():
        project_roots.extend(sorted(run_root.glob("*/projects")))
    for project_root in project_roots:
        if not project_root.is_dir() or project_root.is_symlink():
            continue
        for candidate in sorted(project_root.glob("*.yaml")):
            if candidate.is_symlink() or not candidate.is_file():
                continue
            digest = file_sha256(candidate)
            if digest in required:
                # Approval binds project content, not a presentation filename.
                # A prior agent session may legitimately promote the same
                # validated YAML under another role-specific name.  Keep the
                # first sorted regular-file match so an identical copy does
                # not make an approved workflow impossible to replay.
                matches.setdefault(digest, candidate.resolve())
    missing = sorted(required.difference(matches))
    if missing:
        raise ContractError(
            "workflow approval project artifact is unavailable"
        )
    return tuple(
        TrustedArtifactRefV1(
            artifact_id=f"project-approved-{digest[:16]}",
            kind="project_yaml",
            sha256=digest,
            size_bytes=matches[digest].stat().st_size,
            path=str(matches[digest]),
            cli_value=str(matches[digest]),
        )
        for digest in sorted(required)
    )


def _public_workflow_approval(
    approval: WorkflowExecutionApprovalV1,
    *,
    approved_plan: ScientificWorkflowPlanV2 | None = None,
) -> dict[str, Any]:
    """Expose the reviewed plan while preserving exact digest matching."""

    if approved_plan is not None:
        public_plan = canonical_data(approved_plan)
        for node in public_plan["nodes"]:
            if node.get("charge") is None:
                # The model-facing schema accepts an omitted optional state,
                # not JSON null.  The plan digest projection likewise omits
                # both fields when absent, so this remains an exact and
                # sufficient reproduction of a legacy/inherited-state node.
                node.pop("charge", None)
                node.pop("multiplicity", None)
        return {
            **_public_workflow_approval(approval),
            "approved_scientific_plan": public_plan,
            "approved_plan_sha256": approved_plan.plan_sha256,
            "plan_reproduction_rule": (
                "plan_scientific_workflow must reproduce "
                "approved_scientific_plan exactly; execution refuses any plan "
                "whose digest differs from approved_plan_sha256"
            ),
        }
    return {
        "schema_version": "chemsmart.public-approved-workflow.v1",
        "workflow_id": approval.workflow_id,
        "node_order": tuple(item.node_id for item in approval.node_bindings),
        "nodes": tuple(
            {
                "node_id": item.node_id,
                "program": item.program,
                "engine": item.engine,
                "jobtype": item.jobtype,
                "input_mode": item.input_mode,
            }
            for item in approval.node_bindings
        ),
        "producer_edges": tuple(
            {
                "producer_node_id": item.producer_node_id,
                "consumer_node_id": item.consumer_node_id,
                "artifact_kind": item.artifact_kind,
                "selection_rule": item.selection_rule,
            }
            for item in approval.producer_edges
        ),
        "status": approval.status,
    }


def _parse_frozen_workflow_approval(
    value: Mapping[str, Any],
) -> FrozenWorkflowApprovalV1:
    """Parse the frozen approval and preserve its independently reviewed DAG."""

    raw = dict(value)
    try:
        raw["environment_identity_sha256s"] = tuple(
            raw.get("environment_identity_sha256s", ())
        )
        raw["approved_node_ids"] = tuple(raw.get("approved_node_ids", ()))
        raw["producer_edge_sha256s"] = tuple(
            raw.get("producer_edge_sha256s", ())
        )
        previews = []
        for item in raw.get("materialized_preview_bindings", ()):
            preview = dict(item)
            preview["auxiliary_input_bindings"] = tuple(
                AuxiliaryArtifactBindingV1(**dict(binding))
                for binding in preview.get("auxiliary_input_bindings", ())
            )
            previews.append(FrozenMaterializedNodePreviewV1(**preview))
        raw["materialized_preview_bindings"] = tuple(previews)
        raw["producer_edge_rules"] = tuple(
            FrozenProducerEdgeRuleV1(**dict(item))
            for item in raw.get("producer_edge_rules", ())
        )
        return FrozenWorkflowApprovalV1(**raw)
    except (KeyError, TypeError, ValueError) as exc:
        raise ContractError(
            "frozen workflow approval does not match the v1 schema"
        ) from exc


def _parse_workflow_approval(
    value: Mapping[str, Any],
) -> WorkflowExecutionApprovalV1:
    raw = dict(value)
    try:
        nodes = []
        for item in raw.get("node_bindings", ()):
            node = dict(item)
            node["auxiliary_input_bindings"] = tuple(
                AuxiliaryArtifactBindingV1(**dict(binding))
                for binding in node.get("auxiliary_input_bindings", ())
            )
            nodes.append(ApprovedNodeBindingV1(**node))
        raw["node_bindings"] = tuple(nodes)
        raw["producer_edges"] = tuple(
            ProducerEdgeRuleV1(**dict(item))
            for item in raw.get("producer_edges", ())
        )
        return WorkflowExecutionApprovalV1(**raw)
    except (KeyError, TypeError, ValueError) as exc:
        raise ContractError(
            "workflow approval does not match the v1 schema"
        ) from exc


def _parse_workflow_approval_request(
    value: Mapping[str, Any],
) -> WorkflowApprovalRequestV1:
    raw = dict(value)
    try:
        nodes = []
        for item in raw.get("node_bindings", ()):
            node = dict(item)
            node["auxiliary_input_bindings"] = tuple(
                AuxiliaryArtifactBindingV1(**dict(binding))
                for binding in node.get("auxiliary_input_bindings", ())
            )
            nodes.append(ApprovedNodeBindingV1(**node))
        raw["node_bindings"] = tuple(nodes)
        raw["producer_edges"] = tuple(
            ProducerEdgeRuleV1(**dict(item))
            for item in raw.get("producer_edges", ())
        )
        return WorkflowApprovalRequestV1(**raw)
    except (KeyError, TypeError, ValueError) as exc:
        raise ContractError(
            "workflow approval request does not match the v1 schema"
        ) from exc


def _parse_scientific_toolchain_plan(
    value: Any,
) -> "ScientificToolchainPlanV1 | None":
    """Rebuild the digest-bound analysis chain from a review or bundle."""

    if value is None:
        return None
    if not isinstance(value, Mapping):
        raise ContractError(
            "scientific toolchain plan must be an object or null"
        )
    from chemsmart.agent.scientific_toolchain import (
        AnalysisInputIntentV1,
        AnalysisNodeIntentV1,
        AnalysisOutputIntentV1,
        AnalysisSelectorIntentV1,
        AnalysisValidationRuleIntentV1,
        RegisteredResultInputIntentV1,
        ScientificToolchainPlanV1,
    )

    def _input(item: Mapping[str, Any]) -> Any:
        record = dict(item)
        if record.get("source_kind") == "registered_result":
            return RegisteredResultInputIntentV1(**record)
        return AnalysisInputIntentV1(**record)

    try:
        raw = dict(value)
        raw["analysis_nodes"] = tuple(
            AnalysisNodeIntentV1(
                **{
                    **dict(node),
                    "inputs": tuple(
                        _input(item) for item in node.get("inputs", ())
                    ),
                    "selectors": tuple(
                        AnalysisSelectorIntentV1(**dict(item))
                        for item in node.get("selectors", ())
                    ),
                    "outputs": tuple(
                        AnalysisOutputIntentV1(**dict(item))
                        for item in node.get("outputs", ())
                    ),
                    "expression_nodes": tuple(
                        canonical_data(dict(item))
                        for item in node.get("expression_nodes", ())
                    ),
                    "validation_rules": tuple(
                        AnalysisValidationRuleIntentV1(**dict(item))
                        for item in node.get("validation_rules", ())
                    ),
                }
            )
            for node in raw.get("analysis_nodes", ())
        )
        for key in (
            "calculation_node_ids",
            "node_order",
            "required_output_ids",
        ):
            raw[key] = tuple(raw.get(key, ()))
        raw["calculation_observables"] = tuple(
            (str(node_id), tuple(observables))
            for node_id, observables in raw.get("calculation_observables", ())
        )
        return ScientificToolchainPlanV1(**raw)
    except ContractError:
        raise
    except (KeyError, TypeError, ValueError) as exc:
        raise ContractError(
            "scientific toolchain plan does not match the v1 schema"
        ) from exc


def _parse_stationary_point_policy(
    value: Any,
) -> StationaryPointValidationPolicyV1 | None:
    if value is None:
        return None
    if not isinstance(value, Mapping):
        raise ContractError(
            "stationary-point policy must be an object or null"
        )
    try:
        return StationaryPointValidationPolicyV1(**dict(value))
    except (KeyError, TypeError, ValueError) as exc:
        raise ContractError(
            "stationary-point policy does not match the v1 schema"
        ) from exc


def _parse_bounded_execution_envelope_record(
    value: Any,
    *,
    resources: ExecutionResourceSpecV1,
) -> BoundedExecutionEnvelopeV1:
    if not isinstance(value, Mapping):
        raise ContractError("execution envelope record must be an object")
    raw = dict(value)
    try:
        raw_resources = raw.get("resources")
        if not isinstance(raw_resources, Mapping):
            raise ContractError("execution envelope requires resources")
        parsed_resources = _parse_execution_resources(raw_resources)
        if parsed_resources != resources:
            raise ContractError(
                "execution envelope resources differ from reviewed resources"
            )
        raw_allowlist = raw.get("allowed_program_engines")
        if not isinstance(raw_allowlist, (list, tuple)):
            raise ContractError(
                "execution envelope program allowlist must be a sequence"
            )
        allowlist = tuple(
            (str(item[0]), tuple(str(engine) for engine in item[1]))
            for item in raw_allowlist
        )
        parsed = BoundedExecutionEnvelopeV1(
            schema_version=str(raw["schema_version"]),
            mode=str(raw["mode"]),
            allowed_program_engines=allowlist,
            resources=parsed_resources,
            episode_wall_time_seconds=int(raw["episode_wall_time_seconds"]),
            postprocess_reserve_seconds=int(
                raw["postprocess_reserve_seconds"]
            ),
            max_engine_calls=int(raw["max_engine_calls"]),
            scratch_root=str(raw["scratch_root"]),
        )
    except ContractError:
        raise
    except (IndexError, KeyError, TypeError, ValueError) as exc:
        raise ContractError(
            "execution envelope record does not match the v1 schema"
        ) from exc
    if canonical_data(parsed) != canonical_data(raw):
        raise ContractError(
            "execution envelope record contains unknown fields"
        )
    return parsed


def _read_exact_json_document(
    path: str | Path, *, label: str
) -> Mapping[str, Any]:
    source = Path(path)
    if not source.is_file() or source.is_symlink():
        raise ContractError(f"{label} must be a current regular file")
    try:
        payload = json.loads(source.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ContractError(f"{label} is not valid UTF-8 JSON") from exc
    if not isinstance(payload, Mapping):
        raise ContractError(f"{label} root must be an object")
    return payload


def load_workflow_execution_review(
    path: str | Path,
) -> WorkflowExecutionReviewV1:
    """Load and fully re-hash one inert review packet."""

    payload = _read_exact_json_document(path, label="execution review")
    raw_value = payload.get("workflow_execution_review", payload)
    if not isinstance(raw_value, Mapping):
        raise ContractError("workflow execution review must be an object")
    raw = dict(raw_value)
    try:
        request = raw.get("request")
        plan = raw.get("scientific_plan")
        materialized = raw.get("materialized_workflow")
        resources = raw.get("execution_resources")
        if not all(
            isinstance(item, Mapping)
            for item in (request, plan, materialized, resources)
        ):
            raise ContractError(
                "execution review has incomplete typed records"
            )
        raw["request"] = _parse_workflow_approval_request(request)
        raw["scientific_plan"] = _parse_scientific_workflow_plan(plan)
        raw["materialized_workflow"] = _parse_materialized_workflow(
            materialized
        )
        raw["execution_resources"] = _parse_execution_resources(resources)
        raw["execution_envelope"] = canonical_data(
            _parse_bounded_execution_envelope_record(
                raw.get("execution_envelope"),
                resources=raw["execution_resources"],
            )
        )
        raw["environment_bindings"] = tuple(
            WorkflowEnvironmentBindingV1(**dict(item))
            for item in raw.get("environment_bindings", ())
        )
        raw["node_reviews"] = tuple(
            WorkflowExecutionNodeReviewV1(**dict(item))
            for item in raw.get("node_reviews", ())
        )
        raw["stationary_point_policy"] = _parse_stationary_point_policy(
            raw.get("stationary_point_policy")
        )
        if raw.get("scientific_toolchain_plan") is not None:
            raw["scientific_toolchain_plan"] = (
                _parse_scientific_toolchain_plan(
                    raw.get("scientific_toolchain_plan")
                )
            )
        return WorkflowExecutionReviewV1(**raw)
    except ContractError:
        raise
    except (KeyError, TypeError, ValueError) as exc:
        raise ContractError(
            "workflow execution review does not match the v1 schema"
        ) from exc


def load_workflow_execution_approval_bundle(
    path: str | Path,
    *,
    expected_file_sha256: str = "",
) -> WorkflowExecutionApprovalBundleV1:
    """Load one exact compound bundle from the bytes that were confirmed."""

    source = Path(path)
    if not source.is_file() or source.is_symlink():
        raise ContractError(
            "execution approval bundle must be a current regular file"
        )
    try:
        document = source.read_bytes()
    except OSError as exc:
        raise ContractError(
            "execution approval bundle cannot be read"
        ) from exc
    if expected_file_sha256 and (
        hashlib.sha256(document).hexdigest()
        != str(expected_file_sha256).strip().lower()
    ):
        raise ContractError(
            "execution approval bytes differ from the human-confirmed file"
        )
    try:
        payload = json.loads(document.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ContractError(
            "execution approval bundle is not valid UTF-8 JSON"
        ) from exc
    if not isinstance(payload, Mapping):
        raise ContractError("execution approval bundle root must be an object")
    raw_value = payload.get("workflow_execution_approval_bundle", payload)
    if not isinstance(raw_value, Mapping):
        raise ContractError(
            "workflow execution approval bundle must be an object"
        )
    raw = dict(raw_value)
    try:
        for field in (
            "resolution",
            "workflow_approval",
            "frozen_workflow_approval",
            "approved_scientific_plan",
            "approved_materialized_workflow",
            "execution_resources",
        ):
            if not isinstance(raw.get(field), Mapping):
                raise ContractError(f"execution bundle requires {field}")
        raw["resolution"] = WorkflowReviewResolutionV1(
            **dict(raw["resolution"])
        )
        raw["workflow_approval"] = _parse_workflow_approval(
            raw["workflow_approval"]
        )
        raw["frozen_workflow_approval"] = _parse_frozen_workflow_approval(
            raw["frozen_workflow_approval"]
        )
        raw["approved_scientific_plan"] = _parse_scientific_workflow_plan(
            raw["approved_scientific_plan"]
        )
        raw["approved_materialized_workflow"] = _parse_materialized_workflow(
            raw["approved_materialized_workflow"]
        )
        raw["execution_resources"] = _parse_execution_resources(
            raw["execution_resources"]
        )
        raw["execution_envelope"] = canonical_data(
            _parse_bounded_execution_envelope_record(
                raw.get("execution_envelope"),
                resources=raw["execution_resources"],
            )
        )
        raw["approved_environment_identities"] = tuple(
            raw.get("approved_environment_identities", ())
        )
        raw["node_reviews"] = tuple(
            WorkflowExecutionNodeReviewV1(**dict(item))
            for item in raw.get("node_reviews", ())
        )
        raw["stationary_point_policy"] = _parse_stationary_point_policy(
            raw.get("stationary_point_policy")
        )
        if raw.get("scientific_toolchain_plan") is not None:
            raw["scientific_toolchain_plan"] = (
                _parse_scientific_toolchain_plan(
                    raw.get("scientific_toolchain_plan")
                )
            )
        return WorkflowExecutionApprovalBundleV1(**raw)
    except ContractError:
        raise
    except (KeyError, TypeError, ValueError) as exc:
        raise ContractError(
            "workflow execution approval bundle does not match the v1 schema"
        ) from exc


def write_workflow_execution_review(
    review: WorkflowExecutionReviewV1,
    path: str | Path,
) -> Path:
    """Create the inert review artifact once; never overwrite another body."""

    target = Path(path).expanduser()
    if not target.is_absolute():
        raise ContractError("execution review output path must be absolute")
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.parent.is_symlink():
        raise ContractError(
            "execution review output parent cannot be a symlink"
        )
    _write_private_exact(
        target,
        workflow_execution_review_json(review).encode("utf-8"),
    )
    return target.resolve()


def write_workflow_execution_approval_bundle(
    bundle: WorkflowExecutionApprovalBundleV1,
    path: str | Path,
) -> Path:
    """Create one executable bundle once, refusing every overwrite."""

    target = Path(path).expanduser()
    if not target.is_absolute():
        raise ContractError("execution approval output path must be absolute")
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.parent.is_symlink():
        raise ContractError(
            "execution approval output parent cannot be a symlink"
        )
    _write_private_exact(
        target,
        workflow_execution_approval_bundle_json(bundle).encode("utf-8"),
    )
    return target.resolve()


def resolve_workflow_execution_review(
    *,
    review_file: str | Path,
    reviewed_sha256: str,
    decision: str,
    actor: str,
    output_file: str | Path | None = None,
    decision_log: str | Path | None = None,
    approval_id: str = "",
) -> tuple[
    WorkflowReviewResolutionV1, WorkflowExecutionApprovalBundleV1 | None
]:
    """Human-only library entry point shared by CLI and TUI.

    The caller must supply the complete digest displayed to the human.  This
    function never prompts, consults a model, or starts an engine; it records
    the resolution and, only for ``approve``, creates one exact bundle.
    """

    from chemsmart.agent.execution import (
        approve_workflow_execution_review,
        build_workflow_review_resolution,
    )

    normalized = str(decision).strip().lower()
    reviewed = require_sha256(reviewed_sha256, "reviewed_sha256")
    review_path = Path(review_file).expanduser().resolve()
    review = load_workflow_execution_review(review_path)
    if reviewed != review.review_sha256:
        raise ContractError(
            "reviewed digest does not match the exact review packet"
        )
    if normalized == "approve" and output_file is None:
        raise ContractError("approval requires an output file")
    if normalized != "approve" and output_file is not None:
        raise ContractError("a non-approval decision cannot create a bundle")
    chosen_approval_id = (
        str(approval_id).strip() or "approval-" + review.review_sha256[:16]
    )
    resolution_id = (
        "resolution-" + review.review_sha256[:16] + "-" + normalized
    )
    if normalized == "approve":
        bundle = approve_workflow_execution_review(
            review,
            approval_id=chosen_approval_id,
            approved_review_sha256=reviewed,
            actor=actor,
            resolution_id=resolution_id,
        )
        resolution = bundle.resolution
    else:
        bundle = None
        resolution = build_workflow_review_resolution(
            resolution_id=resolution_id,
            review_sha256=reviewed,
            decision=normalized,
            actor=actor,
        )
    log_path = (
        Path(decision_log).expanduser().resolve()
        if decision_log is not None
        else review_path.with_name(review_path.name + ".decisions.jsonl")
    )
    RuntimeEventStore(
        log_path,
        session_id="review-" + review.review_sha256[:16],
    ).record_workflow_review_resolution(
        turn_id="human-resolution",
        resolution=resolution,
    )
    if bundle is not None:
        write_workflow_execution_approval_bundle(bundle, output_file)
    return resolution, bundle


def inspect_workflow_execution_replay(
    *,
    review_file: str | Path,
    workspace: str | Path,
    task_spec_sha256: str = "",
) -> dict[str, Any]:
    """Report what re-running a stored review would put in front of a human.

    A review is inert and re-loadable, and neither its digest nor an approval
    bundle's contains a run directory, a session id, or a clock, so the same
    stored workflow re-presents identically on an undrifted host.  What was
    missing was a way to *see* that before deciding, and a decision identity
    distinct from the one already spent.

    Drift is reported rather than refused, because the human decides over the
    displayed scientific state.  Approved bytes that are no longer anywhere
    under the workspace are reported as ``missing_approved_artifacts``: they
    make execution certain to fail while resolving inputs, before anything is
    dispatched, so a caller must not offer a decision on them -- but deciding
    that is the command's job, not this function's.

    This reports; it never approves, consumes, or launches.  The argv gate that
    runs immediately before dispatch is untouched, so nothing decided here can
    let drift through to an engine.
    """

    root = _validated_workspace(workspace)
    review = load_workflow_execution_review(
        Path(review_file).expanduser().resolve()
    )
    if Path(review.request.workspace).resolve() != root:
        raise ContractError("execution review targets another workspace")
    expected_task = str(task_spec_sha256).strip()
    if (
        expected_task
        and review.scientific_plan.task_spec_sha256 != expected_task
    ):
        raise ContractError("execution review targets another task")

    # Index the workspace once.  Resolving each binding through a fresh search
    # rehashes the same files for every digest, and replay exists precisely to
    # be run again over a workspace that keeps accumulating outputs -- so the
    # naive form gets slower exactly as the feature gets used.
    wanted = {
        digest
        for binding in review.request.node_bindings
        for digest in (
            binding.initial_artifact_sha256,
            binding.project_artifact_sha256,
        )
        if digest
    }
    found: set[str] = set()
    for candidate in sorted(root.rglob("*")):
        if not candidate.is_file() or candidate.is_symlink():
            continue
        try:
            digest = file_sha256(candidate)
        except OSError:
            continue
        if digest in wanted:
            found.add(digest)
            if found == wanted:
                break

    missing: list[str] = []
    present: list[str] = []
    for binding in review.request.node_bindings:
        for label, digest in (
            ("input", binding.initial_artifact_sha256),
            ("project", binding.project_artifact_sha256),
        ):
            if not digest:
                continue
            record = f"{binding.node_id}:{label}:{digest[:16]}"
            if digest not in found:
                missing.append(record)
            else:
                present.append(record)

    spent = _spent_approval_ids(root, review.review_sha256)
    return {
        "review_sha256": review.review_sha256,
        "workflow_id": review.request.workflow_id,
        "workspace": str(root),
        "task_spec_sha256": review.scientific_plan.task_spec_sha256,
        "node_count": len(review.request.node_bindings),
        "non_executable_node_ids": list(review.non_executable_node_ids),
        "approved_artifacts_present": sorted(set(present)),
        "missing_approved_artifacts": sorted(set(missing)),
        "previously_consumed_approval_ids": spent,
        "canonical_review": canonical_json(
            {"workflow_execution_review": canonical_data(review)}
        ),
    }


def spent_workflow_approval_ids(
    workspace: Path, review_sha256: str
) -> list[str]:
    """Public name for the consumed-decision scan the replay path uses."""

    return _spent_approval_ids(workspace, review_sha256)


def _spent_approval_ids(workspace: Path, review_sha256: str) -> list[str]:
    """Which approval ids for this review have already been burned here.

    The ledger is per approval id, so a second genuine decision needs a
    different one.  Reporting the spent ids is what turns "already consumed"
    from a dead end into an instruction.
    """

    ledger_root = workspace / _PRIVATE_ROOT_NAME / "approvals"
    if not ledger_root.is_dir():
        return []
    spent: list[str] = []
    for entry in sorted(ledger_root.iterdir()):
        events = entry / "consumption-events.jsonl"
        if not entry.is_dir() or not events.is_file():
            continue
        try:
            text = events.read_text(encoding="utf-8")
        except OSError:
            continue
        if review_sha256 in text or entry.name.endswith(review_sha256[:16]):
            spent.append(entry.name)
    return spent


def replay_approval_id() -> str:
    """A decision identity that cannot collide with one already spent.

    The default approval id is derived from the review digest alone, so a
    second honest decision on an identical review reproduces the first one's
    id and its resolution event, and the append is refused as a conflicting
    idempotency key before any bundle is written.  A replay therefore has to
    carry its own identity, the way the interactive path already mints one per
    execution.
    """

    return "replay-" + uuid.uuid4().hex[:16]


def claim_workflow_execution_approval_bundle(
    bundle: WorkflowExecutionApprovalBundleV1,
    *,
    workspace: str | Path,
):
    """Durably consume a bundle once in its exact approved workspace."""

    root = _validated_workspace(workspace)
    if Path(bundle.workflow_approval.workspace).resolve() != root:
        raise ContractError("execution bundle targets another workspace")
    approval_id = bundle.workflow_approval.approval_id
    ledger_root = root / _PRIVATE_ROOT_NAME / "approvals" / approval_id
    if ledger_root.is_symlink():
        raise ContractError("approval ledger cannot be a symbolic link")
    ledger_root.mkdir(parents=True, exist_ok=True, mode=0o700)
    ledger_root.chmod(0o700)
    store = RuntimeEventStore(
        ledger_root / "consumption-events.jsonl",
        session_id="approval-" + approval_id,
    )
    return store.claim_execution_bundle(
        turn_id="execute-claim",
        bundle=bundle,
    )


def _parse_execution_resources(value: Any) -> ExecutionResourceSpecV1:
    """Rehydrate the user-approved resource allocation for execution."""

    if not isinstance(value, Mapping):
        raise ContractError("execution resources must be an object")
    try:
        return ExecutionResourceSpecV1(**dict(value))
    except (KeyError, TypeError, ValueError) as exc:
        raise ContractError(
            "execution resources do not match the v1 schema"
        ) from exc


def _write_execution_server_profile(
    run_directory: Path,
    resources: ExecutionResourceSpecV1,
    *,
    scratch_root: Path | None = None,
) -> Path:
    """Write the local CPU profile from the user-approved allocation."""

    profile = run_directory / "execution-server.yaml"
    num_hours = max(1, math.ceil(resources.node_timeout_seconds / 3600))
    text = (
        "SERVER:\n"
        "  SCHEDULER: null\n"
        f"  NUM_HOURS: {num_hours}\n"
        f"  MEM_GB: {resources.memory_gb:g}\n"
        f"  NUM_CORES: {resources.cores}\n"
        f"  NUM_GPUS: {resources.gpu_count}\n"
        f"  NUM_THREADS: {resources.cores}\n"
        + (
            f"  SCRATCH_DIR: {str(scratch_root)!r}\n"
            if scratch_root is not None
            else "  SCRATCH_DIR: null\n"
        )
        + _local_program_server_blocks(
            preserve_active=scratch_root is not None,
            scratch_root=scratch_root,
        )
    )
    _write_private_exact(profile, text.encode("utf-8"))
    return profile


def _local_result(
    *,
    session_id: str,
    task_spec_sha256: str,
    terminal_state: str,
    execution_requested: bool,
    execution_profile_status: str,
    final_text: str,
) -> LiveAgentSessionResultV1:
    body = {
        "schema_version": "chemsmart.live-agent-session-result.v1",
        "session_id": session_id,
        "task_spec_sha256": task_spec_sha256,
        "terminal_state": terminal_state,
        "execution_requested": bool(execution_requested),
        "execution_profile_status": execution_profile_status,
        "final_text": final_text,
        "artifact_records": (),
        "conformance_records": (),
        "public_transcript": (),
        "successful_tool_calls": 0,
        "failed_tool_calls": 0,
        "execution_review": {},
        "event_stream_head_sha256": "",
    }
    return LiveAgentSessionResultV1(
        **body, result_sha256=canonical_sha256(body)
    )


__all__ = [
    "LiveAgentSessionResultV1",
    "load_workflow_execution_approval_bundle",
    "load_workflow_execution_review",
    "claim_workflow_execution_approval_bundle",
    "run_live_agent_session",
    "resolve_workflow_execution_review",
    "write_workflow_execution_review",
    "write_workflow_execution_approval_bundle",
]
