"""Plan-first live ChemSmart Agent sessions.

This module is the small composition root behind ``chemsmart agent plan`` and
``chemsmart agent run``.  It binds a provider session to exact, pre-existing
XYZ artifacts, the live Click schema, observed program conformance, and a
private Runtime V2 event stream.  It never generates coordinates, engine
inputs, shell commands, or scientific readiness decisions.

The execution profile is deliberately progressive.  Until the tool host has
an approval-bound execution composition API, ``agent run`` uses the complete
planning/preview path and reports ``waiting_for_approval`` rather than
pretending that an engine was invoked.
"""

from __future__ import annotations

import inspect
import json
import math
import os
import shutil
import sys
import uuid
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping

from chemsmart.agent._contracts import (
    ContractError,
    TrustedArtifactRefV1,
    canonical_data,
    canonical_json,
    canonical_sha256,
    file_sha256,
    require_sha256,
)
from chemsmart.agent.adaptive_api_campaign import (
    AdaptiveHypothesisV1,
    AdaptiveNetworkBudgetV1,
)
from chemsmart.agent.analysis_completion import load_analysis_completion_policy
from chemsmart.agent.api_access import load_secret_lease
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
from chemsmart.agent.cli_schema import (
    LiveClickSchemaV1,
    build_live_click_schema,
)
from chemsmart.agent.dependency_context import (
    ContextSelectionReceiptV1,
    TaskDependencyContextV1,
    build_dependency_context_public_projection,
)
from chemsmart.agent.execution import (
    ApprovedNodeBindingV1,
    ExecutionResourceSpecV1,
    FrozenMaterializedNodePreviewV1,
    FrozenProducerEdgeRuleV1,
    FrozenWorkflowApprovalV1,
    ProducerEdgeRuleV1,
    WorkflowExecutionApprovalV1,
    build_execution_resource_spec,
)
from chemsmart.agent.experiments.host_oracle import (
    HostOracleInputBundleV1,
    build_host_oracle_input_bundle,
)
from chemsmart.agent.experiments.qwen_pyscf_dfc import (
    QwenDfcArmV1,
    QwenExperimentPreparationV1,
    QwenPyscfCaseSpecV1,
    build_qwen_experiment_preparation,
)
from chemsmart.agent.identity import (
    ApprovedMolecularIdentityV1,
    validate_identity_for_geometry,
)
from chemsmart.agent.knowledge_packs import (
    BUILTIN_PROGRAM_PACKS,
    activate_program_knowledge,
    skills_for_activation,
)
from chemsmart.agent.live_specialists import (
    LiveSpecialistCampaignV1,
    build_experiment_seed_plan,
    build_f_invariant_critic_candidate,
)
from chemsmart.agent.projects import project_document, render_project_yaml
from chemsmart.agent.provider_config import (
    ALIBABA_TOKEN_PLAN_MODEL,
    ALIBABA_TOKEN_PLAN_PROVIDER,
    AgentProviderProfileV1,
    AgentProviderSelectionV1,
    load_agent_provider_selection,
)
from chemsmart.agent.runtime.contracts import (
    ResourceBudgetV1,
    TaskEnvelopeV1,
    TaskPhase,
)
from chemsmart.agent.runtime.event_store import RuntimeEventStore
from chemsmart.agent.services.unified_session import UnifiedSessionRunner
from chemsmart.agent.skills import (
    SkillDocumentV1,
    resolve_skills,
    skills_enabled,
)
from chemsmart.agent.specialists import READ_ONLY_CRITIC
from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1
from chemsmart.agent.tool_specs import (
    AgentToolSurfaceV1,
    build_approved_execution_tool_surface,
    build_command_compiled_tool_surface,
)
from chemsmart.agent.workflows import (
    HarnessExperimentConfigV1,
    MaterializedNodeV1,
    MaterializedWorkflowV1,
    ScientificWorkflowEdgeV2,
    ScientificWorkflowNodeV2,
    ScientificWorkflowPlanV2,
)
from chemsmart.analysis.result_quantities import (
    QuantityExtractionError,
    validate_pyscf_analysis_artifact,
)

_SESSION_WALL_TIME_SECONDS = 90 * 60
_CHEMISTRY_NODE_TIMEOUT_SECONDS = 10 * 60
_MAX_TOOL_CALLS = 256
_PYSCF_INTERPRETER = (
    Path(os.environ.get("CHEMSMART_PYSCF_INTERPRETER", sys.executable))
    .expanduser()
    .resolve()
)
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
    experiment_observations: dict[str, Any]
    successful_tool_calls: int
    failed_tool_calls: int
    event_stream_head_sha256: str
    result_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.live-agent-session-result.v1":
            raise ContractError("unsupported live session result schema")
        if self.terminal_state not in {
            "complete",
            "planned",
            "failed",
            "blocked",
            "waiting_for_approval",
        }:
            raise ContractError("invalid live session terminal state")
        _validate_path_free_experiment_observations(
            self.experiment_observations
        )
        body = self._body()
        if self.result_sha256 != canonical_sha256(body):
            raise ContractError("live session result digest mismatch")

    def _body(self) -> dict[str, Any]:
        return {
            key: value
            for key, value in asdict(self).items()
            if key != "result_sha256"
        }

    def public_summary_json(self) -> str:
        """Return the exact visible result; no local path is included."""

        return canonical_json(
            {**self._body(), "result_sha256": self.result_sha256}
        )


@dataclass(frozen=True)
class CampaignPreparationHostSnapshotV1:
    """Provider-free, reusable host observations for one campaign artifact.

    Runtime objects are held explicitly by the campaign controller rather than
    hidden in a module cache.  The public digest covers their stable identities
    and every path-free conformance/environment observation.  Per-episode
    workspaces may use the snapshot only when their exact artifact, approved
    identity, provider selection, live schema, registry, and tool surface agree.
    """

    schema_version: str
    provider_profile_sha256: str
    provider_selection_sha256: str
    artifact_sha256s: tuple[str, ...]
    approved_identity_sha256: str
    registry_sha256: str
    live_cli_schema_sha256: str
    tool_schema_sha256: str
    component_conformance_receipts: tuple[
        ProgramComponentConformanceReceiptV1, ...
    ]
    environment_targets: tuple[EnvironmentTargetV1, ...]
    compute_environment_receipts: tuple[
        TrustedComputeEnvironmentReceiptV1, ...
    ]
    conformance_records: tuple[dict[str, Any], ...]
    provider_calls: int
    engine_calls: int
    approval_files: int
    snapshot_sha256: str
    registry: ProgramCapabilityRegistryV1 = field(repr=False, compare=False)
    live_schema: LiveClickSchemaV1 = field(repr=False, compare=False)
    tool_surface: AgentToolSurfaceV1 = field(repr=False, compare=False)

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.campaign-host-snapshot.v1":
            raise ContractError("unsupported campaign host snapshot")
        for name in (
            "provider_profile_sha256",
            "provider_selection_sha256",
            "registry_sha256",
            "live_cli_schema_sha256",
            "tool_schema_sha256",
        ):
            require_sha256(getattr(self, name), name)
        if self.approved_identity_sha256:
            require_sha256(
                self.approved_identity_sha256,
                "approved_identity_sha256",
            )
        if (
            self.artifact_sha256s != tuple(sorted(set(self.artifact_sha256s)))
            or not self.artifact_sha256s
        ):
            raise ContractError("campaign host artifacts must be canonical")
        for digest in self.artifact_sha256s:
            require_sha256(digest, "artifact_sha256")
        if self.provider_calls or self.engine_calls or self.approval_files:
            raise ContractError(
                "campaign host snapshot must remain provider-free"
            )
        if self.registry.registry_sha256 != self.registry_sha256:
            raise ContractError("campaign host registry object mismatch")
        if self.live_schema.schema_sha256 != self.live_cli_schema_sha256:
            raise ContractError("campaign host live schema object mismatch")
        if self.tool_surface.tool_schema_sha256 != self.tool_schema_sha256:
            raise ContractError("campaign host tool surface mismatch")
        for receipt in self.component_conformance_receipts:
            if receipt.registry_sha256 != self.registry_sha256:
                raise ContractError("campaign conformance registry mismatch")
            if receipt.live_cli_schema_sha256 != self.live_cli_schema_sha256:
                raise ContractError("campaign conformance schema mismatch")
        ordered_records = tuple(
            sorted(self.conformance_records, key=_record_sort_key)
        )
        if self.conformance_records != ordered_records:
            raise ContractError("campaign host observations are not canonical")
        _validate_path_free_experiment_observations(
            {"host_observations": self.conformance_records}
        )
        if self.snapshot_sha256 != canonical_sha256(self._body()):
            raise ContractError("campaign host snapshot digest mismatch")

    def _body(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "provider_profile_sha256": self.provider_profile_sha256,
            "provider_selection_sha256": self.provider_selection_sha256,
            "artifact_sha256s": self.artifact_sha256s,
            "approved_identity_sha256": self.approved_identity_sha256,
            "registry_sha256": self.registry_sha256,
            "live_cli_schema_sha256": self.live_cli_schema_sha256,
            "tool_schema_sha256": self.tool_schema_sha256,
            "component_conformance_receipt_sha256s": tuple(
                item.receipt_sha256
                for item in self.component_conformance_receipts
            ),
            "environment_target_sha256s": tuple(
                canonical_sha256(item) for item in self.environment_targets
            ),
            "compute_environment_identity_sha256s": tuple(
                item.evidence_sha256
                for item in self.compute_environment_receipts
            ),
            "conformance_records": self.conformance_records,
            "provider_calls": self.provider_calls,
            "engine_calls": self.engine_calls,
            "approval_files": self.approval_files,
        }

    def public_record(self) -> dict[str, Any]:
        """Return the path-free, provider-free campaign observation."""

        return {**self._body(), "snapshot_sha256": self.snapshot_sha256}


def build_campaign_preparation_host_snapshot(
    *,
    provider: str | None,
    provider_config_file: str | Path | None,
    workspace: str | Path,
    approved_molecular_identity: ApprovedMolecularIdentityV1 | None = None,
) -> CampaignPreparationHostSnapshotV1:
    """Observe schema, fake preview, and environments exactly once.

    This function never loads a credential or creates a provider session.  Its
    only executable work is ChemSmart's non-running conformance/fake-preview
    path and registered environment observation.
    """

    selection = load_agent_provider_selection(
        provider_config_file,
        requested_profile=str(provider).strip() if provider else None,
    )
    profile = selection.active_profile
    _validate_campaign_provider_selection(selection)
    workspace_path = _validated_workspace(workspace)
    observations = _scan_xyz_artifacts(workspace_path)
    if len(observations) != 1:
        raise ContractError(
            "campaign host snapshot requires exactly one approved XYZ"
        )
    _validated_identity_records(observations, approved_molecular_identity)
    probe_directory = _private_preparation_directory(
        workspace=workspace_path,
        episode_id=(
            "campaign-host-snapshot."
            f"{profile.profile_sha256[:16]}."
            f"{observations[0].artifact.sha256[:16]}"
        ),
    )
    registry = load_program_capabilities()
    live_schema = build_live_click_schema()
    conformance, conformance_records = _bootstrap_conformance(
        run_directory=probe_directory,
        input_artifact=observations[0].artifact,
        registry_sha256=registry.registry_sha256,
        live_schema=live_schema,
    )
    environment_targets, compute_receipts, environment_records = (
        _observe_environments()
    )
    records = tuple(
        sorted(
            (*conformance_records, *environment_records), key=_record_sort_key
        )
    )
    surface = build_command_compiled_tool_surface(registry)
    body = {
        "schema_version": "chemsmart.campaign-host-snapshot.v1",
        "provider_profile_sha256": profile.profile_sha256,
        "provider_selection_sha256": selection.selection_sha256,
        "artifact_sha256s": tuple(
            sorted(item.artifact.sha256 for item in observations)
        ),
        "approved_identity_sha256": (
            approved_molecular_identity.identity_sha256
            if approved_molecular_identity is not None
            else ""
        ),
        "registry_sha256": registry.registry_sha256,
        "live_cli_schema_sha256": live_schema.schema_sha256,
        "tool_schema_sha256": surface.tool_schema_sha256,
        "component_conformance_receipt_sha256s": tuple(
            item.receipt_sha256 for item in conformance
        ),
        "environment_target_sha256s": tuple(
            canonical_sha256(item) for item in environment_targets
        ),
        "compute_environment_identity_sha256s": tuple(
            item.evidence_sha256 for item in compute_receipts
        ),
        "conformance_records": records,
        "provider_calls": 0,
        "engine_calls": 0,
        "approval_files": 0,
    }
    return CampaignPreparationHostSnapshotV1(
        schema_version=body["schema_version"],
        provider_profile_sha256=body["provider_profile_sha256"],
        provider_selection_sha256=body["provider_selection_sha256"],
        artifact_sha256s=body["artifact_sha256s"],
        approved_identity_sha256=body["approved_identity_sha256"],
        registry_sha256=body["registry_sha256"],
        live_cli_schema_sha256=body["live_cli_schema_sha256"],
        tool_schema_sha256=body["tool_schema_sha256"],
        component_conformance_receipts=conformance,
        environment_targets=environment_targets,
        compute_environment_receipts=compute_receipts,
        conformance_records=records,
        provider_calls=0,
        engine_calls=0,
        approval_files=0,
        snapshot_sha256=canonical_sha256(body),
        registry=registry,
        live_schema=live_schema,
        tool_surface=surface,
    )


def probe_live_experiment_preparation(
    *,
    task: str,
    provider: str | None,
    provider_config_file: str | Path | None,
    workspace: str | Path,
    experiment_arm: QwenDfcArmV1,
    experiment_case: QwenPyscfCaseSpecV1,
    experiment_repeat_index: int,
    approved_molecular_identity: ApprovedMolecularIdentityV1 | None = None,
    campaign_preparation_snapshot: (
        CampaignPreparationHostSnapshotV1 | None
    ) = None,
) -> QwenExperimentPreparationV1:
    """Observe the exact coordinator base contract without provider access.

    Without an explicit campaign snapshot this performs one local snapshot
    build for backward-compatible standalone use.  Campaign callers pass the
    same snapshot to every case/arm/repeat, avoiding repeated fake previews and
    environment probes.  No credential, provider session, approval, or engine
    launch is involved.
    """

    selection = load_agent_provider_selection(
        provider_config_file,
        requested_profile=str(provider).strip() if provider else None,
    )
    profile = selection.active_profile
    _validate_experiment_request(
        task=task,
        selection=selection,
        arm=experiment_arm,
        case=experiment_case,
        repeat_index=experiment_repeat_index,
        execution_enabled=False,
        approval_file=None,
    )
    workspace_path = _validated_workspace(workspace)
    observations = _scan_xyz_artifacts(workspace_path)
    result_observations = _scan_pyscf_result_artifacts(workspace_path)
    if not observations:
        raise ContractError("experiment preparation requires an exact XYZ")
    identity_records = _validated_identity_records(
        observations, approved_molecular_identity
    )
    task_spec_sha256 = _task_spec_sha256(
        task,
        observations,
        approved_molecular_identity,
        result_observations=result_observations,
    )
    snapshot = campaign_preparation_snapshot
    if snapshot is None:
        snapshot = build_campaign_preparation_host_snapshot(
            provider=provider,
            provider_config_file=provider_config_file,
            workspace=workspace_path,
            approved_molecular_identity=approved_molecular_identity,
        )
    _validate_campaign_snapshot_reuse(
        snapshot=snapshot,
        selection=selection,
        observations=observations,
        approved_molecular_identity=approved_molecular_identity,
    )
    context = _public_context(
        task=task,
        task_spec_sha256=task_spec_sha256,
        observations=observations,
        result_observations=result_observations,
        conformance_records=snapshot.conformance_records,
        registry_sha256=snapshot.registry_sha256,
        live_schema_sha256=snapshot.live_cli_schema_sha256,
        execution_requested=False,
        execution_available=False,
        provider_record=_provider_public_record(
            profile=profile,
            fallback_profiles=selection.fallback_profiles,
            experiment=True,
        ),
        experiment_record=_experiment_public_record(
            arm=experiment_arm,
            case=experiment_case,
            repeat_index=experiment_repeat_index,
        ),
        approved_identity_records=identity_records,
    )
    base_messages = _coordinator_base_messages(
        context=context,
        approved_workflow={},
        experiment_arm=experiment_arm,
        provider_profile=profile,
        task=task,
    )
    return _observed_experiment_preparation(
        arm=experiment_arm,
        case=experiment_case,
        repeat_index=experiment_repeat_index,
        task_spec_sha256=task_spec_sha256,
        artifact_sha256s=tuple(item.artifact.sha256 for item in observations),
        provider_profile=profile,
        base_messages=base_messages,
        tool_schema_sha256=snapshot.tool_schema_sha256,
        host_snapshot_sha256=snapshot.snapshot_sha256,
    )


def run_live_agent_session(
    *,
    task: str,
    provider: str | None,
    provider_config_file: str | Path | None = None,
    secret_file: str | Path,
    workspace: str | Path,
    execution_enabled: bool,
    approval_file: str | Path | None,
    analysis_completion_file: str | Path | None = None,
    experiment_arm: QwenDfcArmV1 | None = None,
    experiment_case: QwenPyscfCaseSpecV1 | None = None,
    experiment_repeat_index: int = 0,
    experiment_config: HarnessExperimentConfigV1 | None = None,
    approved_molecular_identity: ApprovedMolecularIdentityV1 | None = None,
    approved_molecular_identities: Iterable[ApprovedMolecularIdentityV1] = (),
    campaign_preparation_snapshot: (
        CampaignPreparationHostSnapshotV1 | None
    ) = None,
    dependency_context: TaskDependencyContextV1 | None = None,
    dependency_context_selection_receipt: (
        ContextSelectionReceiptV1 | None
    ) = None,
    dependency_public_records: Mapping[str, Mapping[str, Any]] | None = None,
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
    if (experiment_arm is None) != (experiment_case is None):
        raise ContractError(
            "experiment configuration and case must be supplied together"
        )
    if experiment_repeat_index < 0:
        raise ContractError("experiment repeat index must be non-negative")
    if experiment_config is not None and experiment_arm is None:
        raise ContractError(
            "frozen experiment config requires an experiment arm and case"
        )
    if campaign_preparation_snapshot is not None and experiment_arm is None:
        raise ContractError(
            "campaign host snapshot requires an experiment arm and case"
        )
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
    if experiment_arm is not None and experiment_case is not None:
        _validate_experiment_request(
            task=task,
            selection=selection,
            arm=experiment_arm,
            case=experiment_case,
            repeat_index=experiment_repeat_index,
            execution_enabled=execution_enabled,
            approval_file=approval_file,
        )
    workspace_path = _validated_workspace(workspace)
    observations = _scan_xyz_artifacts(workspace_path)
    result_observations = _scan_pyscf_result_artifacts(workspace_path)
    identities = _coerce_approved_identities(
        approved_molecular_identity, approved_molecular_identities
    )
    identity_records = _validated_identity_records(observations, identities)
    task_spec_sha256 = _task_spec_sha256(
        task,
        observations,
        identities,
        result_observations=result_observations,
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
        not observations
        and result_observations
        and analysis_completion_policy is not None
        and not execution_enabled
        and experiment_arm is None
        and campaign_preparation_snapshot is None
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
                "No exact XYZ artifact is present in the approved workspace. "
                "Add a user-approved coordinate file; coordinates were not generated."
            ),
        )

    if campaign_preparation_snapshot is not None:
        if len(identities) > 1:
            raise ContractError(
                "campaign host snapshots currently bind one molecular identity; "
                "run multi-structure paper sessions without a reused snapshot"
            )
        _validate_campaign_snapshot_reuse(
            snapshot=campaign_preparation_snapshot,
            selection=selection,
            observations=observations,
            approved_molecular_identity=(
                identities[0] if identities else None
            ),
        )
        registry = campaign_preparation_snapshot.registry
        live_schema = campaign_preparation_snapshot.live_schema
        conformance = (
            campaign_preparation_snapshot.component_conformance_receipts
        )
        environment_targets = campaign_preparation_snapshot.environment_targets
        compute_receipts = (
            campaign_preparation_snapshot.compute_environment_receipts
        )
        conformance_records = campaign_preparation_snapshot.conformance_records
    elif analysis_only_session:
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

    execution_profile_ready = _execution_composition_available()
    use_execution_surface = bool(
        execution_enabled and approval_file and execution_profile_ready
    )
    approved_project_records: tuple[dict[str, Any], ...] = ()
    approved_workflow_record: dict[str, Any] = {}
    surface = (
        build_approved_execution_tool_surface(registry)
        if use_execution_surface
        else (
            campaign_preparation_snapshot.tool_surface
            if campaign_preparation_snapshot is not None
            else build_command_compiled_tool_surface(registry)
        )
    )

    event_store = RuntimeEventStore(
        run_directory / "events.jsonl", session_id=session_id
    )
    preview_server_path = _ensure_preview_server(run_directory)
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
        # Preview candidates stay inside this session's private workspace.
        # The execution composition below replaces this with the exact
        # user-approved workflow workspace.
        "approved_workspace": run_directory,
        "preview_server": str(preview_server_path),
        "result_functional_evidence": {
            item.validation_receipt_sha256: item.public_record()
            for item in result_observations
        },
        "analysis_completion_policy": analysis_completion_policy,
        "dependency_context": dependency_context,
        "dependency_context_selection_receipt": (
            dependency_context_selection_receipt
        ),
    }
    if use_execution_surface:
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
            # ``_execution_scientific_plan`` resolves the DAG to execute by
            # looking up the frozen approval's ``plan_sha256`` in the host's
            # plan registry, which only holds plans this session produced. The
            # approved plan came from an earlier session, so the lookup failed
            # with "frozen workflow approval has no registered scientific
            # plan" -- observed twelve times in one run. Registering it here
            # means execution can resolve the approved DAG directly instead of
            # waiting for the model to re-derive it byte for byte.
            host_kwargs["scientific_workflow_plan"] = approved_plan
        approved_materialization = execution_inputs.pop(
            "approved_materialized_workflow"
        )
        if approved_materialization is not None:
            # Same story one layer down: execution resolves the frozen
            # approval's ``materialized_workflow_sha256`` from a registry that
            # only this session fills, and refuses with "frozen workflow
            # approval has no canonical materialization" when it cannot.
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
        approved_project_records=approved_project_records,
        approved_workflow_record=approved_workflow_record,
        provider_record=_provider_public_record(
            profile=profile,
            fallback_profiles=selection.fallback_profiles,
            experiment=experiment_arm is not None,
        ),
        experiment_record=(
            _experiment_public_record(
                arm=experiment_arm,
                case=experiment_case,
                repeat_index=experiment_repeat_index,
            )
            if experiment_arm is not None and experiment_case is not None
            else None
        ),
        approved_identity_records=identity_records,
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
        experiment_arm=experiment_arm,
        provider_profile=profile,
        task=task,
    )
    observed_preparation = (
        _observed_experiment_preparation(
            arm=experiment_arm,
            case=experiment_case,
            repeat_index=experiment_repeat_index,
            task_spec_sha256=task_spec_sha256,
            artifact_sha256s=tuple(
                item.artifact.sha256 for item in observations
            ),
            provider_profile=profile,
            base_messages=base_messages,
            tool_schema_sha256=surface.tool_schema_sha256,
            host_snapshot_sha256=(
                campaign_preparation_snapshot.snapshot_sha256
                if campaign_preparation_snapshot is not None
                else ""
            ),
        )
        if experiment_arm is not None and experiment_case is not None
        else None
    )
    if experiment_config is not None and observed_preparation is not None:
        _require_exact_experiment_config(
            frozen=experiment_config,
            observed=observed_preparation.experiment_config,
        )
    experiment_config = (
        observed_preparation.experiment_config
        if observed_preparation is not None
        else None
    )
    specialist_campaign: LiveSpecialistCampaignV1 | None = None
    if (
        experiment_arm is not None
        and experiment_case is not None
        and experiment_config is not None
    ):
        seed_plan = build_experiment_seed_plan(
            case=experiment_case,
            task_spec_sha256=task_spec_sha256,
            artifact_sha256s=tuple(
                item.artifact.sha256 for item in observations
            ),
        )
        specialist_campaign = LiveSpecialistCampaignV1.start(
            arm=experiment_arm,
            experiment_config=experiment_config,
            plan=seed_plan,
            coordinator_session_id=session_id,
            public_context=context,
            source_sha256s=experiment_case.source_sha256s,
            artifact_sha256s=tuple(
                item.artifact.sha256 for item in observations
            ),
            base_tool_surface=surface,
            provider_profile=profile,
            secret_file=secret_file,
            run_directory=run_directory,
            host_builder=_live_specialist_host_builder(
                base_host_kwargs=host_kwargs,
                coordinator_host=host,
                coordinator_surface=surface,
            ),
        )
        context = {
            **context,
            "specialist_advisory": (
                specialist_campaign.coordinator_advisory_record()
            ),
        }
    messages = [
        base_messages[0],
        {"role": "user", "content": canonical_json(context)},
    ]
    network_budget = _network_budget(
        profile,
        max_concurrency=(
            experiment_arm.max_concurrency if experiment_arm is not None else 1
        ),
    )
    hypothesis = _hypothesis(
        session_id=session_id,
        messages=messages,
        task_spec_sha256=task_spec_sha256,
        artifact_sha256s=tuple(
            item.artifact.sha256
            for item in (*observations, *result_observations)
        ),
        tool_schema_sha256=surface.tool_schema_sha256,
        network_budget=network_budget,
        execution_requested=execution_enabled,
        experiment_config=experiment_config,
        experiment_case=experiment_case,
        experiment_repeat_index=experiment_repeat_index,
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
    )
    lease = load_secret_lease(
        provider=normalized_provider,
        path=secret_file,
        # The profile already says which key it bills; selecting a different
        # one here is how a second key for the same provider gets charged by
        # accident.
        label=profile.api_key_env,
        ttl_seconds=_SESSION_WALL_TIME_SECONDS + 60,
    )
    loop_result = UnifiedSessionRunner(
        host=host,
        event_store=event_store,
        credential_lease=lease,
        provider_config=profile.runtime_config(),
    ).run(
        messages=messages,
        envelope=envelope,
        hypothesis=hypothesis,
        network_budget=network_budget,
        feedback_projection=(
            experiment_config.feedback_projection
            if experiment_config is not None
            else "full-v1"
        ),
    )

    experiment_observations: dict[str, Any] = {}
    if specialist_campaign is not None and experiment_case is not None:
        # The grading oracle is reconstructed from the complete canonical host
        # results in Runtime V2, never from the provider-visible F projection.
        # This must happen after the coordinator loop has terminated and before
        # its evidence is combined with detached critic observations.
        coordinator_events = event_store.read_events()
        host_oracle_bundle = build_host_oracle_input_bundle(
            events=coordinator_events,
            session_id=session_id,
            event_stream_head_sha256=loop_result.event_stream_head_sha256,
            successful_tool_calls=loop_result.successful_tool_calls,
            failed_tool_calls=loop_result.failed_tool_calls,
        )
        critic_candidate = build_f_invariant_critic_candidate(
            candidate_id=f"candidate.{session_id}",
            task_spec_sha256=task_spec_sha256,
            host_oracle_input_bundle=host_oracle_bundle,
            coordinator_public_decisions=tuple(
                sorted(
                    host.scientific_decisions.values(),
                    key=lambda item: (
                        item.decision_id,
                        item.record_sha256,
                    ),
                )
            ),
        )
        specialist_campaign.run_critic(
            coordinator_session_id=session_id,
            candidate=critic_candidate,
            public_context={
                "task": task,
                "case_id": experiment_case.case_id,
                "authority": (
                    "Review only; do not repair, approve, execute, or set "
                    "scientific readiness or terminal state."
                ),
            },
            source_sha256s=experiment_case.source_sha256s,
            artifact_sha256s=tuple(
                item.artifact.sha256 for item in observations
            ),
        )
        raw_experiment_observations = (
            specialist_campaign.public_observation_record(
                coordinator_usage=_coordinator_usage_record(
                    event_store=event_store,
                    successful_tool_calls=loop_result.successful_tool_calls,
                    failed_tool_calls=loop_result.failed_tool_calls,
                )
            )
        )
        raw_experiment_observations = _bind_feedback_receipts(
            observations=raw_experiment_observations,
            events=coordinator_events,
        )
        raw_experiment_observations = _bind_host_oracle_bundle(
            observations=raw_experiment_observations,
            bundle=host_oracle_bundle,
        )
        if observed_preparation is None:
            raise ContractError("experiment preparation observation is absent")
        experiment_observations = _bind_preparation_observation(
            observations=raw_experiment_observations,
            preparation=observed_preparation,
        )

    execution_status = (
        "approved_profile_active"
        if use_execution_surface
        else (
            "not_requested"
            if not execution_enabled
            else "planning_complete_execution_profile_unavailable"
        )
    )
    terminal_state = loop_result.terminal_state
    final_text = loop_result.final_text
    if execution_enabled and not use_execution_surface:
        terminal_state = "waiting_for_approval"
        suffix = (
            " Planning and preview evidence were retained, but the current "
            "host does not expose a complete approval-bound execution composition. "
            "No chemistry engine was launched."
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
        "experiment_observations": experiment_observations,
        "successful_tool_calls": loop_result.successful_tool_calls,
        "failed_tool_calls": loop_result.failed_tool_calls,
        "event_stream_head_sha256": loop_result.event_stream_head_sha256,
    }
    return LiveAgentSessionResultV1(
        **body, result_sha256=canonical_sha256(body)
    )


def _validate_campaign_provider_selection(
    selection: AgentProviderSelectionV1,
) -> None:
    profile = selection.active_profile
    if selection.fallback_profiles:
        raise ContractError(
            "campaign host snapshot forbids provider fallbacks"
        )
    runtime_config = profile.runtime_config()
    if (
        runtime_config.provider != profile.provider
        or runtime_config.model != profile.model
        or runtime_config.reasoning_effort != profile.reasoning_effort
    ):
        raise ContractError(
            "campaign provider profile differs from its runtime adapter"
        )


def _validate_campaign_snapshot_reuse(
    *,
    snapshot: CampaignPreparationHostSnapshotV1,
    selection: AgentProviderSelectionV1,
    observations: tuple[_XyzObservation, ...],
    approved_molecular_identity: ApprovedMolecularIdentityV1 | None,
) -> None:
    """Reject cross-artifact, provider, identity, or schema reuse."""

    if not isinstance(snapshot, CampaignPreparationHostSnapshotV1):
        raise ContractError("campaign host snapshot is not typed")
    _validate_campaign_snapshot_integrity(snapshot)
    _validate_campaign_provider_selection(selection)
    if snapshot.provider_profile_sha256 != (
        selection.active_profile.profile_sha256
    ):
        raise ContractError("campaign host snapshot provider profile mismatch")
    if snapshot.provider_selection_sha256 != selection.selection_sha256:
        raise ContractError(
            "campaign host snapshot provider selection mismatch"
        )
    artifact_sha256s = tuple(
        sorted(item.artifact.sha256 for item in observations)
    )
    if snapshot.artifact_sha256s != artifact_sha256s:
        raise ContractError("campaign host snapshot artifact mismatch")
    identity_sha256 = (
        approved_molecular_identity.identity_sha256
        if approved_molecular_identity is not None
        else ""
    )
    if snapshot.approved_identity_sha256 != identity_sha256:
        raise ContractError("campaign host snapshot identity mismatch")


def _validate_campaign_snapshot_integrity(
    snapshot: CampaignPreparationHostSnapshotV1,
) -> None:
    if snapshot.registry.registry_sha256 != snapshot.registry_sha256:
        raise ContractError("campaign host snapshot registry mismatch")
    if snapshot.live_schema.schema_sha256 != snapshot.live_cli_schema_sha256:
        raise ContractError("campaign host snapshot live schema mismatch")
    if snapshot.tool_surface.tool_schema_sha256 != snapshot.tool_schema_sha256:
        raise ContractError("campaign host snapshot tool schema mismatch")
    if snapshot.snapshot_sha256 != canonical_sha256(snapshot._body()):
        raise ContractError("campaign host snapshot digest mismatch")


def validate_campaign_snapshot_binding(
    *,
    snapshot: CampaignPreparationHostSnapshotV1,
    provider_config_file: str | Path,
    provider: str | None,
    artifact_sha256: str,
    approved_identity_sha256: str = "",
) -> None:
    """Validate a controller's path-free source/provider binding.

    The controller has no episode workspace yet, so this entry point checks the
    exact approved coordinate digest and provider selection without probing the
    CLI, environment, or filesystem artifact again.
    """

    if not isinstance(snapshot, CampaignPreparationHostSnapshotV1):
        raise ContractError("campaign host snapshot is not typed")
    _validate_campaign_snapshot_integrity(snapshot)
    selection = load_agent_provider_selection(
        provider_config_file,
        requested_profile=str(provider).strip() if provider else None,
    )
    _validate_campaign_provider_selection(selection)
    require_sha256(artifact_sha256, "artifact_sha256")
    if snapshot.artifact_sha256s != (artifact_sha256,):
        raise ContractError("campaign host snapshot artifact mismatch")
    if (
        snapshot.provider_profile_sha256
        != (selection.active_profile.profile_sha256)
        or snapshot.provider_selection_sha256 != selection.selection_sha256
    ):
        raise ContractError("campaign host snapshot provider mismatch")
    if approved_identity_sha256:
        require_sha256(approved_identity_sha256, "approved_identity_sha256")
    if snapshot.approved_identity_sha256 != approved_identity_sha256:
        raise ContractError("campaign host snapshot identity mismatch")


def _validate_experiment_request(
    *,
    task: str,
    selection: AgentProviderSelectionV1,
    arm: QwenDfcArmV1,
    case: QwenPyscfCaseSpecV1,
    repeat_index: int,
    execution_enabled: bool,
    approval_file: str | Path | None,
) -> None:
    if repeat_index < 0:
        raise ContractError("experiment repeat index must be non-negative")
    if execution_enabled or approval_file is not None:
        raise ContractError("harness experiments are preview-only")
    _validate_campaign_provider_selection(selection)
    if str(task).strip() != case.task:
        raise ContractError(
            "experiment task differs from its preregistered case"
        )
    if not arm.arm_id:
        raise ContractError("experiment arm identity is required")


def _provider_public_record(
    *,
    profile: AgentProviderProfileV1,
    fallback_profiles: Iterable[AgentProviderProfileV1],
    experiment: bool,
) -> dict[str, Any]:
    return {
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
        "fallback_policy": (
            (
                "disabled_for_qwen_experiment"
                if profile.provider == ALIBABA_TOKEN_PLAN_PROVIDER
                and profile.model == ALIBABA_TOKEN_PLAN_MODEL
                else "disabled_for_provider_factorial_experiment"
            )
            if experiment
            else "explicit_attributed_provider_unavailability_only"
        ),
    }


def _experiment_public_record(
    *,
    arm: QwenDfcArmV1,
    case: QwenPyscfCaseSpecV1,
    repeat_index: int,
) -> dict[str, Any]:
    del arm, repeat_index
    return {
        "schema_version": "chemsmart.factor-blind-experiment-context.v1",
        "case_id": case.case_id,
        "case_sha256": case.case_sha256,
        "split": case.split,
        "preview_only": True,
        "engine_calls": 0,
    }


def activated_skill_documents(
    task: str,
) -> tuple[tuple[str, ...], tuple[SkillDocumentV1, ...]]:
    """Resolve the advisory skills a task activates.

    Returns the activated pack digests and the resolved documents.  Activation
    is text-gated by each pack's own ``activation_terms``/``exclusions``; the
    packs are evaluated for every builtin program/engine target so a
    cross-program conventions skill is reachable from any request.
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
    skill_ids: set[str] = set()
    for program, engine in targets:
        receipt = activate_program_knowledge(
            request=task, program=program, engine=engine
        )
        pack_sha256s.update(receipt.activated_pack_sha256s)
        skill_ids.update(skills_for_activation(receipt))
    return tuple(sorted(pack_sha256s)), resolve_skills(
        tuple(sorted(skill_ids))
    )


def _coordinator_base_messages(
    *,
    context: Mapping[str, Any],
    approved_workflow: Mapping[str, Any] | None,
    experiment_arm: QwenDfcArmV1 | None,
    provider_profile: AgentProviderProfileV1 | None = None,
    task: str = "",
) -> list[dict[str, str]]:
    _, documents = activated_skill_documents(task)
    return [
        {
            "role": "system",
            "content": _system_prompt(
                approved_workflow,
                experiment_arm=experiment_arm,
                experiment_provider=(
                    provider_profile.provider
                    if provider_profile is not None
                    else ""
                ),
                skill_index=tuple(item.index_entry() for item in documents),
            ),
        },
        {"role": "user", "content": canonical_json(context)},
    ]


def _observed_experiment_preparation(
    *,
    arm: QwenDfcArmV1,
    case: QwenPyscfCaseSpecV1,
    repeat_index: int,
    task_spec_sha256: str,
    artifact_sha256s: tuple[str, ...],
    provider_profile: AgentProviderProfileV1,
    base_messages: list[dict[str, str]],
    tool_schema_sha256: str,
    host_snapshot_sha256: str,
) -> QwenExperimentPreparationV1:
    return build_qwen_experiment_preparation(
        case=case,
        arm=arm,
        repeat_index=repeat_index,
        task_spec_sha256=task_spec_sha256,
        artifact_sha256s=artifact_sha256s,
        provider_profile_sha256=provider_profile.profile_sha256,
        host_snapshot_sha256=host_snapshot_sha256,
        prompt_sha256=canonical_sha256(base_messages),
        tool_schema_sha256=tool_schema_sha256,
        task_order_sha256=canonical_sha256(
            {
                "case_sha256": case.case_sha256,
                "repeat_index": int(repeat_index),
            }
        ),
        token_budget=provider_profile.context_tokens,
        tool_call_budget=_MAX_TOOL_CALLS,
        wall_time_seconds=_SESSION_WALL_TIME_SECONDS,
        provider_profile=provider_profile,
        episode_namespace=(
            ""
            if provider_profile.provider == ALIBABA_TOKEN_PLAN_PROVIDER
            and provider_profile.model == ALIBABA_TOKEN_PLAN_MODEL
            else provider_profile.provider
        ),
    )


def _require_exact_experiment_config(
    *,
    frozen: HarnessExperimentConfigV1,
    observed: HarnessExperimentConfigV1,
) -> None:
    fields = (
        "schema_version",
        "experiment_id",
        "decomposition",
        "feedback_projection",
        "critic",
        "provider_id",
        "model_id",
        "reasoning_effort",
        "max_concurrency",
        "prompt_sha256",
        "tool_schema_sha256",
        "task_order_sha256",
        "token_budget",
        "tool_call_budget",
        "wall_time_seconds",
        "config_sha256",
    )
    for field in fields:
        if getattr(frozen, field) != getattr(observed, field):
            raise ContractError(f"frozen experiment config mismatch: {field}")


def _bind_preparation_observation(
    *,
    observations: Mapping[str, Any],
    preparation: QwenExperimentPreparationV1,
) -> dict[str, Any]:
    record = dict(observations)
    if record.get("experiment_config_sha256") != (
        preparation.experiment_config.config_sha256
    ):
        raise ContractError(
            "experiment observations contain another configuration"
        )
    record.pop("record_sha256", None)
    record["preparation_sha256"] = preparation.preparation_sha256
    if preparation.host_snapshot_sha256:
        record["host_snapshot_sha256"] = preparation.host_snapshot_sha256
    record["observed_experiment_config_sha256"] = (
        preparation.experiment_config.config_sha256
    )
    return {**record, "record_sha256": canonical_sha256(record)}


def _bind_host_oracle_bundle(
    *,
    observations: Mapping[str, Any],
    bundle: HostOracleInputBundleV1,
) -> dict[str, Any]:
    """Add the coordinator's F-invariant host evidence to an experiment row."""

    record = dict(observations)
    record.pop("record_sha256", None)
    record["host_oracle_input_bundle"] = bundle.public_record()
    return {**record, "record_sha256": canonical_sha256(record)}


def _bind_feedback_receipts(
    *, observations: Mapping[str, Any], events: Iterable[Any]
) -> dict[str, Any]:
    """Persist path-free F-factor receipts in the public experiment record.

    Runtime V2 remains the canonical event source.  The copied receipts let
    the filesystem-independent campaign analyzer measure full versus causal
    feedback without reopening a private run directory.  They contain no
    model reasoning, paths, prompts, or provider-private payloads.
    """

    receipts: list[dict[str, Any]] = []
    existing = observations.get("feedback_receipts", ())
    if not isinstance(existing, (tuple, list)):
        raise ContractError("worker feedback receipts are not a sequence")
    for raw in existing:
        receipts.append(_validated_feedback_receipt(raw))
    for event in events:
        if getattr(event, "kind", "") not in {
            "tool_succeeded",
            "tool_failed",
        }:
            continue
        payload = getattr(event, "payload", None)
        if not isinstance(payload, Mapping):
            raise ContractError("tool event payload is not a mapping")
        raw = payload.get("feedback_equivalence_receipt")
        if not isinstance(raw, Mapping):
            raise ContractError("tool event lacks a feedback receipt")
        receipts.append(_validated_feedback_receipt(raw))
    record = dict(observations)
    record.pop("record_sha256", None)
    record["feedback_receipts"] = tuple(receipts)
    return {**record, "record_sha256": canonical_sha256(record)}


def _validated_feedback_receipt(raw: Any) -> dict[str, Any]:
    if not isinstance(raw, Mapping):
        raise ContractError("tool feedback receipt is not a mapping")
    receipt = canonical_data(dict(raw))
    if receipt.get("schema_version") != (
        "chemsmart.tool-feedback-projection-receipt.v2"
    ):
        raise ContractError("tool event has an unsupported feedback receipt")
    observed_sha256 = str(receipt.get("receipt_sha256") or "")
    body = {
        key: value for key, value in receipt.items() if key != "receipt_sha256"
    }
    if observed_sha256 != canonical_sha256(body):
        raise ContractError("tool feedback receipt digest mismatch")
    return receipt


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


def _private_preparation_directory(
    *, workspace: Path, episode_id: str
) -> Path:
    private_root = workspace / _PRIVATE_ROOT_NAME
    if private_root.is_symlink():
        raise ContractError(
            "private agent directory cannot be a symbolic link"
        )
    private_root.mkdir(exist_ok=True, mode=0o700)
    private_root.chmod(0o700)
    root = private_root / "preparations"
    if root.is_symlink():
        raise ContractError("preparation root cannot be a symbolic link")
    root.mkdir(exist_ok=True, mode=0o700)
    root.chmod(0o700)
    identity = canonical_sha256(
        {
            "schema_version": "chemsmart.live-preparation-directory.v1",
            "episode_id": str(episode_id).lower(),
        }
    )
    episode_root = root / identity
    if episode_root.is_symlink():
        raise ContractError(
            "preparation episode root cannot be a symbolic link"
        )
    episode_root.mkdir(exist_ok=True, mode=0o700)
    target = episode_root / uuid.uuid4().hex
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
    host_node_root = workspace / "nodes"
    for candidate in sorted(workspace.rglob("*.h5")):
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
    result_observations: Iterable[_PySCFResultObservation] = (),
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
                "program": "pyscf",
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
    return canonical_sha256(body)


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
    """Return the core stages ``program`` previews on ``engine``."""

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
    else:
        declared = set(capability.jobtypes)
    return tuple(sorted(declared & _CONFORMANCE_CORE_STAGES))


#: Stages of a route-building program that read their own project section
#: instead of ``gas``/``solv``.  The section carries only ordinary project
#: keys: stage-specific parameters such as the excited-state count are
#: CLI-owned, and putting one here makes the whole project fail to load.
_ROUTE_PROGRAM_STAGE_SECTIONS = frozenset({"td"})


def _conformance_project_sections(
    program: str,
) -> dict[str, dict[str, Any]] | None:
    """Return locked project sections, or ``None`` when unneeded."""

    shape = _CONFORMANCE_PROJECT_SHAPES.get(program)
    if shape == "pyscf_stage_keyed":
        return _locked_pyscf_sections()
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
                sections[stage] = dict(common)
        return sections
    return None


def _bootstrap_conformance(
    *,
    run_directory: Path,
    input_artifact: TrustedArtifactRefV1,
    registry_sha256: str,
    live_schema: Any,
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
        _preview_server_profile().encode("utf-8"),
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


def _ensure_preview_server(run_directory: Path) -> Path:
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
        _preview_server_profile().encode("utf-8"),
    )
    return server_path


def _conformance_engines(program: str) -> tuple[str, ...]:
    """Return the execution engines ChemSmart declares for ``program``."""

    from chemsmart.settings.capabilities import PROGRAM_EXECUTION_ENGINES

    return tuple(PROGRAM_EXECUTION_ENGINES.get(program, ("cpu",)))


def _conformance_programs() -> tuple[str, ...]:
    """Return the programs ChemSmart declares as executable.

    Deriving this from the declaration rather than a hand-written list is what
    makes a newly added program visible to the agent without touching the
    session code.
    """

    from chemsmart.settings.capabilities import EXECUTABLE_PROGRAMS

    return tuple(
        program
        for program in sorted(EXECUTABLE_PROGRAMS)
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
                item.verifier_status,
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
        verifier_receipt_sha256=aggregate("verifier_receipt_sha256"),
        compiler_status=status,
        preview_status=status,
        preflight_status=status,
        verifier_status=status,
    )


def _locked_pyscf_sections() -> dict[str, dict[str, Any]]:
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


def _preview_server_profile() -> str:
    """Scheduler-shaped profile for non-submitting run/sub conformance."""

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
    return (
        "SERVER:\n"
        "  SCHEDULER: PBS\n"
        "  QUEUE_NAME: preview\n"
        "  NUM_HOURS: 1\n"
        "  MEM_GB: 4\n"
        "  NUM_CORES: 4\n"
        "  NUM_GPUS: 0\n"
        "  NUM_THREADS: 4\n"
        "  SUBMIT_COMMAND: true\n"
        "  SCRATCH_DIR: null\n"
        "  PROJECT: preview\n"
        "  USE_HOSTS: false\n" + "".join(blocks)
    )


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
        records.append(
            {
                "record_kind": "program_environment",
                "program": program,
                "engine": "cpu",
                "status": "available" if located else "missing",
                "declared_folder": folder,
                "observation_method": "declared_server_exefolder",
                # A stub satisfies discovery and planning but cannot compute.
                # Recording that keeps "a binary was found" distinct from "a
                # scientific result is obtainable" in the evidence chain.
                "is_discovery_stub": _executable_is_discovery_stub(located),
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
    result_observations: tuple[_PySCFResultObservation, ...] = (),
    conformance_records: tuple[dict[str, Any], ...],
    registry_sha256: str,
    live_schema_sha256: str,
    execution_requested: bool,
    execution_available: bool,
    approved_project_records: tuple[dict[str, Any], ...] = (),
    approved_workflow_record: Mapping[str, Any] | None = None,
    provider_record: Mapping[str, Any] | None = None,
    experiment_record: Mapping[str, Any] | None = None,
    approved_identity_records: tuple[dict[str, Any], ...] = (),
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
        "approved_project_artifacts": approved_project_records,
        "approved_workflow": public_workflow,
        "provider": dict(provider_record or {}),
        "harness_experiment": dict(experiment_record or {}),
        "analysis_completion_policy": dict(analysis_completion_record or {}),
        "approved_execution_contract": _approved_execution_context(
            public_workflow
        ),
        "conformance_observations": conformance_records,
        "execution_requested": bool(execution_requested),
        "execution_tool_available": bool(execution_available),
        "authority": (
            "Use artifact IDs and typed tool fields only. The host owns paths, "
            "CLI argv, project materialization, validation, approval, and execution."
        ),
    }


def _live_specialist_host_builder(
    *,
    base_host_kwargs: Mapping[str, Any],
    coordinator_host: CommandCompiledToolHostV1,
    coordinator_surface: Any,
):
    """Create fresh read-only-session hosts without sharing mutable state."""

    template = dict(base_host_kwargs)

    def build(
        event_store: RuntimeEventStore,
        request: Any,
    ) -> CommandCompiledToolHostV1:
        values = dict(template)
        values["event_store"] = event_store
        # The live specialist facade exposes the narrow read-only projection;
        # the underlying normal host still validates the canonical full profile.
        values["tool_surface"] = coordinator_surface
        if request.role == READ_ONLY_CRITIC:
            # The post-coordinator critic receives detached copies of observable
            # coordinator state.  Its own host and event stream remain fresh.
            values.update(
                {
                    "artifacts": dict(coordinator_host.artifacts),
                    "scientific_identities": dict(
                        coordinator_host.scientific_identities
                    ),
                    "settings_objects": dict(
                        coordinator_host.settings_objects
                    ),
                    "run_receipts": dict(coordinator_host.run_receipts),
                    "scientific_claim_evidence": dict(
                        coordinator_host.scientific_claim_evidence
                    ),
                    "functional_equivalence_receipts": dict(
                        coordinator_host.functional_equivalence_receipts
                    ),
                    "substitution_approvals": dict(
                        coordinator_host.substitution_approvals
                    ),
                    "capability_receipts": dict(coordinator_host.capabilities),
                    "environment_receipts": dict(
                        coordinator_host.environments
                    ),
                    "program_binding_receipts": dict(
                        coordinator_host.program_bindings
                    ),
                    "engine_binding_receipts": dict(
                        coordinator_host.engine_bindings
                    ),
                    "project_validation_receipts": dict(
                        coordinator_host.project_validations
                    ),
                }
            )
        return CommandCompiledToolHostV1(**values)

    return build


def _validate_path_free_experiment_observations(value: Any) -> None:
    """Keep local event/artifact locations outside the public result."""

    forbidden_keys = {
        "cwd",
        "directory",
        "event_store_path",
        "path",
        "run_directory",
        "workspace",
    }

    def visit(item: Any) -> None:
        if isinstance(item, Mapping):
            for key, child in item.items():
                normalized = str(key).strip().lower()
                if normalized in forbidden_keys:
                    raise ContractError(
                        "experiment observation contains a local path field"
                    )
                visit(child)
        elif isinstance(item, (tuple, list)):
            for child in item:
                visit(child)
        elif isinstance(item, str) and (
            item.startswith(("/", "~/"))
            or any(
                marker in item
                for marker in (
                    "/Users/",
                    "/home/",
                    "/private/",
                    "/tmp/",
                    "/var/",
                )
            )
            or (
                len(item) >= 3
                and item[0].isalpha()
                and item[1] == ":"
                and item[2] in {"/", "\\"}
            )
        ):
            raise ContractError(
                "experiment observation contains an absolute local path"
            )

    if not isinstance(value, Mapping):
        raise ContractError("experiment observations must be a mapping")
    visit(value)


def _coordinator_usage_record(
    *,
    event_store: RuntimeEventStore,
    successful_tool_calls: int,
    failed_tool_calls: int,
) -> dict[str, int]:
    attempts = tuple(
        event
        for event in event_store.read_events()
        if event.kind == "api_attempt_observed"
    )
    return {
        "provider_attempts": len(attempts),
        "successful_tool_calls": int(successful_tool_calls),
        "failed_tool_calls": int(failed_tool_calls),
        "input_tokens": sum(
            int(event.payload.get("input_tokens") or 0) for event in attempts
        ),
        "output_tokens": sum(
            int(event.payload.get("output_tokens") or 0) for event in attempts
        ),
        "reasoning_tokens": sum(
            int(event.payload.get("reasoning_tokens") or 0)
            for event in attempts
        ),
        "wall_time_millis": sum(
            int(event.payload.get("latency_ms") or 0) for event in attempts
        ),
    }


def _approved_execution_context(
    public_workflow: Mapping[str, Any],
) -> dict[str, Any]:
    """Describe execution mechanics without inventing scientific settings."""

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
        "resources": {
            "cores": 4,
            "memory_gb": 4,
            "gpu_count": 0,
            "scratch": "none",
        },
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
    experiment_arm: QwenDfcArmV1 | None = None,
    experiment_provider: str = "",
    skill_index: tuple[str, ...] = (),
) -> str:
    execution_available = bool(approved_workflow)
    execution_sentence = (
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
    skill_sentence = ""
    if skill_index:
        listing = " ".join(f"({item})" for item in skill_index)
        skill_sentence = (
            " Domain-knowledge skills are available through "
            "consult_domain_skill(skill_id). Available now: "
            f"{listing}. Consult the relevant skill before you commit to a "
            "reporting convention, an electronic-state assignment, or a "
            "workflow shape it covers, and say when a stated fact came from a "
            "skill. A skill is advisory knowledge only: it never establishes "
            "readiness, approval, terminal state, or an accuracy claim, and it "
            "never replaces a typed host receipt."
        )
    experiment_sentence = ""
    if experiment_arm is not None:
        experiment_label = (
            "Qwen/PySCF"
            if experiment_provider == ALIBABA_TOKEN_PLAN_PROVIDER
            else "provider/PySCF"
        )
        experiment_sentence = (
            f" This is a preregistered {experiment_label} preview episode. Do not "
            "execute a chemistry engine. When specialist_advisory is present, "
            "treat it only as detached typed proposals and resolve it through "
            "normal host tools. Perform the same coordinator duties whether or "
            "not an advisory is present."
        )
    return (
        "You are a professional computational-chemistry planning agent operating "
        "ChemSmart 3.1.4. Work plan-first through typed tools. Inspect program "
        "capability and environment, bind exact artifact identity, render and promote "
        "stage-specific project YAML, validate it, build a scientific tool-chain DAG, compile safe "
        "commands, and preview every currently resolvable node. Keep every future "
        "producer input unresolved until its validated upstream artifact exists. "
        "For every request that ends in a calculated or derived value, use "
        "plan_scientific_workflow to record calculations, result extraction, "
        "validation, mathematics, and claim rendering in one connected DAG. "
        "Its analysis inputs name future producer node/output pairs, so do not wait "
        "for artifact hashes before planning postprocessing. Preserve an unavailable "
        "parser or external analysis as blocked_unsupported instead of deleting the "
        "requested observable. Use plan_command_workflow only for a calculation-only "
        "compatibility task, and use inspect_workflow_frontier for host-derived "
        "next actions. "
        "When the task names a program, plan that program. If its preview is "
        "refused, repair it from the findings preview_command returns, which "
        "name the field, the expected value and the observed one. Only when a "
        "named program still cannot preview green should you plan the "
        "workflow again using a program that can. Every node in the plan "
        "needs a green preview before the plan can be approved, so one node "
        "you have given up on blocks every other node too: plan again without "
        "it rather than leaving it in place beside its replacement. Read "
        "approval_readiness on the planning result for which nodes still "
        "block approval. "
        f"{_workflow_context_sentence()}"
        "Never author native "
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
        "establishes charge or multiplicity. "
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
        "sp, opt, and hess. "
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
        "carry the same source quantity ID, give them distinct task-semantic roles "
        "such as reactant, product, conformer-a, or basis-cardinal-3; never use a "
        "receipt hash as a semantic role. When a numerical condition already exists "
        "as a quantity on a typed receipt, reference that receipt quantity instead "
        "of duplicating it as a literal; use a literal only when no typed source "
        "quantity exists. Use record_analysis_claims to bind each requested reported "
        "number and display unit to an exact receipt quantity; the host, not the "
        "model, supplies the value. The host renders the authoritative final numeric "
        "section from that claim record. Report only those host-rendered claim values. "
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
        + experiment_sentence
    )


def _network_budget(
    provider_profile: AgentProviderProfileV1,
    *,
    max_concurrency: int = 1,
) -> AdaptiveNetworkBudgetV1:
    body = {
        "schema_version": "chemsmart.adaptive-network-budget.v1",
        "allowed_provider": provider_profile.provider,
        "endpoint_origin": provider_profile.endpoint,
        "purpose": "live-command-compiled-computational-chemistry-planning",
        "max_concurrency": int(max_concurrency),
        "max_input_tokens_per_request": provider_profile.context_tokens,
        "max_output_tokens_per_request": provider_profile.max_output_tokens,
        "task_wall_time_seconds": float(_SESSION_WALL_TIME_SECONDS),
        "quota_scope": "current_user_account",
        "top_up_allowed": False,
        "engine_calls": 0,
        "hpc_calls": 0,
    }
    return AdaptiveNetworkBudgetV1(
        **body, budget_sha256=canonical_sha256(body)
    )


def _hypothesis(
    *,
    session_id: str,
    messages: list[dict[str, str]],
    task_spec_sha256: str,
    artifact_sha256s: tuple[str, ...],
    tool_schema_sha256: str,
    network_budget: AdaptiveNetworkBudgetV1,
    execution_requested: bool,
    experiment_config: HarnessExperimentConfigV1 | None = None,
    experiment_case: QwenPyscfCaseSpecV1 | None = None,
    experiment_repeat_index: int = 0,
    dependency_context: TaskDependencyContextV1 | None = None,
) -> AdaptiveHypothesisV1:
    if experiment_config is not None and experiment_case is not None:
        hypothesis_id = experiment_config.experiment_id
        changed_factor = (
            "D/F/C="
            f"d{int(experiment_config.decomposition)}-"
            "f"
            + (
                "causal"
                if experiment_config.feedback_projection == "causal-v1"
                else "full"
            )
            + "-"
            f"c{int(experiment_config.critic)}"
        )
        comparator_id = (
            f"{experiment_case.case_id}.d0-ffull-c0.r{experiment_repeat_index}"
        )
        expected_outcome = experiment_case.expected_observation
        oracle_id = experiment_case.deterministic_oracle_id
        distinct = (
            "Unique preregistered case, D/F/C arm, and repetition tuple."
        )
    elif dependency_context is not None:
        hypothesis_id = session_id
        changed_factor = "task_dependency_context:" + dependency_context.mode
        comparator_id = (
            f"{dependency_context.workflow_id}:no-predecessor-baseline"
        )
        expected_outcome = (
            "The selected predecessor projection preserves every producer "
            "needed by the target while avoiding unrelated workflow branches."
        )
        oracle_id = "paper-workflow-dependency-context-v1"
        distinct = (
            "Unique paper task, dependency-context arm, and live session."
        )
    else:
        hypothesis_id = session_id
        changed_factor = (
            "approval_bound_execution_profile"
            if execution_requested
            else "command_compiled_preview_profile"
        )
        comparator_id = "single-agent-natural-language-baseline"
        expected_outcome = (
            "The model uses typed ChemSmart tools, preserves the exact task, "
            "and leaves every future producer artifact unresolved until host "
            "evidence exists."
        )
        oracle_id = "live-project-command-preview-gates"
        distinct = "Unique live session over an exact task and coordinate-artifact set."
    body = {
        "schema_version": "chemsmart.adaptive-hypothesis.v1",
        "hypothesis_id": hypothesis_id,
        "changed_factor": changed_factor,
        "comparator_id": comparator_id,
        "expected_outcome": expected_outcome,
        "deterministic_oracle_id": oracle_id,
        "source_sha256s": tuple(
            sorted(
                {
                    task_spec_sha256,
                    *artifact_sha256s,
                    *(
                        (experiment_case.case_sha256,)
                        if experiment_case is not None
                        else ()
                    ),
                    *(
                        (dependency_context.plan_sha256,)
                        if dependency_context is not None
                        else ()
                    ),
                }
            )
        ),
        "prompt_sha256": canonical_sha256(messages),
        "tool_schema_sha256": tool_schema_sha256,
        "configuration_sha256": (
            experiment_config.config_sha256
            if experiment_config is not None
            else canonical_sha256(
                {
                    "network_budget": network_budget,
                    "task_spec_sha256": task_spec_sha256,
                    "execution_requested": bool(execution_requested),
                    "dependency_context_sha256": (
                        dependency_context.context_sha256
                        if dependency_context is not None
                        else ""
                    ),
                }
            )
        ),
        "distinct_from_prior": distinct,
    }
    return AdaptiveHypothesisV1(
        **body, hypothesis_sha256=canonical_sha256(body)
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
) -> TaskEnvelopeV1:
    resource = ResourceBudgetV1(
        max_input_tokens_per_request=provider_profile.context_tokens,
        max_output_tokens_per_request=provider_profile.max_output_tokens,
        max_tool_calls=_MAX_TOOL_CALLS,
        wall_time_seconds=float(_SESSION_WALL_TIME_SECONDS),
        chemistry_engine_calls=4 if execution_enabled else 0,
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


def _execution_composition_inputs(
    *,
    host_type: type[CommandCompiledToolHostV1],
    workspace: Path,
    run_directory: Path,
    approval_file: Path,
    task_spec_sha256: str,
) -> dict[str, Any]:
    """Load one exact workflow approval and fixed local resource profile."""

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
    raw_approval = payload.get("workflow_approval", payload)
    if not isinstance(raw_approval, Mapping):
        raise ContractError("workflow approval must be an object")
    resources = _locked_execution_resources()
    approval = _parse_workflow_approval(raw_approval)
    if Path(approval.workspace).resolve() != workspace:
        raise ContractError("workflow approval targets another workspace")
    if approval.task_spec_sha256 != task_spec_sha256:
        raise ContractError("workflow approval targets another task spec")
    if approval.resource_sha256 != resources.resource_sha256:
        raise ContractError(
            "workflow approval resources differ from locked resources"
        )
    server_profile = _write_execution_server_profile(run_directory)
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
        raw["nodes"] = tuple(
            MaterializedNodeV1(**dict(item)) for item in raw.get("nodes", ())
        )
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
    project_root = workspace / "projects"
    if project_root.is_dir() and not project_root.is_symlink():
        for candidate in sorted(project_root.glob("*.yaml")):
            if candidate.is_symlink() or not candidate.is_file():
                continue
            digest = file_sha256(candidate)
            if digest in required:
                if digest in matches:
                    raise ContractError(
                        "approved project bytes have multiple workspace identities"
                    )
                matches[digest] = candidate.resolve()
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
    """Expose the typed execution sequence, and the plan that must be matched.

    Execution refuses unless the session's own plan hashes to the approved
    ``plan_sha256``.  That digest covers nine top-level plan fields, ten fields
    per node and seven per edge; this projection used to disclose one, four and
    two of them.  A session was therefore required to re-derive
    ``project_role``, ``required_observables``, ``complexity_factors`` and
    ``edge_id`` -- free-form strings it has no way to know -- and was told only
    "planned workflow differs from frozen execution approval" when it did not.
    Reproducing an approved plan was not merely hard, it was unspecified.

    Disclosing the approved plan grants nothing: it is the plan the user
    already signed, execution still requires the digest to match, and the
    session cannot widen the approval by reading it.  What it removes is a
    guessing game that no model could win.
    """

    if approved_plan is not None:
        return {
            **_public_workflow_approval(approval),
            "approved_scientific_plan": canonical_data(approved_plan),
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


def _locked_execution_resources() -> ExecutionResourceSpecV1:
    return build_execution_resource_spec(
        execution_target="run",
        cores=4,
        memory_gb=4,
        gpu_count=0,
        scratch_policy="none",
        node_timeout_seconds=_CHEMISTRY_NODE_TIMEOUT_SECONDS,
    )


def _parse_frozen_workflow_approval(
    value: Mapping[str, Any],
) -> FrozenWorkflowApprovalV1:
    """Parse the Runtime V2 body that the execution tool actually requires.

    ``execute_program_node`` refuses a session whose frozen approval is
    ``None`` -- "legacy V1 approval is preview-only" -- and nothing constructed
    one, so every approval file was inert and no node could ever run.

    The frozen body belongs in the approval file next to the V1 approval rather
    than being derived here from the session's own plan.  Deriving it would
    make its ``plan_sha256`` check compare the plan against itself, which is
    the check that stops an approved plan from being swapped for another one.
    A reviewer signs both bodies, and the host cross-checks the pair.
    """

    raw = dict(value)
    try:
        raw["environment_identity_sha256s"] = tuple(
            raw.get("environment_identity_sha256s", ())
        )
        raw["approved_node_ids"] = tuple(raw.get("approved_node_ids", ()))
        raw["producer_edge_sha256s"] = tuple(
            raw.get("producer_edge_sha256s", ())
        )
        raw["materialized_preview_bindings"] = tuple(
            FrozenMaterializedNodePreviewV1(**dict(item))
            for item in raw.get("materialized_preview_bindings", ())
        )
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
        raw["node_bindings"] = tuple(
            ApprovedNodeBindingV1(**dict(item))
            for item in raw.get("node_bindings", ())
        )
        raw["producer_edges"] = tuple(
            ProducerEdgeRuleV1(**dict(item))
            for item in raw.get("producer_edges", ())
        )
        return WorkflowExecutionApprovalV1(**raw)
    except (KeyError, TypeError, ValueError) as exc:
        raise ContractError(
            "workflow approval does not match the v1 schema"
        ) from exc


def _write_execution_server_profile(run_directory: Path) -> Path:
    """Write the deterministic local CPU profile used by approved nodes."""

    xtb = shutil.which("xtb")
    profile = run_directory / "execution-server.yaml"
    text = (
        "SERVER:\n"
        "  SCHEDULER: null\n"
        "  NUM_HOURS: 2\n"
        "  MEM_GB: 4\n"
        "  NUM_CORES: 4\n"
        "  NUM_GPUS: 0\n"
        "  NUM_THREADS: 4\n"
        "  SCRATCH_DIR: null\n"
        "PYSCF:\n"
        f"  EXEFOLDER: {str(_PYSCF_INTERPRETER.parent)!r}\n"
        "  LOCAL_RUN: true\n"
        "  SCRATCH: false\n"
    )
    if xtb:
        text += (
            "XTB:\n"
            f"  EXEFOLDER: {str(Path(xtb).resolve().parent)!r}\n"
            "  LOCAL_RUN: true\n"
            "  SCRATCH: false\n"
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
        "experiment_observations": {},
        "successful_tool_calls": 0,
        "failed_tool_calls": 0,
        "event_stream_head_sha256": "",
    }
    return LiveAgentSessionResultV1(
        **body, result_sha256=canonical_sha256(body)
    )


__all__ = [
    "CampaignPreparationHostSnapshotV1",
    "LiveAgentSessionResultV1",
    "build_campaign_preparation_host_snapshot",
    "probe_live_experiment_preparation",
    "run_live_agent_session",
    "validate_campaign_snapshot_binding",
]
