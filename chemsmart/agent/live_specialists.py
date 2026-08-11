"""Live provider specialist composition for preregistered D/F/C episodes.

The generic contracts and deterministic merge live in ``specialists.py``.
This module supplies the provider-facing composition boundary: every worker
gets a fresh provider-native Runtime V2 session, a private event stream, an
immutable packet, and a read-only projection of the normal ChemSmart host.

Worker results are advisory.  They can add typed proposals to coordinator
context, but they cannot mutate coordinator state or transfer approval,
execution, readiness, filesystem, shell, or terminal authority.
"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import threading
import time
from typing import Any, Callable, Mapping, Sequence

from chemsmart.agent._contracts import (
    ContractError,
    canonical_data,
    canonical_sha256,
    require_identifier,
    require_sha256,
)
from chemsmart.agent.adaptive_api_campaign import (
    AdaptiveHypothesisV1,
    AdaptiveNetworkBudgetV1,
)
from chemsmart.agent.api_access import load_secret_lease
from chemsmart.agent.execution import ScientificDecisionRecordV1
from chemsmart.agent.experiments.host_oracle import HostOracleInputBundleV1
from chemsmart.agent.experiments.qwen_pyscf_dfc import (
    ComplexityGateInputV1,
    ComplexityGateReceiptV1,
    QwenDfcArmV1,
    QwenPyscfCaseSpecV1,
    evaluate_complexity_gate,
)
from chemsmart.agent.provider_config import AgentProviderProfileV1
from chemsmart.agent.runtime.contracts import (
    ResourceBudgetV1,
    TaskEnvelopeV1,
    TaskPhase,
)
from chemsmart.agent.runtime.event_store import RuntimeEventStore
from chemsmart.agent.services.unified_session import UnifiedSessionRunner
from chemsmart.agent.specialists import (
    BoundedSpecialistOrchestratorV1,
    CoordinatorMergeReceiptV1,
    CriticSessionOutcomeV1,
    READ_ONLY_CRITIC,
    SpecialistAdvisoryValidationError,
    SpecialistBudgetV1,
    SpecialistSessionRequestV1,
    SpecialistSessionResponseV1,
    SpecialistSessionOutcomeV1,
    build_specialist_session_response,
    validate_specialist_advisory_value,
)
from chemsmart.agent.tool_specs import AgentToolSurfaceV1
from chemsmart.agent.workflows import (
    HarnessExperimentConfigV1,
    ScientificWorkflowEdgeV2,
    ScientificWorkflowNodeV2,
    ScientificWorkflowPlanV2,
    build_scientific_workflow_plan,
)


_WORKER_ROOT = "specialist-sessions"
_FROZEN_SPECIALIST_BUDGET_SLOTS = 4


@dataclass(frozen=True)
class FInvariantCriticCandidateV1:
    """Projection-independent coordinator record shown to a fresh critic."""

    schema_version: str
    candidate_id: str
    task_spec_sha256: str
    host_oracle_input_bundle: HostOracleInputBundleV1
    coordinator_public_decisions: tuple[ScientificDecisionRecordV1, ...]
    candidate_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.f-invariant-critic-candidate.v1":
            raise ContractError("unsupported F-invariant critic candidate")
        require_identifier(self.candidate_id, "candidate_id")
        require_sha256(self.task_spec_sha256, "task_spec_sha256")
        HostOracleInputBundleV1.from_record(
            self.host_oracle_input_bundle.public_record()
        )
        order = tuple(
            (item.decision_id, item.record_sha256)
            for item in self.coordinator_public_decisions
        )
        if order != tuple(sorted(order)):
            raise ContractError(
                "public scientific decisions are not canonical"
            )
        if len(order) != len(set(order)):
            raise ContractError("public scientific decisions are duplicated")
        if any(
            item.task_spec_sha256 != self.task_spec_sha256
            for item in self.coordinator_public_decisions
        ):
            raise ContractError(
                "public scientific decision targets another task spec"
            )
        if self.candidate_sha256 != canonical_sha256(self.review_record()):
            raise ContractError("F-invariant critic candidate digest mismatch")

    def review_record(self) -> dict[str, Any]:
        """Return the immutable path-free record supplied to the critic."""

        return {
            "schema_version": self.schema_version,
            "candidate_id": self.candidate_id,
            "task_spec_sha256": self.task_spec_sha256,
            "host_oracle_input_bundle": (
                self.host_oracle_input_bundle.public_record()
            ),
            "coordinator_public_decisions": tuple(
                canonical_data(item)
                for item in self.coordinator_public_decisions
            ),
        }


def build_f_invariant_critic_candidate(
    *,
    candidate_id: str,
    task_spec_sha256: str,
    host_oracle_input_bundle: HostOracleInputBundleV1,
    coordinator_public_decisions: Sequence[ScientificDecisionRecordV1],
) -> FInvariantCriticCandidateV1:
    """Build critic input without transcript or feedback-projection content."""

    decisions = tuple(
        sorted(
            coordinator_public_decisions,
            key=lambda item: (item.decision_id, item.record_sha256),
        )
    )
    values = {
        "schema_version": "chemsmart.f-invariant-critic-candidate.v1",
        "candidate_id": str(candidate_id).strip().lower(),
        "task_spec_sha256": str(task_spec_sha256).strip().lower(),
        "host_oracle_input_bundle": host_oracle_input_bundle,
        "coordinator_public_decisions": decisions,
    }
    body = {
        "schema_version": values["schema_version"],
        "candidate_id": values["candidate_id"],
        "task_spec_sha256": values["task_spec_sha256"],
        "host_oracle_input_bundle": host_oracle_input_bundle.public_record(),
        "coordinator_public_decisions": tuple(
            canonical_data(item) for item in decisions
        ),
    }
    return FInvariantCriticCandidateV1(
        **values, candidate_sha256=canonical_sha256(body)
    )


@dataclass(frozen=True)
class LiveSpecialistProviderObservationV1:
    """Path-free observable facts from one disposable provider session."""

    session_id: str
    role: str
    feedback_projection: str
    request_sha256: str
    response_sha256: str
    provider_terminal_state: str
    provider_attempts: int
    successful_tool_calls: int
    failed_tool_calls: int
    input_tokens: int
    output_tokens: int
    reasoning_tokens: int
    wall_time_millis: int
    public_transcript_sha256: str
    event_stream_head_sha256: str
    feedback_receipts: tuple[Mapping[str, Any], ...] = ()
    output_envelope_receipt: Mapping[str, Any] | None = None
    nonsecret_error_class: str = ""

    def __post_init__(self) -> None:
        if self.output_envelope_receipt is not None:
            SpecialistOutputEnvelopeReceiptV1(
                **dict(self.output_envelope_receipt)
            )

    def public_record(self) -> dict[str, Any]:
        return canonical_data(self)


@dataclass(frozen=True)
class SpecialistOutputEnvelopeReceiptV1:
    """Distinguish raw JSON compliance from bounded wire normalization."""

    schema_version: str
    mode: str
    raw_text_sha256: str
    normalized_object_sha256: str
    ignored_prefix_bytes: int
    ignored_suffix_bytes: int
    receipt_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.specialist-output-envelope.v1":
            raise ContractError("unsupported specialist output envelope")
        if self.mode not in {"strict_json", "single_json_object_extracted"}:
            raise ContractError("unsupported specialist output envelope mode")
        require_sha256(self.raw_text_sha256, "raw_text_sha256")
        require_sha256(
            self.normalized_object_sha256, "normalized_object_sha256"
        )
        if min(self.ignored_prefix_bytes, self.ignored_suffix_bytes) < 0:
            raise ContractError("ignored specialist envelope bytes cannot be negative")
        if self.mode == "strict_json" and (
            self.ignored_prefix_bytes or self.ignored_suffix_bytes
        ):
            raise ContractError("strict specialist JSON cannot ignore envelope bytes")
        if self.receipt_sha256 != canonical_sha256(self.receipt_record()):
            raise ContractError("specialist output envelope digest mismatch")

    def receipt_record(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "mode": self.mode,
            "raw_text_sha256": self.raw_text_sha256,
            "normalized_object_sha256": self.normalized_object_sha256,
            "ignored_prefix_bytes": self.ignored_prefix_bytes,
            "ignored_suffix_bytes": self.ignored_suffix_bytes,
        }

    def public_record(self) -> dict[str, Any]:
        return {**self.receipt_record(), "receipt_sha256": self.receipt_sha256}


@dataclass(frozen=True)
class _SpecialistRoleDispatchResultV1:
    """One role-local result, including a path-free failure observation."""

    role: str
    outcomes: tuple[SpecialistSessionOutcomeV1, ...] = ()
    error_class: str = ""
    validation_findings: tuple[Mapping[str, Any], ...] = ()


def _dispatch_pre_coordinator_specialists(
    *,
    roles: Sequence[str],
    max_concurrency: int,
    orchestrator: BoundedSpecialistOrchestratorV1,
    base_tool_surface: AgentToolSurfaceV1,
    session_factory: "LiveProviderSpecialistSessionFactoryV1",
    plan: ScientificWorkflowPlanV2,
    coordinator_session_id: str,
    public_context: Mapping[str, Any],
    source_sha256s: Sequence[str],
    artifact_sha256s: Sequence[str],
    input_record_sha256s: Sequence[str],
    budget: SpecialistBudgetV1,
) -> tuple[_SpecialistRoleDispatchResultV1, ...]:
    """Run independent advisory roles with a bounded provider fan-out.

    The one-worker branch intentionally uses the original shared orchestrator
    and exact role order.  Parallel workers receive role-local orchestrators so
    request construction remains isolated; the shared session factory supplies
    a distinct directory, event store, and credential lease per immutable
    request.  Completion order never affects the returned or merged order.
    """

    normalized_roles = tuple(str(role) for role in roles)
    if not normalized_roles:
        return ()
    if not 1 <= int(max_concurrency) <= 4:
        raise ContractError(
            "specialist provider concurrency must be within [1, 4]"
        )

    common = {
        "plan": plan,
        "coordinator_session_id": coordinator_session_id,
        "public_context": public_context,
        "source_sha256s": source_sha256s,
        "artifact_sha256s": artifact_sha256s,
        "input_record_sha256s": input_record_sha256s,
        "budget": budget,
    }

    def run_role(
        role: str,
        role_orchestrator: BoundedSpecialistOrchestratorV1,
    ) -> _SpecialistRoleDispatchResultV1:
        try:
            outcomes = role_orchestrator.run_before_coordinator(
                **common, roles=(role,)
            )
            return _SpecialistRoleDispatchResultV1(
                role=role,
                outcomes=tuple(outcomes),
            )
        except Exception as exc:
            findings: tuple[Mapping[str, Any], ...] = ()
            if isinstance(exc, SpecialistAdvisoryValidationError):
                findings = (exc.public_finding(),)
            return _SpecialistRoleDispatchResultV1(
                role=role,
                error_class=type(exc).__name__,
                validation_findings=findings,
            )

    if max_concurrency == 1:
        return tuple(run_role(role, orchestrator) for role in normalized_roles)

    # BoundedSpecialistOrchestratorV1 deliberately owns mutable request/session
    # uniqueness state.  Give every concurrent role its own instance rather
    # than sharing that state across threads.
    role_orchestrators = {
        role: BoundedSpecialistOrchestratorV1(
            base_tool_surface=base_tool_surface,
            session_factory=session_factory,
        )
        for role in normalized_roles
    }
    with ThreadPoolExecutor(
        max_workers=min(max_concurrency, len(normalized_roles)),
        thread_name_prefix="chemsmart-specialist",
    ) as executor:
        futures = {
            role: executor.submit(run_role, role, role_orchestrators[role])
            for role in normalized_roles
        }
        # Resolve in role order, not provider completion order.
        return tuple(futures[role].result() for role in sorted(futures))


@dataclass
class LiveSpecialistCampaignV1:
    """One experiment-local bounded specialist and critic composition."""

    arm: QwenDfcArmV1
    experiment_config: HarnessExperimentConfigV1
    plan: ScientificWorkflowPlanV2
    gate: ComplexityGateReceiptV1
    budget: SpecialistBudgetV1
    orchestrator: BoundedSpecialistOrchestratorV1
    session_factory: "LiveProviderSpecialistSessionFactoryV1"
    outcomes: tuple[SpecialistSessionOutcomeV1, ...] = ()
    merge: CoordinatorMergeReceiptV1 | None = None
    critic_candidate: FInvariantCriticCandidateV1 | None = None
    critic_outcome: CriticSessionOutcomeV1 | None = None
    specialist_error_class: str = ""
    specialist_validation_findings: tuple[Mapping[str, Any], ...] = ()
    critic_error_class: str = ""
    specialist_provider_concurrency_limit: int = 1
    specialist_provider_concurrency_observed: int = 0
    specialist_dispatch_wall_time_millis: int = 0

    @classmethod
    def start(
        cls,
        *,
        arm: QwenDfcArmV1,
        experiment_config: HarnessExperimentConfigV1,
        plan: ScientificWorkflowPlanV2,
        coordinator_session_id: str,
        public_context: Mapping[str, Any],
        source_sha256s: Sequence[str],
        artifact_sha256s: Sequence[str],
        base_tool_surface: AgentToolSurfaceV1,
        provider_profile: AgentProviderProfileV1,
        secret_file: str | Path,
        run_directory: str | Path,
        host_builder: Callable[
            [RuntimeEventStore, SpecialistSessionRequestV1], Any
        ],
        runner_factory: Callable[..., Any] = UnifiedSessionRunner,
        lease_loader: Callable[..., Any] = load_secret_lease,
    ) -> "LiveSpecialistCampaignV1":
        runtime_config = provider_profile.runtime_config()
        if (
            experiment_config.provider_id != provider_profile.provider
            or experiment_config.model_id != provider_profile.model
            or experiment_config.reasoning_effort
            != provider_profile.reasoning_effort
            or runtime_config.provider != provider_profile.provider
            or runtime_config.model != provider_profile.model
            or runtime_config.reasoning_effort
            != provider_profile.reasoning_effort
        ):
            raise ContractError(
                "specialist campaign differs from its provider profile"
            )
        if (
            experiment_config.decomposition != arm.decomposition
            or experiment_config.feedback_projection != arm.feedback_projection
            or experiment_config.critic != arm.critic
            or experiment_config.max_concurrency != arm.max_concurrency
        ):
            raise ContractError("specialist campaign differs from its D/F/C arm")

        gate = evaluate_complexity_gate(
            arm, _complexity_signals_from_plan(plan)
        )
        worker_count = len(gate.requested_roles) if gate.activated else 0
        # D and C must not silently change the allowance of another worker.
        # Freeze each disposable specialist/critic to the same quarter-episode
        # allowance (three possible pre-coordinator roles plus one critic).
        # Unused slots are not reassigned.
        budget = derive_specialist_budget(
            experiment_config=experiment_config,
            participant_slots=_FROZEN_SPECIALIST_BUDGET_SLOTS,
        )
        session_factory = LiveProviderSpecialistSessionFactoryV1(
            provider_profile=provider_profile,
            feedback_projection=arm.feedback_projection,
            secret_file=secret_file,
            run_directory=run_directory,
            host_builder=host_builder,
            runner_factory=runner_factory,
            lease_loader=lease_loader,
        )
        orchestrator = BoundedSpecialistOrchestratorV1(
            base_tool_surface=base_tool_surface,
            session_factory=session_factory,
        )
        campaign = cls(
            arm=arm,
            experiment_config=experiment_config,
            plan=plan,
            gate=gate,
            budget=budget,
            orchestrator=orchestrator,
            session_factory=session_factory,
            specialist_provider_concurrency_limit=min(
                experiment_config.max_concurrency,
                max(1, worker_count),
            ),
        )
        if gate.activated:
            dispatch_started = session_factory.clock()
            role_results = _dispatch_pre_coordinator_specialists(
                roles=gate.requested_roles,
                max_concurrency=campaign.specialist_provider_concurrency_limit,
                orchestrator=orchestrator,
                base_tool_surface=base_tool_surface,
                session_factory=session_factory,
                plan=plan,
                coordinator_session_id=coordinator_session_id,
                public_context=public_context,
                source_sha256s=source_sha256s,
                artifact_sha256s=artifact_sha256s,
                input_record_sha256s=(
                    experiment_config.config_sha256,
                    gate.receipt_sha256,
                ),
                budget=budget,
            )
            campaign.specialist_dispatch_wall_time_millis = max(
                0, int((session_factory.clock() - dispatch_started) * 1000)
            )
            campaign.specialist_provider_concurrency_observed = (
                session_factory.max_provider_concurrency_observed()
            )
            accepted_outcomes: list[SpecialistSessionOutcomeV1] = []
            validation_findings: list[Mapping[str, Any]] = []
            for result in role_results:
                accepted_outcomes.extend(result.outcomes)
                if result.error_class and not campaign.specialist_error_class:
                    campaign.specialist_error_class = result.error_class
                validation_findings.extend(result.validation_findings)
            campaign.outcomes = tuple(
                sorted(
                    accepted_outcomes,
                    key=lambda item: item.result_packet.role,
                )
            )
            campaign.specialist_validation_findings = tuple(
                sorted(
                    validation_findings,
                    key=lambda item: (
                        str(item.get("role") or ""),
                        str(item.get("rule_id") or ""),
                    ),
                )
            )
            if campaign.outcomes:
                try:
                    campaign.merge = orchestrator.merge_before_coordinator(
                        campaign.outcomes
                    )
                except ContractError as exc:
                    if not campaign.specialist_error_class:
                        campaign.specialist_error_class = type(exc).__name__
        return campaign

    def coordinator_advisory_record(self) -> dict[str, Any]:
        """Return proposals as non-authoritative, path-free coordinator input."""

        specialists = tuple(
            {
                "role": outcome.result_packet.role,
                "status": outcome.result_packet.status,
                "public_summary": outcome.result_packet.public_summary,
                "result_sha256": outcome.result_packet.result_sha256,
                "outcome_sha256": outcome.outcome_sha256,
            }
            for outcome in sorted(
                self.outcomes, key=lambda item: item.result_packet.role
            )
        )
        merged_fields = ()
        conflicts = ()
        merge_status = "not_dispatched"
        merge_sha256 = ""
        unresolved_fields = ()
        if self.merge is not None:
            merged_fields = tuple(
                {
                    "field_path": item.field_path,
                    "value": json.loads(item.value_json),
                    "source_roles": item.source_roles,
                    "evidence_sha256s": item.evidence_sha256s,
                    "merged_sha256": item.merged_sha256,
                }
                for item in self.merge.merged_fields
            )
            conflicts = tuple(
                {
                    "field_path": item.field_path,
                    "conflict_sha256": item.conflict_sha256,
                }
                for item in self.merge.conflicts
            )
            merge_status = self.merge.status
            merge_sha256 = self.merge.receipt_sha256
            unresolved_fields = self.merge.unresolved_fields
        elif self.specialist_error_class:
            merge_status = "failed"
        elif self.gate.activated:
            merge_status = "no_candidate"
        return {
            "schema_version": "chemsmart.coordinator-specialist-advisory.v1",
            "authority": (
                "advisory_only; coordinator and deterministic host retain "
                "identity, project, DAG, readiness, and terminal authority"
            ),
            "gate": canonical_data(self.gate),
            "specialists": specialists,
            "merge": {
                "status": merge_status,
                "receipt_sha256": merge_sha256,
                "merged_fields": merged_fields,
                "conflicts": conflicts,
                "unresolved_fields": unresolved_fields,
                "nonsecret_error_class": self.specialist_error_class,
                "validation_findings": tuple(
                    canonical_data(item)
                    for item in self.specialist_validation_findings
                ),
            },
        }

    def run_critic(
        self,
        *,
        coordinator_session_id: str,
        candidate: FInvariantCriticCandidateV1,
        public_context: Mapping[str, Any],
        source_sha256s: Sequence[str],
        artifact_sha256s: Sequence[str],
    ) -> None:
        if not self.arm.critic:
            return
        self.critic_candidate = candidate
        try:
            self.critic_outcome = self.orchestrator.run_after_coordinator_critic(
                plan=self.plan,
                coordinator_session_id=coordinator_session_id,
                candidate_id=candidate.candidate_id,
                candidate_record=candidate.review_record(),
                public_context=public_context,
                source_sha256s=source_sha256s,
                artifact_sha256s=artifact_sha256s,
                input_record_sha256s=(
                    self.experiment_config.config_sha256,
                    self.gate.receipt_sha256,
                ),
                budget=self.budget,
            )
            if self.critic_outcome.candidate_sha256_before != (
                candidate.candidate_sha256
            ):
                raise ContractError("critic outcome targets another candidate")
        except ContractError as exc:
            self.critic_error_class = type(exc).__name__

    def public_observation_record(
        self,
        *,
        coordinator_usage: Mapping[str, int] | None = None,
    ) -> dict[str, Any]:
        """Build the path-free experiment record returned by the live session."""

        observations = self.session_factory.public_observations()
        critic: dict[str, Any] = {
            "status": "not_enabled" if not self.arm.critic else "failed",
            "candidate_sha256": (
                self.critic_candidate.candidate_sha256
                if self.critic_candidate is not None
                else ""
            ),
            "candidate_record": (
                self.critic_candidate.review_record()
                if self.critic_candidate is not None
                else {}
            ),
            "review_sha256": "",
            "review": {},
            "outcome_sha256": "",
            "public_summary": "",
            "findings": (),
            "nonsecret_error_class": self.critic_error_class,
        }
        if self.critic_outcome is not None:
            critic = {
                "status": self.critic_outcome.review.status,
                "candidate_sha256": (
                    self.critic_outcome.review.candidate_sha256
                ),
                "candidate_record": (
                    self.critic_candidate.review_record()
                    if self.critic_candidate is not None
                    else {}
                ),
                "review_sha256": self.critic_outcome.review.review_sha256,
                "review": canonical_data(self.critic_outcome.review),
                "outcome_sha256": self.critic_outcome.outcome_sha256,
                "public_summary": (
                    self.critic_outcome.result_packet.public_summary
                ),
                "findings": tuple(
                    canonical_data(item)
                    for item in self.critic_outcome.review.findings
                ),
                "nonsecret_error_class": "",
            }
        specialist_rows = tuple(
            {
                "role": item.result_packet.role,
                "status": item.result_packet.status,
                "public_summary": item.result_packet.public_summary,
                "result_packet": canonical_data(item.result_packet),
                "candidate": canonical_data(item.candidate),
                "result_sha256": item.result_packet.result_sha256,
                "outcome_sha256": item.outcome_sha256,
            }
            for item in sorted(
                self.outcomes, key=lambda row: row.result_packet.role
            )
        )
        coordinator = {
            str(key): int(value)
            for key, value in dict(coordinator_usage or {}).items()
        }
        specialist_observations = tuple(
            item for item in observations if item.role != READ_ONLY_CRITIC
        )
        critic_observations = tuple(
            item for item in observations if item.role == READ_ONLY_CRITIC
        )
        provider_wall_time_sum = sum(
            item.wall_time_millis for item in observations
        )
        critic_observed_wall_time = sum(
            item.wall_time_millis for item in critic_observations
        )
        worker_observed_wall_time = (
            self.specialist_dispatch_wall_time_millis
            + critic_observed_wall_time
        )
        worker_usage = {
            "provider_sessions": len(observations),
            "provider_attempts": sum(
                item.provider_attempts for item in observations
            ),
            "successful_tool_calls": sum(
                item.successful_tool_calls for item in observations
            ),
            "failed_tool_calls": sum(
                item.failed_tool_calls for item in observations
            ),
            "input_tokens": sum(item.input_tokens for item in observations),
            "output_tokens": sum(item.output_tokens for item in observations),
            "reasoning_tokens": sum(
                item.reasoning_tokens for item in observations
            ),
            # The sum is useful for provider-load accounting, but it is not
            # elapsed critical-path time when specialist calls overlap.
            "provider_wall_time_sum_millis": provider_wall_time_sum,
            "observed_wall_time_millis": worker_observed_wall_time,
            "wall_time_millis": worker_observed_wall_time,
            "specialist_dispatch_wall_time_millis": (
                self.specialist_dispatch_wall_time_millis
            ),
            "critic_wall_time_millis": critic_observed_wall_time,
            "provider_concurrency_limit": (
                self.specialist_provider_concurrency_limit
            ),
            "provider_concurrency_observed": (
                self.specialist_provider_concurrency_observed
            ),
        }
        output_envelope_rows = tuple(
            {
                "role": item.role,
                "mode": str(
                    (item.output_envelope_receipt or {}).get("mode") or ""
                ),
                "receipt_sha256": str(
                    (item.output_envelope_receipt or {}).get(
                        "receipt_sha256"
                    )
                    or ""
                ),
            }
            for item in observations
        )
        usage = {
            "provider_sessions": (
                worker_usage["provider_sessions"]
                + int(bool(coordinator_usage))
            ),
            "provider_attempts": (
                worker_usage["provider_attempts"]
                + coordinator.get("provider_attempts", 0)
            ),
            "successful_tool_calls": (
                worker_usage["successful_tool_calls"]
                + coordinator.get("successful_tool_calls", 0)
            ),
            "failed_tool_calls": (
                worker_usage["failed_tool_calls"]
                + coordinator.get("failed_tool_calls", 0)
            ),
            "input_tokens": (
                worker_usage["input_tokens"]
                + coordinator.get("input_tokens", 0)
            ),
            "output_tokens": (
                worker_usage["output_tokens"]
                + coordinator.get("output_tokens", 0)
            ),
            "reasoning_tokens": (
                worker_usage["reasoning_tokens"]
                + coordinator.get("reasoning_tokens", 0)
            ),
            "wall_time_millis": (
                worker_usage["wall_time_millis"]
                + coordinator.get("wall_time_millis", 0)
            ),
            "provider_wall_time_sum_millis": (
                worker_usage["provider_wall_time_sum_millis"]
                + coordinator.get("wall_time_millis", 0)
            ),
            "worker": worker_usage,
            "coordinator": coordinator,
        }
        body = {
            "schema_version": (
                "chemsmart.live-harness-experiment-observations.v1"
            ),
            "experiment_config_sha256": (
                self.experiment_config.config_sha256
            ),
            "feedback_projection": self.session_factory.feedback_projection,
            "gate": canonical_data(self.gate),
            "specialist_budget": canonical_data(self.budget),
            "specialist_dispatch": {
                "requested_roles": tuple(self.gate.requested_roles),
                "provider_concurrency_limit": (
                    self.specialist_provider_concurrency_limit
                ),
                "provider_concurrency_observed": (
                    self.specialist_provider_concurrency_observed
                ),
                "observed_wall_time_millis": (
                    self.specialist_dispatch_wall_time_millis
                ),
                "provider_session_wall_time_sum_millis": sum(
                    item.wall_time_millis for item in specialist_observations
                ),
            },
            "specialists": specialist_rows,
            "merge": (
                {
                    "status": self.merge.status,
                    "receipt_sha256": self.merge.receipt_sha256,
                    "conflict_sha256s": tuple(
                        item.conflict_sha256 for item in self.merge.conflicts
                    ),
                }
                if self.merge is not None
                else {
                    "status": (
                        "failed"
                        if self.specialist_error_class
                        else "not_dispatched"
                    ),
                    "receipt_sha256": "",
                    "conflict_sha256s": (),
                }
            ),
            "critic": critic,
            "provider_sessions": tuple(
                item.public_record() for item in observations
            ),
            "specialist_output_envelopes": {
                "strict_json_count": sum(
                    item["mode"] == "strict_json"
                    for item in output_envelope_rows
                ),
                "normalized_count": sum(
                    item["mode"] == "single_json_object_extracted"
                    for item in output_envelope_rows
                ),
                "rows": output_envelope_rows,
            },
            "feedback_receipts": tuple(
                receipt
                for item in observations
                for receipt in item.feedback_receipts
            ),
            "usage": usage,
            "nonsecret_specialist_error_class": self.specialist_error_class,
            "specialist_validation_findings": tuple(
                canonical_data(item)
                for item in self.specialist_validation_findings
            ),
        }
        return {**body, "record_sha256": canonical_sha256(body)}


class LiveProviderSpecialistSessionFactoryV1:
    """Create one fresh UnifiedSessionRunner session for every worker."""

    def __init__(
        self,
        *,
        provider_profile: AgentProviderProfileV1,
        feedback_projection: str,
        secret_file: str | Path,
        run_directory: str | Path,
        host_builder: Callable[
            [RuntimeEventStore, SpecialistSessionRequestV1], Any
        ],
        runner_factory: Callable[..., Any] = UnifiedSessionRunner,
        lease_loader: Callable[..., Any] = load_secret_lease,
        clock: Callable[[], float] = time.monotonic,
    ) -> None:
        runtime_config = provider_profile.runtime_config()
        if (
            runtime_config.provider != provider_profile.provider
            or runtime_config.model != provider_profile.model
            or runtime_config.reasoning_effort
            != provider_profile.reasoning_effort
        ):
            raise ContractError(
                "specialist provider differs from its runtime adapter"
            )
        self.provider_profile = provider_profile
        if feedback_projection not in {"full-v1", "causal-v1"}:
            raise ContractError("invalid specialist feedback projection")
        self.feedback_projection = feedback_projection
        self.secret_file = Path(secret_file)
        self.run_directory = Path(run_directory)
        self.host_builder = host_builder
        self.runner_factory = runner_factory
        self.lease_loader = lease_loader
        self.clock = clock
        self._observations: list[LiveSpecialistProviderObservationV1] = []
        self._observation_lock = threading.Lock()
        self._active_provider_sessions = 0
        self._max_provider_concurrency_observed = 0

    def __call__(
        self, request: SpecialistSessionRequestV1
    ) -> "_LiveProviderSpecialistSessionV1":
        return _LiveProviderSpecialistSessionV1(factory=self, request=request)

    def public_observations(
        self,
    ) -> tuple[LiveSpecialistProviderObservationV1, ...]:
        with self._observation_lock:
            observations = tuple(self._observations)
        return tuple(
            sorted(observations, key=lambda item: (item.role, item.session_id))
        )

    def max_provider_concurrency_observed(self) -> int:
        """Return the peak number of overlapping disposable provider calls."""

        with self._observation_lock:
            return self._max_provider_concurrency_observed

    def _begin_provider_session(self) -> None:
        with self._observation_lock:
            self._active_provider_sessions += 1
            self._max_provider_concurrency_observed = max(
                self._max_provider_concurrency_observed,
                self._active_provider_sessions,
            )

    def _end_provider_session(self) -> None:
        with self._observation_lock:
            self._active_provider_sessions -= 1
            if self._active_provider_sessions < 0:
                self._active_provider_sessions = 0
                raise ContractError("provider session concurrency underflow")

    def _record_observation(
        self, observation: LiveSpecialistProviderObservationV1
    ) -> None:
        with self._observation_lock:
            self._observations.append(observation)

    def _run(
        self, request: SpecialistSessionRequestV1
    ) -> SpecialistSessionResponseV1:
        started = self.clock()
        worker_directory = (
            self.run_directory / _WORKER_ROOT / request.session_id
        )
        worker_directory.mkdir(parents=True, mode=0o700)
        worker_directory.chmod(0o700)
        event_store = RuntimeEventStore(
            worker_directory / "events.jsonl", session_id=request.session_id
        )
        base_host = self.host_builder(event_store, request)
        surface = AgentToolSurfaceV1(
            schema_version="chemsmart.agent-tool-surface.v1",
            profile=f"{request.role}-read-only",
            tool_definitions=request.tool_definitions,
            tool_schema_sha256=request.context_manifest.tool_schema_sha256,
        )
        host = _ReadOnlySpecialistHostV1(
            base_host=base_host,
            surface=surface,
        )
        messages = _specialist_messages(request)
        envelope = _specialist_envelope(
            request,
            messages=messages,
            provider_profile=self.provider_profile,
        )
        network_budget = _specialist_network_budget(
            request, profile=self.provider_profile
        )
        hypothesis = _specialist_hypothesis(
            request,
            messages=messages,
            network_budget=network_budget,
        )
        lease = self.lease_loader(
            provider=self.provider_profile.provider,
            path=self.secret_file,
            ttl_seconds=request.context_manifest.wall_time_seconds + 60,
        )
        response_sha256 = ""
        output_envelope_record: Mapping[str, Any] | None = None
        transcript_sha256 = ""
        stream_head_sha256 = ""
        terminal_state = "failed"
        successful_tools = 0
        failed_tools = 0
        try:
            result = self.runner_factory(
                host=host,
                event_store=event_store,
                credential_lease=lease,
                provider_config=self.provider_profile.runtime_config(),
            ).run(
                messages=messages,
                envelope=envelope,
                hypothesis=hypothesis,
                network_budget=network_budget,
                feedback_projection=self.feedback_projection,
            )
            terminal_state = str(result.terminal_state)
            successful_tools = int(result.successful_tool_calls)
            failed_tools = int(result.failed_tool_calls)
            transcript_sha256 = str(
                getattr(result, "public_transcript_sha256", "")
                or canonical_sha256(result.public_transcript)
            )
            stream_head_sha256 = str(result.event_stream_head_sha256)
            output, output_envelope = _bounded_public_json_object(
                result.final_text
            )
            output_envelope_record = output_envelope.public_record()
            tool_calls = _tool_call_names(result.public_transcript)
            input_tokens, output_tokens, reasoning_tokens = _event_usage(
                event_store
            )
            response = build_specialist_session_response(
                session_id=request.session_id,
                public_output=output,
                tool_calls=tool_calls,
                input_tokens=input_tokens,
                output_tokens=output_tokens,
                wall_time_millis=max(
                    0, int((self.clock() - started) * 1000)
                ),
            )
            response_sha256 = response.response_sha256
            self._record_observation(
                LiveSpecialistProviderObservationV1(
                    session_id=request.session_id,
                    role=request.role,
                    feedback_projection=self.feedback_projection,
                    request_sha256=request.request_sha256,
                    response_sha256=response_sha256,
                    provider_terminal_state=terminal_state,
                    provider_attempts=_event_attempt_count(event_store),
                    successful_tool_calls=successful_tools,
                    failed_tool_calls=failed_tools,
                    input_tokens=input_tokens,
                    output_tokens=output_tokens,
                    reasoning_tokens=reasoning_tokens,
                    wall_time_millis=max(
                        0, int((self.clock() - started) * 1000)
                    ),
                    public_transcript_sha256=transcript_sha256,
                    event_stream_head_sha256=stream_head_sha256,
                    feedback_receipts=_event_feedback_receipts(event_store),
                    output_envelope_receipt=output_envelope_record,
                )
            )
            return response
        except Exception as exc:
            input_tokens, output_tokens, reasoning_tokens = _event_usage(
                event_store
            )
            self._record_observation(
                LiveSpecialistProviderObservationV1(
                    session_id=request.session_id,
                    role=request.role,
                    feedback_projection=self.feedback_projection,
                    request_sha256=request.request_sha256,
                    response_sha256=response_sha256,
                    provider_terminal_state=terminal_state,
                    provider_attempts=_event_attempt_count(event_store),
                    successful_tool_calls=successful_tools,
                    failed_tool_calls=failed_tools,
                    input_tokens=input_tokens,
                    output_tokens=output_tokens,
                    reasoning_tokens=reasoning_tokens,
                    wall_time_millis=max(
                        0, int((self.clock() - started) * 1000)
                    ),
                    public_transcript_sha256=transcript_sha256,
                    event_stream_head_sha256=stream_head_sha256,
                    feedback_receipts=_event_feedback_receipts(event_store),
                    output_envelope_receipt=output_envelope_record,
                    nonsecret_error_class=(
                        _event_error_class(event_store)
                        or type(exc).__name__
                    ),
                )
            )
            raise


class _LiveProviderSpecialistSessionV1:
    def __init__(
        self,
        *,
        factory: LiveProviderSpecialistSessionFactoryV1,
        request: SpecialistSessionRequestV1,
    ) -> None:
        self.factory = factory
        self.request = request
        self.session_id = request.session_id
        self.closed = False

    def run(
        self, request: SpecialistSessionRequestV1
    ) -> SpecialistSessionResponseV1:
        if self.closed or request is not self.request:
            raise ContractError("specialist session request identity changed")
        self.factory._begin_provider_session()
        try:
            return self.factory._run(request)
        finally:
            self.factory._end_provider_session()

    def close(self) -> None:
        self.closed = True


class _ReadOnlySpecialistHostV1:
    """Narrow surface facade over a fresh normal command-compiled host."""

    def __init__(self, *, base_host: Any, surface: AgentToolSurfaceV1) -> None:
        self._base_host = base_host
        self.surface = surface
        self._allowed = frozenset(
            item["function"]["name"] for item in surface.tool_definitions
        )

    def record_seeded_evidence(self, turn_id: str) -> None:
        self._base_host.record_seeded_evidence(turn_id)

    def dispatch(
        self, *, turn_id: str, tool_name: str, arguments: Mapping[str, Any]
    ) -> dict[str, Any]:
        if tool_name not in self._allowed:
            raise ContractError("tool is outside the read-only specialist surface")
        return self._base_host.dispatch(
            turn_id=turn_id, tool_name=tool_name, arguments=arguments
        )

    def completion_receipts_for_latest_preflight(self) -> tuple[str, ...]:
        raise ContractError("specialists cannot set scientific readiness")

    def latest_workflow_draft_receipt(self) -> str:
        raise ContractError("specialists cannot own the coordinator workflow")


def derive_specialist_budget(
    *,
    experiment_config: HarnessExperimentConfigV1,
    participant_slots: int,
) -> SpecialistBudgetV1:
    """Derive explicit per-session bounds from the frozen episode budget."""

    if participant_slots < 1:
        raise ContractError("specialist budget requires a participant slot")
    return SpecialistBudgetV1(
        token_budget=max(1, experiment_config.token_budget // participant_slots),
        tool_call_budget=max(
            1, experiment_config.tool_call_budget // participant_slots
        ),
        wall_time_seconds=max(
            1, experiment_config.wall_time_seconds // participant_slots
        ),
    )


def build_experiment_seed_plan(
    *,
    case: QwenPyscfCaseSpecV1,
    task_spec_sha256: str,
    artifact_sha256s: Sequence[str],
) -> ScientificWorkflowPlanV2:
    """Parse only request-visible stages into a non-executable advisory DAG.

    This seed is not a benchmark answer or a command plan.  It gives bounded
    workers a shared, immutable description of explicitly requested stages and
    producer semantics before the coordinator creates authoritative YAML and
    command objects.
    """

    task = case.task.lower()
    stage_candidates: list[str] = []
    stage_markers = (
        ("sp", ("single point", "sp(", " sp request", "energy,")),
        ("opt", ("optimization", "optimized", "opt(", "relaxation")),
        ("hess", ("hess", "curvature")),
        ("td", ("tda", "tddft", "vertical-excitation")),
        ("ts", ("transition-state", "transition state")),
        ("irc", ("irc",)),
    )
    for stage, markers in stage_markers:
        if any(marker in task for marker in markers):
            stage_candidates.append(stage)
    if "ts" in stage_candidates and (
        "transition-state optimization" in task
        or "transition state optimization" in task
    ):
        stage_candidates = [
            stage for stage in stage_candidates if stage != "opt"
        ]
    if not stage_candidates:
        stage_candidates.append("unspecified")
    order = {name: index for index, name in enumerate(
        ("sp", "opt", "hess", "td", "ts", "irc", "unspecified")
    )}
    stages = tuple(sorted(set(stage_candidates), key=order.__getitem__))

    unresolved: set[str] = set()
    if "no electronic-structure method has been selected" in task:
        unresolved.add("scientific.method")
    if "charge and multiplicity have not been established" in task:
        unresolved.update(("identity.charge", "identity.multiplicity"))
    if "gives neither model nor solvent identity" in task:
        unresolved.update(("project.solvent", "project.solvent_model"))
    if "unspecified" in stages:
        unresolved.add("workflow.stage")

    engine = "gpu" if "gpu4pyscf" in task else "cpu"
    nodes = tuple(
        ScientificWorkflowNodeV2(
            node_id=f"candidate-{stage}",
            stage=stage,
            requested_program="pyscf",
            program="pyscf",
            engine=engine,
            project_role=f"candidate.pyscf.{stage}",
            unresolved_fields=tuple(sorted(unresolved)),
        )
        for stage in stages
    )
    edges: list[ScientificWorkflowEdgeV2] = []
    if "opt" in stages:
        for consumer in ("hess", "td"):
            if consumer in stages and any(
                marker in task
                for marker in (
                    "optimized geometry",
                    "optimized-geometry",
                    "after optimization",
                    "after relaxation",
                    "derived after relaxation",
                )
            ):
                edges.append(
                    ScientificWorkflowEdgeV2(
                        edge_id=f"opt-to-{consumer}-geometry",
                        source_node_id="candidate-opt",
                        target_node_id=f"candidate-{consumer}",
                        edge_kind="data",
                        artifact_class="geometry_xyz",
                        producer_output_id="optimized-geometry",
                        consumer_input_id="geometry",
                    )
                )
    if "ts" in stages and "irc" in stages:
        edges.append(
            ScientificWorkflowEdgeV2(
                edge_id="ts-to-irc-geometry",
                source_node_id="candidate-ts",
                target_node_id="candidate-irc",
                edge_kind="data",
                artifact_class="geometry_xyz",
                producer_output_id="transition-state-geometry",
                consumer_input_id="geometry",
            )
        )
    explicit_factors = (
        ("solvent_semantics",)
        if "solvent" in task or "solvated" in task
        else ()
    )
    return build_scientific_workflow_plan(
        workflow_id=f"experiment.{case.case_id.lower()}",
        task_spec_sha256=task_spec_sha256,
        scientific_identity_sha256=canonical_sha256(
            {
                "task_spec_sha256": task_spec_sha256,
                "artifact_sha256s": tuple(sorted(set(artifact_sha256s))),
                "binding_state": "coordinator_unbound",
            }
        ),
        nodes=nodes,
        edges=tuple(sorted(edges, key=lambda item: item.edge_id)),
        explicit_complexity_factors=explicit_factors,
    )


def _complexity_signals_from_plan(
    plan: ScientificWorkflowPlanV2,
) -> ComplexityGateInputV1:
    factors = set(plan.complexity_factors)
    return ComplexityGateInputV1(
        multi_stage="multiple_stages" in factors,
        unresolved_scientific_settings=(
            "unresolved_scientific_setting" in factors
        ),
        program_substitution="program_substitution" in factors,
        solvent_semantics="solvent_semantics" in factors,
        gpu_request="gpu_engine" in factors,
        excited_state="excited_state" in factors,
        conflicting_evidence="conflicting_evidence" in factors,
        producer_artifact_edges="producer_artifact_edge" in factors,
    )


def _specialist_messages(
    request: SpecialistSessionRequestV1,
) -> list[dict[str, str]]:
    return [
        {"role": "system", "content": request.system_instruction},
        {"role": "user", "content": request.public_context_json},
    ]


def _specialist_envelope(
    request: SpecialistSessionRequestV1,
    *,
    messages: Sequence[Mapping[str, Any]],
    provider_profile: AgentProviderProfileV1,
) -> TaskEnvelopeV1:
    context = request.context_manifest
    resource = ResourceBudgetV1(
        max_input_tokens_per_request=context.token_budget,
        # The provider adapter owns its wire maximum; the stricter combined
        # worker budget is checked on the typed response by the orchestrator.
        max_output_tokens_per_request=provider_profile.max_output_tokens,
        max_tool_calls=context.tool_call_budget,
        wall_time_seconds=float(context.wall_time_seconds),
        chemistry_engine_calls=0,
        hpc_calls=0,
    )
    body = {
        "schema_version": "chemsmart.task-envelope.v1",
        "task_id": request.task_packet.packet_id,
        "session_id": request.session_id,
        "turn_id": request.session_id + ".turn-1",
        "request_sha256": canonical_sha256(messages),
        "cwd_sha256": canonical_sha256(
            {
                "context_manifest_sha256": context.manifest_sha256,
                "artifact_sha256s": context.artifact_sha256s,
            }
        ),
        "phase": (
            TaskPhase.REVIEW
            if request.role == READ_ONLY_CRITIC
            else TaskPhase.SPECIFY
        ),
        "budget": resource,
        "tool_schema_sha256": context.tool_schema_sha256,
    }
    return TaskEnvelopeV1(
        **body, envelope_sha256=canonical_sha256(body)
    )


def _specialist_network_budget(
    request: SpecialistSessionRequestV1,
    *,
    profile: AgentProviderProfileV1,
) -> AdaptiveNetworkBudgetV1:
    context = request.context_manifest
    body = {
        "schema_version": "chemsmart.adaptive-network-budget.v1",
        "allowed_provider": profile.provider,
        "endpoint_origin": profile.endpoint,
        "purpose": f"bounded-{request.role}-advisory",
        "max_concurrency": 1,
        "max_input_tokens_per_request": context.token_budget,
        "max_output_tokens_per_request": profile.max_output_tokens,
        "task_wall_time_seconds": float(context.wall_time_seconds),
        "quota_scope": "current_user_account",
        "top_up_allowed": False,
        "engine_calls": 0,
        "hpc_calls": 0,
    }
    return AdaptiveNetworkBudgetV1(
        **body, budget_sha256=canonical_sha256(body)
    )


def _specialist_hypothesis(
    request: SpecialistSessionRequestV1,
    *,
    messages: Sequence[Mapping[str, Any]],
    network_budget: AdaptiveNetworkBudgetV1,
) -> AdaptiveHypothesisV1:
    body = {
        "schema_version": "chemsmart.adaptive-hypothesis.v1",
        "hypothesis_id": request.session_id,
        "changed_factor": f"fresh-{request.role}-session",
        "comparator_id": "single-coordinator-without-this-advisory-session",
        "expected_outcome": (
            "Return a typed, public, read-only scientific advisory record "
            "without host-authority transfer."
        ),
        "deterministic_oracle_id": "specialist.typed-read-only-boundary.v1",
        "source_sha256s": tuple(
            sorted(
                {
                    request.request_sha256,
                    request.context_manifest.manifest_sha256,
                    request.task_packet.packet_sha256,
                }
            )
        ),
        "prompt_sha256": canonical_sha256(messages),
        "tool_schema_sha256": request.context_manifest.tool_schema_sha256,
        "configuration_sha256": network_budget.budget_sha256,
        "distinct_from_prior": (
            "Unique immutable task packet, role, and provider session ID."
        ),
    }
    return AdaptiveHypothesisV1(
        **body, hypothesis_sha256=canonical_sha256(body)
    )


def _bounded_public_json_object(
    text: str,
) -> tuple[dict[str, Any], SpecialistOutputEnvelopeReceiptV1]:
    """Parse strict JSON or one unambiguous object wrapped in plain prose.

    The normalized path is intentionally narrow.  It never interprets the
    ignored text as evidence, rejects Markdown fences and multiple/braced
    payloads, and records exact raw-versus-normalized digests.  This keeps a
    serialization lapse distinct from scientific correctness without silently
    reporting the normalized result as raw model success.
    """

    raw = str(text)
    stripped = raw.strip()
    try:
        value = json.loads(
            stripped, object_pairs_hook=_json_object_without_duplicate_keys
        )
        mode = "strict_json"
        prefix = ""
        suffix = ""
    except json.JSONDecodeError:
        if "```" in raw:
            raise ContractError(
                "specialist final response must not use a Markdown fence"
            )
        start = raw.find("{")
        if start < 0:
            raise ContractError(
                "specialist final response does not contain a JSON object"
            )
        decoder = json.JSONDecoder(
            object_pairs_hook=_json_object_without_duplicate_keys
        )
        try:
            value, consumed = decoder.raw_decode(raw[start:])
        except json.JSONDecodeError as exc:
            raise ContractError(
                "specialist final response has no complete JSON object"
            ) from exc
        end = start + consumed
        prefix = raw[:start]
        suffix = raw[end:]
        if any(token in prefix + suffix for token in ("{", "}")):
            raise ContractError(
                "specialist final response contains multiple or ambiguous objects"
            )
        if not prefix.strip() or suffix.strip():
            raise ContractError(
                "specialist JSON normalization requires one plain-prose prefix"
            )
        if len(prefix.encode("utf-8")) > 4096 or prefix.rstrip()[-1] not in (
            ".",
            ":",
            ";",
            "?",
            "!",
        ):
            raise ContractError(
                "specialist JSON prefix is not bounded plain prose"
            )
        try:
            validate_specialist_advisory_value(prefix)
        except ContractError as exc:
            raise ContractError(
                "specialist JSON prefix contains executable or native-input text"
            ) from exc
        mode = "single_json_object_extracted"
    if not isinstance(value, dict):
        raise ContractError("specialist final response must be a JSON object")
    body = {
        "schema_version": "chemsmart.specialist-output-envelope.v1",
        "mode": mode,
        "raw_text_sha256": hashlib.sha256(raw.encode("utf-8")).hexdigest(),
        "normalized_object_sha256": canonical_sha256(value),
        "ignored_prefix_bytes": len(prefix.encode("utf-8")),
        "ignored_suffix_bytes": len(suffix.encode("utf-8")),
    }
    receipt = SpecialistOutputEnvelopeReceiptV1(
        **body, receipt_sha256=canonical_sha256(body)
    )
    return value, receipt


def _json_object_without_duplicate_keys(
    pairs: Sequence[tuple[str, Any]],
) -> dict[str, Any]:
    value: dict[str, Any] = {}
    for key, item in pairs:
        if key in value:
            raise ContractError(
                f"specialist JSON contains duplicate key {key!r}"
            )
        value[key] = item
    return value


def _tool_call_names(
    transcript: Sequence[Mapping[str, Any]],
) -> tuple[str, ...]:
    names: list[str] = []
    for message in transcript:
        for call in message.get("tool_calls") or ():
            function = call.get("function") if isinstance(call, Mapping) else None
            if isinstance(function, Mapping) and function.get("name"):
                names.append(str(function["name"]))
    return tuple(names)


def _event_usage(
    event_store: RuntimeEventStore,
) -> tuple[int, int, int]:
    input_tokens = 0
    output_tokens = 0
    reasoning_tokens = 0
    for event in event_store.read_events():
        if event.kind != "api_attempt_observed":
            continue
        input_tokens += int(event.payload.get("input_tokens") or 0)
        output_tokens += int(event.payload.get("output_tokens") or 0)
        reasoning_tokens += int(event.payload.get("reasoning_tokens") or 0)
    return input_tokens, output_tokens, reasoning_tokens


def _event_feedback_receipts(
    event_store: RuntimeEventStore,
) -> tuple[dict[str, Any], ...]:
    """Return canonical public feedback receipts from one worker stream."""

    receipts: list[dict[str, Any]] = []
    for event in event_store.read_events():
        if event.kind not in {"tool_succeeded", "tool_failed"}:
            continue
        raw = event.payload.get("feedback_equivalence_receipt")
        if not isinstance(raw, Mapping):
            raise ContractError("worker tool event lacks a feedback receipt")
        receipt = canonical_data(dict(raw))
        if receipt.get("schema_version") != (
            "chemsmart.tool-feedback-projection-receipt.v2"
        ):
            raise ContractError("worker feedback receipt schema is unsupported")
        observed_sha256 = str(receipt.get("receipt_sha256") or "")
        body = {
            key: value
            for key, value in receipt.items()
            if key != "receipt_sha256"
        }
        if observed_sha256 != canonical_sha256(body):
            raise ContractError("worker feedback receipt digest mismatch")
        receipts.append(receipt)
    return tuple(receipts)


def _event_error_class(event_store: RuntimeEventStore) -> str:
    for event in reversed(event_store.read_events()):
        if event.kind != "api_attempt_observed":
            continue
        error_class = str(
            event.payload.get("nonsecret_error_class") or ""
        ).strip()
        if error_class:
            return error_class
    return ""


def _event_attempt_count(event_store: RuntimeEventStore) -> int:
    return sum(
        1
        for event in event_store.read_events()
        if event.kind == "api_attempt_observed"
    )


# Backward-compatible Python alias; historical Qwen schemas remain unchanged.
LiveQwenSpecialistSessionFactoryV1 = LiveProviderSpecialistSessionFactoryV1


__all__ = [
    "FInvariantCriticCandidateV1",
    "LiveProviderSpecialistSessionFactoryV1",
    "LiveQwenSpecialistSessionFactoryV1",
    "LiveSpecialistCampaignV1",
    "LiveSpecialistProviderObservationV1",
    "SpecialistOutputEnvelopeReceiptV1",
    "build_experiment_seed_plan",
    "build_f_invariant_critic_candidate",
    "derive_specialist_budget",
]
