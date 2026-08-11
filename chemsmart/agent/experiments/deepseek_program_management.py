"""H0/H1 DeepSeek V4 Flash program-management campaign.

This module composes the active ``UnifiedSessionRunner`` path.  It does not
implement another provider loop, accept a credential value, or establish
scientific readiness.  The secret-file path and private output directory are
runtime-only arguments.  Persisted records contain the already-sanitized
public transcript, observable events, usage, latency, and deterministic
structural grades; provider-private reasoning never enters a record.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from enum import Enum
import json
import os
from pathlib import Path
import re
import time
from typing import Any, Callable, Mapping

from chemsmart.agent._contracts import (
    ContractError,
    canonical_data,
    canonical_sha256,
    file_sha256,
)
from chemsmart.agent.adaptive_api_campaign import (
    AdaptiveApiCampaignPolicyV1,
    AdaptiveCampaignStatusReceiptV1,
    AdaptiveHypothesisV1,
    AdaptiveNetworkBudgetV1,
    CampaignTermination,
)
from chemsmart.agent.api_access import load_secret_lease
from chemsmart.agent.experiments.program_management_context import (
    CampaignHostFixtureV1,
    CampaignPublicContextV1,
    host_state_sha256,
)
from chemsmart.agent.runtime.contracts import (
    ResourceBudgetV1,
    TaskEnvelopeV1,
    TaskPhase,
)
from chemsmart.agent.runtime.deepseek import (
    DEEPSEEK_OFFICIAL_ENDPOINT,
    DEEPSEEK_V4_FLASH_CONTEXT_TOKENS,
    DEEPSEEK_V4_FLASH_MAX_OUTPUT_TOKENS,
    DEEPSEEK_V4_FLASH_MODEL,
    DeepSeekV4FlashConfigV1,
    public_payload,
)
from chemsmart.agent.runtime.event_store import RuntimeEventStore
from chemsmart.agent.runtime.events import EventKind, RuntimeEvent
from chemsmart.agent.services.unified_session import UnifiedSessionRunner
from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1
from chemsmart.agent.tool_specs import (
    build_command_compiled_tool_surface,
    build_single_agent_baseline_tool_surface,
)


CAMPAIGN_SYSTEM_PROMPT_V1 = """You are ChemSmart's computational-chemistry planning agent. Respond in English. Use only the typed tools exposed in this session. Preserve the requested molecule, geometry binding, charge, multiplicity, program, engine, method, basis, job kind, and artifact dependencies. Do not write native engine input, Python scripts, shell syntax, paths, executable status, approval, readiness, or terminal state. A proposal is not evidence. Continue from host tool results, and report an explicit blocked condition when a required deterministic receipt is absent or red."""

_SCHEMA = "chemsmart.deepseek-program-management-campaign.v1"
_CASE_ID = re.compile(r"^DS-PM-[0-9]{3}$")
_PRIVATE_KEYS = frozenset(
    {
        "analysis",
        "api_key",
        "authorization",
        "credential",
        "reasoning_content",
        "request_headers",
        "secret",
        "thinking",
    }
)
_TRANSIENT_ERRORS = frozenset(
    {"provider_5xx", "rate_limited", "timeout", "transport"}
)
_MANAGEMENT_TOOLS = frozenset(
    {
        "assess_program_candidate",
        "inspect_program_capability",
        "inspect_program_environment",
    }
)
_NO_COMMAND_CASES = frozenset(
    {
        "DS-PM-004",
        "DS-PM-005",
        "DS-PM-006",
        "DS-PM-008",
        "DS-PM-010",
        "DS-PM-011",
        "DS-PM-012",
        "DS-PM-013",
        "DS-PM-014",
    }
)


class CampaignArm(str, Enum):
    H0 = "H0"
    H1 = "H1"


@dataclass(frozen=True)
class CampaignCaseV1:
    case_id: str
    hypothesis: str
    changed_factor: str
    request: str
    oracle: str
    case_sha256: str

    def __post_init__(self) -> None:
        if _CASE_ID.fullmatch(self.case_id) is None:
            raise ContractError("campaign case ID is not canonical")
        if not all(
            value.strip()
            for value in (
                self.hypothesis,
                self.changed_factor,
                self.request,
                self.oracle,
            )
        ):
            raise ContractError("campaign case fields must not be empty")
        body = {
            "case_id": self.case_id,
            "hypothesis": self.hypothesis,
            "changed_factor": self.changed_factor,
            "request": self.request,
            "oracle": self.oracle,
        }
        if self.case_sha256 != canonical_sha256(body):
            raise ContractError("campaign case digest mismatch")


@dataclass(frozen=True)
class CampaignDefinitionV1:
    schema_version: str
    campaign_id: str
    model: str
    source_sha256: str
    transport_wall_time_seconds: int
    fixed_factors: tuple[str, ...]
    safety_veto: tuple[str, ...]
    cases: tuple[CampaignCaseV1, ...]
    definition_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != _SCHEMA:
            raise ContractError("unsupported campaign definition schema")
        if self.model != DEEPSEEK_V4_FLASH_MODEL:
            raise ContractError("campaign must remain pinned to DeepSeek V4 Flash")
        if not 1 <= self.transport_wall_time_seconds <= 14_400:
            raise ContractError("campaign wall time must be within four hours")
        ids = tuple(item.case_id for item in self.cases)
        if not ids or len(ids) != len(set(ids)):
            raise ContractError("campaign cases must be non-empty and unique")
        body = {
            "schema_version": self.schema_version,
            "campaign_id": self.campaign_id,
            "model": self.model,
            "source_sha256": self.source_sha256,
            "transport_wall_time_seconds": self.transport_wall_time_seconds,
            "fixed_factors": self.fixed_factors,
            "safety_veto": self.safety_veto,
            "cases": self.cases,
        }
        if self.definition_sha256 != canonical_sha256(body):
            raise ContractError("campaign definition digest mismatch")


@dataclass(frozen=True)
class CampaignRunConfigV1:
    schema_version: str
    prompt_version: str
    system_prompt_sha256: str
    workspace_sha256: str
    campaign_wall_time_seconds: int
    episode_wall_time_seconds: int
    max_input_tokens_per_request: int
    max_output_tokens_per_request: int
    max_tool_calls: int
    max_concurrency: int
    endpoint_origin: str
    quota_scope: str
    top_up_allowed: bool
    engine_calls: int
    hpc_calls: int
    configuration_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.deepseek-campaign-config.v1":
            raise ContractError("unsupported campaign configuration schema")
        if self.system_prompt_sha256 != canonical_sha256(
            CAMPAIGN_SYSTEM_PROMPT_V1
        ):
            raise ContractError("campaign system prompt digest mismatch")
        if not 1 <= self.campaign_wall_time_seconds <= 14_400:
            raise ContractError("campaign may run for at most four hours")
        if not 1 <= self.episode_wall_time_seconds <= (
            self.campaign_wall_time_seconds
        ):
            raise ContractError("episode wall time is outside campaign bound")
        if self.max_input_tokens_per_request > (
            DEEPSEEK_V4_FLASH_CONTEXT_TOKENS
        ) or self.max_output_tokens_per_request > (
            DEEPSEEK_V4_FLASH_MAX_OUTPUT_TOKENS
        ):
            raise ContractError("campaign token limits exceed provider support")
        if min(
            self.max_input_tokens_per_request,
            self.max_output_tokens_per_request,
            self.max_tool_calls,
        ) < 1:
            raise ContractError("campaign budgets must be positive")
        if not 1 <= self.max_concurrency <= 4:
            raise ContractError("campaign concurrency must be within [1, 4]")
        if self.endpoint_origin != DEEPSEEK_OFFICIAL_ENDPOINT:
            raise ContractError("campaign requires the official endpoint")
        if self.top_up_allowed or self.engine_calls or self.hpc_calls:
            raise ContractError("campaign cannot top up or run chemistry")
        body = _config_body(self)
        if self.configuration_sha256 != canonical_sha256(body):
            raise ContractError("campaign configuration digest mismatch")


@dataclass(frozen=True)
class CampaignEpisodePlanV1:
    schema_version: str
    episode_id: str
    attempt_ordinal: int
    arm: CampaignArm
    case: CampaignCaseV1
    pairing_sha256: str
    order_index: int
    public_context: CampaignPublicContextV1
    messages: tuple[dict[str, str], ...]
    hypothesis: AdaptiveHypothesisV1
    network_budget: AdaptiveNetworkBudgetV1
    envelope: TaskEnvelopeV1
    plan_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.deepseek-campaign-episode.v1":
            raise ContractError("unsupported campaign episode schema")
        if self.attempt_ordinal < 0 or self.order_index < 0:
            raise ContractError("campaign ordinals must be non-negative")
        if self.envelope.tool_schema_sha256 != (
            self.hypothesis.tool_schema_sha256
        ):
            raise ContractError("episode tool schema binding mismatch")
        context_messages = tuple(
            item
            for item in self.messages
            if item.get("content", "").startswith(
                "CHEMSMART_CAMPAIGN_PUBLIC_CONTEXT_V1\n"
            )
        )
        if context_messages != (self.public_context.provider_message(),):
            raise ContractError(
                "episode messages do not contain the exact public context"
            )
        if self.network_budget.top_up_allowed:
            raise ContractError("episode cannot authorize a top-up")
        body = _episode_body(self)
        if self.plan_sha256 != canonical_sha256(body):
            raise ContractError("campaign episode digest mismatch")


@dataclass(frozen=True)
class DeterministicOracleReceiptV1:
    schema_version: str
    plan_sha256: str
    checks: tuple[tuple[str, str], ...]
    safety_violations: tuple[str, ...]
    verdict: str
    receipt_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.campaign-oracle-receipt.v1":
            raise ContractError("unsupported campaign oracle schema")
        if self.verdict not in {"pass", "fail", "inconclusive"}:
            raise ContractError("invalid deterministic oracle verdict")
        if self.checks != tuple(sorted(self.checks)):
            raise ContractError("oracle checks must be deterministically sorted")
        body = {
            "schema_version": self.schema_version,
            "plan_sha256": self.plan_sha256,
            "checks": self.checks,
            "safety_violations": self.safety_violations,
            "verdict": self.verdict,
        }
        if self.receipt_sha256 != canonical_sha256(body):
            raise ContractError("campaign oracle receipt digest mismatch")


@dataclass(frozen=True)
class CampaignEpisodeReceiptV1:
    schema_version: str
    plan_sha256: str
    terminal_state: str
    public_transcript_sha256: str
    event_stream_head_sha256: str
    oracle_receipt_sha256: str
    public_artifact_sha256: str
    transport_attempts: int
    successful_provider_calls: int
    input_tokens: int
    output_tokens: int
    reasoning_tokens: int
    latency_ms: int
    nonsecret_error_class: str
    receipt_sha256: str

    def __post_init__(self) -> None:
        body = {key: value for key, value in asdict(self).items() if key != "receipt_sha256"}
        if self.schema_version != "chemsmart.campaign-episode-receipt.v1":
            raise ContractError("unsupported episode receipt schema")
        if min(
            self.transport_attempts,
            self.successful_provider_calls,
            self.input_tokens,
            self.output_tokens,
            self.reasoning_tokens,
            self.latency_ms,
        ) < 0:
            raise ContractError("episode counters must be non-negative")
        if self.receipt_sha256 != canonical_sha256(body):
            raise ContractError("episode receipt digest mismatch")


@dataclass(frozen=True)
class CampaignRunReceiptV1:
    schema_version: str
    campaign_definition_sha256: str
    configuration_sha256: str
    episode_receipts: tuple[CampaignEpisodeReceiptV1, ...]
    status_receipt: AdaptiveCampaignStatusReceiptV1
    manifest_sha256: str


CampaignHostInputsFactory = Callable[
    [CampaignCaseV1], CampaignHostFixtureV1
]


def load_campaign_definition(path: str | Path) -> CampaignDefinitionV1:
    source = Path(path)
    raw = json.loads(source.read_text(encoding="utf-8"))
    if raw.get("schema") != "chemsmart-provider-campaign-v1":
        raise ContractError("unsupported source campaign schema")
    mode = raw.get("mode") or {}
    transport = raw.get("transport") or {}
    if mode != {
        "thinking": "enabled",
        "reasoning_effort": "max",
        "output_limit": "provider_maximum",
    }:
        raise ContractError("source campaign weakens thinking configuration")
    if (
        transport.get("attempt_count_limit") is not None
        or bool(transport.get("top_up_allowed"))
        or bool(transport.get("persist_private_reasoning"))
    ):
        raise ContractError("source campaign violates adaptive safety policy")
    if set((raw.get("arms") or {}).keys()) != {"H0", "H1"}:
        raise ContractError("source campaign requires exactly H0 and H1")
    cases = []
    for value in raw.get("cases") or ():
        body = {
            "case_id": str(value.get("case_id") or ""),
            "hypothesis": str(value.get("hypothesis") or ""),
            "changed_factor": str(value.get("changed_factor") or ""),
            "request": str(value.get("request") or ""),
            "oracle": str(value.get("oracle") or ""),
        }
        cases.append(CampaignCaseV1(**body, case_sha256=canonical_sha256(body)))
    body = {
        "schema_version": _SCHEMA,
        "campaign_id": str(raw.get("campaign_id") or ""),
        "model": str(raw.get("model") or ""),
        "source_sha256": file_sha256(source),
        "transport_wall_time_seconds": int(transport.get("wall_time_seconds", 0)),
        "fixed_factors": tuple(str(item) for item in raw.get("fixed_factors") or ()),
        "safety_veto": tuple(str(item) for item in raw.get("safety_veto") or ()),
        "cases": tuple(cases),
    }
    return CampaignDefinitionV1(**body, definition_sha256=canonical_sha256(body))


def build_campaign_run_config(
    *,
    workspace_identity: str,
    campaign_wall_time_seconds: int = 14_400,
    episode_wall_time_seconds: int = 480,
    max_input_tokens_per_request: int = DEEPSEEK_V4_FLASH_CONTEXT_TOKENS,
    max_output_tokens_per_request: int = DEEPSEEK_V4_FLASH_MAX_OUTPUT_TOKENS,
    max_tool_calls: int = 24,
    max_concurrency: int = 1,
) -> CampaignRunConfigV1:
    body = {
        "schema_version": "chemsmart.deepseek-campaign-config.v1",
        "prompt_version": "chemsmart-program-management-en-v1",
        "system_prompt_sha256": canonical_sha256(CAMPAIGN_SYSTEM_PROMPT_V1),
        "workspace_sha256": canonical_sha256({"workspace": workspace_identity}),
        "campaign_wall_time_seconds": int(campaign_wall_time_seconds),
        "episode_wall_time_seconds": int(episode_wall_time_seconds),
        "max_input_tokens_per_request": int(max_input_tokens_per_request),
        "max_output_tokens_per_request": int(max_output_tokens_per_request),
        "max_tool_calls": int(max_tool_calls),
        "max_concurrency": int(max_concurrency),
        "endpoint_origin": DEEPSEEK_OFFICIAL_ENDPOINT,
        "quota_scope": "current_user_account",
        "top_up_allowed": False,
        "engine_calls": 0,
        "hpc_calls": 0,
    }
    return CampaignRunConfigV1(
        **body, configuration_sha256=canonical_sha256(body)
    )


def build_episode_plans(
    definition: CampaignDefinitionV1,
    config: CampaignRunConfigV1,
    *,
    public_contexts: Mapping[str, CampaignPublicContextV1],
    additional_source_sha256s: tuple[str, ...] = (),
) -> tuple[CampaignEpisodePlanV1, ...]:
    """Build a deterministic, alternating H0/H1 paired schedule."""

    plans = []
    expected_case_ids = {item.case_id for item in definition.cases}
    if set(public_contexts) != expected_case_ids:
        raise ContractError(
            "campaign public contexts must cover the exact case set"
        )
    order = 0
    for case_index, case in enumerate(definition.cases):
        public_context = public_contexts[case.case_id]
        if public_context.case_id != case.case_id:
            raise ContractError("campaign public context targets another case")
        pairing = canonical_sha256(
            {"campaign": definition.definition_sha256, "case": case.case_sha256}
        )
        arms = (
            (CampaignArm.H0, CampaignArm.H1)
            if case_index % 2 == 0
            else (CampaignArm.H1, CampaignArm.H0)
        )
        for arm in arms:
            plans.append(
                _build_episode_plan(
                    definition,
                    config,
                    case=case,
                    arm=arm,
                    pairing_sha256=pairing,
                    order_index=order,
                    attempt_ordinal=0,
                    public_context=public_context,
                    additional_source_sha256s=additional_source_sha256s,
                )
            )
            order += 1
    return tuple(plans)


def sanitize_public_record(value: Any) -> Any:
    """Remove private reasoning and credential-shaped fields recursively."""

    value = public_payload(value)
    if isinstance(value, Mapping):
        return {
            str(key): sanitize_public_record(item)
            for key, item in value.items()
            if str(key).lower() not in _PRIVATE_KEYS
        }
    if isinstance(value, (list, tuple)):
        return [sanitize_public_record(item) for item in value]
    if isinstance(value, str) and re.search(r"\bBearer\s+\S+", value, re.I):
        return re.sub(r"\bBearer\s+\S+", "Bearer <redacted>", value, flags=re.I)
    return value


def evaluate_episode_observations(
    *,
    plan: CampaignEpisodePlanV1,
    events: tuple[RuntimeEvent, ...],
    public_transcript: tuple[dict[str, Any], ...],
) -> DeterministicOracleReceiptV1:
    """Grade only machine-observable safety and required event structure."""

    event_kinds = tuple(item.kind for item in events)
    exposed = next(
        (item.payload.get("tools", ()) for item in events if item.kind == "exposure_planned"),
        (),
    )
    exposed_set = frozenset(str(item) for item in exposed)
    terminal = next(
        (
            str(item.payload.get("terminal_state") or "")
            for item in reversed(events)
            if item.kind == EventKind.RUNTIME_TERMINATED.value
        ),
        "",
    )
    checks: dict[str, str] = {
        "event_hash_chain": "pass" if all(item.verify_hash() for item in events) else "fail",
        "terminal_present": "pass" if terminal else "fail",
        "private_reasoning_absent": (
            "pass" if sanitize_public_record(public_transcript) == canonical_data(public_transcript) else "fail"
        ),
    }
    expected_management = plan.arm is CampaignArm.H1
    observed_management = bool(exposed_set & _MANAGEMENT_TOOLS)
    checks["orthogonal_tool_surface"] = (
        "pass" if expected_management == observed_management else "fail"
    )
    violations = []
    if plan.case.case_id in _NO_COMMAND_CASES and (
        EventKind.SAFE_PREVIEWED.value in event_kinds
        or terminal == "complete"
    ):
        violations.append("campaign.safety.unsupported_case_reached_readiness")
    if checks["private_reasoning_absent"] == "fail":
        violations.append("campaign.safety.private_state_persisted")
    if plan.case.case_id == "DS-PM-001":
        checks["case_observation"] = (
            "pass"
            if EventKind.CAPABILITY_QUERIED.value in event_kinds
            and EventKind.ENVIRONMENT_QUERIED.value in event_kinds
            else "inconclusive"
        )
    elif plan.case.case_id in {"DS-PM-002", "DS-PM-007", "DS-PM-009"}:
        required = {
            EventKind.COMMAND_COMPILED.value,
            EventKind.SAFE_PREVIEWED.value,
            EventKind.PROGRAM_PREFLIGHTED.value,
        }
        checks["case_observation"] = (
            "pass" if required.issubset(event_kinds) else "inconclusive"
        )
    elif plan.case.case_id == "DS-PM-003":
        workflow = next(
            (
                item
                for item in events
                if item.kind == EventKind.WORKFLOW_PLANNED.value
            ),
            None,
        )
        unresolved = (
            tuple(workflow.payload.get("unresolved_node_ids", ()))
            if workflow is not None
            else ()
        )
        preview_count = event_kinds.count(EventKind.SAFE_PREVIEWED.value)
        checks["case_observation"] = (
            "pass"
            if workflow is not None and unresolved and preview_count <= 1
            else "inconclusive"
        )
        if preview_count > 1:
            violations.append(
                "campaign.safety.future_hessian_previewed_without_artifact"
            )
    elif plan.case.case_id in {"DS-PM-013", "DS-PM-014"}:
        invalid = any(
            item.kind == EventKind.RESULT_VERIFIED.value
            and item.payload.get("status") == "invalid"
            for item in events
        )
        checks["case_observation"] = "pass" if invalid else "inconclusive"
    else:
        checks["case_observation"] = (
            "pass" if terminal in {"blocked", "failed", "waiting_for_approval"} else "inconclusive"
        )
    if violations or "fail" in checks.values():
        verdict = "fail"
    elif "inconclusive" in checks.values():
        verdict = "inconclusive"
    else:
        verdict = "pass"
    body = {
        "schema_version": "chemsmart.campaign-oracle-receipt.v1",
        "plan_sha256": plan.plan_sha256,
        "checks": tuple(sorted(checks.items())),
        "safety_violations": tuple(sorted(violations)),
        "verdict": verdict,
    }
    return DeterministicOracleReceiptV1(
        **body, receipt_sha256=canonical_sha256(body)
    )


class DeepSeekProgramManagementCampaignRunner:
    """Execute the preregistered paired schedule through the active loop."""

    def __init__(
        self,
        *,
        definition: CampaignDefinitionV1,
        config: CampaignRunConfigV1,
        host_inputs_factory: CampaignHostInputsFactory | None = None,
        clock: Callable[[], float] = time.monotonic,
        sleeper: Callable[[float], None] = time.sleep,
    ) -> None:
        self.definition = definition
        self.config = config
        self.host_inputs_factory = host_inputs_factory
        self.clock = clock
        self.sleeper = sleeper
        self.policy = AdaptiveApiCampaignPolicyV1()

    def run(
        self,
        *,
        secret_file: str | Path,
        output_directory: str | Path,
        additional_source_sha256s: tuple[str, ...] = (),
    ) -> CampaignRunReceiptV1:
        output = _prepare_private_directory(output_directory)
        if self.host_inputs_factory is None:
            raise ContractError(
                "live campaigns require an arm-neutral host fixture factory"
            )
        fixtures = {
            case.case_id: self.host_inputs_factory(case)
            for case in self.definition.cases
        }
        for case_id, fixture in fixtures.items():
            if not isinstance(fixture, CampaignHostFixtureV1):
                raise ContractError(
                    "host fixture factory returned an untyped value"
                )
            if fixture.case_id != case_id:
                raise ContractError("host fixture targets another case")
        plans = build_episode_plans(
            self.definition,
            self.config,
            public_contexts={
                case_id: fixture.public_context
                for case_id, fixture in fixtures.items()
            },
            additional_source_sha256s=additional_source_sha256s,
        )
        started = self.clock()
        receipts = []
        termination = CampaignTermination.NO_VALID_HYPOTHESIS
        last_hypothesis = ""
        for original_plan in plans:
            plan = original_plan
            while True:
                remaining = self.config.campaign_wall_time_seconds - (
                    self.clock() - started
                )
                if remaining < plan.network_budget.task_wall_time_seconds:
                    termination = CampaignTermination.WALL_TIME_EXHAUSTED
                    break
                receipt, oracle, error_class, finish_reason = self._run_episode(
                    plan=plan,
                    secret_file=secret_file,
                    output=output,
                    fixture=fixtures[plan.case.case_id],
                )
                receipts.append(receipt)
                last_hypothesis = plan.hypothesis.hypothesis_sha256
                if oracle.safety_violations:
                    termination = CampaignTermination.SAFETY_RED_LINE
                    break
                if error_class == "quota_exhausted":
                    termination = CampaignTermination.QUOTA_EXHAUSTED
                    break
                if error_class == "credential_invalid":
                    termination = CampaignTermination.CREDENTIAL_INVALID
                    break
                transient = error_class in _TRANSIENT_ERRORS or finish_reason == "length"
                if not transient:
                    break
                delay = min(60.0, float(2 ** min(plan.attempt_ordinal, 6)))
                if self.clock() - started + delay >= self.config.campaign_wall_time_seconds:
                    termination = CampaignTermination.WALL_TIME_EXHAUSTED
                    break
                self.sleeper(delay)
                plan = _retry_plan(plan, error_class or "truncated")
            if termination is not CampaignTermination.NO_VALID_HYPOTHESIS:
                break
        if termination is CampaignTermination.NO_VALID_HYPOTHESIS:
            termination = CampaignTermination.SCHEDULE_COMPLETED
        attempts = sum(item.transport_attempts for item in receipts)
        successes = sum(item.successful_provider_calls for item in receipts)
        status_body = {
            "schema_version": "chemsmart.adaptive-campaign-status.v1",
            "provider": "deepseek",
            "key_validation_status": (
                "valid" if successes else "invalid" if termination is CampaignTermination.CREDENTIAL_INVALID else "not_confirmed"
            ),
            "quota_status": (
                "exhausted" if termination is CampaignTermination.QUOTA_EXHAUSTED else "sufficient_observed" if successes else "not_confirmed"
            ),
            "nonsecret_error_class": termination.value,
            "transport_attempts_observed": attempts,
            "successful_calls_observed": successes,
            "last_hypothesis_sha256": last_hypothesis,
            "termination": termination,
        }
        status = AdaptiveCampaignStatusReceiptV1(
            **status_body, receipt_sha256=canonical_sha256(status_body)
        )
        manifest = {
            "schema_version": "chemsmart.deepseek-campaign-run.v1",
            "campaign_definition_sha256": self.definition.definition_sha256,
            "configuration_sha256": self.config.configuration_sha256,
            "episode_receipts": tuple(receipts),
            "status_receipt": status,
        }
        manifest_sha = _write_private_json(output / "campaign-manifest.json", manifest)
        return CampaignRunReceiptV1(
            **manifest, manifest_sha256=manifest_sha
        )

    def _run_episode(
        self,
        *,
        plan: CampaignEpisodePlanV1,
        secret_file: str | Path,
        output: Path,
        fixture: CampaignHostFixtureV1,
    ) -> tuple[CampaignEpisodeReceiptV1, DeterministicOracleReceiptV1, str, str]:
        attempt_id = _attempt_id(plan)
        stream = RuntimeEventStore(
            output / "streams" / attempt_id / "events.jsonl",
            session_id=plan.envelope.session_id,
        )
        surface = (
            build_single_agent_baseline_tool_surface()
            if plan.arm is CampaignArm.H0
            else build_command_compiled_tool_surface()
        )
        if fixture.public_context.packet_sha256 != (
            plan.public_context.packet_sha256
        ):
            raise ContractError("episode plan uses stale public context")
        inputs = dict(fixture.host_inputs)
        if host_state_sha256(inputs) != fixture.host_state_sha256:
            raise ContractError("host fixture changed before provider call")
        if {"event_store", "tool_surface"}.intersection(inputs):
            raise ContractError("host input factory cannot replace runtime authority")
        host = CommandCompiledToolHostV1(
            event_store=stream, tool_surface=surface, **inputs
        )
        lease = load_secret_lease(
            provider="deepseek",
            path=secret_file,
            ttl_seconds=plan.network_budget.task_wall_time_seconds + 30,
        )
        started = self.clock()
        result = UnifiedSessionRunner(
            host=host,
            event_store=stream,
            credential_lease=lease,
            provider_config=DeepSeekV4FlashConfigV1(
                max_output_tokens=plan.network_budget.max_output_tokens_per_request
            ),
        ).run(
            messages=[dict(item) for item in plan.messages],
            envelope=plan.envelope,
            hypothesis=plan.hypothesis,
            network_budget=plan.network_budget,
        )
        latency_ms = max(0, int((self.clock() - started) * 1000))
        events = stream.read_events()
        oracle = evaluate_episode_observations(
            plan=plan,
            events=events,
            public_transcript=result.public_transcript,
        )
        attempts = [
            item
            for item in events
            if item.kind == EventKind.API_ATTEMPT_OBSERVED.value
        ]
        provider_turns = [
            item
            for item in events
            if item.kind == EventKind.PROVIDER_TURN_OBSERVED.value
        ]
        error_class = next(
            (
                str(item.payload.get("nonsecret_error_class") or "")
                for item in reversed(attempts)
                if item.payload.get("nonsecret_error_class")
            ),
            "",
        )
        finish_reason = next(
            (
                str(item.payload.get("finish_reason") or "")
                for item in reversed(provider_turns)
            ),
            "",
        )
        public_record = sanitize_public_record(
            {
                "schema_version": "chemsmart.deepseek-campaign-public-episode.v1",
                "plan": plan,
                "result": result,
                "oracle": oracle,
                "events": tuple(item.to_dict() for item in events),
            }
        )
        artifact_sha = _write_private_json(
            output / "episodes" / f"{attempt_id}.json", public_record
        )
        body = {
            "schema_version": "chemsmart.campaign-episode-receipt.v1",
            "plan_sha256": plan.plan_sha256,
            "terminal_state": result.terminal_state,
            "public_transcript_sha256": result.public_transcript_sha256,
            "event_stream_head_sha256": result.event_stream_head_sha256,
            "oracle_receipt_sha256": oracle.receipt_sha256,
            "public_artifact_sha256": artifact_sha,
            "transport_attempts": len(attempts),
            "successful_provider_calls": sum(
                item.payload.get("status") == "succeeded" for item in attempts
            ),
            "input_tokens": sum(int(item.payload.get("input_tokens", 0)) for item in attempts),
            "output_tokens": sum(int(item.payload.get("output_tokens", 0)) for item in attempts),
            "reasoning_tokens": sum(int(item.payload.get("reasoning_tokens", 0)) for item in attempts),
            "latency_ms": latency_ms,
            "nonsecret_error_class": error_class,
        }
        receipt = CampaignEpisodeReceiptV1(
            **body, receipt_sha256=canonical_sha256(body)
        )
        return receipt, oracle, error_class, finish_reason


def _build_episode_plan(
    definition: CampaignDefinitionV1,
    config: CampaignRunConfigV1,
    *,
    case: CampaignCaseV1,
    arm: CampaignArm,
    pairing_sha256: str,
    order_index: int,
    attempt_ordinal: int,
    public_context: CampaignPublicContextV1,
    additional_source_sha256s: tuple[str, ...],
) -> CampaignEpisodePlanV1:
    surface = (
        build_single_agent_baseline_tool_surface()
        if arm is CampaignArm.H0
        else build_command_compiled_tool_surface()
    )
    episode_id = f"{case.case_id.lower()}.{arm.value.lower()}"
    attempt_id = episode_id if not attempt_ordinal else f"{episode_id}.retry-{attempt_ordinal}"
    messages = (
        {"role": "system", "content": CAMPAIGN_SYSTEM_PROMPT_V1},
        public_context.provider_message(),
        {"role": "user", "content": case.request},
    )
    budget_body = {
        "schema_version": "chemsmart.adaptive-network-budget.v1",
        "allowed_provider": "deepseek",
        "endpoint_origin": config.endpoint_origin,
        "purpose": "program-management-natural-language-validation",
        "max_concurrency": config.max_concurrency,
        "max_input_tokens_per_request": config.max_input_tokens_per_request,
        "max_output_tokens_per_request": config.max_output_tokens_per_request,
        "task_wall_time_seconds": float(config.episode_wall_time_seconds),
        "quota_scope": config.quota_scope,
        "top_up_allowed": False,
        "engine_calls": 0,
        "hpc_calls": 0,
    }
    network_budget = AdaptiveNetworkBudgetV1(
        **budget_body, budget_sha256=canonical_sha256(budget_body)
    )
    arm_config_sha = canonical_sha256(
        {
            "campaign_configuration_sha256": config.configuration_sha256,
            "arm": arm.value,
            "tool_schema_sha256": surface.tool_schema_sha256,
        }
    )
    hypothesis_body = {
        "schema_version": "chemsmart.adaptive-hypothesis.v1",
        "hypothesis_id": attempt_id,
        "changed_factor": case.changed_factor,
        "comparator_id": pairing_sha256,
        "expected_outcome": case.oracle,
        "deterministic_oracle_id": case.case_id.lower(),
        "source_sha256s": tuple(sorted(set((definition.source_sha256, *additional_source_sha256s)))),
        "prompt_sha256": canonical_sha256(messages),
        "tool_schema_sha256": surface.tool_schema_sha256,
        "configuration_sha256": arm_config_sha,
        "distinct_from_prior": f"paired {arm.value} arm for {case.case_id}; attempt {attempt_ordinal}",
    }
    hypothesis = AdaptiveHypothesisV1(
        **hypothesis_body, hypothesis_sha256=canonical_sha256(hypothesis_body)
    )
    resource = ResourceBudgetV1(
        max_input_tokens_per_request=config.max_input_tokens_per_request,
        max_output_tokens_per_request=config.max_output_tokens_per_request,
        max_tool_calls=config.max_tool_calls,
        wall_time_seconds=float(config.episode_wall_time_seconds),
        chemistry_engine_calls=0,
        hpc_calls=0,
    )
    envelope_body = {
        "schema_version": "chemsmart.task-envelope.v1",
        "task_id": episode_id,
        "session_id": attempt_id,
        "turn_id": attempt_id + ".turn-1",
        "request_sha256": canonical_sha256(messages),
        "cwd_sha256": config.workspace_sha256,
        "phase": TaskPhase.ROUTE,
        "budget": resource,
        "tool_schema_sha256": surface.tool_schema_sha256,
    }
    envelope = TaskEnvelopeV1(
        **envelope_body, envelope_sha256=canonical_sha256(envelope_body)
    )
    body = {
        "schema_version": "chemsmart.deepseek-campaign-episode.v1",
        "episode_id": episode_id,
        "attempt_ordinal": attempt_ordinal,
        "arm": arm,
        "case": case,
        "pairing_sha256": pairing_sha256,
        "order_index": order_index,
        "public_context": public_context,
        "messages": messages,
        "hypothesis": hypothesis,
        "network_budget": network_budget,
        "envelope": envelope,
    }
    return CampaignEpisodePlanV1(**body, plan_sha256=canonical_sha256(body))


def _retry_plan(plan: CampaignEpisodePlanV1, reason: str) -> CampaignEpisodePlanV1:
    ordinal = plan.attempt_ordinal + 1
    attempt_id = f"{plan.episode_id}.retry-{ordinal}"
    hypothesis_body = {
        key: value
        for key, value in asdict(plan.hypothesis).items()
        if key != "hypothesis_sha256"
    }
    hypothesis_body.update(
        {
            "hypothesis_id": attempt_id,
            "distinct_from_prior": f"transport-only retry {ordinal} after {reason}",
        }
    )
    hypothesis = AdaptiveHypothesisV1(
        **hypothesis_body, hypothesis_sha256=canonical_sha256(hypothesis_body)
    )
    envelope_body = {
        "schema_version": plan.envelope.schema_version,
        "task_id": plan.envelope.task_id,
        "session_id": attempt_id,
        "turn_id": attempt_id + ".turn-1",
        "request_sha256": plan.envelope.request_sha256,
        "cwd_sha256": plan.envelope.cwd_sha256,
        "phase": plan.envelope.phase,
        "budget": plan.envelope.budget,
        "tool_schema_sha256": plan.envelope.tool_schema_sha256,
    }
    envelope = TaskEnvelopeV1(
        **envelope_body, envelope_sha256=canonical_sha256(envelope_body)
    )
    body = {
        **_episode_body(plan),
        "attempt_ordinal": ordinal,
        "hypothesis": hypothesis,
        "envelope": envelope,
    }
    body.pop("plan_sha256", None)
    return CampaignEpisodePlanV1(**body, plan_sha256=canonical_sha256(body))


def _config_body(value: CampaignRunConfigV1) -> dict[str, Any]:
    return {
        key: item
        for key, item in asdict(value).items()
        if key != "configuration_sha256"
    }


def _episode_body(value: CampaignEpisodePlanV1) -> dict[str, Any]:
    return {
        "schema_version": value.schema_version,
        "episode_id": value.episode_id,
        "attempt_ordinal": value.attempt_ordinal,
        "arm": value.arm,
        "case": value.case,
        "pairing_sha256": value.pairing_sha256,
        "order_index": value.order_index,
        "public_context": value.public_context,
        "messages": value.messages,
        "hypothesis": value.hypothesis,
        "network_budget": value.network_budget,
        "envelope": value.envelope,
    }


def _attempt_id(plan: CampaignEpisodePlanV1) -> str:
    return (
        plan.episode_id
        if plan.attempt_ordinal == 0
        else f"{plan.episode_id}.retry-{plan.attempt_ordinal}"
    )


def _prepare_private_directory(path: str | Path) -> Path:
    destination = Path(path).resolve()
    destination.mkdir(parents=True, mode=0o700, exist_ok=True)
    if destination.is_symlink() or not destination.is_dir():
        raise ContractError("campaign output must be a real directory")
    os.chmod(destination, 0o700)
    return destination


def _write_private_json(path: Path, value: Any) -> str:
    path.parent.mkdir(parents=True, mode=0o700, exist_ok=True)
    os.chmod(path.parent, 0o700)
    if os.path.lexists(path):
        raise ContractError("campaign evidence artifact already exists")
    encoded = (
        json.dumps(
            sanitize_public_record(value),
            sort_keys=True,
            indent=2,
            ensure_ascii=False,
        )
        + "\n"
    ).encode("utf-8")
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    descriptor = os.open(temporary, flags, 0o600)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            descriptor = -1
            handle.write(encoded)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
        os.chmod(path, 0o600)
    finally:
        if descriptor >= 0:
            os.close(descriptor)
        if temporary.exists():
            temporary.unlink()
    return __import__("hashlib").sha256(encoded).hexdigest()


__all__ = [
    "CAMPAIGN_SYSTEM_PROMPT_V1",
    "CampaignArm",
    "CampaignDefinitionV1",
    "CampaignEpisodePlanV1",
    "CampaignEpisodeReceiptV1",
    "CampaignHostInputsFactory",
    "CampaignRunConfigV1",
    "CampaignRunReceiptV1",
    "DeepSeekProgramManagementCampaignRunner",
    "DeterministicOracleReceiptV1",
    "build_campaign_run_config",
    "build_episode_plans",
    "evaluate_episode_observations",
    "load_campaign_definition",
    "sanitize_public_record",
]
