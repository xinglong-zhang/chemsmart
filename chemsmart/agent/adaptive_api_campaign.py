"""Count-unbounded but hypothesis-, quota-, and time-bounded API campaigns."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Any

from chemsmart.agent._contracts import (
    ContractError,
    canonical_sha256,
    require_sha256,
)


class CampaignTermination(str, Enum):
    ACTIVE = "active"
    QUOTA_EXHAUSTED = "quota_exhausted"
    NO_VALID_HYPOTHESIS = "no_valid_hypothesis"
    SCHEDULE_COMPLETED = "schedule_completed"
    CREDENTIAL_INVALID = "credential_invalid"
    SAFETY_RED_LINE = "safety_red_line"
    WALL_TIME_EXHAUSTED = "wall_time_exhausted"


@dataclass(frozen=True)
class AdaptiveApiCampaignPolicyV1:
    schema_version: str = "chemsmart.adaptive-api-campaign-policy.v1"
    transport_attempt_limit: None = None
    quota_source: str = "current_user_account"
    top_up_allowed: bool = False
    unique_hypothesis_required: bool = True
    repeat_prompt_without_measurement: bool = False

    def __post_init__(self) -> None:
        if self.transport_attempt_limit is not None:
            raise ContractError("adaptive campaigns do not use a call-count cap")
        if self.top_up_allowed:
            raise ContractError("campaign policy cannot authorize quota top-up")
        if not self.unique_hypothesis_required:
            raise ContractError("every provider call requires a hypothesis")
        if self.repeat_prompt_without_measurement:
            raise ContractError("unmeasured duplicate calls are not allowed")


@dataclass(frozen=True)
class AdaptiveNetworkBudgetV1:
    schema_version: str
    allowed_provider: str
    endpoint_origin: str
    purpose: str
    max_concurrency: int
    max_input_tokens_per_request: int
    max_output_tokens_per_request: int
    task_wall_time_seconds: float
    quota_scope: str
    top_up_allowed: bool
    engine_calls: int
    hpc_calls: int
    budget_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.adaptive-network-budget.v1":
            raise ContractError("unsupported adaptive network budget schema")
        if not 1 <= self.max_concurrency <= 4:
            raise ContractError("adaptive concurrency must be within [1, 4]")
        if min(
            self.max_input_tokens_per_request,
            self.max_output_tokens_per_request,
        ) < 1 or self.task_wall_time_seconds <= 0:
            raise ContractError("adaptive token and wall-time bounds are required")
        if self.top_up_allowed or self.engine_calls or self.hpc_calls:
            raise ContractError("campaign budget cannot top up or run chemistry")
        body = {
            "schema_version": self.schema_version,
            "allowed_provider": self.allowed_provider,
            "endpoint_origin": self.endpoint_origin,
            "purpose": self.purpose,
            "max_concurrency": self.max_concurrency,
            "max_input_tokens_per_request": self.max_input_tokens_per_request,
            "max_output_tokens_per_request": self.max_output_tokens_per_request,
            "task_wall_time_seconds": self.task_wall_time_seconds,
            "quota_scope": self.quota_scope,
            "top_up_allowed": self.top_up_allowed,
            "engine_calls": self.engine_calls,
            "hpc_calls": self.hpc_calls,
        }
        if self.budget_sha256 != canonical_sha256(body):
            raise ContractError("adaptive network budget digest mismatch")


@dataclass(frozen=True)
class AdaptiveHypothesisV1:
    schema_version: str
    hypothesis_id: str
    changed_factor: str
    comparator_id: str
    expected_outcome: str
    deterministic_oracle_id: str
    source_sha256s: tuple[str, ...]
    prompt_sha256: str
    tool_schema_sha256: str
    configuration_sha256: str
    distinct_from_prior: str
    hypothesis_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.adaptive-hypothesis.v1":
            raise ContractError("unsupported adaptive hypothesis schema")
        required = (
            self.hypothesis_id,
            self.changed_factor,
            self.expected_outcome,
            self.deterministic_oracle_id,
            self.distinct_from_prior,
        )
        if not all(str(item).strip() for item in required):
            raise ContractError("adaptive hypothesis fields must not be empty")
        body = {
            "schema_version": self.schema_version,
            "hypothesis_id": self.hypothesis_id,
            "changed_factor": self.changed_factor,
            "comparator_id": self.comparator_id,
            "expected_outcome": self.expected_outcome,
            "deterministic_oracle_id": self.deterministic_oracle_id,
            "source_sha256s": self.source_sha256s,
            "prompt_sha256": self.prompt_sha256,
            "tool_schema_sha256": self.tool_schema_sha256,
            "configuration_sha256": self.configuration_sha256,
            "distinct_from_prior": self.distinct_from_prior,
        }
        if self.hypothesis_sha256 != canonical_sha256(body):
            raise ContractError("adaptive hypothesis digest mismatch")


@dataclass(frozen=True)
class AdaptiveCampaignStatusReceiptV1:
    schema_version: str
    provider: str
    key_validation_status: str
    quota_status: str
    nonsecret_error_class: str
    transport_attempts_observed: int
    successful_calls_observed: int
    last_hypothesis_sha256: str
    termination: CampaignTermination
    receipt_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.adaptive-campaign-status.v1":
            raise ContractError("unsupported adaptive campaign status schema")
        if min(
            self.transport_attempts_observed,
            self.successful_calls_observed,
        ) < 0:
            raise ContractError("observed call counts must be non-negative")
        if self.successful_calls_observed > self.transport_attempts_observed:
            raise ContractError("successful calls exceed transport attempts")
        body = {
            "schema_version": self.schema_version,
            "provider": self.provider,
            "key_validation_status": self.key_validation_status,
            "quota_status": self.quota_status,
            "nonsecret_error_class": self.nonsecret_error_class,
            "transport_attempts_observed": self.transport_attempts_observed,
            "successful_calls_observed": self.successful_calls_observed,
            "last_hypothesis_sha256": self.last_hypothesis_sha256,
            "termination": self.termination,
        }
        if self.receipt_sha256 != canonical_sha256(body):
            raise ContractError("adaptive campaign status digest mismatch")


@dataclass(frozen=True)
class ApiAttemptReceiptV1:
    """One observable provider attempt bound to hypothesis and budget."""

    schema_version: str
    attempt_id: str
    provider: str
    endpoint_origin: str
    hypothesis_sha256: str
    budget_sha256: str
    request_sha256: str
    response_sha256: str
    status: str
    latency_ms: int
    input_tokens: int
    output_tokens: int
    reasoning_tokens: int
    reported_cost_usd: str
    retry_ordinal: int
    nonsecret_error_class: str
    receipt_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.api-attempt-receipt.v1":
            raise ContractError("unsupported API attempt receipt schema")
        if not self.attempt_id or not self.provider or not self.endpoint_origin:
            raise ContractError("API attempt identity must not be empty")
        for field_name, digest in (
            ("hypothesis_sha256", self.hypothesis_sha256),
            ("budget_sha256", self.budget_sha256),
            ("request_sha256", self.request_sha256),
        ):
            require_sha256(digest, field_name)
        if self.response_sha256:
            require_sha256(self.response_sha256, "response_sha256")
        if self.status not in {
            "succeeded",
            "quota_exhausted",
            "credential_invalid",
            "rate_limited",
            "timeout",
            "transport_failed",
            "protocol_failed",
        }:
            raise ContractError("invalid API attempt status")
        if min(
            self.latency_ms,
            self.input_tokens,
            self.output_tokens,
            self.reasoning_tokens,
            self.retry_ordinal,
        ) < 0:
            raise ContractError("API attempt counters must be non-negative")
        body = {
            "schema_version": self.schema_version,
            "attempt_id": self.attempt_id,
            "provider": self.provider,
            "endpoint_origin": self.endpoint_origin,
            "hypothesis_sha256": self.hypothesis_sha256,
            "budget_sha256": self.budget_sha256,
            "request_sha256": self.request_sha256,
            "response_sha256": self.response_sha256,
            "status": self.status,
            "latency_ms": self.latency_ms,
            "input_tokens": self.input_tokens,
            "output_tokens": self.output_tokens,
            "reasoning_tokens": self.reasoning_tokens,
            "reported_cost_usd": self.reported_cost_usd,
            "retry_ordinal": self.retry_ordinal,
            "nonsecret_error_class": self.nonsecret_error_class,
        }
        if self.receipt_sha256 != canonical_sha256(body):
            raise ContractError("API attempt receipt digest mismatch")


def build_api_attempt_receipt(**values: Any) -> ApiAttemptReceiptV1:
    body = {
        "schema_version": "chemsmart.api-attempt-receipt.v1",
        "attempt_id": str(values["attempt_id"]),
        "provider": str(values["provider"]),
        "endpoint_origin": str(values["endpoint_origin"]),
        "hypothesis_sha256": str(values["hypothesis_sha256"]),
        "budget_sha256": str(values["budget_sha256"]),
        "request_sha256": str(values["request_sha256"]),
        "response_sha256": str(values.get("response_sha256") or ""),
        "status": str(values["status"]),
        "latency_ms": int(values.get("latency_ms", 0)),
        "input_tokens": int(values.get("input_tokens", 0)),
        "output_tokens": int(values.get("output_tokens", 0)),
        "reasoning_tokens": int(values.get("reasoning_tokens", 0)),
        "reported_cost_usd": str(values.get("reported_cost_usd") or ""),
        "retry_ordinal": int(values.get("retry_ordinal", 0)),
        "nonsecret_error_class": str(
            values.get("nonsecret_error_class") or ""
        ),
    }
    return ApiAttemptReceiptV1(
        **body, receipt_sha256=canonical_sha256(body)
    )


__all__ = [
    "AdaptiveApiCampaignPolicyV1",
    "AdaptiveCampaignStatusReceiptV1",
    "AdaptiveHypothesisV1",
    "AdaptiveNetworkBudgetV1",
    "ApiAttemptReceiptV1",
    "CampaignTermination",
    "build_api_attempt_receipt",
]
