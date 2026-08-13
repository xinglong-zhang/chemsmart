"""Host-neutral provenance and bounds for one provider request loop."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from chemsmart.agent._contracts import (
    ContractError,
    canonical_sha256,
    require_sha256,
)


@dataclass(frozen=True)
class ProviderNetworkBudgetV1:
    """Provider, token, concurrency, and elapsed-time bounds."""

    schema_version: str
    allowed_provider: str
    endpoint_origin: str
    max_concurrency: int
    max_input_tokens_per_request: int
    max_output_tokens_per_request: int
    task_wall_time_seconds: float
    budget_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.provider-network-budget.v1":
            raise ContractError("unsupported provider network budget schema")
        if not self.allowed_provider.strip() or not self.endpoint_origin.strip():
            raise ContractError("provider budget identity must not be empty")
        if not 1 <= self.max_concurrency <= 4:
            raise ContractError("provider concurrency must be within [1, 4]")
        if min(
            self.max_input_tokens_per_request,
            self.max_output_tokens_per_request,
        ) < 1 or self.task_wall_time_seconds <= 0:
            raise ContractError("provider token and wall-time bounds are required")
        body = {
            "schema_version": self.schema_version,
            "allowed_provider": self.allowed_provider,
            "endpoint_origin": self.endpoint_origin,
            "max_concurrency": self.max_concurrency,
            "max_input_tokens_per_request": self.max_input_tokens_per_request,
            "max_output_tokens_per_request": self.max_output_tokens_per_request,
            "task_wall_time_seconds": self.task_wall_time_seconds,
        }
        if self.budget_sha256 != canonical_sha256(body):
            raise ContractError("provider network budget digest mismatch")


@dataclass(frozen=True)
class RequestContextProvenanceV1:
    """Digests that identify the exact request without encoding a benchmark."""

    schema_version: str
    task_spec_sha256: str
    prompt_sha256: str
    tool_schema_sha256: str
    configuration_sha256: str
    provider_budget_sha256: str
    provenance_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.request-context-provenance.v1":
            raise ContractError("unsupported request context provenance schema")
        for name, digest in (
            ("task_spec_sha256", self.task_spec_sha256),
            ("prompt_sha256", self.prompt_sha256),
            ("tool_schema_sha256", self.tool_schema_sha256),
            ("configuration_sha256", self.configuration_sha256),
            ("provider_budget_sha256", self.provider_budget_sha256),
        ):
            require_sha256(digest, name)
        body = {
            "schema_version": self.schema_version,
            "task_spec_sha256": self.task_spec_sha256,
            "prompt_sha256": self.prompt_sha256,
            "tool_schema_sha256": self.tool_schema_sha256,
            "configuration_sha256": self.configuration_sha256,
            "provider_budget_sha256": self.provider_budget_sha256,
        }
        if self.provenance_sha256 != canonical_sha256(body):
            raise ContractError("request context provenance digest mismatch")


@dataclass(frozen=True)
class ProviderAttemptReceiptV1:
    """Sanitized observation of one provider transport attempt."""

    schema_version: str
    attempt_id: str
    provider: str
    endpoint_origin: str
    request_context_sha256: str
    provider_budget_sha256: str
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
        if self.schema_version != "chemsmart.provider-attempt-receipt.v1":
            raise ContractError("unsupported provider attempt receipt schema")
        if not self.attempt_id or not self.provider or not self.endpoint_origin:
            raise ContractError("provider attempt identity must not be empty")
        for name, digest in (
            ("request_context_sha256", self.request_context_sha256),
            ("provider_budget_sha256", self.provider_budget_sha256),
            ("request_sha256", self.request_sha256),
        ):
            require_sha256(digest, name)
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
            raise ContractError("invalid provider attempt status")
        if min(
            self.latency_ms,
            self.input_tokens,
            self.output_tokens,
            self.reasoning_tokens,
            self.retry_ordinal,
        ) < 0:
            raise ContractError("provider attempt counters must be non-negative")
        body = {
            key: value
            for key, value in self.__dict__.items()
            if key != "receipt_sha256"
        }
        if self.receipt_sha256 != canonical_sha256(body):
            raise ContractError("provider attempt receipt digest mismatch")


def build_provider_attempt_receipt(**values: Any) -> ProviderAttemptReceiptV1:
    body = {
        "schema_version": "chemsmart.provider-attempt-receipt.v1",
        "attempt_id": str(values["attempt_id"]),
        "provider": str(values["provider"]),
        "endpoint_origin": str(values["endpoint_origin"]),
        "request_context_sha256": str(values["request_context_sha256"]),
        "provider_budget_sha256": str(values["provider_budget_sha256"]),
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
    return ProviderAttemptReceiptV1(
        **body, receipt_sha256=canonical_sha256(body)
    )


def build_request_context_provenance(
    *,
    task_spec_sha256: str,
    prompt_sha256: str,
    tool_schema_sha256: str,
    configuration_sha256: str,
    provider_budget_sha256: str,
) -> RequestContextProvenanceV1:
    """Bind one request to its exact task, prompt, tools, and runtime inputs."""

    body = {
        "schema_version": "chemsmart.request-context-provenance.v1",
        "task_spec_sha256": task_spec_sha256,
        "prompt_sha256": prompt_sha256,
        "tool_schema_sha256": tool_schema_sha256,
        "configuration_sha256": configuration_sha256,
        "provider_budget_sha256": provider_budget_sha256,
    }
    return RequestContextProvenanceV1(
        **body, provenance_sha256=canonical_sha256(body)
    )


__all__ = [
    "ProviderAttemptReceiptV1",
    "ProviderNetworkBudgetV1",
    "RequestContextProvenanceV1",
    "build_provider_attempt_receipt",
    "build_request_context_provenance",
]
