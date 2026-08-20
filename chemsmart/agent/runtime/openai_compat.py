"""The generic OpenAI chat-completions adapter (provider ``openai``).

The DeepSeek session already speaks plain chat completions with tool
continuation; this adapter subclasses it, binds the transport to the
official OpenAI endpoint, and strips the payload of every vendor
extension: no ``thinking`` key, and ``reasoning_effort`` only when the
profile set one (reasoning-series models accept it; chat models reject
unknown keys).
"""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass, field
from typing import Any

from chemsmart.agent._contracts import ContractError
from chemsmart.agent.runtime.deepseek import (
    DeepSeekHttpsTransport,
    DeepSeekV4ToolSession,
    _require_explicit_model_id,
    _validate_token_limits,
)
from chemsmart.agent.runtime.transport import ProviderTurnDeadlinesV1

OPENAI_OFFICIAL_ENDPOINT = "https://api.openai.com/v1"

_EFFORT_VOCABULARY = {"", "low", "medium", "high"}


@dataclass(frozen=True)
class OpenAICompatibleConfigV1:
    """OpenAI adapter configuration bound to an explicit model profile."""

    model: str
    context_tokens: int
    max_output_tokens: int
    provider: str = "openai"
    endpoint: str = OPENAI_OFFICIAL_ENDPOINT
    reasoning_effort: str = ""
    sdk_max_retries: int = 0
    turn_deadlines: ProviderTurnDeadlinesV1 = field(
        default_factory=ProviderTurnDeadlinesV1
    )

    def __post_init__(self) -> None:
        if self.provider != "openai":
            raise ContractError("OpenAI provider identity is immutable")
        _require_explicit_model_id(self.model, provider="OpenAI")
        _validate_token_limits(
            context_tokens=self.context_tokens,
            max_output_tokens=self.max_output_tokens,
        )
        if self.endpoint.rstrip("/") != OPENAI_OFFICIAL_ENDPOINT:
            raise ContractError(
                "the OpenAI adapter requires the official endpoint"
            )
        if self.reasoning_effort not in _EFFORT_VOCABULARY:
            raise ContractError(
                "OpenAI reasoning effort must be low, medium, high, or "
                "omitted"
            )
        if self.sdk_max_retries != 0:
            raise ContractError(
                "provider retries require a separately authorized attempt"
            )


class OpenAIHttpsTransport(DeepSeekHttpsTransport):
    """The shared bounded HTTPS transport, bound to api.openai.com."""

    _ENDPOINT_ERROR = "the OpenAI transport requires the official endpoint"

    @staticmethod
    def _endpoint_is_registered(endpoint: str) -> bool:
        return endpoint.rstrip("/") == OPENAI_OFFICIAL_ENDPOINT


class OpenAICompatibleToolSession(DeepSeekV4ToolSession):
    """Plain chat-completions continuation without vendor extensions."""

    def request_payload(
        self, *, tools: list[dict[str, Any]] | None = None
    ) -> dict[str, Any]:
        payload: dict[str, Any] = {
            "model": self.config.model,
            "messages": deepcopy(self._history),
            "max_tokens": self.config.max_output_tokens,
        }
        if self.config.reasoning_effort:
            payload["reasoning_effort"] = self.config.reasoning_effort
        if tools:
            payload["tools"] = deepcopy(tools)
        return payload


__all__ = [
    "OPENAI_OFFICIAL_ENDPOINT",
    "OpenAICompatibleConfigV1",
    "OpenAICompatibleToolSession",
    "OpenAIHttpsTransport",
]
