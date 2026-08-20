"""Alibaba Token Plan streaming reasoning adapter."""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass, field
import json
import time
from typing import Any, Callable, Mapping
from urllib.parse import urlsplit

from chemsmart.agent._contracts import ContractError
from chemsmart.agent.provider_config import (
    ALIBABA_TOKEN_PLAN_CONTEXT_TOKENS,
    ALIBABA_TOKEN_PLAN_ENDPOINT,
    ALIBABA_TOKEN_PLAN_MAX_OUTPUT_TOKENS,
    ALIBABA_TOKEN_PLAN_MODEL,
    ALIBABA_TOKEN_PLAN_PROVIDER,
)
from chemsmart.agent.runtime.deepseek import (
    DeepSeekHttpsTransport,
    DeepSeekProtocolError,
    DeepSeekTransportError,
    DeepSeekV4ToolSession,
    ProviderCapabilitiesV1,
    _require_explicit_model_id,
    _validate_token_limits,
)
from chemsmart.agent.runtime.transport import (
    ProviderTurnDeadlinesV1,
    iter_bounded_response_lines,
    open_bounded_https_response,
)


@dataclass(frozen=True)
class AlibabaTokenPlanConfigV1:
    """One model selected from the Alibaba Token Plan catalog."""

    model: str
    context_tokens: int
    max_output_tokens: int
    provider: str = ALIBABA_TOKEN_PLAN_PROVIDER
    endpoint: str = ALIBABA_TOKEN_PLAN_ENDPOINT
    reasoning_effort: str = "xhigh"
    preserve_thinking: bool = True
    sdk_max_retries: int = 0
    turn_deadlines: ProviderTurnDeadlinesV1 = field(
        default_factory=ProviderTurnDeadlinesV1
    )

    def __post_init__(self) -> None:
        if self.provider != ALIBABA_TOKEN_PLAN_PROVIDER:
            raise ContractError(
                "Alibaba Token Plan provider identity is immutable"
            )
        _require_explicit_model_id(self.model, provider="Alibaba Token Plan")
        _validate_token_limits(
            context_tokens=self.context_tokens,
            max_output_tokens=self.max_output_tokens,
        )
        if not _is_token_plan_endpoint(self.endpoint):
            raise ContractError(
                "Alibaba adapter requires the Token Plan endpoint"
            )
        if self.reasoning_effort not in {"high", "max", "xhigh"}:
            raise ContractError(
                "Alibaba Token Plan reasoning effort must be high, max or "
                "xhigh"
            )
        if not self.preserve_thinking:
            raise ContractError(
                "Alibaba tool continuation must preserve thinking"
            )
        if self.sdk_max_retries != 0:
            raise ContractError(
                "provider retries require a separately authorized attempt"
            )


@dataclass(frozen=True)
class Qwen38MaxConfigV1(AlibabaTokenPlanConfigV1):
    """Backward-compatible production Qwen 3.8 Max configuration."""

    model: str = ALIBABA_TOKEN_PLAN_MODEL
    context_tokens: int = ALIBABA_TOKEN_PLAN_CONTEXT_TOKENS
    max_output_tokens: int = ALIBABA_TOKEN_PLAN_MAX_OUTPUT_TOKENS

    def __post_init__(self) -> None:
        super().__post_init__()
        if self.model != ALIBABA_TOKEN_PLAN_MODEL:
            raise ContractError("Qwen adapter requires production qwen3.8-max")
        if self.reasoning_effort != "xhigh":
            raise ContractError("Qwen 3.8 Max reasoning effort must be xhigh")


class AlibabaTokenPlanHttpsTransport(DeepSeekHttpsTransport):
    """Streaming Token Plan transport with sanitized failure classes.

    Only the endpoint registration, the Token Plan lease shape, the SSE
    Accept header, and the streamed-body assembly differ from the base
    transport; the deadline discipline and failure ladder are inherited.
    """

    _ENDPOINT_ERROR = "Alibaba transport requires the Token Plan endpoint"
    _KEY_ERROR = "Alibaba transport requires a Token Plan lease"
    _TIMEOUT_ERROR = "Alibaba timeout must be positive"
    _CLOSED_ERROR = "Alibaba credential lease is closed"
    _ACCEPT = "text/event-stream"

    @staticmethod
    def _endpoint_is_registered(endpoint: str) -> bool:
        return _is_token_plan_endpoint(endpoint)

    @staticmethod
    def _api_key_is_leased(api_key: str) -> bool:
        return api_key.startswith("sk-sp-")

    def __init__(
        self,
        *,
        api_key: str,
        endpoint: str = ALIBABA_TOKEN_PLAN_ENDPOINT,
        timeout_seconds: float | None = None,
        turn_deadlines: ProviderTurnDeadlinesV1 | None = None,
        clock: Callable[[], float] = time.monotonic,
    ) -> None:
        super().__init__(
            api_key=api_key,
            endpoint=endpoint,
            timeout_seconds=timeout_seconds,
            turn_deadlines=turn_deadlines,
            clock=clock,
        )

    def _open_response(self, request, *, deadline):
        return open_bounded_https_response(request, deadline=deadline)

    def _read_response(self, response, *, deadline, payload):
        return _assemble_sse_response(
            response,
            expected_model=str(payload.get("model") or ""),
            deadline=deadline,
        )


class AlibabaTokenPlanToolSession(DeepSeekV4ToolSession):
    """Provider-native continuation with private reasoning held in RAM."""

    def __init__(
        self,
        *,
        transport,
        messages: list[dict[str, Any]],
        config: AlibabaTokenPlanConfigV1,
    ) -> None:
        super().__init__(
            transport=transport,
            messages=messages,
            config=config,
        )

    @property
    def capabilities(self) -> ProviderCapabilitiesV1:
        return ProviderCapabilitiesV1(
            provider=self.config.provider,
            model=self.config.model,
            context_tokens=self.config.context_tokens,
            max_output_tokens=self.config.max_output_tokens,
        )

    def request_payload(
        self, *, tools: list[dict[str, Any]] | None = None
    ) -> dict[str, Any]:
        payload: dict[str, Any] = {
            "model": self.config.model,
            "messages": deepcopy(self._history),
            "stream": True,
            "stream_options": {"include_usage": True},
            "max_tokens": self.config.max_output_tokens,
            "reasoning_effort": self.config.reasoning_effort,
            "preserve_thinking": self.config.preserve_thinking,
        }
        if tools:
            payload["tools"] = deepcopy(tools)
        return payload


class Qwen38MaxToolSession(AlibabaTokenPlanToolSession):
    """Backward-compatible production Qwen 3.8 Max session."""

    def __init__(
        self,
        *,
        transport,
        messages: list[dict[str, Any]],
        config: Qwen38MaxConfigV1 | None = None,
    ) -> None:
        super().__init__(
            transport=transport,
            messages=messages,
            config=config or Qwen38MaxConfigV1(),
        )


def _assemble_sse_response(
    response,
    *,
    expected_model: str,
    deadline: ProviderTurnDeadline | None = None,
) -> dict[str, Any]:
    response_id = ""
    observed_model = ""
    finish_reason = ""
    usage: dict[str, Any] = {}
    role = "assistant"
    content_parts: list[str] = []
    reasoning_parts: list[str] = []
    tool_calls: dict[int, dict[str, Any]] = {}
    saw_event = False

    lines = (
        iter_bounded_response_lines(response, deadline=deadline)
        if deadline is not None
        else response
    )
    for raw_line in lines:
        if deadline is not None:
            deadline.check()
        line = raw_line.decode("utf-8").strip()
        if not line or line.startswith(":"):
            continue
        if not line.startswith("data:"):
            raise DeepSeekProtocolError("Alibaba stream contains a non-SSE event")
        data = line[5:].strip()
        if data == "[DONE]":
            if deadline is not None:
                deadline.event_observed()
            break
        chunk = json.loads(data)
        if not isinstance(chunk, Mapping):
            raise DeepSeekProtocolError("Alibaba stream chunk must be an object")
        if deadline is not None:
            deadline.event_observed()
        saw_event = True
        response_id = response_id or str(chunk.get("id") or "")
        observed_model = observed_model or str(chunk.get("model") or "")
        if isinstance(chunk.get("usage"), Mapping):
            usage = deepcopy(dict(chunk["usage"]))
        choices = chunk.get("choices") or []
        if not isinstance(choices, list):
            raise DeepSeekProtocolError("Alibaba stream choices must be a list")
        for choice in choices:
            if not isinstance(choice, Mapping):
                raise DeepSeekProtocolError("Alibaba stream choice must be an object")
            finish_reason = str(choice.get("finish_reason") or finish_reason)
            delta = choice.get("delta") or {}
            if not isinstance(delta, Mapping):
                raise DeepSeekProtocolError("Alibaba stream delta must be an object")
            role = str(delta.get("role") or role)
            if delta.get("content"):
                content_parts.append(str(delta["content"]))
            if delta.get("reasoning_content"):
                reasoning_parts.append(str(delta["reasoning_content"]))
            _merge_tool_call_deltas(tool_calls, delta.get("tool_calls") or [])

    if not saw_event:
        raise DeepSeekProtocolError("Alibaba stream returned no events")
    if observed_model != expected_model:
        raise DeepSeekProtocolError(
            "Alibaba stream did not confirm the requested model"
        )
    message: dict[str, Any] = {
        "role": role,
        "content": "".join(content_parts),
    }
    if reasoning_parts or tool_calls:
        message["reasoning_content"] = "".join(reasoning_parts)
    if tool_calls:
        message["tool_calls"] = [tool_calls[index] for index in sorted(tool_calls)]
    return {
        "id": response_id,
        "model": observed_model,
        "choices": [{"finish_reason": finish_reason, "message": message}],
        "usage": usage,
    }


def _merge_tool_call_deltas(
    assembled: dict[int, dict[str, Any]], raw_calls: Any
) -> None:
    if not isinstance(raw_calls, list):
        raise DeepSeekProtocolError("Alibaba tool-call delta must be a list")
    for raw in raw_calls:
        if not isinstance(raw, Mapping):
            raise DeepSeekProtocolError("Alibaba tool-call fragment must be an object")
        raw_index = raw.get("index", 0)
        if (
            isinstance(raw_index, bool)
            or not isinstance(raw_index, int)
            or raw_index < 0
        ):
            raise DeepSeekProtocolError(
                "Alibaba tool-call index must be a non-negative integer",
                failing_field=(
                    "choices[*].delta.tool_calls[*].index"
                ),
            )
        index = raw_index
        target = assembled.setdefault(
            index,
            {
                "id": "",
                "type": "function",
                "function": {"name": "", "arguments": ""},
            },
        )
        if raw.get("id"):
            target["id"] += str(raw["id"])
        if raw.get("type"):
            target["type"] = str(raw["type"])
        function = raw.get("function") or {}
        if not isinstance(function, Mapping):
            raise DeepSeekProtocolError("Alibaba tool function must be an object")
        if function.get("name"):
            target["function"]["name"] += str(function["name"])
        if function.get("arguments"):
            target["function"]["arguments"] += str(function["arguments"])


def _is_token_plan_endpoint(endpoint: str) -> bool:
    try:
        parsed = urlsplit(str(endpoint).strip())
    except ValueError:
        return False
    return (
        parsed.scheme.lower() == "https"
        and parsed.hostname == "token-plan.ap-southeast-1.maas.aliyuncs.com"
        and parsed.port in {None, 443}
        and parsed.path.rstrip("/") == "/compatible-mode/v1"
        and not parsed.username
        and not parsed.password
        and not parsed.query
        and not parsed.fragment
    )


__all__ = [
    "AlibabaTokenPlanConfigV1",
    "AlibabaTokenPlanHttpsTransport",
    "AlibabaTokenPlanToolSession",
    "Qwen38MaxConfigV1",
    "Qwen38MaxToolSession",
]
