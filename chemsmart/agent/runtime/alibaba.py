"""Alibaba Token Plan streaming reasoning adapter."""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass, field
from http.client import HTTPException
import json
import ssl
import time
from typing import Any, Callable, Mapping
from urllib.error import HTTPError, URLError
from urllib.parse import urlsplit
from urllib.request import Request

from chemsmart.agent._contracts import ContractError
from chemsmart.agent.provider_config import (
    ALIBABA_TOKEN_PLAN_CONTEXT_TOKENS,
    ALIBABA_TOKEN_PLAN_ENDPOINT,
    ALIBABA_TOKEN_PLAN_MAX_OUTPUT_TOKENS,
    ALIBABA_TOKEN_PLAN_MODEL,
    ALIBABA_TOKEN_PLAN_PROVIDER,
)
from chemsmart.agent.runtime.deepseek import (
    DeepSeekProtocolError,
    DeepSeekTransportError,
    DeepSeekV4ToolSession,
    ProviderCapabilitiesV1,
    _deadline_transport_error,
    _require_explicit_model_id,
    _validate_token_limits,
)
from chemsmart.agent.runtime.transport import (
    ProviderDeadlineExceeded,
    ProviderTurnDeadline,
    ProviderTurnDeadlinesV1,
    is_socket_timeout,
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


class AlibabaTokenPlanHttpsTransport:
    """Streaming Token Plan transport with sanitized failure classes."""

    def __init__(
        self,
        *,
        api_key: str,
        endpoint: str = ALIBABA_TOKEN_PLAN_ENDPOINT,
        timeout_seconds: float | None = None,
        turn_deadlines: ProviderTurnDeadlinesV1 | None = None,
        clock: Callable[[], float] = time.monotonic,
    ) -> None:
        if not _is_token_plan_endpoint(endpoint):
            raise ContractError("Alibaba transport requires the Token Plan endpoint")
        if not api_key.startswith("sk-sp-"):
            raise ContractError("Alibaba transport requires a Token Plan lease")
        policy = turn_deadlines or ProviderTurnDeadlinesV1()
        if timeout_seconds is not None:
            value = float(timeout_seconds)
            if value <= 0:
                raise ContractError("Alibaba timeout must be positive")
            policy = ProviderTurnDeadlinesV1(
                connect_seconds=min(policy.connect_seconds, value),
                first_event_seconds=min(policy.first_event_seconds, value),
                inter_event_seconds=min(policy.inter_event_seconds, value),
                absolute_seconds=min(policy.absolute_seconds, value),
            )
        self.endpoint = endpoint.rstrip("/")
        self.turn_deadlines = policy
        self._clock = clock
        self._next_turn_timeout_seconds: float | None = None
        self._last_deadline_record = policy.public_record()
        self._api_key = api_key
        self._closed = False

    def set_timeout_seconds(self, timeout_seconds: float) -> None:
        """Tighten the next Qwen turn to the remaining task allowance."""

        value = float(timeout_seconds)
        if value <= 0:
            raise ContractError("Alibaba timeout must be positive")
        self._next_turn_timeout_seconds = min(
            value, self.turn_deadlines.absolute_seconds
        )

    def public_deadline_record(self) -> dict[str, float]:
        return dict(self._last_deadline_record)

    def __call__(self, payload: dict[str, Any]) -> Mapping[str, Any]:
        if self._closed or not self._api_key:
            raise ContractError("Alibaba credential lease is closed")
        encoded = json.dumps(
            payload, separators=(",", ":"), ensure_ascii=False
        ).encode("utf-8")
        request = Request(
            self.endpoint + "/chat/completions",
            data=encoded,
            method="POST",
            headers={
                "Authorization": "Bearer " + self._api_key,
                "Content-Type": "application/json",
                "Accept": "text/event-stream",
                "User-Agent": "chemsmart-agent/1",
            },
        )
        deadline = ProviderTurnDeadline(
            self.turn_deadlines,
            turn_limit_seconds=self._next_turn_timeout_seconds,
            clock=self._clock,
        )
        self._next_turn_timeout_seconds = None
        self._last_deadline_record = deadline.public_record()
        try:
            with open_bounded_https_response(
                request, deadline=deadline
            ) as response:
                if deadline.connected_at is None:
                    # Focused in-memory transports may inject an already-open
                    # response; the production opener marks this itself.
                    deadline.response_acquired()
                return _assemble_sse_response(
                    response,
                    expected_model=str(payload.get("model") or ""),
                    deadline=deadline,
                )
        except ProviderDeadlineExceeded as exc:
            raise _deadline_transport_error(exc) from None
        except HTTPError as exc:
            raise DeepSeekTransportError(
                _http_error_class(exc), http_status=exc.code
            ) from None
        except TimeoutError:
            raise _deadline_transport_error(deadline.timeout_error()) from None
        except URLError as exc:
            reason = getattr(exc, "reason", None)
            if isinstance(reason, BaseException) and is_socket_timeout(reason):
                raise _deadline_transport_error(
                    deadline.timeout_error()
                ) from None
            raise DeepSeekTransportError("transport") from None
        except (HTTPException, ssl.SSLError, ConnectionError, OSError):
            raise DeepSeekTransportError("transport") from None
        except (UnicodeDecodeError, json.JSONDecodeError) as exc:
            raise DeepSeekTransportError("invalid_json") from exc

    def close(self) -> None:
        self._api_key = ""
        self._closed = True


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


def _http_error_class(exc: HTTPError) -> str:
    if exc.code == 401:
        return "credential_invalid"
    if exc.code == 404:
        return "model_unavailable"
    if exc.code == 429:
        # Never read an untrusted error body here: a peer could drip it past
        # the turn cap, and provider text must not enter Runtime evidence.
        return "rate_limited"
    if 500 <= exc.code < 600:
        return "provider_5xx"
    return "http_error"


__all__ = [
    "AlibabaTokenPlanConfigV1",
    "AlibabaTokenPlanHttpsTransport",
    "AlibabaTokenPlanToolSession",
    "Qwen38MaxConfigV1",
    "Qwen38MaxToolSession",
]
