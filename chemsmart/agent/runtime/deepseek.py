"""DeepSeek V4 Flash thinking/tool continuation with private ephemeral state.

The transport is injected so this module neither owns credentials nor adds an
SDK dependency to ChemSmart v3.1.4.  It preserves the provider's complete
``reasoning_content`` only inside an uninterrupted in-memory tool session,
while every public/persistable projection removes private reasoning fields.
"""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass, field
from http.client import HTTPException, IncompleteRead
import json
import re
import ssl
import time
from typing import Any, Callable, Mapping
from urllib.error import HTTPError, URLError
from urllib.parse import urlsplit
from urllib.request import Request

from chemsmart.agent._contracts import (
    ContractError,
    canonical_json,
    canonical_sha256,
)
from chemsmart.agent.runtime.transport import (
    ProviderDeadlineExceeded,
    ProviderTurnDeadline,
    ProviderTurnDeadlinesV1,
    is_socket_timeout,
    open_bounded_https_response,
    read_bounded_response_body,
)


DEEPSEEK_V4_FLASH_MODEL = "deepseek-v4-flash"
DEEPSEEK_OFFICIAL_ENDPOINT = "https://api.deepseek.com"
DEEPSEEK_V4_FLASH_CONTEXT_TOKENS = 1_000_000
DEEPSEEK_V4_FLASH_MAX_OUTPUT_TOKENS = 384_000
_PRIVATE_REASONING_KEYS = frozenset(
    {"reasoning_content", "thinking", "analysis", "<think>"}
)
_THINK_BLOCK = re.compile(
    r"<think\b[^>]*>.*?</think>", flags=re.IGNORECASE | re.DOTALL
)
_UNCLOSED_THINK = re.compile(r"<think\b", flags=re.IGNORECASE)


class DeepSeekProtocolError(RuntimeError):
    """Raised for an unsafe provider response without retaining its contents.

    A protocol failure can happen after the provider has returned a complete
    envelope but before any tool call is admissible.  Keep only a digest and
    byte count of the public (reasoning-scrubbed) envelope plus the structural
    location of the failure.  The raw response and private continuation never
    become exception attributes and therefore cannot leak into Runtime V2
    events or persisted tracebacks.
    """

    def __init__(
        self,
        message: str,
        *,
        failing_field: str = "",
        json_offset: int | None = None,
        response_envelope_sha256: str = "",
        response_envelope_bytes: int = 0,
    ) -> None:
        super().__init__(message)
        self.error_class = type(self).__name__
        self.failing_field = str(failing_field)
        self.json_offset = (
            None if json_offset is None else max(0, int(json_offset))
        )
        self.response_envelope_sha256 = str(response_envelope_sha256)
        self.response_envelope_bytes = max(0, int(response_envelope_bytes))

    def with_public_response(
        self, response: Mapping[str, Any]
    ) -> DeepSeekProtocolError:
        """Bind nonsecret envelope evidence without retaining the envelope."""

        public = public_provider_response(response)
        encoded = canonical_json(public).encode("utf-8")
        return type(self)(
            str(self),
            failing_field=self.failing_field,
            json_offset=self.json_offset,
            response_envelope_sha256=canonical_sha256(public),
            response_envelope_bytes=len(encoded),
        )

    def public_observation(self) -> dict[str, Any]:
        """Return the closed, nonsecret protocol-failure observation."""

        return {
            "schema_version": "chemsmart.provider-protocol-failure.v1",
            "error_class": self.error_class,
            "failing_field": self.failing_field,
            "json_offset": self.json_offset,
            "response_envelope_sha256": self.response_envelope_sha256,
            "response_envelope_bytes": self.response_envelope_bytes,
            "retry_decision": "not_retried",
            "recovery_requirement": "new_independent_attempt",
        }


class DeepSeekTransportError(RuntimeError):
    """Non-secret transport failure with a stable error class."""

    def __init__(
        self,
        error_class: str,
        *,
        http_status: int = 0,
        timeout_phase: str = "",
        timeout_seconds: float = 0.0,
    ) -> None:
        super().__init__(error_class)
        self.error_class = error_class
        self.http_status = int(http_status)
        self.timeout_phase = str(timeout_phase)
        self.timeout_seconds = max(0.0, float(timeout_seconds))

    def public_observation(self) -> dict[str, Any]:
        """Return only stable, non-secret transport and deadline facts."""

        return {
            "schema_version": "chemsmart.provider-transport-failure.v1",
            "error_class": self.error_class,
            "http_status": self.http_status,
            "timeout_phase": self.timeout_phase,
            "timeout_seconds": self.timeout_seconds,
            "retry_decision": "not_retried",
            "recovery_requirement": "new_independent_attempt",
        }


class DeepSeekHttpsTransport:
    """Minimal official HTTPS transport owned by one credential lease."""

    def __init__(
        self,
        *,
        api_key: str,
        endpoint: str = DEEPSEEK_OFFICIAL_ENDPOINT,
        timeout_seconds: float | None = None,
        turn_deadlines: ProviderTurnDeadlinesV1 | None = None,
        clock: Callable[[], float] = time.monotonic,
    ) -> None:
        if not _is_official_endpoint(endpoint):
            raise ContractError("DeepSeek transport requires the official endpoint")
        if not api_key:
            raise ContractError("DeepSeek transport requires a leased credential")
        policy = turn_deadlines or ProviderTurnDeadlinesV1()
        if timeout_seconds is not None:
            value = float(timeout_seconds)
            if value <= 0:
                raise ContractError("DeepSeek timeout must be positive")
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
        """Tighten the next provider turn to the remaining task allowance."""

        value = float(timeout_seconds)
        if value <= 0:
            raise ContractError("DeepSeek timeout must be positive")
        self._next_turn_timeout_seconds = min(
            value, self.turn_deadlines.absolute_seconds
        )

    def public_deadline_record(self) -> dict[str, float]:
        return dict(self._last_deadline_record)

    def __call__(self, payload: dict[str, Any]) -> Mapping[str, Any]:
        if self._closed or not self._api_key:
            raise ContractError("DeepSeek credential lease is closed")
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
                "Accept": "application/json",
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
                raw = read_bounded_response_body(
                    response, deadline=deadline
                )
        except ProviderDeadlineExceeded as exc:
            raise _deadline_transport_error(exc) from None
        except IncompleteRead:
            # A replay would be a distinct paid provider attempt. Close this
            # attempt visibly; only a separately authorized new conversation
            # may retry it.
            raise DeepSeekTransportError("incomplete_read") from None
        except HTTPError as exc:
            raise DeepSeekTransportError(
                _http_error_class(exc.code), http_status=exc.code
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
        try:
            decoded = json.loads(raw.decode("utf-8"))
        except (UnicodeDecodeError, json.JSONDecodeError) as exc:
            raise DeepSeekTransportError("invalid_json") from exc
        if not isinstance(decoded, Mapping):
            raise DeepSeekTransportError("invalid_response_shape")
        return decoded

    def close(self) -> None:
        self._api_key = ""
        self._closed = True


@dataclass(frozen=True)
class ProviderCapabilitiesV1:
    schema_version: str = "chemsmart.provider-capabilities.v1"
    provider: str = "deepseek"
    wire_protocol: str = "openai-chat-completions"
    model: str = DEEPSEEK_V4_FLASH_MODEL
    context_tokens: int = DEEPSEEK_V4_FLASH_CONTEXT_TOKENS
    max_output_tokens: int = DEEPSEEK_V4_FLASH_MAX_OUTPUT_TOKENS
    thinking_with_tools: bool = True
    continuation_mode: str = "assistant_reasoning_and_tool_results"
    private_reasoning_persisted: bool = False


@dataclass(frozen=True)
class DeepSeekV4FlashConfigV1:
    provider: str = "deepseek"
    model: str = DEEPSEEK_V4_FLASH_MODEL
    endpoint: str = DEEPSEEK_OFFICIAL_ENDPOINT
    thinking_mode: str = "enabled"
    reasoning_effort: str = "max"
    max_output_tokens: int = DEEPSEEK_V4_FLASH_MAX_OUTPUT_TOKENS
    context_tokens: int = DEEPSEEK_V4_FLASH_CONTEXT_TOKENS
    sdk_max_retries: int = 0
    turn_deadlines: ProviderTurnDeadlinesV1 = field(
        default_factory=ProviderTurnDeadlinesV1
    )

    def __post_init__(self) -> None:
        if self.provider != "deepseek":
            raise ContractError("DeepSeek provider identity is immutable")
        if self.model != DEEPSEEK_V4_FLASH_MODEL:
            raise ContractError("this adapter is pinned to deepseek-v4-flash")
        if not _is_official_endpoint(self.endpoint):
            raise ContractError("thinking continuation requires official DeepSeek")
        if self.thinking_mode != "enabled":
            raise ContractError("DeepSeek V4 tool experiments require thinking")
        if self.reasoning_effort not in {"high", "max"}:
            raise ContractError("reasoning_effort must be high or max")
        if not 1 <= self.max_output_tokens <= (
            DEEPSEEK_V4_FLASH_MAX_OUTPUT_TOKENS
        ):
            raise ContractError(
                "max_output_tokens exceeds DeepSeek V4 Flash capability"
            )
        if self.sdk_max_retries != 0:
            raise ContractError(
                "provider retries require a separately authorized attempt"
            )


@dataclass(frozen=True)
class ProviderTurnReceiptV1:
    schema_version: str
    provider: str
    requested_model: str
    observed_model: str
    request_sha256: str
    tool_schema_sha256: str
    input_tokens: int
    output_tokens: int
    reasoning_tokens: int
    finish_reason: str
    tool_calls_present: bool
    reasoning_continuation_present: bool
    private_reasoning_persisted: bool
    receipt_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.provider-turn-receipt.v1":
            raise ContractError("unsupported provider turn receipt schema")
        if self.private_reasoning_persisted:
            raise ContractError("private provider reasoning cannot be persisted")
        expected = canonical_sha256(
            {
                "schema_version": self.schema_version,
                "provider": self.provider,
                "requested_model": self.requested_model,
                "observed_model": self.observed_model,
                "request_sha256": self.request_sha256,
                "tool_schema_sha256": self.tool_schema_sha256,
                "input_tokens": self.input_tokens,
                "output_tokens": self.output_tokens,
                "reasoning_tokens": self.reasoning_tokens,
                "finish_reason": self.finish_reason,
                "tool_calls_present": self.tool_calls_present,
                "reasoning_continuation_present": (
                    self.reasoning_continuation_present
                ),
                "private_reasoning_persisted": self.private_reasoning_persisted,
            }
        )
        if self.receipt_sha256 != expected:
            raise ContractError("provider turn receipt digest mismatch")


class DeepSeekV4ToolSession:
    """One in-memory conversation with full provider-native continuation."""

    def __init__(
        self,
        *,
        transport: Callable[[dict[str, Any]], Mapping[str, Any]],
        messages: list[dict[str, Any]],
        config: DeepSeekV4FlashConfigV1 | None = None,
    ) -> None:
        self.config = config or DeepSeekV4FlashConfigV1()
        self._transport = transport
        self._history = deepcopy(messages)
        self._receipts: list[ProviderTurnReceiptV1] = []
        self._outstanding_tool_call_ids: tuple[str, ...] = ()
        self._seen_tool_call_ids: set[str] = set()

    @property
    def capabilities(self) -> ProviderCapabilitiesV1:
        return ProviderCapabilitiesV1(
            provider=self.config.provider,
            model=self.config.model,
            context_tokens=self.config.context_tokens,
            max_output_tokens=self.config.max_output_tokens,
        )

    @property
    def receipts(self) -> tuple[ProviderTurnReceiptV1, ...]:
        return tuple(self._receipts)

    @property
    def outstanding_tool_call_ids(self) -> tuple[str, ...]:
        return self._outstanding_tool_call_ids

    def request_payload(
        self, *, tools: list[dict[str, Any]] | None = None
    ) -> dict[str, Any]:
        """Return the exact provider request with ephemeral reasoning intact."""

        payload: dict[str, Any] = {
            "model": self.config.model,
            "messages": deepcopy(self._history),
            "max_tokens": self.config.max_output_tokens,
            "reasoning_effort": self.config.reasoning_effort,
            # This is the official raw wire shape. An OpenAI-SDK bridge may
            # place it in ``extra_body`` at its own adapter boundary.
            "thinking": {"type": "enabled"},
        }
        if tools:
            payload["tools"] = deepcopy(tools)
        return payload

    def set_turn_timeout_seconds(self, timeout_seconds: float) -> None:
        """Bind a real transport turn to the loop's remaining wall time."""

        setter = getattr(self._transport, "set_timeout_seconds", None)
        if setter is None:
            return
        setter(float(timeout_seconds))

    def public_transport_deadline_record(self) -> dict[str, float]:
        """Expose only the transport's effective public deadline policy."""

        getter = getattr(self._transport, "public_deadline_record", None)
        return dict(getter()) if getter is not None else {}

    def turn(
        self, *, tools: list[dict[str, Any]] | None = None
    ) -> tuple[dict[str, Any], ProviderTurnReceiptV1]:
        """Call the injected transport and retain private continuation in RAM."""

        if self._outstanding_tool_call_ids:
            raise DeepSeekProtocolError(
                "all tool results must be appended before the next provider turn"
            )
        request = self.request_payload(tools=tools)
        response = self._transport(deepcopy(request))
        if not isinstance(response, Mapping):
            raise DeepSeekProtocolError("provider response must be a mapping")
        payload = deepcopy(dict(response))
        try:
            assistant = _assistant_message(payload)
            tool_call_ids = _validate_tool_calls(assistant)
            if self._seen_tool_call_ids.intersection(tool_call_ids):
                raise DeepSeekProtocolError(
                    "provider reused a tool-call ID within one session",
                    failing_field="choices[0].message.tool_calls[*].id",
                )
        except DeepSeekProtocolError as exc:
            # Validation is intentionally all-or-nothing and precedes every
            # history mutation, receipt, and tool dispatch.  Do not repair or
            # replay malformed provider-authored arguments.
            raise exc.with_public_response(payload) from None
        if tool_call_ids:
            # Short follow-up tool turns may carry an explicitly empty
            # continuation. Replay the exact empty field rather than inventing
            # private reasoning or aborting an otherwise valid exchange.
            assistant.setdefault("reasoning_content", "")
        # Preserve the entire provider message. In particular, do not trim or
        # summarize reasoning_content before the next tool-result subturn.
        self._history.append(deepcopy(assistant))
        receipt = _turn_receipt(
            request=request,
            response=payload,
            assistant=assistant,
            tools=tools or [],
            config=self.config,
        )
        self._receipts.append(receipt)
        self._outstanding_tool_call_ids = tool_call_ids
        self._seen_tool_call_ids.update(tool_call_ids)
        return public_provider_response(payload), receipt

    def append_tool_results(
        self, results: list[dict[str, Any]]
    ) -> None:
        """Append OpenAI-shaped tool results after the preserved assistant."""

        observed_ids = tuple(str(item.get("tool_call_id") or "") for item in results)
        if tuple(sorted(observed_ids)) != tuple(
            sorted(self._outstanding_tool_call_ids)
        ):
            raise ContractError("tool results must exactly match outstanding calls")
        if len(observed_ids) != len(set(observed_ids)):
            raise ContractError("tool results contain duplicate call IDs")
        for result in results:
            if result.get("role") != "tool" or not result.get("tool_call_id"):
                raise ContractError("tool results require role and tool_call_id")
            if not isinstance(result.get("content"), str):
                raise ContractError("tool result content must be a JSON string")
            self._history.append(deepcopy(result))
        self._outstanding_tool_call_ids = ()

    def public_history(self) -> list[dict[str, Any]]:
        """Return the only history shape permitted in events or artifacts."""

        return [
            _public_message(message)
            for message in self._history
        ]

    def close(self) -> None:
        """Discard provider-private continuation state deterministically."""

        self._history = self.public_history()
        self._outstanding_tool_call_ids = ()
        self._seen_tool_call_ids.clear()


def public_provider_response(response: Mapping[str, Any]) -> dict[str, Any]:
    sanitized = _public_payload_at_path(
        deepcopy(dict(response)), (), "provider_response"
    )
    return sanitized if isinstance(sanitized, dict) else {}


def public_provider_request(request: Mapping[str, Any]) -> dict[str, Any]:
    """Return one reasoning-free request without changing tool actions."""

    sanitized = _public_payload_at_path(
        deepcopy(dict(request)), (), "provider_request"
    )
    return sanitized if isinstance(sanitized, dict) else {}


def _public_message(message: Mapping[str, Any]) -> dict[str, Any]:
    context = (
        "assistant_message"
        if str(message.get("role") or "") == "assistant"
        else "generic"
    )
    sanitized = _public_payload_at_path(deepcopy(dict(message)), (), context)
    return sanitized if isinstance(sanitized, dict) else {}


def public_payload(value: Any) -> Any:
    """Recursively remove private reasoning fields from persistable values."""

    return _public_payload_at_path(value, (), "generic")


def _is_tool_function_path(
    path: tuple[Any, ...], context: str
) -> bool:
    direct_message = (
        context == "assistant_message"
        and len(path) == 3
        and path[0] == "tool_calls"
        and isinstance(path[1], int)
        and path[2] == "function"
    )
    provider_response = (
        context == "provider_response"
        and len(path) == 6
        and path[0] == "choices"
        and path[1] == 0
        and path[2] == "message"
        and path[3] == "tool_calls"
        and isinstance(path[4], int)
        and path[5] == "function"
    )
    provider_request = (
        context == "provider_request"
        and len(path) == 5
        and path[0] == "messages"
        and isinstance(path[1], int)
        and path[2] == "tool_calls"
        and isinstance(path[3], int)
        and path[4] == "function"
    )
    return direct_message or provider_response or provider_request


def _public_payload_at_path(
    value: Any, path: tuple[Any, ...], context: str
) -> Any:
    """Sanitize one value while retaining its exact provider-envelope path."""

    if isinstance(value, dict):
        if str(value.get("type", "")).lower() in _PRIVATE_REASONING_KEYS:
            return None
        result = {}
        for key, item in value.items():
            if key in _PRIVATE_REASONING_KEYS:
                continue
            if _is_tool_function_path(path, context) and key == "arguments":
                # Function arguments are a validated public action, not a
                # reasoning channel.  Rewriting substrings inside their JSON
                # text can change the action or make a later call undecodable.
                result[key] = item
                continue
            public = _public_payload_at_path(
                item, (*path, key), context
            )
            if public is not None:
                result[key] = public
        return result
    if isinstance(value, list):
        result = []
        for index, item in enumerate(value):
            public = _public_payload_at_path(
                item, (*path, index), context
            )
            if public is not None:
                result.append(public)
        return result
    if isinstance(value, str):
        public = _THINK_BLOCK.sub("", value)
        unclosed = _UNCLOSED_THINK.search(public)
        if unclosed is not None:
            public = public[: unclosed.start()]
        return public
    return value


def _assistant_message(response: Mapping[str, Any]) -> dict[str, Any]:
    choices = response.get("choices")
    if not isinstance(choices, list) or not choices:
        raise DeepSeekProtocolError(
            "provider response has no choices", failing_field="choices"
        )
    choice = choices[0]
    message = choice.get("message") if isinstance(choice, Mapping) else None
    if not isinstance(message, Mapping):
        raise DeepSeekProtocolError(
            "provider response has no assistant message",
            failing_field="choices[0].message",
        )
    result = deepcopy(dict(message))
    result.setdefault("role", "assistant")
    # Official thinking+tools continuation requires replaying a non-null
    # assistant content field alongside the verbatim reasoning_content.
    if result.get("tool_calls") and result.get("content") is None:
        result["content"] = ""
    return result


def _validate_tool_calls(assistant: Mapping[str, Any]) -> tuple[str, ...]:
    raw = assistant.get("tool_calls") or []
    if not isinstance(raw, list):
        raise DeepSeekProtocolError(
            "assistant tool_calls must be a list",
            failing_field="choices[0].message.tool_calls",
        )
    identifiers = []
    for index, item in enumerate(raw):
        field = f"choices[0].message.tool_calls[{index}]"
        if not isinstance(item, Mapping):
            raise DeepSeekProtocolError(
                "assistant tool call must be an object", failing_field=field
            )
        identifier = str(item.get("id") or "")
        function = item.get("function")
        if not identifier or not isinstance(function, Mapping):
            raise DeepSeekProtocolError(
                "tool call lacks ID or function", failing_field=field
            )
        if not str(function.get("name") or ""):
            raise DeepSeekProtocolError(
                "tool call lacks function name",
                failing_field=field + ".function.name",
            )
        arguments = function.get("arguments")
        if not isinstance(arguments, str):
            raise DeepSeekProtocolError(
                "tool call arguments must be JSON text",
                failing_field=field + ".function.arguments",
            )
        try:
            decoded = json.loads(arguments)
        except json.JSONDecodeError as exc:
            raise DeepSeekProtocolError(
                "tool call arguments are invalid JSON",
                failing_field=field + ".function.arguments",
                json_offset=exc.pos,
            ) from None
        if not isinstance(decoded, dict):
            raise DeepSeekProtocolError(
                "tool call arguments must decode to an object",
                failing_field=field + ".function.arguments",
            )
        identifiers.append(identifier)
    if len(identifiers) != len(set(identifiers)):
        raise DeepSeekProtocolError(
            "assistant reused a tool-call ID",
            failing_field="choices[0].message.tool_calls[*].id",
        )
    return tuple(identifiers)


def _turn_receipt(
    *,
    request: dict[str, Any],
    response: dict[str, Any],
    assistant: dict[str, Any],
    tools: list[dict[str, Any]],
    config: Any,
) -> ProviderTurnReceiptV1:
    choices = response.get("choices") or [{}]
    choice = choices[0] if isinstance(choices[0], Mapping) else {}
    usage = response.get("usage")
    usage = usage if isinstance(usage, Mapping) else {}
    details = usage.get("completion_tokens_details")
    details = details if isinstance(details, Mapping) else {}
    body = {
        "schema_version": "chemsmart.provider-turn-receipt.v1",
        "provider": str(getattr(config, "provider", "deepseek")),
        "requested_model": config.model,
        "observed_model": str(response.get("model") or ""),
        "request_sha256": canonical_sha256(public_provider_request(request)),
        "tool_schema_sha256": canonical_sha256(tools),
        "input_tokens": int(
            usage.get("prompt_tokens", usage.get("input_tokens", 0)) or 0
        ),
        "output_tokens": int(
            usage.get("completion_tokens", usage.get("output_tokens", 0)) or 0
        ),
        "reasoning_tokens": int(details.get("reasoning_tokens", 0) or 0),
        "finish_reason": str(choice.get("finish_reason") or ""),
        "tool_calls_present": bool(assistant.get("tool_calls")),
        "reasoning_continuation_present": isinstance(
            assistant.get("reasoning_content"), str
        ),
        "private_reasoning_persisted": False,
    }
    return ProviderTurnReceiptV1(
        **body, receipt_sha256=canonical_sha256(body)
    )


def _is_official_endpoint(endpoint: str) -> bool:
    try:
        parsed = urlsplit(str(endpoint).strip())
        return (
            parsed.scheme.lower() == "https"
            and parsed.hostname == "api.deepseek.com"
            and parsed.port in {None, 443}
            and parsed.path.rstrip("/") in {"", "/v1"}
            and not parsed.username
            and not parsed.password
            and not parsed.query
            and not parsed.fragment
        )
    except ValueError:
        return False


def _http_error_class(status: int) -> str:
    if status == 401:
        return "credential_invalid"
    if status == 402:
        return "quota_exhausted"
    if status == 429:
        return "rate_limited"
    if 500 <= status < 600:
        return "provider_5xx"
    return "http_error"


def _deadline_transport_error(
    error: ProviderDeadlineExceeded,
) -> DeepSeekTransportError:
    error_class = {
        "connect": "connect_timeout",
        "first_event": "first_event_timeout",
        "inter_event": "inter_event_timeout",
        "absolute": "turn_deadline_exceeded",
    }.get(error.phase, "timeout")
    return DeepSeekTransportError(
        error_class,
        timeout_phase=error.phase,
        timeout_seconds=error.timeout_seconds,
    )


__all__ = [
    "DEEPSEEK_OFFICIAL_ENDPOINT",
    "DEEPSEEK_V4_FLASH_CONTEXT_TOKENS",
    "DEEPSEEK_V4_FLASH_MAX_OUTPUT_TOKENS",
    "DEEPSEEK_V4_FLASH_MODEL",
    "DeepSeekProtocolError",
    "DeepSeekTransportError",
    "DeepSeekHttpsTransport",
    "DeepSeekV4FlashConfigV1",
    "DeepSeekV4ToolSession",
    "ProviderCapabilitiesV1",
    "ProviderTurnReceiptV1",
    "ProviderTurnDeadlinesV1",
    "public_payload",
    "public_provider_request",
    "public_provider_response",
]
