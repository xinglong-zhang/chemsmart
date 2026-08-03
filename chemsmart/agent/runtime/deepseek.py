"""DeepSeek V4 Flash thinking/tool continuation with private ephemeral state.

The transport is injected so this module neither owns credentials nor adds an
SDK dependency to ChemSmart v3.1.4.  It preserves the provider's complete
``reasoning_content`` only inside an uninterrupted in-memory tool session,
while every public/persistable projection removes private reasoning fields.
"""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
import json
import re
from typing import Any, Callable, Mapping
from urllib.error import HTTPError, URLError
from urllib.parse import urlsplit
from urllib.request import Request, urlopen

from chemsmart.agent._contracts import ContractError, canonical_sha256


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
    """Raised for an unsafe or malformed provider response."""


class DeepSeekTransportError(RuntimeError):
    """Non-secret transport failure with a stable error class."""

    def __init__(self, error_class: str, *, http_status: int = 0) -> None:
        super().__init__(error_class)
        self.error_class = error_class
        self.http_status = int(http_status)


class DeepSeekHttpsTransport:
    """Minimal official HTTPS transport owned by one credential lease."""

    def __init__(
        self,
        *,
        api_key: str,
        endpoint: str = DEEPSEEK_OFFICIAL_ENDPOINT,
        timeout_seconds: float = 120.0,
    ) -> None:
        if not _is_official_endpoint(endpoint):
            raise ContractError("DeepSeek transport requires the official endpoint")
        if not api_key:
            raise ContractError("DeepSeek transport requires a leased credential")
        if timeout_seconds <= 0:
            raise ContractError("DeepSeek timeout must be positive")
        self.endpoint = endpoint.rstrip("/")
        self.timeout_seconds = float(timeout_seconds)
        self._api_key = api_key
        self._closed = False

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
        try:
            with urlopen(request, timeout=self.timeout_seconds) as response:
                raw = response.read()
        except HTTPError as exc:
            raise DeepSeekTransportError(
                _http_error_class(exc.code), http_status=exc.code
            ) from None
        except TimeoutError:
            raise DeepSeekTransportError("timeout") from None
        except URLError as exc:
            reason = getattr(exc, "reason", None)
            error_class = "timeout" if isinstance(reason, TimeoutError) else "transport"
            raise DeepSeekTransportError(error_class) from None
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
    model: str = DEEPSEEK_V4_FLASH_MODEL
    endpoint: str = DEEPSEEK_OFFICIAL_ENDPOINT
    thinking_mode: str = "enabled"
    reasoning_effort: str = "max"
    max_output_tokens: int = DEEPSEEK_V4_FLASH_MAX_OUTPUT_TOKENS
    sdk_max_retries: int = 0

    def __post_init__(self) -> None:
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
        if self.sdk_max_retries < 0:
            raise ContractError("sdk_max_retries must be non-negative")


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
            model=self.config.model,
            max_output_tokens=DEEPSEEK_V4_FLASH_MAX_OUTPUT_TOKENS,
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
        assistant = _assistant_message(payload)
        tool_call_ids = _validate_tool_calls(assistant)
        if self._seen_tool_call_ids.intersection(tool_call_ids):
            raise DeepSeekProtocolError(
                "provider reused a tool-call ID within one session"
            )
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
            public_payload(message)
            for message in self._history
        ]

    def close(self) -> None:
        """Discard provider-private continuation state deterministically."""

        self._history = self.public_history()
        self._outstanding_tool_call_ids = ()
        self._seen_tool_call_ids.clear()


def public_provider_response(response: Mapping[str, Any]) -> dict[str, Any]:
    sanitized = public_payload(deepcopy(dict(response)))
    return sanitized if isinstance(sanitized, dict) else {}


def public_payload(value: Any) -> Any:
    """Recursively remove private reasoning fields from persistable values."""

    if isinstance(value, dict):
        if str(value.get("type", "")).lower() in _PRIVATE_REASONING_KEYS:
            return None
        result = {}
        for key, item in value.items():
            if key in _PRIVATE_REASONING_KEYS:
                continue
            public = public_payload(item)
            if public is not None:
                result[key] = public
        return result
    if isinstance(value, list):
        result = []
        for item in value:
            public = public_payload(item)
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
        raise DeepSeekProtocolError("provider response has no choices")
    choice = choices[0]
    message = choice.get("message") if isinstance(choice, Mapping) else None
    if not isinstance(message, Mapping):
        raise DeepSeekProtocolError("provider response has no assistant message")
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
        raise DeepSeekProtocolError("assistant tool_calls must be a list")
    identifiers = []
    for item in raw:
        if not isinstance(item, Mapping):
            raise DeepSeekProtocolError("assistant tool call must be an object")
        identifier = str(item.get("id") or "")
        function = item.get("function")
        if not identifier or not isinstance(function, Mapping):
            raise DeepSeekProtocolError("tool call lacks ID or function")
        if not str(function.get("name") or ""):
            raise DeepSeekProtocolError("tool call lacks function name")
        arguments = function.get("arguments")
        if not isinstance(arguments, str):
            raise DeepSeekProtocolError("tool call arguments must be JSON text")
        try:
            decoded = json.loads(arguments)
        except json.JSONDecodeError as exc:
            raise DeepSeekProtocolError("tool call arguments are invalid JSON") from exc
        if not isinstance(decoded, dict):
            raise DeepSeekProtocolError("tool call arguments must decode to an object")
        identifiers.append(identifier)
    if len(identifiers) != len(set(identifiers)):
        raise DeepSeekProtocolError("assistant reused a tool-call ID")
    return tuple(identifiers)


def _turn_receipt(
    *,
    request: dict[str, Any],
    response: dict[str, Any],
    assistant: dict[str, Any],
    tools: list[dict[str, Any]],
    config: DeepSeekV4FlashConfigV1,
) -> ProviderTurnReceiptV1:
    choices = response.get("choices") or [{}]
    choice = choices[0] if isinstance(choices[0], Mapping) else {}
    usage = response.get("usage")
    usage = usage if isinstance(usage, Mapping) else {}
    details = usage.get("completion_tokens_details")
    details = details if isinstance(details, Mapping) else {}
    body = {
        "schema_version": "chemsmart.provider-turn-receipt.v1",
        "provider": "deepseek",
        "requested_model": config.model,
        "observed_model": str(response.get("model") or ""),
        "request_sha256": canonical_sha256(public_payload(request)),
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
    "public_payload",
    "public_provider_response",
]
