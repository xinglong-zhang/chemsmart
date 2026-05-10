from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Any, Literal

ToolOutcomeStatus = Literal[
    "ok",
    "denied",
    "error",
    "skipped",
    "interrupted",
]


@dataclass(frozen=True)
class ToolRequest:
    request_id: str
    provider: str
    provider_call_id: str
    name: str
    arguments_json: str
    arguments: dict[str, Any]
    raw: dict[str, Any]


@dataclass(frozen=True)
class ToolOutcome:
    request_id: str
    provider_call_id: str
    name: str
    status: ToolOutcomeStatus
    result: Any = None
    display_result: Any = None
    error_type: str | None = None
    error_message: str | None = None


def normalize_response(
    provider: str,
    response: Any,
) -> tuple[str, list[ToolRequest], str | None]:
    payload = _response_payload(response)
    if provider == "anthropic":
        return _normalize_anthropic_response(payload)
    if provider == "openai":
        return _normalize_openai_response(payload)
    raise ValueError(f"Unsupported provider {provider!r}")


def build_tool_result_messages(
    provider: str,
    outcomes: list[ToolOutcome],
) -> list[dict[str, Any]]:
    if not outcomes:
        return []

    if provider == "anthropic":
        return [
            {
                "role": "user",
                "content": [
                    _anthropic_tool_result_block(outcome)
                    for outcome in outcomes
                ],
            }
        ]

    if provider == "openai":
        return [
            {
                "role": "tool",
                "tool_call_id": outcome.provider_call_id,
                "content": json.dumps(
                    _outcome_payload(outcome),
                    sort_keys=True,
                ),
            }
            for outcome in outcomes
        ]

    raise ValueError(f"Unsupported provider {provider!r}")


def _normalize_anthropic_response(
    payload: dict[str, Any],
) -> tuple[str, list[ToolRequest], str | None]:
    text_parts: list[str] = []
    requests: list[ToolRequest] = []
    for index, block in enumerate(payload.get("content") or []):
        if not isinstance(block, dict):
            continue
        block_type = block.get("type")
        if block_type == "text":
            text = block.get("text")
            if isinstance(text, str) and text:
                text_parts.append(text)
            continue
        if block_type != "tool_use":
            continue
        provider_call_id = str(block.get("id") or f"anthropic-{index}")
        arguments = block.get("input") or {}
        if not isinstance(arguments, dict):
            raise ValueError("Anthropic tool_use.input must be a JSON object")
        requests.append(
            ToolRequest(
                request_id=f"anthropic:{provider_call_id}",
                provider="anthropic",
                provider_call_id=provider_call_id,
                name=str(block.get("name") or ""),
                arguments_json=json.dumps(arguments, sort_keys=True),
                arguments=arguments,
                raw=block,
            )
        )
    return "\n".join(text_parts), requests, payload.get("stop_reason")


def _normalize_openai_response(
    payload: dict[str, Any],
) -> tuple[str, list[ToolRequest], str | None]:
    choices = payload.get("choices") or []
    if not choices:
        return "", [], None

    choice = choices[0] or {}
    message = choice.get("message") or {}
    text = _openai_message_text(message)
    requests: list[ToolRequest] = []
    for index, call in enumerate(message.get("tool_calls") or []):
        if not isinstance(call, dict):
            continue
        function = call.get("function") or {}
        arguments_json = function.get("arguments") or "{}"
        if not isinstance(arguments_json, str):
            raise ValueError("OpenAI function arguments must be a JSON string")
        arguments = (
            json.loads(arguments_json) if arguments_json.strip() else {}
        )
        if not isinstance(arguments, dict):
            raise ValueError(
                "OpenAI function arguments must decode to a JSON object"
            )
        provider_call_id = str(call.get("id") or f"openai-{index}")
        requests.append(
            ToolRequest(
                request_id=f"openai:{provider_call_id}",
                provider="openai",
                provider_call_id=provider_call_id,
                name=str(function.get("name") or ""),
                arguments_json=arguments_json,
                arguments=arguments,
                raw=call,
            )
        )
    return text, requests, choice.get("finish_reason")


def _openai_message_text(message: dict[str, Any]) -> str:
    content = message.get("content")
    if isinstance(content, str):
        return content
    if not isinstance(content, list):
        return ""

    text_parts: list[str] = []
    for part in content:
        if isinstance(part, dict) and part.get("type") == "text":
            text = part.get("text")
            if isinstance(text, str) and text:
                text_parts.append(text)
    return "\n".join(text_parts)


def _response_payload(response: Any) -> dict[str, Any]:
    if isinstance(response, dict):
        return response
    if hasattr(response, "model_dump"):
        payload = response.model_dump()
        if isinstance(payload, dict):
            return payload
    raise TypeError("Provider response must be a dict-like payload")


def _anthropic_tool_result_block(outcome: ToolOutcome) -> dict[str, Any]:
    block = {
        "type": "tool_result",
        "tool_use_id": outcome.provider_call_id,
        "content": json.dumps(
            _outcome_payload(outcome),
            sort_keys=True,
        ),
    }
    if outcome.status in {"denied", "error"}:
        block["is_error"] = True
    return block


def _outcome_payload(outcome: ToolOutcome) -> Any:
    if outcome.display_result is not None:
        return outcome.display_result
    if outcome.result is not None:
        return outcome.result
    if outcome.status == "ok":
        return {"ok": True}
    return {
        "ok": False,
        "error": {
            "type": outcome.error_type or _default_error_type(outcome.status),
            "message": outcome.error_message
            or _default_error_message(outcome),
            "tool": outcome.name,
        },
    }


def _default_error_type(status: ToolOutcomeStatus) -> str:
    return {
        "denied": "PermissionDenied",
        "error": "ToolError",
        "skipped": "ToolSkipped",
        "interrupted": "Interrupted",
        "ok": "ToolError",
    }[status]


def _default_error_message(outcome: ToolOutcome) -> str:
    return f"{outcome.name} returned status={outcome.status}"
