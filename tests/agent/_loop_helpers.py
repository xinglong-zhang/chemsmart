from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Any


@dataclass
class DummyMolecule:
    label: str


class ScriptedRegistry:
    def __init__(self, results: dict[str, Any]):
        self._results = dict(results)
        self.calls: list[tuple[str, dict[str, Any]]] = []

    def call(self, name: str, args: dict[str, Any] | None = None) -> Any:
        payload = dict(args or {})
        self.calls.append((name, payload))
        result = self._results[name]
        return result(payload) if callable(result) else result

    def tool_defs_for_provider(
        self, provider_name: str
    ) -> list[dict[str, Any]]:
        return self.openai_tool_defs()

    def openai_tool_defs(self) -> list[dict[str, Any]]:
        tool_names = list(self._results)
        return [openai_tool_def(name) for name in tool_names]


def openai_tool_def(name: str) -> dict[str, Any]:
    return {
        "type": "function",
        "function": {
            "name": name,
            "description": name,
            "parameters": {"type": "object", "properties": {}},
        },
    }


def openai_tool_call_response(
    *tool_calls: dict[str, Any],
    content: str = "",
    finish_reason: str = "tool_calls",
) -> dict[str, Any]:
    return {
        "choices": [
            {
                "message": {
                    "role": "assistant",
                    "content": content,
                    "tool_calls": list(tool_calls),
                },
                "finish_reason": finish_reason,
            }
        ],
        "usage": {"prompt_tokens": 10, "completion_tokens": 5},
    }


def openai_final_response(text: str) -> dict[str, Any]:
    return {
        "choices": [
            {
                "message": {"role": "assistant", "content": text},
                "finish_reason": "stop",
            }
        ],
        "usage": {"prompt_tokens": 8, "completion_tokens": 4},
    }


def anthropic_tool_use_response(
    *tool_uses: dict[str, Any],
    text: str = "",
    stop_reason: str = "tool_use",
) -> dict[str, Any]:
    content = []
    if text:
        content.append({"type": "text", "text": text})
    content.extend(tool_uses)
    return {
        "content": content,
        "stop_reason": stop_reason,
        "usage": {"input_tokens": 10, "output_tokens": 5},
    }


def anthropic_final_response(text: str) -> dict[str, Any]:
    return {
        "content": [{"type": "text", "text": text}],
        "stop_reason": "end_turn",
        "usage": {"input_tokens": 8, "output_tokens": 4},
    }


def tool_call(
    call_id: str,
    name: str,
    arguments: dict[str, Any],
) -> dict[str, Any]:
    return {
        "id": call_id,
        "type": "function",
        "function": {
            "name": name,
            "arguments": json.dumps(arguments, sort_keys=True),
        },
    }


def anthropic_tool_use(
    call_id: str,
    name: str,
    arguments: dict[str, Any],
) -> dict[str, Any]:
    return {
        "type": "tool_use",
        "id": call_id,
        "name": name,
        "input": arguments,
    }
