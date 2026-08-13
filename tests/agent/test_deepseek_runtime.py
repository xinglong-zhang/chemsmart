from __future__ import annotations

import json
from http.client import IncompleteRead

import pytest

from chemsmart.agent._contracts import ContractError
from chemsmart.agent.runtime.deepseek import (
    DEEPSEEK_V4_FLASH_CONTEXT_TOKENS,
    DEEPSEEK_V4_FLASH_MAX_OUTPUT_TOKENS,
    DeepSeekHttpsTransport,
    DeepSeekTransportError,
    DeepSeekV4FlashConfigV1,
    DeepSeekV4ToolSession,
)


def test_deepseek_thinking_tool_continuation_is_full_but_ephemeral():
    requests = []
    responses = iter(
        [
            {
                "model": "deepseek-v4-flash",
                "choices": [
                    {
                        "finish_reason": "tool_calls",
                        "message": {
                            "role": "assistant",
                            "content": "<think>inline private</think>Inspecting.",
                            "reasoning_content": "complete private reasoning",
                            "tool_calls": [
                                {
                                    "id": "call-1",
                                    "type": "function",
                                    "function": {
                                        "name": "query_program_capability",
                                        "arguments": '{"program":"pyscf"}',
                                    },
                                }
                            ],
                        },
                    }
                ],
                "usage": {
                    "prompt_tokens": 20,
                    "completion_tokens": 10,
                    "completion_tokens_details": {"reasoning_tokens": 8},
                },
            },
            {
                "model": "deepseek-v4-flash",
                "choices": [
                    {
                        "finish_reason": "stop",
                        "message": {
                            "role": "assistant",
                            "content": "Done.",
                            "reasoning_content": "second private reasoning",
                        },
                    }
                ],
                "usage": {"prompt_tokens": 40, "completion_tokens": 5},
            },
        ]
    )

    def transport(payload):
        requests.append(payload)
        return next(responses)

    session = DeepSeekV4ToolSession(
        transport=transport,
        messages=[{"role": "user", "content": "Plan a PySCF node."}],
    )
    public_first, first_receipt = session.turn(tools=[{"type": "function"}])
    session.append_tool_results(
        [
            {
                "role": "tool",
                "tool_call_id": "call-1",
                "content": '{"status":"supported"}',
            }
        ]
    )
    public_second, _ = session.turn(tools=[{"type": "function"}])

    assert requests[0]["thinking"] == {"type": "enabled"}
    assert "extra_body" not in requests[0]
    assert requests[0]["max_tokens"] == 384_000
    assert requests[1]["messages"][1]["reasoning_content"] == (
        "complete private reasoning"
    )
    assert "reasoning_content" not in json.dumps(public_first)
    assert "inline private" not in json.dumps(public_first)
    assert public_first["choices"][0]["message"]["content"] == "Inspecting."
    assert "second private reasoning" not in json.dumps(public_second)
    assert first_receipt.reasoning_continuation_present is True
    assert first_receipt.private_reasoning_persisted is False
    assert DEEPSEEK_V4_FLASH_CONTEXT_TOKENS == 1_000_000
    assert DEEPSEEK_V4_FLASH_MAX_OUTPUT_TOKENS == 384_000

    session.close()
    assert "private reasoning" not in json.dumps(session.public_history())


def test_incomplete_response_is_not_silently_retried(monkeypatch):
    import chemsmart.agent.runtime.deepseek as deepseek

    calls = []

    class _IncompleteResponse:
        def __enter__(self):
            return self

        def __exit__(self, *_args):
            return False

        def read(self):
            raise IncompleteRead(b"", 1)

    def _urlopen(*args, **kwargs):
        calls.append((args, kwargs))
        return _IncompleteResponse()

    monkeypatch.setattr(deepseek, "urlopen", _urlopen)
    transport = DeepSeekHttpsTransport(api_key="not-persisted")

    with pytest.raises(DeepSeekTransportError, match="incomplete_read"):
        transport({"model": "deepseek-v4-flash", "messages": []})

    assert len(calls) == 1


def test_deepseek_config_rejects_implicit_provider_retries():
    with pytest.raises(ContractError, match="separately authorized attempt"):
        DeepSeekV4FlashConfigV1(sdk_max_retries=1)
