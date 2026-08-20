from __future__ import annotations

import json
from http.client import IncompleteRead

import pytest

from chemsmart.agent._contracts import ContractError, canonical_sha256
from chemsmart.agent.runtime.deepseek import (
    DEEPSEEK_V4_FLASH_CONTEXT_TOKENS,
    DEEPSEEK_V4_FLASH_MAX_OUTPUT_TOKENS,
    DeepSeekHttpsTransport,
    DeepSeekTransportError,
    DeepSeekV4FlashConfigV1,
    DeepSeekV4ToolSession,
    ProviderTurnReceiptV1,
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
        config=DeepSeekV4FlashConfigV1(
            model="deepseek-v4-flash",
            context_tokens=1_000_000,
            max_output_tokens=384_000,
        ),
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

    def _open_response(*args, **kwargs):
        calls.append((args, kwargs))
        return _IncompleteResponse()

    monkeypatch.setattr(
        deepseek, "open_bounded_https_response", _open_response
    )
    transport = DeepSeekHttpsTransport(api_key="not-persisted")

    with pytest.raises(DeepSeekTransportError, match="incomplete_read"):
        transport({"model": "deepseek-v4-flash", "messages": []})

    assert len(calls) == 1


def test_deepseek_config_rejects_implicit_provider_retries():
    with pytest.raises(ContractError, match="separately authorized attempt"):
        DeepSeekV4FlashConfigV1(
            model="deepseek-v4-flash",
            context_tokens=131072,
            max_output_tokens=8192,
            sdk_max_retries=1,
        )


@pytest.mark.parametrize(
    "endpoint",
    (
        "https://user:SUPERSECRET@api.deepseek.com/v1",
        "https://api.deepseek.com/v1?api_key=SUPERSECRET",
        "https://api.deepseek.com/v1#SUPERSECRET",
    ),
)
def test_deepseek_endpoint_rejects_non_authority_data_without_reflection(
    endpoint,
):
    with pytest.raises(ContractError) as config_error:
        DeepSeekV4FlashConfigV1(
            model="deepseek-v4-flash",
            context_tokens=131072,
            max_output_tokens=8192,
            endpoint=endpoint,
        )
    with pytest.raises(ContractError) as transport_error:
        DeepSeekHttpsTransport(api_key="not-persisted", endpoint=endpoint)

    assert "SUPERSECRET" not in str(config_error.value)
    assert "SUPERSECRET" not in str(transport_error.value)


def test_historical_v1_provider_turn_receipt_digest_remains_unchanged():
    body = {
        "schema_version": "chemsmart.provider-turn-receipt.v1",
        "provider": "deepseek",
        "requested_model": "deepseek-v4-flash",
        "observed_model": "deepseek-v4-flash",
        "request_sha256": "a" * 64,
        "tool_schema_sha256": "b" * 64,
        "input_tokens": 10,
        "output_tokens": 4,
        "reasoning_tokens": 3,
        "finish_reason": "stop",
        "tool_calls_present": False,
        "reasoning_continuation_present": True,
        "private_reasoning_persisted": False,
    }

    receipt = ProviderTurnReceiptV1(
        **body, receipt_sha256=canonical_sha256(body)
    )

    assert receipt.receipt_sha256 == canonical_sha256(body)
    assert "transport_deadlines" not in receipt.__dict__
