"""The OpenAI adapter speaks plain chat completions -- no vendor keys.

The shared session machinery (history, receipts, tool continuation,
reasoning scrubbing) is inherited; what this adapter owns is the payload
shape and the endpoint binding.
"""

from __future__ import annotations

import json

import pytest

from chemsmart.agent._contracts import ContractError
from chemsmart.agent.runtime.openai_compat import (
    OPENAI_OFFICIAL_ENDPOINT,
    OpenAICompatibleConfigV1,
    OpenAICompatibleToolSession,
    OpenAIHttpsTransport,
)

_CONFIG = dict(
    model="gpt-5.2",
    context_tokens=400000,
    max_output_tokens=64000,
)


def test_the_payload_carries_no_vendor_extensions():
    session = OpenAICompatibleToolSession(
        transport=lambda payload: payload,
        messages=[{"role": "user", "content": "Plan."}],
        config=OpenAICompatibleConfigV1(**_CONFIG),
    )

    payload = session.request_payload(tools=[{"type": "function"}])

    assert payload["model"] == "gpt-5.2"
    assert payload["max_tokens"] == 64000
    assert "thinking" not in payload
    assert "preserve_thinking" not in payload
    assert "reasoning_effort" not in payload, "omitted when unset"


def test_a_set_effort_is_forwarded_and_vocabulary_is_enforced():
    session = OpenAICompatibleToolSession(
        transport=lambda payload: payload,
        messages=[{"role": "user", "content": "Plan."}],
        config=OpenAICompatibleConfigV1(**_CONFIG, reasoning_effort="high"),
    )

    assert session.request_payload()["reasoning_effort"] == "high"

    with pytest.raises(ContractError, match="low, medium, high"):
        OpenAICompatibleConfigV1(**_CONFIG, reasoning_effort="xhigh")


def test_the_transport_binds_the_official_endpoint_only():
    transport = OpenAIHttpsTransport(
        api_key="sk-test", endpoint=OPENAI_OFFICIAL_ENDPOINT
    )
    assert transport.endpoint == OPENAI_OFFICIAL_ENDPOINT

    with pytest.raises(ContractError, match="official endpoint"):
        OpenAIHttpsTransport(
            api_key="sk-test", endpoint="https://proxy.example/v1"
        )
    with pytest.raises(ContractError, match="official endpoint"):
        OpenAICompatibleConfigV1(
            **_CONFIG, endpoint="https://proxy.example/v1"
        )


def test_a_canned_turn_round_trips_through_the_shared_machinery():
    def transport(payload):
        assert "thinking" not in payload
        return {
            "id": "turn-1",
            "model": "gpt-5.2",
            "choices": [
                {
                    "finish_reason": "stop",
                    "message": {
                        "role": "assistant",
                        "content": "The plan stands.",
                    },
                }
            ],
            "usage": {"prompt_tokens": 4, "completion_tokens": 3},
        }

    session = OpenAICompatibleToolSession(
        transport=transport,
        messages=[{"role": "user", "content": "Plan."}],
        config=OpenAICompatibleConfigV1(**_CONFIG),
    )
    public, receipt = session.turn(tools=[])

    assert public["choices"][0]["message"]["content"] == "The plan stands."
    assert receipt.observed_model == "gpt-5.2"
    assert (
        json.loads(json.dumps(session.public_history()))[-1]["content"]
        == "The plan stands."
    )
