from __future__ import annotations

import io
import json
from types import SimpleNamespace

import pytest

from chemsmart.agent._contracts import ContractError
from chemsmart.agent.api_access import load_secret_lease
from chemsmart.agent.runtime.alibaba import (
    AlibabaTokenPlanConfigV1,
    AlibabaTokenPlanToolSession,
    Qwen38MaxConfigV1,
    Qwen38MaxToolSession,
    _assemble_sse_response,
)
from chemsmart.agent.services.unified_session import UnifiedSessionRunner


def _sse(*chunks):
    lines = []
    for chunk in chunks:
        lines.append(b"data: " + json.dumps(chunk).encode("utf-8") + b"\n\n")
    lines.append(b"data: [DONE]\n\n")
    return io.BytesIO(b"".join(lines))


def test_qwen_payload_uses_xhigh_and_omits_conflicting_limits():
    session = Qwen38MaxToolSession(
        transport=lambda payload: payload,
        messages=[{"role": "user", "content": "Plan."}],
    )

    payload = session.request_payload(tools=[{"type": "function"}])

    assert payload["model"] == "qwen3.8-max"
    assert payload["reasoning_effort"] == "xhigh"
    assert payload["preserve_thinking"] is True
    assert payload["stream"] is True
    assert payload["stream_options"] == {"include_usage": True}
    assert "thinking_budget" not in payload
    assert "thinking" not in payload
    assert "max_tokens" not in payload
    assert "max_completion_tokens" not in payload


def test_token_plan_payload_uses_profile_selected_deepseek_model():
    session = AlibabaTokenPlanToolSession(
        transport=lambda payload: payload,
        messages=[{"role": "user", "content": "Plan."}],
        config=AlibabaTokenPlanConfigV1(
            model="deepseek-v4-flash-0731",
            reasoning_effort="max",
        ),
    )

    payload = session.request_payload(tools=[{"type": "function"}])

    assert payload["model"] == "deepseek-v4-flash-0731"
    assert payload["reasoning_effort"] == "max"
    assert payload["preserve_thinking"] is True


def test_qwen_config_rejects_preview_and_reduced_reasoning():
    with pytest.raises(ContractError, match="production qwen3.8-max"):
        Qwen38MaxConfigV1(model="qwen3.8-max-preview")
    with pytest.raises(ContractError, match="must be xhigh"):
        Qwen38MaxConfigV1(reasoning_effort="high")
    with pytest.raises(ContractError, match="separately authorized attempt"):
        Qwen38MaxConfigV1(sdk_max_retries=1)


def test_qwen_stream_rejects_an_unexpected_observed_model():
    with pytest.raises(ContractError, match="production qwen3.8-max"):
        # Config validation catches a requested preview alias before transport.
        Qwen38MaxConfigV1(model="qwen3.8-max-preview")
    from chemsmart.agent.runtime.deepseek import DeepSeekProtocolError

    with pytest.raises(DeepSeekProtocolError, match="did not confirm"):
        _assemble_sse_response(
            _sse(
                {
                    "id": "chat-preview",
                    "model": "qwen3.8-max-preview",
                    "choices": [
                        {
                            "delta": {"role": "assistant", "content": "x"},
                            "finish_reason": "stop",
                        }
                    ],
                }
            )
        )


def test_qwen_stream_rejects_malformed_tool_index_without_echoing_value():
    from chemsmart.agent.runtime.deepseek import DeepSeekProtocolError

    secret_index = "sk-sp-provider-authored-secret"
    with pytest.raises(
        DeepSeekProtocolError, match="non-negative integer"
    ) as observed:
        _assemble_sse_response(
            _sse(
                {
                    "id": "chat-malformed-index",
                    "model": "qwen3.8-max",
                    "choices": [
                        {
                            "delta": {
                                "role": "assistant",
                                "tool_calls": [
                                    {
                                        "index": secret_index,
                                        "id": "call-1",
                                        "type": "function",
                                        "function": {
                                            "name": "inspect_program_capability",
                                            "arguments": "{}",
                                        },
                                    }
                                ],
                            },
                            "finish_reason": "tool_calls",
                        }
                    ],
                }
            )
        )

    assert secret_index not in str(observed.value)
    assert observed.value.failing_field == (
        "choices[*].delta.tool_calls[*].index"
    )


def test_qwen_sse_assembles_reasoning_and_fragmented_tool_call():
    response = _assemble_sse_response(
        _sse(
            {
                "id": "chat-1",
                "model": "qwen3.8-max",
                "choices": [
                    {
                        "delta": {
                            "role": "assistant",
                            "reasoning_content": "private ",
                            "tool_calls": [
                                {
                                    "index": 0,
                                    "id": "call-1",
                                    "type": "function",
                                    "function": {
                                        "name": "inspect_program_",
                                        "arguments": '{"program":',
                                    },
                                }
                            ],
                        },
                        "finish_reason": None,
                    }
                ],
            },
            {
                "id": "chat-1",
                "model": "qwen3.8-max",
                "choices": [
                    {
                        "delta": {
                            "reasoning_content": "reasoning",
                            "tool_calls": [
                                {
                                    "index": 0,
                                    "function": {
                                        "name": "capability",
                                        "arguments": '"pyscf"}',
                                    },
                                }
                            ],
                        },
                        "finish_reason": "tool_calls",
                    }
                ],
            },
            {
                "id": "chat-1",
                "model": "qwen3.8-max",
                "choices": [],
                "usage": {
                    "prompt_tokens": 10,
                    "completion_tokens": 7,
                    "completion_tokens_details": {"reasoning_tokens": 5},
                },
            },
        )
    )

    message = response["choices"][0]["message"]
    assert message["reasoning_content"] == "private reasoning"
    assert message["tool_calls"][0]["id"] == "call-1"
    assert message["tool_calls"][0]["function"] == {
        "name": "inspect_program_capability",
        "arguments": '{"program":"pyscf"}',
    }
    assert response["usage"]["completion_tokens"] == 7


def test_qwen_private_reasoning_replays_in_ram_but_not_publicly():
    requests = []
    responses = iter(
        [
            {
                "model": "qwen3.8-max",
                "choices": [
                    {
                        "finish_reason": "tool_calls",
                        "message": {
                            "role": "assistant",
                            "content": "Inspecting.",
                            "reasoning_content": "private continuation",
                            "tool_calls": [
                                {
                                    "id": "call-1",
                                    "type": "function",
                                    "function": {
                                        "name": "inspect_program_capability",
                                        "arguments": '{"program":"pyscf"}',
                                    },
                                }
                            ],
                        },
                    }
                ],
                "usage": {
                    "prompt_tokens": 10,
                    "completion_tokens": 7,
                    "completion_tokens_details": {"reasoning_tokens": 5},
                },
            },
            {
                "model": "qwen3.8-max",
                "choices": [
                    {
                        "finish_reason": "stop",
                        "message": {
                            "role": "assistant",
                            "content": "Done.",
                            "reasoning_content": "private final",
                        },
                    }
                ],
                "usage": {"prompt_tokens": 20, "completion_tokens": 3},
            },
        ]
    )

    def transport(payload):
        requests.append(payload)
        return next(responses)

    session = Qwen38MaxToolSession(
        transport=transport,
        messages=[{"role": "user", "content": "Inspect."}],
    )
    first, receipt = session.turn(tools=[{"type": "function"}])
    session.append_tool_results(
        [
            {
                "role": "tool",
                "tool_call_id": "call-1",
                "content": '{"status":"supported"}',
            }
        ]
    )
    second, _ = session.turn(tools=[{"type": "function"}])

    assert requests[1]["messages"][1]["reasoning_content"] == (
        "private continuation"
    )
    assert "reasoning_content" not in json.dumps(first)
    assert "reasoning_content" not in json.dumps(second)
    assert receipt.provider == "alibaba-token-plan"
    assert receipt.requested_model == "qwen3.8-max"
    session.close()
    assert "private" not in json.dumps(session.public_history())


def test_unified_runner_selects_qwen_and_passes_feedback_projection(
    tmp_path, monkeypatch
):
    import chemsmart.agent.runtime.alibaba as alibaba
    import chemsmart.agent.services.unified_session as unified

    secret_file = tmp_path / "api.env"
    secret_file.write_text("ALIBABA_TOKEN_PLAN_KEY=sk-sp-test\n")
    observed = {}

    class FakeTransport:
        def __init__(self, **kwargs):
            observed["transport"] = kwargs

        def close(self):
            observed["transport_closed"] = True

    class FakeSession:
        def __init__(self, **kwargs):
            self.config = kwargs["config"]
            observed["session"] = kwargs

        def close(self):
            observed["session_closed"] = True

    class FakeLoop:
        def __init__(self, **kwargs):
            observed["loop_init"] = kwargs

        def run(self, **kwargs):
            observed["loop_run"] = kwargs
            return "normal-session-result"

    monkeypatch.setattr(alibaba, "AlibabaTokenPlanHttpsTransport", FakeTransport)
    monkeypatch.setattr(
        alibaba, "AlibabaTokenPlanToolSession", FakeSession
    )
    monkeypatch.setattr(unified, "ToolLoopRunner", FakeLoop)
    lease = load_secret_lease(
        provider="alibaba-token-plan", path=secret_file
    )
    runner = UnifiedSessionRunner(
        host=object(),
        event_store=object(),
        credential_lease=lease,
        provider_config=Qwen38MaxConfigV1(),
    )

    result = runner.run(
        messages=[{"role": "user", "content": "Plan."}],
        envelope=SimpleNamespace(
            budget=SimpleNamespace(max_output_tokens_per_request=262_144)
        ),
        hypothesis=object(),
        network_budget=SimpleNamespace(
            max_output_tokens_per_request=262_144,
            task_wall_time_seconds=60.0,
        ),
        feedback_projection="causal-v1",
    )

    assert result == "normal-session-result"
    assert observed["session"]["config"].model == "qwen3.8-max"
    assert observed["loop_init"]["feedback_projection"] == "causal-v1"
    assert observed["transport_closed"] is True
    assert observed["session_closed"] is True


def test_qwen_session_propagates_remaining_turn_timeout():
    observed = []

    class TimeoutAwareTransport:
        def set_timeout_seconds(self, value):
            observed.append(value)

        def __call__(self, payload):
            return {
                "id": "response-timeout",
                "model": "qwen3.8-max",
                "choices": [
                    {
                        "finish_reason": "stop",
                        "message": {
                            "role": "assistant",
                            "content": "Observed.",
                        },
                    }
                ],
                "usage": {},
            }

    session = Qwen38MaxToolSession(
        transport=TimeoutAwareTransport(),
        messages=[{"role": "user", "content": "Inspect."}],
    )
    session.set_turn_timeout_seconds(17.5)
    session.turn()

    assert observed == [17.5]
