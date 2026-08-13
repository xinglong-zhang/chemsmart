import pytest

from chemsmart.agent._contracts import ContractError
from chemsmart.agent.runtime.deepseek import (
    DeepSeekProtocolError,
    DeepSeekV4ToolSession,
)
from chemsmart.agent.runtime.event_store import RuntimeEventStore
from chemsmart.agent.runtime.events import EventKind
from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1


def _tool_response(*, call_id="call-1", arguments="{}", content=None):
    return {
        "model": "deepseek-v4-flash",
        "choices": [
            {
                "finish_reason": "tool_calls",
                "message": {
                    "role": "assistant",
                    "content": content,
                    "reasoning_content": "private continuation",
                    "tool_calls": [
                        {
                            "id": call_id,
                            "type": "function",
                            "function": {
                                "name": "inspect_program_capability",
                                "arguments": arguments,
                            },
                        }
                    ],
                },
            }
        ],
        "usage": {"prompt_tokens": 1, "completion_tokens": 1},
    }


def test_null_content_is_normalized_only_for_private_tool_continuation():
    observed_requests = []

    def transport(payload):
        observed_requests.append(payload)
        return _tool_response()

    session = DeepSeekV4ToolSession(
        transport=transport,
        messages=[{"role": "user", "content": "Inspect capability."}],
    )
    session.turn(tools=[])
    assert session.public_history()[-1]["content"] == ""
    assert "reasoning_content" not in session.public_history()[-1]
    assert "tool_choice" not in observed_requests[0]


def test_mismatched_tool_result_id_is_rejected():
    session = DeepSeekV4ToolSession(
        transport=lambda _: _tool_response(),
        messages=[{"role": "user", "content": "Inspect capability."}],
    )
    session.turn(tools=[])
    with pytest.raises(ContractError, match="exactly match"):
        session.append_tool_results(
            [{"role": "tool", "tool_call_id": "wrong", "content": "{}"}]
        )


def test_duplicate_tool_call_id_in_one_turn_is_rejected():
    response = _tool_response()
    duplicate = dict(response["choices"][0]["message"]["tool_calls"][0])
    response["choices"][0]["message"]["tool_calls"].append(duplicate)
    session = DeepSeekV4ToolSession(
        transport=lambda _: response,
        messages=[{"role": "user", "content": "Inspect capability."}],
    )
    with pytest.raises(DeepSeekProtocolError, match="reused"):
        session.turn(tools=[])


def test_conflicting_tool_call_id_reuse_across_turns_is_rejected():
    responses = iter([_tool_response(), _tool_response()])
    session = DeepSeekV4ToolSession(
        transport=lambda _: next(responses),
        messages=[{"role": "user", "content": "Inspect capability."}],
    )
    session.turn(tools=[])
    session.append_tool_results(
        [{"role": "tool", "tool_call_id": "call-1", "content": "{}"}]
    )
    with pytest.raises(DeepSeekProtocolError, match="reused"):
        session.turn(tools=[])


def test_malformed_tool_json_is_rejected_before_dispatch():
    session = DeepSeekV4ToolSession(
        transport=lambda _: _tool_response(arguments="{"),
        messages=[{"role": "user", "content": "Inspect capability."}],
    )
    with pytest.raises(DeepSeekProtocolError, match="invalid JSON") as observed:
        session.turn(tools=[])
    error = observed.value
    assert error.failing_field == (
        "choices[0].message.tool_calls[0].function.arguments"
    )
    assert error.json_offset == 1
    assert len(error.response_envelope_sha256) == 64
    assert error.response_envelope_bytes > 0
    assert error.public_observation()["retry_decision"] == "not_retried"
    assert error.public_observation()["recovery_requirement"] == (
        "new_independent_attempt"
    )
    assert session.receipts == ()
    assert session.outstanding_tool_call_ids == ()


def test_unknown_tool_is_absent_from_host_dispatch(tmp_path):
    store = RuntimeEventStore(tmp_path / "events" / "runtime.jsonl", session_id="s1")
    host = CommandCompiledToolHostV1(event_store=store)
    with pytest.raises(ContractError, match="not exposed"):
        host.dispatch(turn_id="t1", tool_name="run_local", arguments={})


def test_terminal_complete_rejects_red_preflight(tmp_path):
    store = RuntimeEventStore(tmp_path / "events" / "runtime.jsonl", session_id="s1")
    receipt = "a" * 64
    store.append(
        turn_id="t1",
        kind=EventKind.PROGRAM_PREFLIGHTED.value,
        payload={
            "receipt_sha256": receipt,
            "plan_state": "blocked",
            "critical_finding_count": 1,
            "safe_preview_receipt_sha256": "",
        },
    )
    with pytest.raises(ContractError, match="gate is red"):
        store.terminate(
            turn_id="t1",
            terminal_state="complete",
            reason="model claimed completion",
            required_receipt_sha256s=(receipt,),
        )
