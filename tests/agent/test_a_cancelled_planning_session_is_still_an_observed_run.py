"""Double-esc in the terminal must yield an observed run, not a corpse.

Cancellation is checked only at the provider-turn boundary, so every turn
that happened is fully accounted: attempts and receipts are complete, the
public transcript artifact is persisted, and the event stream terminates
with an explicit ``cancelled`` state. A cancelled session must also never
mint pending authority -- no review packet, no ``waiting_for_approval``.
"""

from __future__ import annotations

import inspect

from chemsmart.agent.loop import ToolLoopRunner
from chemsmart.agent.runtime.alibaba import (
    Qwen38MaxConfigV1,
    Qwen38MaxToolSession,
)
from chemsmart.agent.runtime.event_store import RuntimeEventStore

from tests.agent.test_provider_protocol_failure_evidence import (
    _DispatchSpyHost,
    _run_contracts,
)

_TOOL_TURN = {
    "id": "first-turn",
    "model": "qwen3.8-max",
    "choices": [
        {
            "finish_reason": "tool_calls",
            "message": {
                "role": "assistant",
                "content": "",
                "reasoning_content": "",
                "tool_calls": [
                    {
                        "id": "call-1",
                        "type": "function",
                        "function": {
                            "name": "inspect_program_capability",
                            "arguments": '{"program": "orca"}',
                        },
                    }
                ],
            },
        }
    ],
    "usage": {"prompt_tokens": 10, "completion_tokens": 5},
}


def _terminal_event(store):
    return next(
        event.payload
        for event in reversed(store.read_events())
        if event.kind == "runtime_terminated"
    )


def test_cancellation_before_any_turn_makes_no_provider_request(tmp_path):
    def transport(_payload):
        raise AssertionError("a cancelled session must not call the provider")

    config = Qwen38MaxConfigV1()
    session = Qwen38MaxToolSession(
        transport=transport,
        messages=[{"role": "user", "content": "Plan the workflow."}],
        config=config,
    )
    store = RuntimeEventStore(
        tmp_path / "events" / "runtime.jsonl", session_id="protocol-session"
    )
    host = _DispatchSpyHost()
    envelope, request_context, network = _run_contracts(host, config)

    result = ToolLoopRunner(host=host, event_store=store).run(
        session=session,
        envelope=envelope,
        request_context=request_context,
        provider_budget=network,
        should_stop=lambda: True,
    )

    assert result.terminal_state == "cancelled"
    assert "cancelled by the human" in result.final_text
    assert result.public_transcript_artifact_sha256
    assert store.state().terminal_state == "cancelled"
    assert _terminal_event(store)["terminal_state"] == "cancelled"


def test_cancellation_lands_between_turns_with_the_first_turn_accounted(
    tmp_path,
):
    calls = {"turns": 0}

    def transport(_payload):
        calls["turns"] += 1
        return _TOOL_TURN

    config = Qwen38MaxConfigV1()
    session = Qwen38MaxToolSession(
        transport=transport,
        messages=[{"role": "user", "content": "Plan the workflow."}],
        config=config,
    )
    store = RuntimeEventStore(
        tmp_path / "events" / "runtime.jsonl", session_id="protocol-session"
    )
    host = _DispatchSpyHost()
    envelope, request_context, network = _run_contracts(host, config)

    result = ToolLoopRunner(host=host, event_store=store).run(
        session=session,
        envelope=envelope,
        request_context=request_context,
        provider_budget=network,
        should_stop=lambda: calls["turns"] >= 1,
    )

    assert calls["turns"] == 1
    assert result.terminal_state == "cancelled"
    assert result.successful_tool_calls == 1
    assert store.state().terminal_state == "cancelled"


def test_a_cancelled_session_cannot_become_an_approval_request():
    """The live-session composition refuses to build a review when cancelled.

    Running the full composition needs a provider profile and a bootstrapped
    workspace, so this pins the seam the way the executor gates are pinned:
    by the composition source itself.
    """

    from chemsmart.agent import live_session

    source = inspect.getsource(live_session.run_live_agent_session)
    assert 'loop_result.terminal_state != "cancelled"' in source
    assert source.count("on_run_directory(run_directory)") == 1
    signature = inspect.signature(live_session.run_live_agent_session)
    assert "on_run_directory" in signature.parameters
    assert "should_stop" in signature.parameters
