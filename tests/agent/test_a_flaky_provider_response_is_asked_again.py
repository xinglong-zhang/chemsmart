"""One malformed provider response must not end an otherwise-live session.

A real cycle-040 session previewed its entire four-node workflow -- both
platform arms green -- and then died one turn before its decision record
because the provider emitted malformed tool-call JSON once. The adapter's
own classification of that failure says `new_independent_attempt`, and the
adapter validates all-or-nothing before any history append, so the failed
attempt left the session exactly as it was. The loop now implements the
recovery it was already naming: bounded independent re-asks, every attempt
kept in the event trail, termination only when failure is persistent.
"""

from __future__ import annotations

from chemsmart.agent.loop import ToolLoopRunner
from chemsmart.agent.runtime.alibaba import (
    Qwen38MaxConfigV1,
    Qwen38MaxToolSession,
)
from chemsmart.agent.runtime.event_store import RuntimeEventStore
from chemsmart.agent.runtime.events import EventKind
from tests.agent.test_provider_protocol_failure_evidence import (
    _DispatchSpyHost,
    _run_contracts,
)

_MALFORMED = {
    "id": "flaky-attempt",
    "model": "qwen3.8-max",
    "choices": [
        {
            "finish_reason": "tool_calls",
            "message": {
                "role": "assistant",
                "content": "",
                "reasoning_content": "PRIVATE",
                "tool_calls": [
                    {
                        "id": "broken-call",
                        "type": "function",
                        "function": {
                            "name": "inspect_program_capability",
                            "arguments": '{"program":',
                        },
                    }
                ],
            },
        }
    ],
    "usage": {"prompt_tokens": 10, "completion_tokens": 5},
}

_FINAL = {
    "id": "recovered-attempt",
    "model": "qwen3.8-max",
    "choices": [
        {
            "finish_reason": "stop",
            "message": {
                "role": "assistant",
                "content": "The planning summary stands.",
                "reasoning_content": "",
            },
        }
    ],
    "usage": {"prompt_tokens": 12, "completion_tokens": 6},
}


def _flaky_once_transport():
    responses = iter((_MALFORMED, _FINAL))
    return lambda _payload: next(responses)


def test_one_malformed_response_is_retried_and_the_session_completes(
    tmp_path,
):
    config = Qwen38MaxConfigV1()
    session = Qwen38MaxToolSession(
        transport=_flaky_once_transport(),
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
    )

    assert result.terminal_state not in ("failed",)
    assert result.final_text == "The planning summary stands."

    attempts = [
        event
        for event in store.read_events()
        if event.kind == EventKind.API_ATTEMPT_OBSERVED.value
    ]
    assert [a.payload["status"] for a in attempts] == [
        "protocol_failed",
        "succeeded",
    ]
    failure = attempts[0].payload["protocol_failure"]
    assert failure["retry_decision"] == "retried_independently"
    assert failure["recovery_requirement"] == "new_independent_attempt"
