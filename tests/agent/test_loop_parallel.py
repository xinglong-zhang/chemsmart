from __future__ import annotations

from chemsmart.agent.core import DecisionLog
from chemsmart.agent.handles import HandleStore
from chemsmart.agent.loop import ToolLoop

from ._agent_session_helpers import FakeProvider
from ._loop_helpers import (
    ScriptedRegistry,
    openai_final_response,
    openai_tool_call_response,
    tool_call,
)


def test_tool_loop_executes_parallel_tool_requests_in_one_step(tmp_path):
    provider = FakeProvider(
        [
            {
                "__raw_response__": openai_tool_call_response(
                    tool_call("call_1", "recommend_method", {"task": "opt"}),
                    tool_call(
                        "call_2", "validate_runtime", {"job": "job_target"}
                    ),
                )
            },
            {"__raw_response__": openai_final_response("Ready.")},
        ]
    )
    registry = ScriptedRegistry(
        {
            "recommend_method": {"method": "wb97xd"},
            "validate_runtime": {"ok": "ok", "local_issues": []},
        }
    )
    loop = ToolLoop(
        provider=provider,
        registry=registry,
        handle_store=HandleStore(tmp_path),
        decision_log=DecisionLog(tmp_path / "decision_log.jsonl"),
    )

    result = loop.run_turn(
        messages=[{"role": "user", "content": "Prepare and validate."}],
        tool_defs=registry.openai_tool_defs(),
    )

    assert result["assistant_text"] == "Ready."
    assert [
        request.provider_call_id for request in result["tool_requests"]
    ] == [
        "call_1",
        "call_2",
    ]
    assert [outcome.status for outcome in result["tool_outcomes"]] == [
        "ok",
        "ok",
    ]
    assert registry.calls == [
        ("recommend_method", {"task": "opt"}),
        ("validate_runtime", {"job": "job_target"}),
    ]
    second_call_messages = provider.calls[1]["messages"]
    assert [message["role"] for message in second_call_messages].count(
        "tool"
    ) == 2
