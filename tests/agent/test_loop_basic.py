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


def test_tool_loop_runs_single_tool_then_final_assistant(tmp_path):
    provider = FakeProvider(
        [
            {
                "__raw_response__": openai_tool_call_response(
                    tool_call("call_1", "recommend_method", {"task": "opt"}),
                    content="I'll recommend a method first.",
                )
            },
            {"__raw_response__": openai_final_response("Use B3LYP/6-31G*.")},
        ]
    )
    registry = ScriptedRegistry(
        {
            "recommend_method": {
                "method": "B3LYP/6-31G*",
                "why": "Balanced default",
            }
        }
    )
    loop = ToolLoop(
        provider=provider,
        registry=registry,
        handle_store=HandleStore(tmp_path),
        decision_log=DecisionLog(tmp_path / "decision_log.jsonl"),
    )

    result = loop.run_turn(
        messages=[{"role": "user", "content": "Recommend a method."}],
        tool_defs=registry.openai_tool_defs(),
    )

    assert result["assistant_text"] == "Use B3LYP/6-31G*."
    assert result["stop_reason"] == "stop"
    assert result["model_steps"] == 2
    assert result["limit_reason"] is None
    assert [request.name for request in result["tool_requests"]] == [
        "recommend_method"
    ]
    assert [outcome.status for outcome in result["tool_outcomes"]] == ["ok"]
    assert registry.calls == [
        ("recommend_method", {"task": "opt"}),
    ]
    assert provider.calls[0]["tools"] == registry.openai_tool_defs()
    assert any(
        message["role"] == "tool" for message in provider.calls[1]["messages"]
    )
