from __future__ import annotations

from chemsmart.agent.core import DecisionLog
from chemsmart.agent.handles import HandleStore
from chemsmart.agent.loop import ToolLoop, ToolLoopBudgets

from ._agent_session_helpers import FakeProvider
from ._loop_helpers import (
    ScriptedRegistry,
    openai_final_response,
    openai_tool_call_response,
    tool_call,
)


def test_tool_loop_writes_new_decision_log_kinds_in_order(tmp_path):
    provider = FakeProvider(
        [
            {
                "__raw_response__": openai_tool_call_response(
                    tool_call("call_1", "recommend_method", {"task": "opt"})
                )
            },
            {"__raw_response__": openai_final_response("Done.")},
        ]
    )
    decision_log = DecisionLog(tmp_path / "decision_log.jsonl")
    loop = ToolLoop(
        provider=provider,
        registry=ScriptedRegistry({"recommend_method": {"method": "b3lyp"}}),
        handle_store=HandleStore(tmp_path),
        decision_log=decision_log,
        budgets=ToolLoopBudgets(log_provider_turn_raw=True),
    )

    loop.run_turn(
        messages=[{"role": "user", "content": "Recommend."}],
        tool_defs=[
            {
                "type": "function",
                "function": {
                    "name": "recommend_method",
                    "description": "recommend_method",
                    "parameters": {"type": "object", "properties": {}},
                },
            }
        ],
    )

    kinds = [entry["kind"] for entry in decision_log.read_all()]
    assert kinds == [
        "provider_turn_raw",
        "assistant_turn",
        "tool_use_request",
        "tool_use_result",
        "provider_turn_raw",
        "assistant_turn",
    ]
