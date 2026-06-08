from __future__ import annotations

import pytest

from chemsmart.agent.core import DecisionLog
from chemsmart.agent.handles import HandleStore
from chemsmart.agent.loop import ToolLoop, ToolLoopBudgets

from ._agent_session_helpers import FakeProvider
from ._loop_helpers import (
    ScriptedRegistry,
    openai_tool_call_response,
    tool_call,
)


@pytest.mark.parametrize(
    (
        "budgets",
        "responses",
        "registry",
        "expected_reason",
        "expected_statuses",
    ),
    [
        (
            ToolLoopBudgets(max_model_steps_per_turn=1),
            [
                {
                    "__raw_response__": openai_tool_call_response(
                        tool_call(
                            "call_1", "recommend_method", {"task": "opt"}
                        )
                    )
                },
                {
                    "__raw_response__": openai_tool_call_response(
                        tool_call(
                            "call_2", "recommend_method", {"task": "freq"}
                        )
                    )
                },
            ],
            ScriptedRegistry(
                {
                    "recommend_method": {"method": "b3lyp"},
                }
            ),
            "max_model_steps",
            ["ok"],
        ),
        (
            ToolLoopBudgets(max_total_tool_calls_per_turn=1),
            [
                {
                    "__raw_response__": openai_tool_call_response(
                        tool_call(
                            "call_1", "recommend_method", {"task": "opt"}
                        ),
                        tool_call(
                            "call_2", "recommend_method", {"task": "freq"}
                        ),
                    )
                },
            ],
            ScriptedRegistry(
                {
                    "recommend_method": {"method": "b3lyp"},
                }
            ),
            "max_tool_calls",
            ["ok", "skipped"],
        ),
        (
            ToolLoopBudgets(max_consecutive_tool_errors=2),
            [
                {
                    "__raw_response__": openai_tool_call_response(
                        tool_call(
                            "call_1", "recommend_method", {"task": "opt"}
                        )
                    )
                },
                {
                    "__raw_response__": openai_tool_call_response(
                        tool_call(
                            "call_2", "recommend_method", {"task": "freq"}
                        )
                    )
                },
            ],
            ScriptedRegistry(
                {
                    "recommend_method": {
                        "ok": False,
                        "error": {
                            "type": "RuntimeError",
                            "message": "boom",
                            "tool": "recommend_method",
                        },
                    },
                }
            ),
            "max_consecutive_errors",
            ["error", "error"],
        ),
        (
            ToolLoopBudgets(max_same_signature_retries=2),
            [
                {
                    "__raw_response__": openai_tool_call_response(
                        tool_call(
                            "call_1", "recommend_method", {"task": "opt"}
                        )
                    )
                },
                {
                    "__raw_response__": openai_tool_call_response(
                        tool_call(
                            "call_2", "recommend_method", {"task": "opt"}
                        )
                    )
                },
                {
                    "__raw_response__": openai_tool_call_response(
                        tool_call(
                            "call_3", "recommend_method", {"task": "opt"}
                        )
                    )
                },
            ],
            ScriptedRegistry(
                {
                    "recommend_method": {"method": "b3lyp"},
                }
            ),
            "repeat_signature",
            ["ok", "ok", "skipped"],
        ),
    ],
)
def test_tool_loop_stops_cleanly_on_budget_limits(
    tmp_path,
    budgets,
    responses,
    registry,
    expected_reason,
    expected_statuses,
):
    provider = FakeProvider(responses)
    loop = ToolLoop(
        provider=provider,
        registry=registry,
        handle_store=HandleStore(tmp_path),
        decision_log=DecisionLog(tmp_path / "decision_log.jsonl"),
        budgets=budgets,
    )

    result = loop.run_turn(
        messages=[{"role": "user", "content": "Run the loop."}],
        tool_defs=registry.openai_tool_defs(),
    )

    assert result["limit_reason"] == expected_reason
    assert [
        outcome.status for outcome in result["tool_outcomes"]
    ] == expected_statuses
    entries = loop.decision_log.read_all()
    assert entries[-1]["kind"] == "loop_limit_exceeded"
    assert entries[-1]["payload"]["limit_reason"] == expected_reason
