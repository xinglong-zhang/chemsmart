from __future__ import annotations

import pytest

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


class FlakyProvider(FakeProvider):
    def __init__(self, failures, responses):
        super().__init__(responses)
        self.failures = failures

    def chat(self, messages, tools=None, timeout_s=30):
        if self.failures:
            self.failures -= 1
            raise TimeoutError("provider timeout")
        return super().chat(messages, tools=tools, timeout_s=timeout_s)


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


def test_tool_loop_retries_one_provider_timeout_then_recovers(tmp_path):
    provider = FlakyProvider(
        1,
        [{"__raw_response__": openai_final_response("Recovered.")}],
    )
    registry = ScriptedRegistry({"recommend_method": {"method": "b3lyp"}})
    loop = ToolLoop(
        provider=provider,
        registry=registry,
        handle_store=HandleStore(tmp_path),
        decision_log=DecisionLog(tmp_path / "decision_log.jsonl"),
    )

    result = loop.run_turn(
        messages=[{"role": "user", "content": "Try once more."}],
        tool_defs=registry.openai_tool_defs(),
    )

    assert result["assistant_text"] == "Recovered."
    assert result["provider_errors"] == 1
    assert result["limit_reason"] is None
    assert any(
        entry["kind"] == "provider_turn_error"
        for entry in loop.decision_log.read_all()
    )


def test_tool_loop_stops_after_provider_error_budget(tmp_path):
    provider = FlakyProvider(3, [])
    registry = ScriptedRegistry({"recommend_method": {"method": "b3lyp"}})
    loop = ToolLoop(
        provider=provider,
        registry=registry,
        handle_store=HandleStore(tmp_path),
        decision_log=DecisionLog(tmp_path / "decision_log.jsonl"),
        budgets=ToolLoopBudgets(max_provider_errors_per_turn=2),
    )

    result = loop.run_turn(
        messages=[{"role": "user", "content": "Try twice."}],
        tool_defs=registry.openai_tool_defs(),
    )

    assert result["provider_errors"] == 2
    assert result["limit_reason"] == "provider_errors"
