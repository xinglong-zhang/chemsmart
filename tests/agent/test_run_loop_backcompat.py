from __future__ import annotations

from chemsmart.agent.core import AgentSession

from ._agent_session_helpers import FakeProvider
from ._loop_helpers import (
    ScriptedRegistry,
    openai_final_response,
    openai_tool_call_response,
    tool_call,
)


def test_run_loop_returns_old_and_new_keys_while_run_stays_unchanged(tmp_path):
    registry = ScriptedRegistry(
        {
            "recommend_method": {"method": "b3lyp"},
        }
    )
    loop_provider = FakeProvider(
        [
            {
                "__raw_response__": openai_tool_call_response(
                    tool_call("call_1", "recommend_method", {"task": "opt"})
                )
            },
            {"__raw_response__": openai_final_response("Use B3LYP.")},
        ]
    )
    loop_session = AgentSession(
        provider=loop_provider,
        registry=registry,
        session_root=tmp_path,
    )

    loop_result = loop_session.run_loop("Recommend a method.")

    for key in (
        "session_id",
        "session_dir",
        "completed_steps",
        "blocked",
        "results",
        "assistant_output",
        "tool_requests",
        "tool_outcomes",
        "loop_state",
        "final_message",
        "limit_reason",
        "plan",
        "plan_text",
    ):
        assert key in loop_result
    assert loop_result["assistant_output"] == "Use B3LYP."
    assert loop_result["final_message"] == "Use B3LYP."
    assert loop_result["completed_steps"] == 1
    assert len(loop_result["tool_requests"]) == 1
    assert len(loop_result["tool_outcomes"]) == 1
    assert loop_result["blocked"] is False
    assert loop_result["limit_reason"] is None
    assert loop_result["results"] == [{"method": "b3lyp"}]

    run_provider = FakeProvider(
        [
            {
                "steps": [],
                "rationale": "Answer directly.",
                "estimated_cost": "none",
                "intent": "advisory",
            }
        ]
    )
    run_session = AgentSession(
        provider=run_provider,
        registry=registry,
        session_root=tmp_path / "run_only",
    )

    run_result = run_session.run("Hello")

    assert "assistant_output" not in run_result
    assert "tool_requests" not in run_result
    assert run_result["blocked"] is False
    assert run_result["advisory_only"] is True
