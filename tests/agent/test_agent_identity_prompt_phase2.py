from __future__ import annotations

from typing import Any

from chemsmart.agent.core import AgentSession
from chemsmart.agent.prompts.identity import build_system_prompt
from chemsmart.agent.registry import ToolRegistry

from ._loop_helpers import openai_tool_call_response, tool_call


class Phase2IdentityAwareProvider:
    def __init__(self) -> None:
        self.calls: list[dict[str, Any]] = []
        self.name = "openai"
        self.default_model = "gpt-5.4-mock"

    def chat(self, messages, tools=None, timeout_s=30):
        self.calls.append(
            {"messages": messages, "tools": tools, "timeout_s": timeout_s}
        )
        system_prompt = messages[0]["content"]
        if (
            "STRUCTURED SLOT" in system_prompt
            and "ask_user" in system_prompt
            and "FORBIDDEN" in system_prompt
        ):
            return openai_tool_call_response(
                tool_call(
                    "call_ask",
                    "ask_user",
                    {
                        "question": "어느 서버?",
                        "options": ["chemnode1", "chemnode2"],
                    },
                )
            )
        return {
            "choices": [
                {
                    "message": {
                        "role": "assistant",
                        "content": "Which server do you mean?",
                    },
                    "finish_reason": "stop",
                }
            ],
            "usage": {"prompt_tokens": 8, "completion_tokens": 4},
        }


def test_build_system_prompt_contains_structured_slot_ask_user_guardrails():
    prompt = build_system_prompt(
        registry=ToolRegistry.default(),
        stage_instructions="Answer briefly.",
    )

    assert "ask_user tool" in prompt
    assert "FORBIDDEN" in prompt


def test_build_system_prompt_contains_advisory_guardrails():
    prompt = build_system_prompt(
        registry=ToolRegistry.default(),
        stage_instructions="Answer briefly.",
    )

    assert "ADVISORY" in prompt
    assert "walltime" in prompt


def test_build_system_prompt_contains_remote_path_precedence_guardrails():
    prompt = build_system_prompt(
        registry=ToolRegistry.default(),
        stage_instructions="Answer briefly.",
    )

    assert "prefer `log_tail(server=last_server" in prompt


def test_run_loop_uses_ask_user_for_structured_slot_ambiguity(tmp_path):
    provider = Phase2IdentityAwareProvider()
    session = AgentSession(
        provider=provider,
        registry=ToolRegistry.default(),
        session_root=tmp_path,
    )

    result = session.run_loop("서버 상태 봐줘")

    assert result["ask_user_question"] == {
        "question": "어느 서버?",
        "options": ["chemnode1", "chemnode2"],
    }
    assert result["assistant_output"] == ""
