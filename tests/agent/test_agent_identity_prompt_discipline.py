from __future__ import annotations

from typing import Any

from chemsmart.agent.core import AgentSession
from chemsmart.agent.prompts.identity import build_system_prompt
from chemsmart.agent.registry import ToolRegistry

from ._loop_helpers import openai_final_response


class ToolResultDisciplineProvider:
    def __init__(self) -> None:
        self.calls: list[dict[str, Any]] = []
        self.name = "openai"
        self.default_model = "gpt-5.4-mock"

    def chat(self, messages, tools=None, timeout_s=30):
        self.calls.append(
            {"messages": messages, "tools": tools, "timeout_s": timeout_s}
        )
        system_prompt = messages[0]["content"]
        if "Never summarize a tool result beyond its literal content." in (
            system_prompt
        ):
            text = "0 queued, 2 running"
        else:
            text = "queue almost full"
        return openai_final_response(text)


class ScopeDisciplineProvider:
    def __init__(self) -> None:
        self.calls: list[dict[str, Any]] = []
        self.name = "openai"
        self.default_model = "gpt-5.4-mock"

    def chat(self, messages, tools=None, timeout_s=30):
        self.calls.append(
            {"messages": messages, "tools": tools, "timeout_s": timeout_s}
        )
        system_prompt = messages[0]["content"]
        if "That is outside my scope as the chemsmart agent." in system_prompt:
            text = (
                "That is outside my scope as the chemsmart agent. "
                "I can help with computational chemistry and HPC workflows."
            )
        else:
            text = "try kimchi jjigae"
        return openai_final_response(text)


def test_build_system_prompt_contains_tool_result_discipline_literal():
    prompt = build_system_prompt(
        registry=ToolRegistry.default(),
        stage_instructions="Answer briefly.",
    )

    assert "summarize a tool result beyond" in prompt


def test_build_system_prompt_contains_scope_discipline_refusal_literal():
    prompt = build_system_prompt(
        registry=ToolRegistry.default(),
        stage_instructions="Answer briefly.",
    )

    assert "outside my scope as the chemsmart agent" in prompt


def test_run_loop_tool_result_discipline_uses_literal_queue_summary(tmp_path):
    provider = ToolResultDisciplineProvider()
    session = AgentSession(
        provider=provider,
        registry=ToolRegistry.default(),
        session_root=tmp_path,
    )

    result = session.run_loop(
        "Summarize this tool result: Queued:0, Running:2"
    )

    assert result["assistant_output"] == "0 queued, 2 running"


def test_run_loop_scope_discipline_refuses_off_topic_chitchat(tmp_path):
    provider = ScopeDisciplineProvider()
    session = AgentSession(
        provider=provider,
        registry=ToolRegistry.default(),
        session_root=tmp_path,
    )

    result = session.run_loop("What should I eat for lunch tomorrow?")

    assert result["assistant_output"] == (
        "That is outside my scope as the chemsmart agent. "
        "I can help with computational chemistry and HPC workflows."
    )
