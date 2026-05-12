from __future__ import annotations

from typing import Any

from chemsmart.agent.core import AgentSession
from chemsmart.agent.prompts.identity import build_system_prompt
from chemsmart.agent.registry import ToolRegistry

from ._loop_helpers import openai_final_response


class InstallCommandPolicyProvider:
    def __init__(self) -> None:
        self.calls: list[dict[str, Any]] = []
        self.name = "openai"
        self.default_model = "gpt-5.4-mock"

    def chat(self, messages, tools=None, timeout_s=30) -> dict[str, Any]:
        self.calls.append(
            {"messages": messages, "tools": tools, "timeout_s": timeout_s}
        )
        system_prompt = messages[0]["content"]
        if "Never print shell installation commands in prose" in system_prompt:
            text = "Please follow the project README for installation."
        else:
            text = "pip install chemsmart-agent"
        return openai_final_response(text)


def test_build_system_prompt_contains_install_command_policy_literal():
    prompt = build_system_prompt(
        registry=ToolRegistry.default(),
        stage_instructions="Answer briefly.",
    )

    assert "Never print shell installation commands" in prompt


def test_build_system_prompt_contains_read_only_tool_initiative_literal():
    prompt = build_system_prompt(
        registry=ToolRegistry.default(),
        stage_instructions="Answer briefly.",
    )

    assert "invoke the appropriate read-only tool directly" in prompt


def test_build_system_prompt_contains_entity_slot_type_discipline_literal():
    prompt = build_system_prompt(
        registry=ToolRegistry.default(),
        stage_instructions="Answer briefly.",
    )

    assert "Never pass `last_job_id` as a `path` argument" in prompt


def test_run_loop_install_command_policy_refuses_shell_install_text(tmp_path):
    provider = InstallCommandPolicyProvider()
    session = AgentSession(
        provider=provider,
        registry=ToolRegistry.default(),
        session_root=tmp_path,
    )

    result = session.run_loop("pip install chemsmart-agent 해줘")

    assert (
        result["assistant_output"]
        == "Please follow the project README for installation."
    )
    assert "pip install" not in result["assistant_output"]
