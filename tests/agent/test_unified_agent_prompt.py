from __future__ import annotations

from pathlib import Path

from chemsmart.agent.prompts.identity import build_system_prompt
from chemsmart.agent.registry import ToolRegistry
from chemsmart.agent.tools_command import execute_chemsmart_command


def test_explicit_submit_policy_requires_execute_not_test_mode() -> None:
    prompt = (
        Path(__file__).parents[2]
        / "chemsmart"
        / "agent"
        / "prompts"
        / "unified_agent.md"
    ).read_text(encoding="utf-8")

    assert "configured mock\n  scheduler" in prompt
    assert "execute_chemsmart_command(test=false)" in prompt
    assert "Test mode is not a successful submission" in prompt
    docstring = execute_chemsmart_command.__doc__ or ""
    assert "test=False" in docstring
    assert "explicit run or submission" in docstring


def test_unified_outer_prompt_fits_4096_character_budget() -> None:
    stage = (
        Path(__file__).parents[2]
        / "chemsmart"
        / "agent"
        / "prompts"
        / "unified_agent.md"
    ).read_text(encoding="utf-8")
    prompt = build_system_prompt(
        registry=ToolRegistry.default(),
        stage_instructions=stage,
        session_meta={
            "stage": "tool_loop",
            "approval_mode": "driving",
            "yolo": False,
            "session_id": "20260715T170000Z-12345678",
            "turn_index": 12,
            "request_intent": "workflow",
        },
        conversation_context={
            "entities": {"last_project": "demo", "last_server": "cluster"},
            "recent_turns": [{"request": "x" * 2000}],
        },
        request="create a project YAML then submit an ORCA TS on cluster",
        max_chars=4096,
    )
    assert len(prompt) <= 4096
    assert "Compact workflow state" in prompt
