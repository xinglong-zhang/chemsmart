from __future__ import annotations

from pathlib import Path

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
