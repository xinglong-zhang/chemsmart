"""CHEMSMART.md rules: loading, precedence, policy parsing, injection."""

from __future__ import annotations

from pathlib import Path

from chemsmart.agent.behavior_rules import (
    MAX_RULES_CHARS_PER_FILE,
    load_behavior_rules,
    parse_policy_block,
    read_behavior_rules,
    write_behavior_rules,
)
from chemsmart.agent.command_answerer import _system_prompt_with_rules
from chemsmart.agent.prompts.identity import build_system_prompt
from chemsmart.agent.registry import ToolRegistry

_WORKSPACE_RULES = """ChemSmart agent — reproducible workflows.

## Style
Answer in Korean. Keep answers concise.

## Policy
- xtb_real_runs: auto
- verbosity: concise

## Notes
Prefer xTB pre-optimization before DFT.
"""

_USER_RULES = """## Style
Prefer detailed English answers.

## Policy
- xtb_real_runs: never
- answer_language: en
"""


def _isolate(monkeypatch, tmp_path: Path) -> tuple[Path, Path]:
    home = tmp_path / "home"
    workspace = tmp_path / "workspace"
    (home / ".chemsmart").mkdir(parents=True)
    workspace.mkdir()
    monkeypatch.setattr(Path, "home", staticmethod(lambda: home))
    monkeypatch.chdir(workspace)
    return home, workspace


def test_policy_block_parses_only_policy_section() -> None:
    policy = parse_policy_block(_WORKSPACE_RULES)
    assert policy == {"xtb_real_runs": "auto", "verbosity": "concise"}


def test_policy_block_accepts_plain_and_bullet_lines() -> None:
    text = "## Policy\nxtb_real_runs: ask\n* creativity: exploratory\n"
    assert parse_policy_block(text) == {
        "xtb_real_runs": "ask",
        "creativity": "exploratory",
    }


def test_load_merges_scopes_and_workspace_policy_wins(
    monkeypatch, tmp_path
) -> None:
    home, workspace = _isolate(monkeypatch, tmp_path)
    (home / ".chemsmart" / "CHEMSMART.md").write_text(_USER_RULES)
    (workspace / "CHEMSMART.md").write_text(_WORKSPACE_RULES)

    rules = load_behavior_rules()

    assert rules.loaded
    # Workspace overrides user for shared keys; user-only keys survive.
    assert rules.policy["xtb_real_runs"] == "auto"
    assert rules.policy["answer_language"] == "en"
    # User section rides first so the workspace section dominates attention.
    assert rules.text.index("[user rules") < rules.text.index(
        "[workspace rules"
    )


def test_load_without_any_rules_files_is_empty(monkeypatch, tmp_path) -> None:
    _isolate(monkeypatch, tmp_path)
    rules = load_behavior_rules()
    assert not rules.loaded
    assert rules.policy == {}


def test_oversized_rules_are_clamped_for_prompt_use(
    monkeypatch, tmp_path
) -> None:
    _, workspace = _isolate(monkeypatch, tmp_path)
    (workspace / "CHEMSMART.md").write_text(
        "R" * (MAX_RULES_CHARS_PER_FILE * 3)
    )
    rules = load_behavior_rules()
    assert rules.truncated
    assert len(rules.text) < MAX_RULES_CHARS_PER_FILE * 2


def test_write_refuses_overwrite_empty_and_unknown_scope(
    monkeypatch, tmp_path
) -> None:
    _isolate(monkeypatch, tmp_path)
    assert write_behavior_rules("", scope="workspace")["error"] == (
        "empty_content"
    )
    assert write_behavior_rules("x", scope="cluster")["error"] == (
        "unknown_scope"
    )
    first = write_behavior_rules(_WORKSPACE_RULES, scope="workspace")
    assert first["policy"]["xtb_real_runs"] == "auto"
    assert first["prompt_bounded"] is True
    again = write_behavior_rules("y", scope="workspace")
    assert again["error"] == "file_exists"
    replaced = write_behavior_rules(
        "# new\n", scope="workspace", overwrite=True
    )
    assert replaced["path"] == first["path"]


def test_write_user_scope_creates_parent_and_read_reports_both(
    monkeypatch, tmp_path
) -> None:
    home, workspace = _isolate(monkeypatch, tmp_path)
    (home / ".chemsmart").rmdir()
    write_behavior_rules(_USER_RULES, scope="user")
    (workspace / "CHEMSMART.md").write_text(_WORKSPACE_RULES)

    report = read_behavior_rules()

    assert report["user_path"] == str(home / ".chemsmart" / "CHEMSMART.md")
    assert report["workspace_path"] == str(workspace / "CHEMSMART.md")
    assert report["policy"]["xtb_real_runs"] == "auto"
    assert "[workspace rules" in report["content"]


def test_system_prompt_injects_bounded_rules_block(
    monkeypatch, tmp_path
) -> None:
    _, workspace = _isolate(monkeypatch, tmp_path)
    (workspace / "CHEMSMART.md").write_text(_WORKSPACE_RULES)
    rules = load_behavior_rules()

    prompt = build_system_prompt(
        registry=ToolRegistry.default(),
        stage_instructions="Stage instructions.",
        behavior_rules=rules.text,
        max_chars=4096,
    )

    assert "User rules (CHEMSMART.md)" in prompt
    assert "Answer in Korean" in prompt

    huge = build_system_prompt(
        registry=ToolRegistry.default(),
        stage_instructions="Stage instructions.",
        behavior_rules="R" * 50_000,
        max_chars=4096,
    )
    # Budget expands by exactly the clamped rules block, never unbounded.
    assert len(huge) <= 4096 + 1536 + 2


def test_answerer_prompt_gains_style_rules_but_keeps_contract(
    monkeypatch, tmp_path
) -> None:
    _, workspace = _isolate(monkeypatch, tmp_path)
    base = "BASE CONTRACT"
    assert _system_prompt_with_rules(base) == base

    (workspace / "CHEMSMART.md").write_text(_WORKSPACE_RULES)
    combined = _system_prompt_with_rules(base)
    assert combined.startswith(base)
    assert "never override grounding" in combined
    assert "Answer in Korean" in combined
