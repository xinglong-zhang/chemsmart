"""`chemsmart config agent` writes a selection the live loader accepts.

The credential goes to the managed store (0600) and never into agent.yaml;
the written yaml is round-trip validated and rolled back on refusal; the
legacy api.env migrates verbatim, unknown labels included.
"""

from __future__ import annotations

import yaml
from click.testing import CliRunner

from chemsmart.cli.config import config


def _home(monkeypatch, tmp_path):
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.delenv("CHEMSMART_AGENT_CONFIG", raising=False)
    monkeypatch.delenv("CHEMSMART_AGENT_KEYS", raising=False)
    return home


def test_flags_write_a_loader_valid_profile_and_store_the_key(
    tmp_path, monkeypatch
):
    home = _home(monkeypatch, tmp_path)

    result = CliRunner().invoke(
        config,
        [
            "agent",
            "--provider",
            "openai",
            "--model",
            "gpt-5.2",
            "--reasoning-effort",
            "high",
            "--context-tokens",
            "400000",
            "--max-output-tokens",
            "64000",
            "--api-key-stdin",
        ],
        input="sk-test-key\n",
    )

    assert result.exit_code == 0, result.output
    agent_yaml = home / ".chemsmart" / "agent" / "agent.yaml"
    payload = yaml.safe_load(agent_yaml.read_text(encoding="utf-8"))
    assert payload["active"] == "openai"
    assert payload["fallback"] == []
    profile = payload["providers"]["openai"]
    assert profile["model"] == "gpt-5.2"
    assert profile["base_url"] == "https://api.openai.com/v1"
    assert "api_key" not in profile, "keys never live in agent.yaml"
    assert "sk-test-key" not in agent_yaml.read_text(encoding="utf-8")
    keys = home / ".chemsmart" / "agent" / "keys.env"
    assert keys.read_text(encoding="utf-8") == "OPENAI_API_KEY=sk-test-key\n"
    assert (keys.stat().st_mode & 0o777) == 0o600
    assert "sk-test-key" not in result.output, "never echoed"

    from chemsmart.agent.provider_config import (
        load_agent_provider_selection,
    )

    profile = load_agent_provider_selection(agent_yaml).active_profile
    assert profile.provider == "openai"


def test_prompts_walk_a_new_user_through_the_choices(tmp_path, monkeypatch):
    home = _home(monkeypatch, tmp_path)

    result = CliRunner().invoke(
        config,
        ["agent", "--skip-key"],
        input="deepseek\ndeepseek-v4-flash\nmax\n1000000\n384000\n",
    )

    assert result.exit_code == 0, result.output
    payload = yaml.safe_load(
        (home / ".chemsmart" / "agent" / "agent.yaml").read_text(
            encoding="utf-8"
        )
    )
    assert payload["providers"]["deepseek"]["model"] == "deepseek-v4-flash"
    assert "export DEEPSEEK" in result.output


def test_an_invalid_profile_is_rolled_back(tmp_path, monkeypatch):
    home = _home(monkeypatch, tmp_path)
    agent_dir = home / ".chemsmart" / "agent"
    agent_dir.mkdir(parents=True)
    original = "active: alibaba-token-plan\nfallback: []\nproviders: {}\n"
    (agent_dir / "agent.yaml").write_text(original, encoding="utf-8")

    result = CliRunner().invoke(
        config,
        [
            "agent",
            "--provider",
            "openai",
            "--model",
            "gpt-5.2",
            "--context-tokens",
            "1000",
            "--max-output-tokens",
            "64000",  # exceeds context -> loader refuses
            "--skip-key",
        ],
    )

    assert result.exit_code != 0
    assert "rolled back" in result.output
    assert (agent_dir / "agent.yaml").read_text(encoding="utf-8") == original


def test_the_api_env_migration_carries_unknown_labels_verbatim(
    tmp_path, monkeypatch
):
    home = _home(monkeypatch, tmp_path)
    agent_dir = home / ".chemsmart" / "agent"
    agent_dir.mkdir(parents=True)
    (agent_dir / "api.env").write_text(
        "ALIBABA-TOKEN_PLAN_KEY=sk-sp-legacy\n"
        "SEMANTIC-SCHOLAR-API-KEY=ss-legacy\n",
        encoding="utf-8",
    )

    result = CliRunner().invoke(
        config,
        [
            "agent",
            "--provider",
            "deepseek",
            "--model",
            "deepseek-v4-flash",
            "--reasoning-effort",
            "max",
            "--context-tokens",
            "1000000",
            "--max-output-tokens",
            "384000",
            "--skip-key",
        ],
        input="y\n",
    )

    assert result.exit_code == 0, result.output
    keys = (agent_dir / "keys.env").read_text(encoding="utf-8")
    assert "ALIBABA-TOKEN_PLAN_KEY=sk-sp-legacy" in keys
    assert "SEMANTIC-SCHOLAR-API-KEY=ss-legacy" in keys
    assert (agent_dir / "api.env").exists(), "the legacy file is untouched"
    assert "imported 2 credential label(s)" in result.output
