"""
Tests for `chemsmart agent doctor`.

Test cases:
  a. valid api.env + AI_PROVIDER=anthropic -> green output.
  b. Missing AI_PROVIDER -> exit != 0, message names env var.
  c. Empty/missing ai_api_key -> exit != 0, message names key.
  d. Unsupported AI_PROVIDER -> exit != 0, lists supported providers.
"""

import sys
from unittest.mock import MagicMock

import pytest
from click.testing import CliRunner

from chemsmart.agent.cli import agent


@pytest.fixture
def api_env_file(tmp_path):
    env_file = tmp_path / "api.env"
    env_file.write_text("ai_api_key=testkey12345678901234567890\n")
    return str(env_file)


def test_doctor_valid_anthropic(monkeypatch, api_env_file):
    """a. valid api.env + AI_PROVIDER=anthropic -> green output."""
    monkeypatch.setenv("AI_PROVIDER", "anthropic")
    monkeypatch.setattr(
        "chemsmart.agent.providers._API_ENV_PATH", api_env_file
    )

    mock_anthropic = MagicMock()
    monkeypatch.setitem(sys.modules, "anthropic", mock_anthropic)

    runner = CliRunner()
    result = runner.invoke(agent, ["doctor"])

    assert result.exit_code == 0, result.output
    assert "AI_PROVIDER=anthropic OK" in result.output
    assert "api.env: OK (key length=" in result.output
    assert "tools registered: 0" in result.output


def test_doctor_missing_ai_provider(monkeypatch, api_env_file):
    """b. Missing AI_PROVIDER -> exit != 0, message names env var."""
    monkeypatch.delenv("AI_PROVIDER", raising=False)
    monkeypatch.setattr(
        "chemsmart.agent.providers._API_ENV_PATH", api_env_file
    )

    runner = CliRunner()
    result = runner.invoke(agent, ["doctor"])

    assert result.exit_code != 0
    assert "AI_PROVIDER" in result.output


def test_doctor_empty_api_key(monkeypatch, tmp_path):
    """c. Empty/missing ai_api_key -> exit != 0, message names key."""
    env_file = tmp_path / "api.env"
    env_file.write_text("ai_api_key=\n")

    monkeypatch.setenv("AI_PROVIDER", "anthropic")
    monkeypatch.setattr(
        "chemsmart.agent.providers._API_ENV_PATH", str(env_file)
    )

    runner = CliRunner()
    result = runner.invoke(agent, ["doctor"])

    assert result.exit_code != 0
    assert "ai_api_key" in result.output


def test_doctor_unsupported_provider(monkeypatch, api_env_file):
    """d. Unsupported AI_PROVIDER -> exit != 0, lists supported providers."""
    monkeypatch.setenv("AI_PROVIDER", "azure")
    monkeypatch.setattr(
        "chemsmart.agent.providers._API_ENV_PATH", api_env_file
    )

    runner = CliRunner()
    result = runner.invoke(agent, ["doctor"])

    assert result.exit_code != 0
    assert "azure" in result.output
    assert "anthropic" in result.output
