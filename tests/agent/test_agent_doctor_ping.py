"""Ping-specific tests for `chemsmart agent doctor`."""

import sys
from unittest.mock import MagicMock

from click.testing import CliRunner

from chemsmart.agent.cli import agent
from chemsmart.agent.providers import ProviderError


def test_doctor_no_ping_skips_provider_ping(monkeypatch, tmp_path):
    """--no-ping should not call provider.ping()."""
    env_file = tmp_path / "api.env"
    env_file.write_text("ai_api_key=testkey12345678901234567890\n")

    monkeypatch.setenv("AI_PROVIDER", "openai")
    monkeypatch.setattr(
        "chemsmart.agent.providers._API_ENV_PATH", str(env_file)
    )

    ping_calls = {"count": 0}

    def _count_ping(self):
        ping_calls["count"] += 1
        return {
            "ok": True,
            "resolved_model": "gpt-5.4-2026-03-05",
            "latency_ms": 42,
        }

    monkeypatch.setattr(
        "chemsmart.agent.providers.OpenAIProvider.ping",
        _count_ping,
    )

    mock_openai = MagicMock()
    monkeypatch.setitem(sys.modules, "openai", mock_openai)

    runner = CliRunner()
    result = runner.invoke(agent, ["doctor", "--no-ping"])

    assert result.exit_code == 0, result.output
    assert ping_calls["count"] == 0
    assert "ping: skipped (--no-ping)" in result.output


def test_doctor_ping_failure_exits_nonzero(monkeypatch, tmp_path):
    """A ping failure should make doctor exit non-zero."""
    env_file = tmp_path / "api.env"
    env_file.write_text("ai_api_key=testkey12345678901234567890\n")

    monkeypatch.setenv("AI_PROVIDER", "openai")
    monkeypatch.setattr(
        "chemsmart.agent.providers._API_ENV_PATH", str(env_file)
    )
    monkeypatch.setattr(
        "chemsmart.agent.providers.OpenAIProvider.ping",
        lambda self: (_ for _ in ()).throw(ProviderError("boom")),
    )

    mock_openai = MagicMock()
    monkeypatch.setitem(sys.modules, "openai", mock_openai)

    runner = CliRunner()
    result = runner.invoke(agent, ["doctor"])

    assert result.exit_code != 0
    assert "gateway: https://factchat-cloud.mindlogic.ai/v1/gateway" in (
        result.output
    )
    assert "ping: FAILED boom" in result.output
    assert "tools registered: 9" in result.output
