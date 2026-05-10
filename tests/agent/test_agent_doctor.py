"""Tests for `chemsmart agent doctor`."""

import logging
import sys
from unittest.mock import MagicMock

import pytest
from click.testing import CliRunner

from chemsmart.agent.cli import agent
from chemsmart.agent.registry import ToolRegistry

_PING_RESULT = {
    "ok": True,
    "resolved_model": "gpt-5.4-2026-03-05",
    "latency_ms": 42,
}
_REGISTERED_TOOLS_LINE = (
    f"tools registered: {len(ToolRegistry.default().list_tools())}"
)


@pytest.fixture
def api_env_file(tmp_path):
    env_file = tmp_path / "api.env"
    env_file.write_text("ai_api_key=testkey12345678901234567890\n")
    return str(env_file)


def test_doctor_valid_anthropic(monkeypatch, api_env_file):
    """valid api.env + AI_PROVIDER=anthropic -> green output."""
    monkeypatch.setenv("AI_PROVIDER", "anthropic")
    monkeypatch.setenv("CHEMSMART_API_ENV", api_env_file)
    monkeypatch.setattr(
        "chemsmart.agent.providers.AnthropicProvider.ping",
        lambda self: _PING_RESULT,
    )

    mock_anthropic = MagicMock()
    monkeypatch.setitem(sys.modules, "anthropic", mock_anthropic)

    runner = CliRunner()
    result = runner.invoke(agent, ["doctor"])

    assert result.exit_code == 0, result.output
    assert "AI_PROVIDER=anthropic OK" in result.output
    assert "api.env: OK (key length=" in result.output
    assert (
        "gateway: https://factchat-cloud.mindlogic.ai/v1/gateway/claude"
        in result.output
    )
    assert "ping: ok (model=gpt-5.4-2026-03-05, latency=42ms)" in (
        result.output
    )
    assert _REGISTERED_TOOLS_LINE in result.output
    assert "WARN: this gateway tenant may not have Anthropic access" in (
        result.output
    )


def test_doctor_valid_openai(monkeypatch, api_env_file):
    """AI_PROVIDER=openai + valid api.env -> green output."""
    monkeypatch.setenv("AI_PROVIDER", "openai")
    monkeypatch.setenv("CHEMSMART_API_ENV", api_env_file)
    monkeypatch.setattr(
        "chemsmart.agent.providers.OpenAIProvider.ping",
        lambda self: _PING_RESULT,
    )

    mock_openai = MagicMock()
    monkeypatch.setitem(sys.modules, "openai", mock_openai)

    runner = CliRunner()
    result = runner.invoke(agent, ["doctor"])

    assert result.exit_code == 0, result.output
    assert "AI_PROVIDER=openai OK" in result.output
    assert "api.env: OK (key length=" in result.output
    assert "gateway: https://factchat-cloud.mindlogic.ai/v1/gateway" in (
        result.output
    )
    assert "ping: ok (model=gpt-5.4-2026-03-05, latency=42ms)" in (
        result.output
    )
    assert _REGISTERED_TOOLS_LINE in result.output


def test_doctor_default_output_has_no_debug_or_warnings(
    monkeypatch, api_env_file, tmp_path
):
    """Default doctor output should be quiet while still capturing logs."""
    monkeypatch.setenv("AI_PROVIDER", "openai")
    monkeypatch.setenv("CHEMSMART_API_ENV", api_env_file)
    monkeypatch.setattr(
        "chemsmart.agent.cli._default_session_root",
        lambda: str(tmp_path / "sessions"),
    )

    def noisy_ping(self):
        logging.getLogger("chemsmart.agent.registry").warning(
            "registry warning should stay out of stderr"
        )
        logging.getLogger("numexpr.utils").info(
            "numexpr info should stay out of stderr"
        )
        logging.getLogger("openai.http").debug(
            "openai debug should stay out of stderr"
        )
        return _PING_RESULT

    monkeypatch.setattr(
        "chemsmart.agent.providers.OpenAIProvider.ping",
        noisy_ping,
    )

    mock_openai = MagicMock()
    monkeypatch.setitem(sys.modules, "openai", mock_openai)

    runner = CliRunner()
    result = runner.invoke(agent, ["doctor"], catch_exceptions=False)

    assert result.exit_code == 0, result.output
    assert "registry warning should stay out of stderr" not in result.output
    assert "numexpr info should stay out of stderr" not in result.output
    assert "openai debug should stay out of stderr" not in result.output

    log_path = tmp_path / "_tui.log"
    assert log_path.exists()
    log_text = log_path.read_text(encoding="utf-8")
    assert "registry warning should stay out of stderr" in log_text
    assert "numexpr info should stay out of stderr" in log_text
    assert "openai debug should stay out of stderr" in log_text


def test_doctor_verbose_emits_debug(monkeypatch, api_env_file, caplog):
    """--verbose should keep DEBUG output on the console."""
    monkeypatch.setenv("AI_PROVIDER", "openai")
    monkeypatch.setenv("CHEMSMART_API_ENV", api_env_file)

    def verbose_ping(self):
        logging.getLogger("openai.http").debug("verbose doctor debug")
        return _PING_RESULT

    monkeypatch.setattr(
        "chemsmart.agent.providers.OpenAIProvider.ping",
        verbose_ping,
    )

    mock_openai = MagicMock()
    monkeypatch.setitem(sys.modules, "openai", mock_openai)

    caplog.set_level(logging.DEBUG)
    runner = CliRunner()
    result = runner.invoke(agent, ["--verbose", "doctor"])

    assert result.exit_code == 0, result.output
    assert "verbose doctor debug" in caplog.text


def test_doctor_no_ping_reports_skip(monkeypatch, api_env_file):
    """doctor --no-ping prints skip and still succeeds."""
    monkeypatch.setenv("AI_PROVIDER", "openai")
    monkeypatch.setenv("CHEMSMART_API_ENV", api_env_file)

    def _ping_should_not_run(self):  # pragma: no cover - defensive
        raise AssertionError("ping should not run when --no-ping is used")

    monkeypatch.setattr(
        "chemsmart.agent.providers.OpenAIProvider.ping",
        _ping_should_not_run,
    )

    mock_openai = MagicMock()
    monkeypatch.setitem(sys.modules, "openai", mock_openai)

    runner = CliRunner()
    result = runner.invoke(agent, ["doctor", "--no-ping"])

    assert result.exit_code == 0, result.output
    assert "ping: skipped (--no-ping)" in result.output
    assert _REGISTERED_TOOLS_LINE in result.output


def test_doctor_missing_ai_provider(monkeypatch, api_env_file):
    """Missing AI_PROVIDER -> exit != 0, message names env var."""
    monkeypatch.delenv("AI_PROVIDER", raising=False)
    monkeypatch.setenv("CHEMSMART_API_ENV", api_env_file)

    runner = CliRunner()
    result = runner.invoke(agent, ["doctor"])

    assert result.exit_code != 0
    assert "AI_PROVIDER" in result.output


def test_doctor_empty_api_key(monkeypatch, tmp_path):
    """Empty/missing ai_api_key -> exit != 0, message names key."""
    env_file = tmp_path / "api.env"
    env_file.write_text("ai_api_key=\n")

    monkeypatch.setenv("AI_PROVIDER", "anthropic")
    monkeypatch.setenv("CHEMSMART_API_ENV", str(env_file))

    runner = CliRunner()
    result = runner.invoke(agent, ["doctor"])

    assert result.exit_code != 0
    assert "ai_api_key" in result.output


def test_doctor_unsupported_provider(monkeypatch, api_env_file):
    """Unsupported AI_PROVIDER -> exit != 0, lists supported providers."""
    monkeypatch.setenv("AI_PROVIDER", "azure")
    monkeypatch.setenv("CHEMSMART_API_ENV", api_env_file)

    runner = CliRunner()
    result = runner.invoke(agent, ["doctor"])

    assert result.exit_code != 0
    assert "azure" in result.output
    assert "anthropic" in result.output
    assert "openai" in result.output


def test_doctor_unsupported_now_lists_both(monkeypatch, api_env_file):
    """Unsupported AI_PROVIDER now lists Anthropic and OpenAI."""
    monkeypatch.setenv("AI_PROVIDER", "bedrock")
    monkeypatch.setenv("CHEMSMART_API_ENV", api_env_file)

    runner = CliRunner()
    result = runner.invoke(agent, ["doctor"])

    assert result.exit_code != 0
    assert "anthropic" in result.output
    assert "openai" in result.output
