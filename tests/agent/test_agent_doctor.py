"""Tests for `chemsmart agent doctor` (agent.yaml provider-config contract).

`doctor` reads the active provider from ``~/.chemsmart/agent/agent.yaml`` via
``load_active_provider_config`` and pings it. These tests drive that contract by
monkeypatching the config loader + ``providers.get_provider`` (both imported
inside the command), so they never touch the real config or the network.
"""

import logging
import sys
from unittest.mock import MagicMock

from click.testing import CliRunner

from chemsmart.agent.cli import agent
from chemsmart.agent.provider_config import (
    AgentProviderConfig,
    AgentProviderConfigError,
)
from chemsmart.agent.providers import ProviderError
from chemsmart.agent.registry import ToolRegistry

_PING_RESULT = {
    "ok": True,
    "resolved_model": "gpt-5.4-2026-03-05",
    "latency_ms": 42,
}
_REGISTERED_TOOLS_LINE = (
    f"tools registered: {len(ToolRegistry.default().list_tools())}"
)


def _config(provider_type, *, name="main", model="chemsmart-test"):
    return AgentProviderConfig(
        name=name,
        type=provider_type,
        api_key="testkey12345678901234567890",
        model=model,
        base_url="https://gateway.test/v1",
        extra_headers={},
    )


class _FakeProvider:
    def __init__(self, *, ping_result=None, ping_exc=None):
        self._ping_result = ping_result
        self._ping_exc = ping_exc

    def ping(self):
        if self._ping_exc is not None:
            raise self._ping_exc
        return self._ping_result


def _install(
    monkeypatch,
    *,
    config="__unset__",
    config_exc=None,
    provider=None,
    provider_exc=None,
):
    """Drive doctor: control load_active_provider_config + get_provider."""
    def _load(*_args, **_kwargs):
        if config_exc is not None:
            raise config_exc
        return None if config == "__unset__" else config

    def _get(*_args, **_kwargs):
        if provider_exc is not None:
            raise provider_exc
        return provider

    monkeypatch.setattr(
        "chemsmart.agent.provider_config.load_active_provider_config", _load
    )
    monkeypatch.setattr("chemsmart.agent.providers.get_provider", _get)


def test_doctor_valid_anthropic(monkeypatch):
    """active anthropic provider in agent.yaml + ping ok -> green output."""
    _install(
        monkeypatch,
        config=_config("anthropic"),
        provider=_FakeProvider(ping_result=_PING_RESULT),
    )
    monkeypatch.setitem(sys.modules, "anthropic", MagicMock())

    result = CliRunner().invoke(agent, ["doctor"])

    assert result.exit_code == 0, result.output
    assert "agent.yaml:" in result.output
    assert "active: main" in result.output
    assert "type: anthropic" in result.output
    assert "model: chemsmart-test" in result.output
    assert "ping: ok (model=gpt-5.4-2026-03-05, latency=42ms)" in result.output
    assert _REGISTERED_TOOLS_LINE in result.output
    assert "WARN: this gateway tenant may not have Anthropic access" in (
        result.output
    )


def test_doctor_valid_openai(monkeypatch):
    """active openai provider + ping ok -> green output, no anthropic warn."""
    _install(
        monkeypatch,
        config=_config("openai"),
        provider=_FakeProvider(ping_result=_PING_RESULT),
    )
    monkeypatch.setitem(sys.modules, "openai", MagicMock())

    result = CliRunner().invoke(agent, ["doctor"])

    assert result.exit_code == 0, result.output
    assert "type: openai" in result.output
    assert "ping: ok (model=gpt-5.4-2026-03-05, latency=42ms)" in result.output
    assert _REGISTERED_TOOLS_LINE in result.output
    assert "Anthropic access" not in result.output


def test_doctor_valid_local(monkeypatch):
    """active local v13.1 provider + ping ok -> green output."""
    _install(
        monkeypatch,
        config=_config(
            "local",
            name="local_chemsmart_v13_1",
            model="chemsmart-qwen2.5-coder-3b-instruct-v13_1",
        ),
        provider=_FakeProvider(ping_result=_PING_RESULT),
    )

    result = CliRunner().invoke(agent, ["doctor"])

    assert result.exit_code == 0, result.output
    assert "active: local_chemsmart_v13_1" in result.output
    assert "type: local" in result.output
    assert "model: chemsmart-qwen2.5-coder-3b-instruct-v13_1" in result.output
    assert "ping: ok" in result.output


def test_doctor_default_output_has_no_debug_or_warnings(monkeypatch, tmp_path):
    """Default doctor output stays quiet; logs are captured to the tui log."""
    monkeypatch.setattr(
        "chemsmart.agent.cli_commands._default_session_root",
        lambda: str(tmp_path / "sessions"),
    )

    def noisy_ping():
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

    provider = _FakeProvider()
    monkeypatch.setattr(provider, "ping", noisy_ping)
    _install(monkeypatch, config=_config("openai"), provider=provider)

    result = CliRunner().invoke(agent, ["doctor"], catch_exceptions=False)

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


def test_doctor_verbose_emits_debug(monkeypatch, caplog):
    """--verbose keeps DEBUG output on the console."""
    def verbose_ping():
        logging.getLogger("openai.http").debug("verbose doctor debug")
        return _PING_RESULT

    provider = _FakeProvider()
    monkeypatch.setattr(provider, "ping", verbose_ping)
    _install(monkeypatch, config=_config("openai"), provider=provider)

    caplog.set_level(logging.DEBUG)
    result = CliRunner().invoke(agent, ["--verbose", "doctor"])

    assert result.exit_code == 0, result.output
    assert "verbose doctor debug" in caplog.text


def test_doctor_no_ping_reports_skip(monkeypatch):
    """doctor --no-ping prints skip and still succeeds (ping not called)."""
    def _ping_should_not_run():  # pragma: no cover - defensive
        raise AssertionError("ping should not run when --no-ping is used")

    provider = _FakeProvider()
    monkeypatch.setattr(provider, "ping", _ping_should_not_run)
    _install(monkeypatch, config=_config("openai"), provider=provider)

    result = CliRunner().invoke(agent, ["doctor", "--no-ping"])

    assert result.exit_code == 0, result.output
    assert "ping: skipped (--no-ping)" in result.output
    assert _REGISTERED_TOOLS_LINE in result.output


def test_doctor_malformed_config_exits_nonzero(monkeypatch):
    """A malformed agent.yaml -> AgentProviderConfigError -> exit != 0."""
    _install(
        monkeypatch,
        config_exc=AgentProviderConfigError(
            "agent.yaml: 'active' does not name a configured provider"
        ),
    )

    result = CliRunner().invoke(agent, ["doctor"])

    assert result.exit_code != 0
    assert "active" in result.output


def test_doctor_unsupported_provider_lists_supported(monkeypatch):
    """An unsupported provider type -> exit != 0, names the supported set."""
    _install(
        monkeypatch,
        config_exc=AgentProviderConfigError(
            "unsupported provider type 'azure'; supported: "
            "openai, anthropic, local"
        ),
    )

    result = CliRunner().invoke(agent, ["doctor"])

    assert result.exit_code != 0
    assert "azure" in result.output
    assert "openai" in result.output
    assert "anthropic" in result.output
    assert "local" in result.output


def test_doctor_provider_build_error_exits_nonzero(monkeypatch):
    """A configured provider that fails to construct -> exit != 0."""
    _install(
        monkeypatch,
        config=_config("openai"),
        provider_exc=ProviderError("missing api key for provider 'main'"),
    )

    result = CliRunner().invoke(agent, ["doctor"])

    assert result.exit_code != 0
    assert "api key" in result.output


def test_doctor_no_config_and_no_env_exits_nonzero(monkeypatch):
    """No agent.yaml and no usable env provider -> exit != 0."""
    monkeypatch.delenv("AI_PROVIDER", raising=False)
    _install(
        monkeypatch,
        config="__unset__",  # load_active_provider_config returns None
        provider_exc=ProviderError("no provider configured"),
    )

    result = CliRunner().invoke(agent, ["doctor"])

    assert result.exit_code != 0
    assert "provider" in result.output.lower()
