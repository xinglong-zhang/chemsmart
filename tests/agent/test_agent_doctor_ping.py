"""Ping-specific tests for `chemsmart agent doctor` (agent.yaml contract)."""

from click.testing import CliRunner

from chemsmart.agent.cli import agent
from chemsmart.agent.provider_config import AgentProviderConfig
from chemsmart.agent.providers import ProviderError
from chemsmart.agent.registry import ToolRegistry

_REGISTERED_TOOLS_LINE = (
    f"tools registered: {len(ToolRegistry.default().list_tools())}"
)


def _config(provider_type="openai"):
    return AgentProviderConfig(
        name="main",
        type=provider_type,
        api_key="testkey12345678901234567890",
        model="chemsmart-test",
        base_url="https://gateway.test/v1",
        extra_headers={},
    )


def _install(monkeypatch, *, provider):
    monkeypatch.setattr(
        "chemsmart.agent.provider_config.load_active_provider_config",
        lambda *_a, **_k: _config(),
    )
    monkeypatch.setattr(
        "chemsmart.agent.providers.get_provider",
        lambda *_a, **_k: provider,
    )


class _Provider:
    def __init__(self, *, ping_exc=None):
        self.ping_calls = 0
        self._ping_exc = ping_exc

    def ping(self):
        self.ping_calls += 1
        if self._ping_exc is not None:
            raise self._ping_exc
        return {
            "ok": True,
            "resolved_model": "gpt-5.4-2026-03-05",
            "latency_ms": 42,
        }


def test_doctor_no_ping_skips_provider_ping(monkeypatch):
    """--no-ping should not call provider.ping()."""
    provider = _Provider()
    _install(monkeypatch, provider=provider)

    result = CliRunner().invoke(agent, ["doctor", "--no-ping"])

    assert result.exit_code == 0, result.output
    assert provider.ping_calls == 0
    assert "ping: skipped (--no-ping)" in result.output


def test_doctor_ping_failure_exits_nonzero(monkeypatch):
    """A ping failure should make doctor exit non-zero and report it."""
    provider = _Provider(ping_exc=ProviderError("boom"))
    _install(monkeypatch, provider=provider)

    result = CliRunner().invoke(agent, ["doctor"])

    assert result.exit_code != 0
    assert "ping: FAILED boom" in result.output
    assert _REGISTERED_TOOLS_LINE in result.output
