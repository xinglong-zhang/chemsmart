"""Tests for provider selection."""

from __future__ import annotations

import pytest

from chemsmart.agent import providers
from chemsmart.agent.provider_config import AgentProviderConfig
from chemsmart.agent.providers import ProviderError


class DummyOpenAIProvider:
    def __init__(self, api_key: str, **kwargs) -> None:
        self.api_key = api_key
        self.kwargs = kwargs


class DummyAnthropicProvider:
    def __init__(self, api_key: str, **kwargs) -> None:
        self.api_key = api_key
        self.kwargs = kwargs


def test_get_provider_prefers_yaml_config(monkeypatch):
    config = AgentProviderConfig(
        name="main",
        type="openai",
        api_key="yaml-key",
        model="gpt-yaml",
        base_url="https://example.test/v1",
        extra_headers={"X-Test": "value"},
    )
    monkeypatch.setattr(
        providers, "load_active_provider_config", lambda: config
    )
    monkeypatch.setattr(providers, "OpenAIProvider", DummyOpenAIProvider)
    monkeypatch.setenv("AI_PROVIDER", "anthropic")
    monkeypatch.setenv("ai_api_key", "legacy-key")

    provider = providers.get_provider()

    assert isinstance(provider, DummyOpenAIProvider)
    assert provider.api_key == "yaml-key"
    assert provider.kwargs == {
        "model": "gpt-yaml",
        "base_url": "https://example.test/v1",
        "extra_headers": {"X-Test": "value"},
    }


def test_get_provider_legacy_fallback_emits_deprecation_warning(
    monkeypatch, tmp_path
):
    env_file = tmp_path / "api.env"
    env_file.write_text("AI_PROVIDER=anthropic\nai_api_key=legacy-key\n")
    monkeypatch.setattr(providers, "load_active_provider_config", lambda: None)
    monkeypatch.setattr(providers, "AnthropicProvider", DummyAnthropicProvider)
    monkeypatch.delenv("AI_PROVIDER", raising=False)
    monkeypatch.delenv("ai_api_key", raising=False)

    with pytest.warns(DeprecationWarning, match="agent.yaml missing"):
        provider = providers.get_provider(env_path=str(env_file))

    assert isinstance(provider, DummyAnthropicProvider)
    assert provider.api_key == "legacy-key"


def test_get_provider_both_yaml_and_legacy_missing_raises(
    monkeypatch, tmp_path
):
    monkeypatch.setattr(providers, "load_active_provider_config", lambda: None)
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    monkeypatch.chdir(tmp_path)
    monkeypatch.delenv("AI_PROVIDER", raising=False)
    monkeypatch.delenv("ai_api_key", raising=False)
    monkeypatch.delenv("CHEMSMART_API_ENV", raising=False)

    with pytest.warns(DeprecationWarning, match="agent.yaml missing"):
        with pytest.raises(ProviderError, match="AI_PROVIDER"):
            providers.get_provider()


@pytest.mark.parametrize(
    ("model", "expected"),
    [
        ("gpt-5.4", True),
        ("gpt-5o", True),
        ("o1-preview", True),
        ("o3-mini", True),
        ("o4-mini", True),
        ("gpt-4o", False),
        ("gpt-4.1-turbo", False),
        ("gpt-4-turbo", False),
    ],
)
def test_openai_uses_max_completion_tokens_predicate(
    model: str, expected: bool
) -> None:
    assert providers._openai_uses_max_completion_tokens(model) is expected
