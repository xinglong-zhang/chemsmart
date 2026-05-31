"""Tests for agent provider YAML configuration."""

from __future__ import annotations

import pytest

from chemsmart.agent.provider_config import (
    AgentProviderConfig,
    AgentProviderConfigError,
    load_active_provider_config,
)


def _write_agent_yaml(path, text: str):
    path.write_text(text, encoding="utf-8")
    return path


def test_load_active_provider_config_valid_yaml(tmp_path):
    yaml_path = _write_agent_yaml(
        tmp_path / "agent.yaml",
        """
active: main
providers:
  main:
    type: openai
    api_key_env: OPENAI_API_KEY
    api_key: literal-key
    model: gpt-test
    base_url: https://example.test/v1
    extra_headers:
      X-Test: value
""",
    )

    config = load_active_provider_config(yaml_path)

    assert config == AgentProviderConfig(
        name="main",
        type="openai",
        api_key="literal-key",
        model="gpt-test",
        base_url="https://example.test/v1",
        extra_headers={"X-Test": "value"},
    )


def test_load_active_provider_config_missing_yaml_returns_none(tmp_path):
    assert load_active_provider_config(tmp_path / "missing.yaml") is None


@pytest.mark.parametrize(
    "text",
    [
        "active: [broken\n",
        "providers: {}\n",
    ],
)
def test_load_active_provider_config_rejects_malformed_or_missing_active(
    tmp_path, text
):
    yaml_path = _write_agent_yaml(tmp_path / "agent.yaml", text)

    with pytest.raises(AgentProviderConfigError):
        load_active_provider_config(yaml_path)


def test_load_active_provider_config_rejects_unknown_type(tmp_path):
    yaml_path = _write_agent_yaml(
        tmp_path / "agent.yaml",
        """
active: localish
providers:
  localish:
    type: local
    api_key: key
    model: model
    base_url: ""
    extra_headers: {}
""",
    )

    with pytest.raises(AgentProviderConfigError):
        load_active_provider_config(yaml_path)


def test_load_active_provider_config_resolves_api_key_env(
    monkeypatch, tmp_path
):
    monkeypatch.setenv("ANTHROPIC_API_KEY", "env-key")
    yaml_path = _write_agent_yaml(
        tmp_path / "agent.yaml",
        """
active: claude
providers:
  claude:
    type: anthropic
    api_key_env: ANTHROPIC_API_KEY
    api_key: ""
    model: claude-test
    base_url: ""
    extra_headers: {}
""",
    )

    config = load_active_provider_config(yaml_path)

    assert config is not None
    assert config.api_key == "env-key"
    assert config.type == "anthropic"
