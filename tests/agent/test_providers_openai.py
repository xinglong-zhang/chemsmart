"""Tests for the OpenAI provider adapter."""

import sys
from unittest.mock import MagicMock

from chemsmart.agent.providers import OpenAIProvider


def test_openai_provider_chat_returns_dict_and_forwards_tools(monkeypatch):
    """OpenAIProvider.chat returns a dict and forwards tools."""
    messages = [{"role": "user", "content": "Hello"}]
    tools = [{"type": "function", "function": {"name": "demo"}}]

    response = MagicMock()
    response.model_dump.return_value = {"id": "chatcmpl_test"}

    completions = MagicMock()
    completions.create.return_value = response

    client = MagicMock()
    client.chat.completions = completions

    openai_module = MagicMock()
    openai_module.OpenAI.return_value = client
    monkeypatch.setitem(sys.modules, "openai", openai_module)

    provider = OpenAIProvider("test-key")
    result = provider.chat(messages, tools=tools)

    assert result == {"id": "chatcmpl_test"}
    openai_module.OpenAI.assert_called_once_with(
        api_key="test-key",
        base_url="https://factchat-cloud.mindlogic.ai/v1/gateway",
    )
    completions.create.assert_called_once_with(
        model="gpt-5.4",
        messages=messages,
        tools=tools,
        timeout=30,
    )


def test_openai_provider_ping_returns_resolved_model(monkeypatch):
    """OpenAIProvider.ping returns model metadata from the response."""
    response = MagicMock()
    response.model = "gpt-5.4-2026-03-05"

    completions = MagicMock()
    completions.create.return_value = response

    client = MagicMock()
    client.chat.completions = completions

    openai_module = MagicMock()
    openai_module.OpenAI.return_value = client
    monkeypatch.setitem(sys.modules, "openai", openai_module)

    provider = OpenAIProvider("test-key")
    result = provider.ping()

    assert result["ok"] is True
    assert result["resolved_model"] == "gpt-5.4-2026-03-05"
    assert result["latency_ms"] >= 0
    completions.create.assert_called_once_with(
        model="gpt-5.4",
        messages=[{"role": "user", "content": "ping"}],
        max_completion_tokens=5,
        timeout=30,
    )


def test_openai_provider_ping_uses_max_tokens_for_legacy_models(monkeypatch):
    """OpenAIProvider.ping keeps max_tokens for pre GPT-5 models."""
    response = MagicMock()
    response.model = "gpt-4o-2024-08-06"

    completions = MagicMock()
    completions.create.return_value = response

    client = MagicMock()
    client.chat.completions = completions

    openai_module = MagicMock()
    openai_module.OpenAI.return_value = client
    monkeypatch.setitem(sys.modules, "openai", openai_module)

    provider = OpenAIProvider("test-key", model="gpt-4o")
    result = provider.ping()

    assert result["ok"] is True
    completions.create.assert_called_once_with(
        model="gpt-4o",
        messages=[{"role": "user", "content": "ping"}],
        max_tokens=5,
        timeout=30,
    )


def test_openai_provider_tool_probe_forces_and_validates_call(monkeypatch):
    response = MagicMock()
    response.model = "deepseek-v4-pro"
    response.model_dump.return_value = {
        "model": "deepseek-v4-pro",
        "choices": [
            {
                "message": {
                    "role": "assistant",
                    "content": None,
                    "tool_calls": [
                        {
                            "id": "probe-1",
                            "type": "function",
                            "function": {
                                "name": "chemsmart_doctor_probe",
                                "arguments": ('{"nonce": "chemsmart-doctor"}'),
                            },
                        }
                    ],
                },
                "finish_reason": "tool_calls",
            }
        ],
    }
    completions = MagicMock()
    completions.create.return_value = response
    client = MagicMock()
    client.chat.completions = completions
    openai_module = MagicMock()
    openai_module.OpenAI.return_value = client
    monkeypatch.setitem(sys.modules, "openai", openai_module)

    result = OpenAIProvider(
        "test-key",
        model="deepseek-v4-pro",
    ).tool_probe()

    assert result["ok"] is True
    assert result["resolved_model"] == "deepseek-v4-pro"
    assert result["tool"] == "chemsmart_doctor_probe"
    kwargs = completions.create.call_args.kwargs
    assert kwargs["tool_choice"]["function"]["name"] == (
        "chemsmart_doctor_probe"
    )
    assert (
        kwargs["tools"][0]["function"]["parameters"]["additionalProperties"]
        is False
    )
    assert kwargs["timeout"] == 30


def test_anthropic_provider_tool_probe_forces_and_validates_call(monkeypatch):
    from chemsmart.agent.providers import AnthropicProvider

    response = MagicMock()
    response.model = "claude-sonnet-4-6"
    response.model_dump.return_value = {
        "model": "claude-sonnet-4-6",
        "content": [
            {
                "type": "tool_use",
                "id": "probe-1",
                "name": "chemsmart_doctor_probe",
                "input": {"nonce": "chemsmart-doctor"},
            }
        ],
        "stop_reason": "tool_use",
    }
    client = MagicMock()
    client.messages.create.return_value = response
    anthropic_module = MagicMock()
    anthropic_module.Anthropic.return_value = client
    monkeypatch.setitem(sys.modules, "anthropic", anthropic_module)

    result = AnthropicProvider("test-key").tool_probe()

    assert result["ok"] is True
    assert result["resolved_model"] == "claude-sonnet-4-6"
    assert result["tool"] == "chemsmart_doctor_probe"
    kwargs = client.messages.create.call_args.kwargs
    assert kwargs["tool_choice"] == {
        "type": "tool",
        "name": "chemsmart_doctor_probe",
    }
    assert kwargs["tools"][0]["input_schema"]["additionalProperties"] is False
    assert kwargs["timeout"] == 30
