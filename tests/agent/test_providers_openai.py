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
        max_tokens=5,
    )
