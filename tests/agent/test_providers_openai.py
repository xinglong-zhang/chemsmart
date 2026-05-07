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
    openai_module.OpenAI.assert_called_once_with(api_key="test-key")
    completions.create.assert_called_once_with(
        model="gpt-4o-mini",
        messages=messages,
        tools=tools,
    )
