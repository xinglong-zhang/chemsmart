"""
Provider adapters for chemsmart agent.

Reads api.env via python-dotenv (key: ai_api_key).
Dispatches on AI_PROVIDER env var; v1 supports Anthropic only.
"""

import os
from typing import Optional

from dotenv import load_dotenv

_API_ENV_PATH = "/Users/hongjiseung/developer/chemsmart/api.env"
_SUPPORTED = frozenset({"anthropic"})


class ProviderError(Exception):
    pass


class AnthropicProvider:
    name = "anthropic"

    def __init__(self, api_key: str) -> None:
        import anthropic

        self._client = anthropic.Anthropic(api_key=api_key)

    def chat(
        self,
        messages: list,
        tools: Optional[list] = None,
    ) -> dict:
        kwargs: dict = {
            "model": "claude-opus-4-5",
            "max_tokens": 4096,
            "messages": messages,
        }
        if tools:
            kwargs["tools"] = tools
        response = self._client.messages.create(**kwargs)
        return response.model_dump()


def get_provider(
    env_path: Optional[str] = None,
) -> AnthropicProvider:
    """Return a configured provider instance; raises ProviderError on failure.

    Validates AI_PROVIDER and ai_api_key before constructing the provider.
    """
    if env_path is None:
        env_path = _API_ENV_PATH

    provider_name = os.environ.get("AI_PROVIDER")
    if not provider_name or not provider_name.strip():
        raise ProviderError("AI_PROVIDER env var is not set")
    provider_name = provider_name.strip()

    if provider_name not in _SUPPORTED:
        raise ProviderError(
            f"AI_PROVIDER={provider_name!r} is not supported; "
            f"supported: {sorted(_SUPPORTED)}"
        )

    load_dotenv(env_path, override=True)
    api_key = os.environ.get("ai_api_key", "").strip()
    if not api_key:
        raise ProviderError("api.env: ai_api_key is empty or missing")

    if provider_name == "anthropic":
        return AnthropicProvider(api_key)

    raise ProviderError(f"AI_PROVIDER={provider_name!r} is not supported")
