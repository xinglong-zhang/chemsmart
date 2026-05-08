"""
Provider adapters for chemsmart agent.

Reads api.env via python-dotenv (key: ai_api_key).
Dispatches on AI_PROVIDER env var; v1 supports Anthropic and OpenAI.
"""

from __future__ import annotations

import os
import time
from typing import Any, Optional

from dotenv import load_dotenv

_API_ENV_PATH = "/Users/hongjiseung/developer/chemsmart/api.env"
_GATEWAY_URL_OPENAI = "https://factchat-cloud.mindlogic.ai/v1/gateway"
_GATEWAY_URL_ANTHROPIC = (
    "https://factchat-cloud.mindlogic.ai/v1/gateway/claude"
)
_AVAILABLE_MODELS = {
    "openai": [
        "gpt-5.4",
        "gpt-5.4-mini",
        "gpt-5.4-nano",
        "gemini-3.1-pro-preview",
        "grok-4",
        "sonar-pro",
    ],
    "anthropic": ["claude-sonnet-4-6"],
}
_SUPPORTED = frozenset(_AVAILABLE_MODELS)
_PING_MESSAGES = [{"role": "user", "content": "ping"}]
DEFAULT_TIMEOUT_S = 30


class ProviderError(Exception):
    pass


class AnthropicProvider:
    name = "anthropic"
    default_model = "claude-sonnet-4-6"
    gateway_url = _GATEWAY_URL_ANTHROPIC

    def __init__(self, api_key: str) -> None:
        import anthropic

        self._client = anthropic.Anthropic(
            api_key=api_key,
            base_url=self.gateway_url,
        )

    def chat(
        self,
        messages: list,
        tools: Optional[list] = None,
        timeout_s: float = DEFAULT_TIMEOUT_S,
    ) -> dict:
        kwargs: dict[str, Any] = {
            "model": self.default_model,
            "max_tokens": 4096,
            "messages": messages,
            "timeout": timeout_s,
        }
        if tools:
            kwargs["tools"] = tools
        response = self._client.messages.create(**kwargs)
        return response.model_dump()

    def ping(self) -> dict[str, Any]:
        started = time.perf_counter()
        try:
            response = self._client.messages.create(
                model=self.default_model,
                max_tokens=5,
                messages=_PING_MESSAGES,
                timeout=DEFAULT_TIMEOUT_S,
            )
        except Exception as exc:
            raise ProviderError(f"ping failed: {exc}") from exc

        return {
            "ok": True,
            "resolved_model": _resolve_model(response, self.default_model),
            "latency_ms": _latency_ms(started),
        }


class OpenAIProvider:
    name = "openai"
    default_model = "gpt-5.4"
    gateway_url = _GATEWAY_URL_OPENAI

    def __init__(self, api_key: str) -> None:
        import openai

        self._client = openai.OpenAI(
            api_key=api_key,
            base_url=self.gateway_url,
        )

    def chat(
        self,
        messages: list,
        tools: Optional[list] = None,
        timeout_s: float = DEFAULT_TIMEOUT_S,
    ) -> dict:
        kwargs: dict[str, Any] = {
            "model": self.default_model,
            "messages": messages,
            "timeout": timeout_s,
        }
        if tools:
            kwargs["tools"] = tools
        response = self._client.chat.completions.create(**kwargs)
        return response.model_dump()

    def ping(self) -> dict[str, Any]:
        started = time.perf_counter()
        try:
            response = self._client.chat.completions.create(
                model=self.default_model,
                messages=_PING_MESSAGES,
                max_tokens=5,
                timeout=DEFAULT_TIMEOUT_S,
            )
        except Exception as exc:
            raise ProviderError(f"ping failed: {exc}") from exc

        return {
            "ok": True,
            "resolved_model": _resolve_model(response, self.default_model),
            "latency_ms": _latency_ms(started),
        }


def get_provider(
    env_path: Optional[str] = None,
) -> AnthropicProvider | OpenAIProvider:
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
    if provider_name == "openai":
        return OpenAIProvider(api_key)

    raise ProviderError(f"AI_PROVIDER={provider_name!r} is not supported")


def _resolve_model(response: Any, fallback: str) -> str:
    model = getattr(response, "model", None)
    if isinstance(model, str) and model.strip():
        return model

    if hasattr(response, "model_dump"):
        payload = response.model_dump()
        dumped_model = payload.get("model")
        if isinstance(dumped_model, str) and dumped_model.strip():
            return dumped_model

    return fallback


def _latency_ms(started: float) -> int:
    return max(0, int((time.perf_counter() - started) * 1000))


def extract_response_usage(response: Any) -> dict[str, int | None]:
    usage = None
    if isinstance(response, dict):
        usage = response.get("usage")
    elif hasattr(response, "model_dump"):
        usage = response.model_dump().get("usage")

    if not isinstance(usage, dict):
        return {"input_tokens": None, "output_tokens": None}

    input_tokens = usage.get("input_tokens", usage.get("prompt_tokens"))
    output_tokens = usage.get("output_tokens", usage.get("completion_tokens"))
    return {
        "input_tokens": (
            int(input_tokens) if isinstance(input_tokens, int) else None
        ),
        "output_tokens": (
            int(output_tokens) if isinstance(output_tokens, int) else None
        ),
    }
