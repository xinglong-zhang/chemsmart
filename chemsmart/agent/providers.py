"""
Provider adapters for chemsmart agent.

Reads api.env via python-dotenv (key: ai_api_key).
Dispatches on AI_PROVIDER env var; v1 supports Anthropic and OpenAI.
"""

from __future__ import annotations

import json
import os
import re
import time
import warnings
from pathlib import Path
from typing import Any, Optional

from dotenv import load_dotenv

from chemsmart.agent.provider_config import (
    AgentProviderConfigError,
    load_active_provider_config,
)

_GATEWAY_URL_OPENAI = "https://factchat-cloud.mindlogic.ai/v1/gateway"
_GATEWAY_URL_ANTHROPIC = (
    "https://factchat-cloud.mindlogic.ai/v1/gateway/claude"
)
_SUPPORTED = frozenset({"openai", "anthropic"})
_PING_MESSAGES: Any = [{"role": "user", "content": "ping"}]
DEFAULT_TIMEOUT_S = 30
_OPENAI_USES_MCT = re.compile(r"^(gpt-5|o1|o3|o4)")


def _openai_uses_max_completion_tokens(model: str | None) -> bool:
    """Return whether an OpenAI-family model requires max_completion_tokens."""
    return bool(_OPENAI_USES_MCT.match((model or "").strip().lower()))


class ProviderError(Exception):
    pass


class AnthropicProvider:
    name = "anthropic"
    wire_protocol = "anthropic"
    default_model = "claude-sonnet-4-6"
    gateway_url = _GATEWAY_URL_ANTHROPIC

    def __init__(
        self,
        api_key: str,
        model: str | None = None,
        base_url: str | None = None,
        extra_headers: dict[str, str] | None = None,
    ) -> None:
        import anthropic

        self.default_model = model or type(self).default_model
        client_kwargs: dict[str, Any] = {
            "api_key": api_key,
            "base_url": base_url or self.gateway_url,
        }
        if extra_headers:
            client_kwargs["default_headers"] = extra_headers
        self._client = anthropic.Anthropic(**client_kwargs)

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
    wire_protocol = "openai"
    default_model = "gpt-5.4"
    gateway_url = _GATEWAY_URL_OPENAI

    def __init__(
        self,
        api_key: str,
        model: str | None = None,
        base_url: str | None = None,
        extra_headers: dict[str, str] | None = None,
        provider_name: str | None = None,
    ) -> None:
        import openai

        self.name = (provider_name or type(self).name).strip()
        self.default_model = model or type(self).default_model
        client_kwargs: dict[str, Any] = {
            "api_key": api_key,
            "base_url": base_url or self.gateway_url,
        }
        if extra_headers:
            client_kwargs["default_headers"] = extra_headers
        self._client = openai.OpenAI(**client_kwargs)

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
            kwargs: dict[str, Any] = {
                "model": self.default_model,
                "messages": _PING_MESSAGES,
                "timeout": DEFAULT_TIMEOUT_S,
            }
            if _openai_uses_max_completion_tokens(self.default_model):
                kwargs["max_completion_tokens"] = 5
            else:
                kwargs["max_tokens"] = 5
            response = self._client.chat.completions.create(**kwargs)
        except Exception as exc:
            raise ProviderError(f"ping failed: {exc}") from exc

        return {
            "ok": True,
            "resolved_model": _resolve_model(response, self.default_model),
            "latency_ms": _latency_ms(started),
        }


class LocalProvider:
    """In-process v13.1 local provider for ``type: local`` agent.yaml entries.

    Loads :mod:`chemsmart.agent.local` lazily on the first ``chat`` call so
    ``chemsmart agent doctor`` / ``chemsmart config agent`` stay fast and do
    not require torch/peft at import time.
    """

    name = "local"
    wire_protocol = "openai"
    default_model = "chemsmart-qwen2.5-coder-3b-instruct-v13_1"

    def __init__(
        self,
        api_key: str = "",
        model: str | None = None,
        base_url: str | None = None,
        extra_headers: dict[str, str] | None = None,
        base_model_id: str = "",
        adapter_repo_id: str = "",
        runtime: str = "",
        project: str = "",
    ) -> None:
        del base_url, extra_headers
        self.default_model = model or type(self).default_model
        self._hf_token = api_key or os.environ.get("HF_TOKEN", "")
        self._base_model_id = base_model_id
        self._adapter_repo_id = adapter_repo_id
        self._runtime = runtime or None
        self._project = project
        self._bundle: Any = None

    def _ensure_loaded(self) -> Any:
        if self._bundle is not None:
            return self._bundle
        if self._runtime == "mlx":
            from chemsmart.agent.local.mlx_loader import (
                MODEL_REPO_ID,
                load_mlx_model,
            )

            self._bundle = load_mlx_model(
                model_id=self._base_model_id or MODEL_REPO_ID,
                hf_token=self._hf_token or None,
            )
            return self._bundle

        from chemsmart.agent.local.loader import (
            ADAPTER_REPO_ID,
            BASE_MODEL_ID,
            load_lora_model,
        )

        self._bundle = load_lora_model(
            base_model_id=self._base_model_id or BASE_MODEL_ID,
            adapter_repo_id=self._adapter_repo_id or ADAPTER_REPO_ID,
            hf_token=self._hf_token or None,
            runtime=self._runtime,
        )
        return self._bundle

    def chat(
        self,
        messages: list,
        tools: Optional[list] = None,
        timeout_s: float = DEFAULT_TIMEOUT_S,
    ) -> dict:
        del tools, timeout_s
        from chemsmart.agent.local.adapter import plan_to_synthesis_result
        from chemsmart.agent.local.generator import generate_plan

        user_query = ""
        history: list[dict[str, str]] = []
        for msg in messages:
            if not isinstance(msg, dict):
                continue
            role = str(msg.get("role", ""))
            content = msg.get("content", "")
            if not isinstance(content, str):
                content = str(content)
            if role == "user":
                if user_query:
                    history.append({"role": "user", "content": user_query})
                user_query = content
            elif role == "assistant":
                history.append({"role": "assistant", "content": content})
        if not user_query:
            raise ProviderError("local provider requires a user message")

        bundle = self._ensure_loaded()
        try:
            plan = generate_plan(bundle, user_query, history=history)
        except ValueError as exc:
            raise ProviderError(f"local provider decode failed: {exc}") from exc
        from chemsmart.agent.kind_disambiguator import disambiguate

        plan, _changed = disambiguate(user_query, plan)
        result = plan_to_synthesis_result(
            plan,
            user_query,
            default_project=self._project or None,
        )
        return {
            "id": "local-completion",
            "model": self.default_model,
            "choices": [
                {
                    "index": 0,
                    "message": {
                        "role": "assistant",
                        "content": json.dumps(result, ensure_ascii=False),
                    },
                    "finish_reason": "stop",
                }
            ],
            # The model's own SPEC (not the adapted status/command result) is what it
            # must see as assistant history for multi-turn edits. The session replays
            # this verbatim so follow-ups ("change the server", "make it a ts") carry
            # context, matching the format the model was prompted/generalizes with.
            "raw_plan": json.dumps(plan, ensure_ascii=False),
            "usage": {"prompt_tokens": 0, "completion_tokens": 0},
        }

    def ping(self) -> dict[str, Any]:
        started = time.perf_counter()
        if self._runtime == "mlx":
            return self._ping_mlx(started)

        try:
            from chemsmart.agent.local.loader import BASE_MODEL_ID
        except Exception as exc:  # pragma: no cover - import guard
            raise ProviderError(f"local provider ping failed: {exc}") from exc

        try:
            from transformers import AutoTokenizer  # type: ignore[import-not-found]
        except ImportError as exc:
            raise ProviderError(
                "local provider requires transformers; PEFT/bitsandbytes are "
                "only needed when adapter_repo_id is set. "
                "Install with: pip install 'huggingface_hub>=0.34.0,<1.0' "
                "'transformers==4.56.2' 'accelerate==1.10.0'"
            ) from exc

        try:
            AutoTokenizer.from_pretrained(
                self._base_model_id or BASE_MODEL_ID,
                token=(self._hf_token or None),
            )
        except Exception as exc:
            raise ProviderError(
                f"local provider ping failed (cache missing? run "
                f"`chemsmart config agent` to pre-fetch). underlying: {exc}"
            ) from exc
        return {
            "ok": True,
            "resolved_model": self.default_model,
            "latency_ms": _latency_ms(started),
        }

    def _ping_mlx(self, started: float) -> dict[str, Any]:
        try:
            from chemsmart.agent.local.mlx_loader import MODEL_REPO_ID
            from huggingface_hub import HfApi
        except Exception as exc:  # pragma: no cover - optional runtime guard
            raise ProviderError(f"MLX local provider ping failed: {exc}") from exc

        repo_id = self._base_model_id or MODEL_REPO_ID
        token = self._hf_token or None
        try:
            HfApi(token=token).model_info(repo_id)
        except Exception as exc:
            raise ProviderError(
                f"MLX local provider ping failed for {repo_id!r}: {exc}"
            ) from exc
        return {
            "ok": True,
            "resolved_model": self.default_model,
            "latency_ms": _latency_ms(started),
        }


def _resolve_api_env_path(explicit: str | None) -> str | None:
    if isinstance(explicit, str) and explicit.strip():
        return explicit.strip()

    env_path = os.environ.get("CHEMSMART_API_ENV")
    if isinstance(env_path, str) and env_path.strip():
        return env_path.strip()

    candidates = (
        Path.home() / ".chemsmart" / "api.env",
        Path.cwd() / "api.env",
    )
    for candidate in candidates:
        if candidate.is_file():
            return str(candidate)

    return None


def get_provider(
    env_path: Optional[str] = None,
) -> AnthropicProvider | OpenAIProvider | LocalProvider:
    """Return a configured provider instance; raises ProviderError on failure.

    Prefer ``~/.chemsmart/agent/agent.yaml``. If it is absent, fall back to the
    legacy ``api.env``/``AI_PROVIDER`` path for backwards compatibility.
    """
    if isinstance(env_path, str) and env_path.strip():
        load_dotenv(env_path.strip(), override=False)

    try:
        config = load_active_provider_config()
    except AgentProviderConfigError as exc:
        raise ProviderError(str(exc)) from exc

    if config is not None:
        if config.type == "anthropic":
            return AnthropicProvider(
                config.api_key,
                model=config.model,
                base_url=config.base_url or None,
                extra_headers=config.extra_headers,
            )
        if config.type == "openai":
            return OpenAIProvider(
                config.api_key,
                model=config.model,
                base_url=config.base_url or None,
                extra_headers=config.extra_headers,
                provider_name=config.name,
            )
        if config.type == "local":
            return LocalProvider(
                api_key=config.hf_token or config.api_key,
                model=config.model,
                base_model_id=config.base_model_id,
                adapter_repo_id=config.adapter_repo_id,
                runtime=config.runtime,
                project=config.project,
            )
        raise ProviderError(f"provider type {config.type!r} is not supported")

    warnings.warn(
        "agent.yaml missing; legacy api.env will be removed. "
        "Run `make configure`.",
        DeprecationWarning,
        stacklevel=2,
    )

    env_path = _resolve_api_env_path(env_path)

    # Load api.env first so AI_PROVIDER and ai_api_key can both come from the
    # file — the user should not need to export AI_PROVIDER in their shell.
    if env_path is not None:
        load_dotenv(env_path, override=False)

    provider_name = os.environ.get("AI_PROVIDER")
    if not provider_name or not provider_name.strip():
        raise ProviderError(
            "AI_PROVIDER is not set. "
            "Add AI_PROVIDER=openai to api.env or export it in your shell."
        )
    provider_name = provider_name.strip()

    if provider_name not in _SUPPORTED:
        raise ProviderError(
            f"AI_PROVIDER={provider_name!r} is not supported; "
            f"supported: {sorted(_SUPPORTED)}"
        )

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
