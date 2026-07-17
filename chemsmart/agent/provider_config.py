"""Configuration loader for chemsmart agent providers."""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

from chemsmart.io.yaml import YAMLFile

_SUPPORTED_PROVIDER_TYPES = frozenset({"openai", "anthropic", "local"})


@dataclass
class AgentProviderConfig:
    """Resolved agent provider configuration."""

    name: str
    type: str
    api_key: str
    model: str
    base_url: str
    extra_headers: dict[str, str]
    base_model_id: str = ""
    adapter_repo_id: str = ""
    hf_token: str = ""
    runtime: str = ""
    project: str = ""


class AgentProviderConfigError(Exception):
    """Raised when ``~/.chemsmart/agent/agent.yaml`` is invalid."""


def _load_legacy_env(path: Path) -> None:
    """Load the one-release ``api.env`` compatibility path lazily."""
    try:
        from dotenv import load_dotenv
    except ImportError as exc:
        raise AgentProviderConfigError(
            "Legacy api.env loading requires the agent extra. Install with "
            "`pip install 'chemsmart[agent]'`, or move the provider settings "
            "to ~/.chemsmart/agent/agent.yaml."
        ) from exc
    load_dotenv(path, override=False)


def _default_yaml_path() -> Path:
    """Return the default agent provider YAML path."""
    return Path.home() / ".chemsmart" / "agent" / "agent.yaml"


def _load_provider_environment() -> Path | None:
    """Load the first configured agent environment file, if one exists."""
    explicit = os.environ.get("CHEMSMART_API_ENV", "").strip()
    candidates = (
        Path(explicit) if explicit else None,
        Path.home() / ".chemsmart" / "api.env",
        Path.cwd() / "api.env",
    )
    for candidate in candidates:
        if candidate is not None and candidate.is_file():
            _load_legacy_env(candidate)
            return candidate
    return None


def load_active_provider_config(
    yaml_path: str | os.PathLike[str] | None = None,
) -> AgentProviderConfig | None:
    """Load and resolve the active agent provider configuration.

    Args:
        yaml_path: Optional path to an ``agent.yaml`` file. Defaults to
            ``~/.chemsmart/agent/agent.yaml``.

    Returns:
        The resolved active provider config, or ``None`` when the YAML file is
        absent.

    Raises:
        AgentProviderConfigError: If the YAML is malformed, lacks an active
            provider, uses an unknown provider type, or does not resolve an API
            key.
    """
    _load_provider_environment()
    path = Path(yaml_path) if yaml_path is not None else _default_yaml_path()
    if not path.is_file():
        return None

    try:
        payload = YAMLFile(str(path)).yaml_contents_dict
    except (OSError, yaml.YAMLError) as exc:
        raise AgentProviderConfigError(
            f"Failed to read agent provider config {path}: {exc}"
        ) from exc

    if not isinstance(payload, dict):
        raise AgentProviderConfigError(
            f"agent provider config {path} must be a YAML mapping"
        )

    active = payload.get("active")
    if not isinstance(active, str) or not active.strip():
        raise AgentProviderConfigError(
            f"agent provider config {path} is missing non-empty 'active'"
        )
    active = active.strip()

    providers = payload.get("providers")
    if not isinstance(providers, dict):
        raise AgentProviderConfigError(
            f"agent provider config {path} is missing 'providers' mapping"
        )

    provider_entry = providers.get(active)
    if not isinstance(provider_entry, dict):
        raise AgentProviderConfigError(
            f"active provider {active!r} is not defined in {path}"
        )

    provider_type = _string_field(provider_entry, "type").lower()
    if provider_type not in _SUPPORTED_PROVIDER_TYPES:
        raise AgentProviderConfigError(
            f"provider {active!r} has unsupported type {provider_type!r}; "
            f"supported: {sorted(_SUPPORTED_PROVIDER_TYPES)}"
        )

    model = _string_field(provider_entry, "model")
    base_url = _optional_string_field(provider_entry, "base_url")
    extra_headers = _extra_headers(provider_entry, active)
    base_model_id = _optional_string_field(provider_entry, "base_model_id")
    adapter_repo_id = _optional_string_field(provider_entry, "adapter_repo_id")
    runtime = _optional_string_field(provider_entry, "runtime")
    project = _optional_string_field(provider_entry, "project")

    if provider_type == "local":
        api_key = _resolve_local_hf_token(provider_entry)
    else:
        api_key = _resolve_api_key(provider_entry, active)

    return AgentProviderConfig(
        name=active,
        type=provider_type,
        api_key=api_key,
        model=model,
        base_url=base_url,
        extra_headers=extra_headers,
        base_model_id=base_model_id,
        adapter_repo_id=adapter_repo_id,
        hf_token=api_key if provider_type == "local" else "",
        runtime=runtime,
        project=project,
    )


def _resolve_local_hf_token(entry: dict[str, Any]) -> str:
    """Return the HF token for a ``type: local`` provider.

    Order: literal ``hf_token`` → ``hf_token_env`` → ``HF_TOKEN`` env var →
    empty string (loader will raise if the model is not already cached).
    """
    literal = _optional_string_field(entry, "hf_token")
    if literal:
        return literal
    token_env = _optional_string_field(entry, "hf_token_env") or "HF_TOKEN"
    return os.environ.get(token_env, "").strip()


def _string_field(entry: dict[str, Any], key: str) -> str:
    value = entry.get(key)
    if not isinstance(value, str) or not value.strip():
        raise AgentProviderConfigError(
            f"provider entry is missing non-empty {key!r}"
        )
    return value.strip()


def _optional_string_field(entry: dict[str, Any], key: str) -> str:
    value = entry.get(key, "")
    if value is None:
        return ""
    if not isinstance(value, str):
        raise AgentProviderConfigError(f"provider field {key!r} must be text")
    return value.strip()


def _resolve_api_key(entry: dict[str, Any], provider_name: str) -> str:
    literal_api_key = _optional_string_field(entry, "api_key")
    if literal_api_key:
        return literal_api_key

    api_key_env = _optional_string_field(entry, "api_key_env")
    if api_key_env:
        env_api_key = os.environ.get(api_key_env, "").strip()
        if env_api_key:
            return env_api_key

    raise AgentProviderConfigError(
        f"provider {provider_name!r} must set api_key or api_key_env"
    )


def _extra_headers(
    entry: dict[str, Any], provider_name: str
) -> dict[str, str]:
    raw_headers = entry.get("extra_headers", {})
    if raw_headers is None:
        return {}
    if not isinstance(raw_headers, dict):
        raise AgentProviderConfigError(
            f"provider {provider_name!r} extra_headers must be a mapping"
        )
    headers: dict[str, str] = {}
    for key, value in raw_headers.items():
        if not isinstance(key, str) or not isinstance(value, str):
            raise AgentProviderConfigError(
                f"provider {provider_name!r} extra_headers keys and values "
                "must be strings"
            )
        headers[key] = value
    return headers
