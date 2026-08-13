"""Provider-neutral ``agent.yaml`` configuration for Runtime V2."""

from __future__ import annotations

from dataclasses import dataclass
import os
from pathlib import Path
from typing import Any, Mapping

from chemsmart.agent._contracts import ContractError, canonical_sha256
from chemsmart.agent.api_access import require_provider_key_label
from chemsmart.agent.runtime.deepseek import (
    DEEPSEEK_OFFICIAL_ENDPOINT,
    DEEPSEEK_V4_FLASH_CONTEXT_TOKENS,
    DEEPSEEK_V4_FLASH_MAX_OUTPUT_TOKENS,
    DEEPSEEK_V4_FLASH_MODEL,
    DeepSeekV4FlashConfigV1,
)
from chemsmart.agent.runtime.transport import ProviderTurnDeadlinesV1
from chemsmart.io.yaml import YAMLFile

ALIBABA_TOKEN_PLAN_PROVIDER = "alibaba-token-plan"
ALIBABA_TOKEN_PLAN_ENDPOINT = (
    "https://token-plan.ap-southeast-1.maas.aliyuncs.com/compatible-mode/v1"
)
ALIBABA_TOKEN_PLAN_MODEL = "qwen3.8-max"
ALIBABA_TOKEN_PLAN_CONTEXT_TOKENS = 1_000_000
ALIBABA_TOKEN_PLAN_MAX_OUTPUT_TOKENS = 262_144


@dataclass(frozen=True)
class AgentProviderProfileV1:
    """One public provider profile; credentials remain in the secret file.

    V1 records remain byte-identical and receive the immutable default
    deadline policy at runtime.  Explicit deadline profiles use the V2
    subclass below so historical dataclass projections also remain identical.
    """

    schema_version: str
    profile_name: str
    provider: str
    wire_protocol: str
    api_key_env: str
    model: str
    endpoint: str
    reasoning_effort: str
    preserve_thinking: bool
    context_tokens: int
    max_output_tokens: int
    profile_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.agent-provider-profile.v1":
            raise ContractError("unsupported agent provider profile schema")
        body = {
            key: value
            for key, value in self.__dict__.items()
            if key != "profile_sha256"
        }
        if self.profile_sha256 != canonical_sha256(body):
            raise ContractError("agent provider profile digest mismatch")

    def runtime_config(self):
        """Build the provider-native ephemeral continuation configuration."""

        turn_deadlines = getattr(
            self, "transport_deadlines", ProviderTurnDeadlinesV1()
        )

        if self.provider == ALIBABA_TOKEN_PLAN_PROVIDER:
            from chemsmart.agent.runtime.alibaba import (
                AlibabaTokenPlanConfigV1,
            )

            return AlibabaTokenPlanConfigV1(
                model=self.model,
                endpoint=self.endpoint,
                reasoning_effort=self.reasoning_effort,
                preserve_thinking=self.preserve_thinking,
                max_output_tokens=self.max_output_tokens,
                turn_deadlines=turn_deadlines,
            )
        if self.provider == "deepseek":
            return DeepSeekV4FlashConfigV1(
                model=self.model,
                endpoint=self.endpoint,
                reasoning_effort=self.reasoning_effort,
                max_output_tokens=self.max_output_tokens,
                turn_deadlines=turn_deadlines,
            )
        raise ContractError("provider profile has no runtime adapter")


@dataclass(frozen=True)
class AgentProviderProfileV2(AgentProviderProfileV1):
    """Provider profile whose explicit transport policy is digest-bound."""

    transport_deadlines: ProviderTurnDeadlinesV1

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.agent-provider-profile.v2":
            raise ContractError("unsupported agent provider profile schema")
        body = {
            key: value
            for key, value in self.__dict__.items()
            if key not in {"profile_sha256", "transport_deadlines"}
        }
        body["transport_deadlines"] = (
            self.transport_deadlines.configuration_record()
        )
        if self.profile_sha256 != canonical_sha256(body):
            raise ContractError("agent provider profile digest mismatch")


@dataclass(frozen=True)
class AgentProviderSelectionV1:
    """Active profile plus ordered, explicitly attributed fallbacks."""

    schema_version: str
    active_profile: AgentProviderProfileV1
    fallback_profiles: tuple[AgentProviderProfileV1, ...]
    selection_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.agent-provider-selection.v1":
            raise ContractError("unsupported agent provider selection schema")
        names = (
            self.active_profile.profile_name,
            *(item.profile_name for item in self.fallback_profiles),
        )
        if len(names) != len(set(names)):
            raise ContractError(
                "provider selection contains duplicate profiles"
            )
        body = {
            "schema_version": self.schema_version,
            "active_profile_sha256": self.active_profile.profile_sha256,
            "fallback_profile_sha256s": tuple(
                item.profile_sha256 for item in self.fallback_profiles
            ),
        }
        if self.selection_sha256 != canonical_sha256(body):
            raise ContractError("agent provider selection digest mismatch")


def default_agent_config_path() -> Path:
    explicit = os.environ.get("CHEMSMART_AGENT_CONFIG", "").strip()
    return (
        Path(explicit).expanduser()
        if explicit
        else Path.home() / ".chemsmart" / "agent" / "agent.yaml"
    )


def load_agent_provider_selection(
    path: str | Path | None = None,
    *,
    requested_profile: str | None = None,
) -> AgentProviderSelectionV1:
    """Load one active profile without reading or embedding its credential."""

    config_path = (
        Path(path).expanduser() if path else default_agent_config_path()
    )
    if not config_path.is_file() or config_path.is_symlink():
        raise ContractError(
            "agent provider config must be a regular YAML file"
        )
    payload = YAMLFile(str(config_path)).yaml_contents_dict
    if not isinstance(payload, Mapping):
        raise ContractError("agent provider config must be a YAML mapping")
    providers = payload.get("providers")
    if not isinstance(providers, Mapping) or not providers:
        raise ContractError("agent provider config requires providers")

    selected_name = str(
        requested_profile or payload.get("active") or ""
    ).strip()
    if not selected_name:
        raise ContractError("agent provider config requires an active profile")
    if selected_name not in providers or not isinstance(
        providers[selected_name], Mapping
    ):
        raise ContractError(
            "requested agent provider profile is not configured"
        )

    raw_fallbacks = payload.get("fallback", payload.get("fallbacks", ()))
    if raw_fallbacks is None:
        fallback_names: tuple[str, ...] = ()
    elif isinstance(raw_fallbacks, str):
        fallback_names = (
            (raw_fallbacks.strip(),) if raw_fallbacks.strip() else ()
        )
    elif isinstance(raw_fallbacks, (list, tuple)) and all(
        isinstance(item, str) and item.strip() for item in raw_fallbacks
    ):
        fallback_names = tuple(item.strip() for item in raw_fallbacks)
    else:
        raise ContractError("agent provider fallbacks must be profile names")
    missing = tuple(
        name
        for name in fallback_names
        if name not in providers or not isinstance(providers[name], Mapping)
    )
    if missing:
        raise ContractError(
            "agent provider fallback profile is not configured"
        )
    # Only profiles that can enter this selection are part of the Runtime V2
    # trust boundary.  Idle profiles may belong to another client or legacy
    # workflow; inspecting them made a valid explicitly selected profile fail
    # for unrelated configuration.  Selected and fallback profiles still fail
    # closed on literal credentials.
    selected_entries = {
        name: providers[name]
        for name in (selected_name, *fallback_names)
        if name in providers
    }
    _reject_literal_secrets(selected_entries)
    selected_profile = _build_profile(selected_name, providers[selected_name])
    fallback_profiles = tuple(
        _build_profile(name, providers[name])
        for name in fallback_names
        if name != selected_name
    )
    body = {
        "schema_version": "chemsmart.agent-provider-selection.v1",
        "active_profile_sha256": selected_profile.profile_sha256,
        "fallback_profile_sha256s": tuple(
            item.profile_sha256 for item in fallback_profiles
        ),
    }
    return AgentProviderSelectionV1(
        schema_version=body["schema_version"],
        active_profile=selected_profile,
        fallback_profiles=fallback_profiles,
        selection_sha256=canonical_sha256(body),
    )


def _reject_literal_secrets(providers: Mapping[str, Any]) -> None:
    """Reject credentials in profiles selected for the active fallback chain.

    The check intentionally inspects field names and presence only; a secret
    value is never copied into an exception or receipt.  Unselected profiles
    are outside this loader's authority and are not parsed as Runtime V2 data.
    """

    secret_fields = {
        "api_key",
        "access_token",
        "token",
        "secret",
        "password",
        "hf_token",
    }
    secret_headers = {"authorization", "api-key", "x-api-key"}
    for profile_name, raw_profile in providers.items():
        if not isinstance(raw_profile, Mapping):
            continue
        if any(
            str(raw_profile.get(field) or "").strip()
            for field in secret_fields
        ):
            raise ContractError(
                "agent.yaml must not contain literal API keys or other "
                "secrets; use an environment label in every provider profile"
            )
        headers = raw_profile.get("extra_headers")
        if isinstance(headers, Mapping) and any(
            str(name).strip().lower() in secret_headers
            and str(value or "").strip()
            for name, value in headers.items()
        ):
            raise ContractError(
                "agent.yaml contains a credential header; use a short-lived "
                "provider lease"
            )


def _build_profile(
    profile_name: str, entry: Mapping[str, Any]
) -> AgentProviderProfileV1:
    provider_type = str(entry.get("type") or "").strip().lower()
    if provider_type not in {"openai", "openai-chat-completions"}:
        raise ContractError("active Runtime V2 profiles require OpenAI chat")
    if str(entry.get("api_key") or "").strip():
        raise ContractError("agent.yaml must not contain literal API keys")
    api_key_env = str(entry.get("api_key_env") or "").strip()
    model = str(entry.get("model") or "").strip()
    endpoint = str(
        entry.get("base_url") or entry.get("endpoint") or ""
    ).rstrip("/")
    reasoning_effort = str(entry.get("reasoning_effort") or "").strip().lower()
    preserve_thinking = entry.get("preserve_thinking", True)
    if not isinstance(preserve_thinking, bool):
        raise ContractError("preserve_thinking must be boolean")
    raw_deadlines = entry.get("transport_deadlines")
    if raw_deadlines is None:
        transport_deadlines = None
    elif not isinstance(raw_deadlines, Mapping):
        raise ContractError("transport_deadlines must be a mapping")
    else:
        required_deadlines = {
            "connect_seconds",
            "first_event_seconds",
            "inter_event_seconds",
            "absolute_turn_seconds",
        }
        if set(raw_deadlines) != required_deadlines:
            raise ContractError(
                "transport_deadlines requires exactly connect_seconds, "
                "first_event_seconds, inter_event_seconds, and "
                "absolute_turn_seconds"
            )
        try:
            transport_deadlines = ProviderTurnDeadlinesV1(
                connect_seconds=float(raw_deadlines["connect_seconds"]),
                first_event_seconds=float(
                    raw_deadlines["first_event_seconds"]
                ),
                inter_event_seconds=float(
                    raw_deadlines["inter_event_seconds"]
                ),
                absolute_seconds=float(
                    raw_deadlines["absolute_turn_seconds"]
                ),
            )
        except (TypeError, ValueError) as exc:
            raise ContractError(
                "transport deadline values must be numbers"
            ) from exc

    if endpoint == ALIBABA_TOKEN_PLAN_ENDPOINT:
        provider = ALIBABA_TOKEN_PLAN_PROVIDER
        if not model:
            raise ContractError("Alibaba Token Plan profile requires a model")
        if model == "qwen3.8-max-preview":
            raise ContractError(
                "Alibaba Token Plan profile requires production qwen3.8-max"
            )
        require_provider_key_label(api_key_env, provider=provider)
        reasoning_effort = reasoning_effort or (
            "xhigh" if model == ALIBABA_TOKEN_PLAN_MODEL else "max"
        )
        if model == ALIBABA_TOKEN_PLAN_MODEL and reasoning_effort != "xhigh":
            raise ContractError(
                "qwen3.8-max maximum reasoning requires xhigh"
            )
        if reasoning_effort not in {"high", "max", "xhigh"}:
            raise ContractError(
                "Alibaba Token Plan reasoning effort must be high, max or "
                "xhigh"
            )
        if not preserve_thinking:
            raise ContractError(
                "Alibaba Token Plan tool continuation must be preserved"
            )
        context_tokens = ALIBABA_TOKEN_PLAN_CONTEXT_TOKENS
        max_output_tokens = ALIBABA_TOKEN_PLAN_MAX_OUTPUT_TOKENS
    elif endpoint == DEEPSEEK_OFFICIAL_ENDPOINT:
        provider = "deepseek"
        if model != DEEPSEEK_V4_FLASH_MODEL:
            raise ContractError("DeepSeek fallback requires deepseek-v4-flash")
        require_provider_key_label(api_key_env, provider=provider)
        reasoning_effort = reasoning_effort or "max"
        if reasoning_effort not in {"high", "max"}:
            raise ContractError(
                "DeepSeek reasoning effort must be high or max"
            )
        context_tokens = DEEPSEEK_V4_FLASH_CONTEXT_TOKENS
        max_output_tokens = DEEPSEEK_V4_FLASH_MAX_OUTPUT_TOKENS
    else:
        raise ContractError("agent provider endpoint is not registered")

    body = {
        "schema_version": (
            "chemsmart.agent-provider-profile.v2"
            if transport_deadlines is not None
            else "chemsmart.agent-provider-profile.v1"
        ),
        "profile_name": profile_name,
        "provider": provider,
        "wire_protocol": "openai-chat-completions",
        "api_key_env": api_key_env,
        "model": model,
        "endpoint": endpoint,
        "reasoning_effort": reasoning_effort,
        "preserve_thinking": preserve_thinking,
        "context_tokens": context_tokens,
        "max_output_tokens": max_output_tokens,
    }
    if transport_deadlines is not None:
        body["transport_deadlines"] = (
            transport_deadlines.configuration_record()
        )
    constructor = dict(body)
    constructor.pop("transport_deadlines", None)
    profile_type = (
        AgentProviderProfileV2
        if transport_deadlines is not None
        else AgentProviderProfileV1
    )
    values = {
        **constructor,
        "profile_sha256": canonical_sha256(body),
    }
    if transport_deadlines is not None:
        values["transport_deadlines"] = transport_deadlines
    return profile_type(**values)


__all__ = [
    "ALIBABA_TOKEN_PLAN_CONTEXT_TOKENS",
    "ALIBABA_TOKEN_PLAN_ENDPOINT",
    "ALIBABA_TOKEN_PLAN_MAX_OUTPUT_TOKENS",
    "ALIBABA_TOKEN_PLAN_MODEL",
    "ALIBABA_TOKEN_PLAN_PROVIDER",
    "AgentProviderProfileV1",
    "AgentProviderProfileV2",
    "AgentProviderSelectionV1",
    "default_agent_config_path",
    "load_agent_provider_selection",
]
