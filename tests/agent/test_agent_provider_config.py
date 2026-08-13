from __future__ import annotations

from pathlib import Path

import pytest

from chemsmart.agent._contracts import ContractError
from chemsmart.agent._contracts import canonical_sha256
from chemsmart.agent.provider_config import (
    AgentProviderProfileV1,
    load_agent_provider_selection,
)
from chemsmart.agent.runtime.deepseek import DeepSeekV4FlashConfigV1


def _write_config(
    path: Path,
    *,
    model: str = "qwen3.8-max",
    reasoning_effort: str | None = None,
) -> Path:
    reasoning_effort = reasoning_effort or (
        "xhigh" if model == "qwen3.8-max" else "max"
    )
    path.write_text(
        f"""
active: alibaba-token-plan
fallback: []
providers:
  unrelated-local:
    type: local
    model: ignored
  unrelated-openrouter:
    type: openai
    api_key_env: OPENROUTER_API_KEY
    model: ignored
    base_url: https://openrouter.ai/api/v1
  alibaba-token-plan:
    type: openai
    api_key_env: ALIBABA_TOKEN_PLAN_KEY
    model: {model}
    context_tokens: 1000000
    max_output_tokens: 262144
    base_url: https://token-plan.ap-southeast-1.maas.aliyuncs.com/compatible-mode/v1
    reasoning_effort: {reasoning_effort}
    preserve_thinking: true
  deepseek:
    type: openai
    api_key_env: DEEPSEEK-api-key
    model: deepseek-v4-flash
    context_tokens: 1000000
    max_output_tokens: 384000
    base_url: https://api.deepseek.com
    reasoning_effort: max
    preserve_thinking: true
""".lstrip(),
        encoding="utf-8",
    )
    return path


def test_agent_yaml_selects_explicit_model_and_limits(tmp_path):
    selection = load_agent_provider_selection(_write_config(tmp_path / "agent.yaml"))

    assert selection.active_profile.profile_name == "alibaba-token-plan"
    assert selection.active_profile.provider == "alibaba-token-plan"
    assert selection.active_profile.model == "qwen3.8-max"
    assert selection.active_profile.reasoning_effort == "xhigh"
    assert selection.active_profile.preserve_thinking is True
    assert selection.active_profile.context_tokens == 1_000_000
    assert selection.active_profile.max_output_tokens == 262_144
    assert selection.fallback_profiles == ()
    runtime = selection.active_profile.runtime_config()
    assert runtime.model == "qwen3.8-max"
    assert runtime.context_tokens == 1_000_000
    assert runtime.max_output_tokens == 262_144


def test_historical_provider_profile_v1_digest_is_byte_compatible():
    body = {
        "schema_version": "chemsmart.agent-provider-profile.v1",
        "profile_name": "alibaba-token-plan",
        "provider": "alibaba-token-plan",
        "wire_protocol": "openai-chat-completions",
        "api_key_env": "ALIBABA_TOKEN_PLAN_KEY",
        "model": "qwen3.8-max",
        "endpoint": (
            "https://token-plan.ap-southeast-1.maas.aliyuncs.com/"
            "compatible-mode/v1"
        ),
        "reasoning_effort": "xhigh",
        "preserve_thinking": True,
        "context_tokens": 1_000_000,
        "max_output_tokens": 262_144,
    }

    profile = AgentProviderProfileV1(
        **body, profile_sha256=canonical_sha256(body)
    )

    assert profile.__dict__ == {
        **body,
        "profile_sha256": canonical_sha256(body),
    }


def test_agent_yaml_profile_override_selects_one_explicit_attempt(tmp_path):
    selection = load_agent_provider_selection(
        _write_config(tmp_path / "agent.yaml"), requested_profile="deepseek"
    )

    assert selection.active_profile.provider == "deepseek"
    assert selection.active_profile.context_tokens == 1_000_000
    assert selection.active_profile.max_output_tokens == 384_000
    assert selection.fallback_profiles == ()


def test_deepseek_adapter_uses_the_explicit_profile_model(tmp_path):
    path = _write_config(tmp_path / "agent.yaml")
    content = path.read_text(encoding="utf-8").replace(
        "model: deepseek-v4-flash\n    base_url: https://api.deepseek.com",
        "model: deepseek-reasoner\n    base_url: https://api.deepseek.com",
    )
    path.write_text(content, encoding="utf-8")

    profile = load_agent_provider_selection(
        path, requested_profile="deepseek"
    ).active_profile

    assert profile.model == "deepseek-reasoner"
    assert profile.runtime_config().model == "deepseek-reasoner"
    assert profile.runtime_config().context_tokens == 1_000_000
    assert profile.runtime_config().max_output_tokens == 384_000
    with pytest.raises(TypeError):
        DeepSeekV4FlashConfigV1()


def test_agent_yaml_accepts_arbitrary_explicit_token_plan_model(tmp_path):
    selection = load_agent_provider_selection(
        _write_config(
            tmp_path / "agent.yaml",
            model="vendor/deployment:model-2026-08",
            reasoning_effort="max",
        )
    )

    assert selection.active_profile.model == "vendor/deployment:model-2026-08"


def test_agent_yaml_accepts_catalog_model_on_token_plan(tmp_path):
    selection = load_agent_provider_selection(
        _write_config(
            tmp_path / "agent.yaml",
            model="deepseek-v4-flash-0731",
            reasoning_effort="max",
        )
    )

    profile = selection.active_profile
    assert profile.provider == "alibaba-token-plan"
    assert profile.model == "deepseek-v4-flash-0731"
    assert profile.reasoning_effort == "max"
    runtime = profile.runtime_config()
    assert runtime.model == profile.model
    assert runtime.reasoning_effort == profile.reasoning_effort
    assert runtime.context_tokens == profile.context_tokens
    assert runtime.max_output_tokens == profile.max_output_tokens


@pytest.mark.parametrize(
    ("field_name", "replacement", "message"),
    (
        ("context_tokens", "context_tokens: 0", "positive context_tokens"),
        ("context_tokens", "context_tokens: true", "positive context_tokens"),
        (
            "max_output_tokens",
            "max_output_tokens: -1",
            "positive max_output_tokens",
        ),
        (
            "max_output_tokens",
            "max_output_tokens: 1000001",
            "cannot exceed context_tokens",
        ),
    ),
)
def test_agent_yaml_rejects_invalid_explicit_token_limits(
    tmp_path, field_name, replacement, message
):
    path = _write_config(tmp_path / "agent.yaml")
    original = (
        "context_tokens: 1000000"
        if field_name == "context_tokens"
        else "max_output_tokens: 262144"
    )
    path.write_text(
        path.read_text(encoding="utf-8").replace(
            original, replacement, 1
        ),
        encoding="utf-8",
    )

    with pytest.raises(ContractError, match=message):
        load_agent_provider_selection(path)


def test_agent_yaml_requires_model_and_both_token_limits(tmp_path):
    for omitted_line, message in (
        ("    model: qwen3.8-max\n", "explicit model"),
        ("    context_tokens: 1000000\n", "positive context_tokens"),
        ("    max_output_tokens: 262144\n", "positive max_output_tokens"),
    ):
        path = _write_config(tmp_path / (omitted_line.strip().split(":")[0] + ".yaml"))
        path.write_text(
            path.read_text(encoding="utf-8").replace(omitted_line, "", 1),
            encoding="utf-8",
        )
        with pytest.raises(ContractError, match=message):
            load_agent_provider_selection(path)


def test_provider_profile_digest_binds_token_limits():
    body = {
        "schema_version": "chemsmart.agent-provider-profile.v1",
        "profile_name": "deepseek",
        "provider": "deepseek",
        "wire_protocol": "openai-chat-completions",
        "api_key_env": "DEEPSEEK_API_KEY",
        "model": "deepseek-reasoner",
        "endpoint": "https://api.deepseek.com",
        "reasoning_effort": "max",
        "preserve_thinking": True,
        "context_tokens": 128_000,
        "max_output_tokens": 8_192,
    }
    digest = canonical_sha256(body)
    body["max_output_tokens"] = 8_193

    with pytest.raises(ContractError, match="digest mismatch"):
        AgentProviderProfileV1(**body, profile_sha256=digest)


def test_agent_yaml_rejects_nonempty_fallback(tmp_path):
    path = _write_config(tmp_path / "agent.yaml")
    path.write_text(
        path.read_text(encoding="utf-8").replace(
            "fallback: []", "fallback: [deepseek]"
        ),
        encoding="utf-8",
    )

    with pytest.raises(ContractError, match="does not execute fallback"):
        load_agent_provider_selection(path)


def test_agent_yaml_binds_explicit_provider_turn_deadlines(tmp_path):
    path = _write_config(
        tmp_path / "agent.yaml",
        model="deepseek-v4-flash-0731",
        reasoning_effort="max",
    )
    content = path.read_text(encoding="utf-8").replace(
        "    preserve_thinking: true\n",
        "    preserve_thinking: true\n"
        "    transport_deadlines:\n"
        "      connect_seconds: 15\n"
        "      first_event_seconds: 90\n"
        "      inter_event_seconds: 90\n"
        "      absolute_turn_seconds: 300\n",
        1,
    )
    path.write_text(content, encoding="utf-8")

    profile = load_agent_provider_selection(path).active_profile

    assert profile.schema_version == "chemsmart.agent-provider-profile.v2"
    assert profile.transport_deadlines.configuration_record() == {
        "connect_seconds": 15.0,
        "first_event_seconds": 90.0,
        "inter_event_seconds": 90.0,
        "absolute_turn_seconds": 300.0,
    }
    assert profile.runtime_config().turn_deadlines == profile.transport_deadlines


def test_agent_yaml_rejects_partial_provider_turn_deadlines(tmp_path):
    path = _write_config(tmp_path / "agent.yaml")
    content = path.read_text(encoding="utf-8").replace(
        "    preserve_thinking: true\n",
        "    preserve_thinking: true\n"
        "    transport_deadlines:\n"
        "      connect_seconds: 15\n",
        1,
    )
    path.write_text(content, encoding="utf-8")

    with pytest.raises(ContractError, match="requires exactly"):
        load_agent_provider_selection(path)


def test_selected_agent_profile_rejects_literal_secret(tmp_path):
    path = _write_config(tmp_path / "agent.yaml")
    content = path.read_text(encoding="utf-8").replace(
        "    api_key_env: ALIBABA_TOKEN_PLAN_KEY\n",
        "    api_key_env: ALIBABA_TOKEN_PLAN_KEY\n    api_key: forbidden\n",
    )
    path.write_text(content, encoding="utf-8")

    with pytest.raises(ContractError, match="literal API keys"):
        load_agent_provider_selection(path)


def test_unselected_agent_profile_is_outside_the_runtime_selection(tmp_path):
    path = _write_config(tmp_path / "agent.yaml")
    content = path.read_text(encoding="utf-8").replace(
        "    model: ignored\n",
        "    model: ignored\n    api_key: forbidden\n",
        1,
    )
    path.write_text(content, encoding="utf-8")

    selection = load_agent_provider_selection(path)

    assert selection.active_profile.profile_name == "alibaba-token-plan"


def test_unselected_agent_profile_header_is_not_loaded(tmp_path):
    path = _write_config(tmp_path / "agent.yaml")
    content = path.read_text(encoding="utf-8").replace(
        "    base_url: https://openrouter.ai/api/v1\n",
        "    base_url: https://openrouter.ai/api/v1\n"
        "    extra_headers:\n      Authorization: Bearer forbidden\n",
    )
    path.write_text(content, encoding="utf-8")

    selection = load_agent_provider_selection(path)

    assert selection.active_profile.profile_name == "alibaba-token-plan"


def test_nonempty_fallback_is_rejected_without_reading_that_profile(tmp_path):
    path = _write_config(tmp_path / "agent.yaml")
    content = path.read_text(encoding="utf-8").replace(
        "    model: deepseek-v4-flash\n",
        "    model: deepseek-v4-flash\n    api_key: forbidden\n",
    )
    content = content.replace("fallback: []", "fallback: [deepseek]")
    path.write_text(content, encoding="utf-8")

    with pytest.raises(ContractError, match="does not execute fallback"):
        load_agent_provider_selection(path)
