from __future__ import annotations

from pathlib import Path

import pytest

from chemsmart.agent._contracts import ContractError
from chemsmart.agent.provider_config import load_agent_provider_selection


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
fallback:
  - deepseek
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
    base_url: https://token-plan.ap-southeast-1.maas.aliyuncs.com/compatible-mode/v1
    reasoning_effort: {reasoning_effort}
    preserve_thinking: true
  deepseek:
    type: openai
    api_key_env: DEEPSEEK-api-key
    model: deepseek-v4-flash
    base_url: https://api.deepseek.com
    reasoning_effort: max
    preserve_thinking: true
""".lstrip(),
        encoding="utf-8",
    )
    return path


def test_agent_yaml_selects_qwen_production_and_explicit_fallback(tmp_path):
    selection = load_agent_provider_selection(_write_config(tmp_path / "agent.yaml"))

    assert selection.active_profile.profile_name == "alibaba-token-plan"
    assert selection.active_profile.provider == "alibaba-token-plan"
    assert selection.active_profile.model == "qwen3.8-max"
    assert selection.active_profile.reasoning_effort == "xhigh"
    assert selection.active_profile.preserve_thinking is True
    assert tuple(item.profile_name for item in selection.fallback_profiles) == (
        "deepseek",
    )
    assert selection.fallback_profiles[0].provider == "deepseek"


def test_agent_yaml_profile_override_keeps_fallback_attributed(tmp_path):
    selection = load_agent_provider_selection(
        _write_config(tmp_path / "agent.yaml"), requested_profile="deepseek"
    )

    assert selection.active_profile.provider == "deepseek"
    assert selection.fallback_profiles == ()


def test_agent_yaml_rejects_preview_alias(tmp_path):
    with pytest.raises(ContractError, match="production qwen3.8-max"):
        load_agent_provider_selection(
            _write_config(tmp_path / "agent.yaml", model="qwen3.8-max-preview")
        )


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


def test_selected_fallback_rejects_a_literal_secret(tmp_path):
    path = _write_config(tmp_path / "agent.yaml")
    content = path.read_text(encoding="utf-8").replace(
        "    model: deepseek-v4-flash\n",
        "    model: deepseek-v4-flash\n    api_key: forbidden\n",
    )
    path.write_text(content, encoding="utf-8")

    with pytest.raises(ContractError, match="literal API keys"):
        load_agent_provider_selection(path)
