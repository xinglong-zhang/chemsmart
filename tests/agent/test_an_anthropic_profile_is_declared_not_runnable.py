"""An anthropic profile is valid configuration that refuses to run.

The agent.yaml structure allows all four providers; selecting anthropic
fails closed with the adapter message rather than at parse time, so a
written profile survives loading and an inactive block sits untouched.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from chemsmart.agent._contracts import ContractError
from chemsmart.agent.provider_config import load_agent_provider_selection

_YAML = """\
active: anthropic
fallback: []
providers:
  anthropic:
    type: openai
    api_key_env: ANTHROPIC_API_KEY
    model: claude-fable-5
    context_tokens: 1000000
    max_output_tokens: 128000
    base_url: https://api.anthropic.com
    reasoning_effort: high
"""


def _config(tmp_path: Path) -> Path:
    path = tmp_path / "agent.yaml"
    path.write_text(_YAML, encoding="utf-8")
    return path


def test_the_profile_loads_and_derives_its_provider(tmp_path):
    profile = load_agent_provider_selection(_config(tmp_path)).active_profile

    assert profile.provider == "anthropic"
    assert profile.model == "claude-fable-5"
    assert profile.api_key_env == "ANTHROPIC_API_KEY"


def test_selecting_it_for_execution_refuses_with_the_adapter_message(
    tmp_path,
):
    profile = load_agent_provider_selection(_config(tmp_path)).active_profile

    with pytest.raises(
        ContractError, match="anthropic adapter is not registered"
    ):
        profile.runtime_config()


def test_an_inactive_anthropic_block_never_blocks_the_active_profile(
    tmp_path,
):
    path = tmp_path / "agent.yaml"
    path.write_text(
        _YAML.replace("active: anthropic", "active: openai") + """\
  openai:
    type: openai
    api_key_env: OPENAI_API_KEY
    model: gpt-5.2
    context_tokens: 400000
    max_output_tokens: 64000
    base_url: https://api.openai.com/v1
""",
        encoding="utf-8",
    )

    profile = load_agent_provider_selection(path).active_profile

    assert profile.provider == "openai"
    assert profile.runtime_config().model == "gpt-5.2"
    assert profile.runtime_config().reasoning_effort == ""
