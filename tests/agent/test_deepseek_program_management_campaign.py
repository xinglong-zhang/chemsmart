from __future__ import annotations

from pathlib import Path

from chemsmart.agent.experiments.deepseek_program_management import (
    CampaignArm,
    build_campaign_run_config,
    build_episode_plans,
    load_campaign_definition,
    sanitize_public_record,
)
from chemsmart.agent.runtime.deepseek import (
    DEEPSEEK_V4_FLASH_CONTEXT_TOKENS,
    DEEPSEEK_V4_FLASH_MAX_OUTPUT_TOKENS,
)


def _definition_path() -> Path:
    return (
        Path(__file__).resolve().parents[2]
        / "docs"
        / "maintenance"
        / "deepseek-v4-flash-validation-v1.json"
    )


def test_campaign_builds_counterbalanced_unique_h0_h1_plans():
    definition = load_campaign_definition(_definition_path())
    config = build_campaign_run_config(workspace_identity="integration-tree")
    plans = build_episode_plans(definition, config)

    assert len(plans) == len(definition.cases) * 2
    assert len({plan.plan_sha256 for plan in plans}) == len(plans)
    assert len({plan.hypothesis.hypothesis_id for plan in plans}) == len(plans)
    assert [plan.arm for plan in plans[:4]] == [
        CampaignArm.H0,
        CampaignArm.H1,
        CampaignArm.H1,
        CampaignArm.H0,
    ]
    assert all(plan.network_budget.top_up_allowed is False for plan in plans)
    assert all(plan.network_budget.engine_calls == 0 for plan in plans)
    assert all(plan.network_budget.hpc_calls == 0 for plan in plans)


def test_campaign_uses_provider_maxima_without_attempt_count_cap():
    definition = load_campaign_definition(_definition_path())
    config = build_campaign_run_config(workspace_identity="integration-tree")
    plans = build_episode_plans(definition, config)

    assert config.max_input_tokens_per_request == DEEPSEEK_V4_FLASH_CONTEXT_TOKENS
    assert config.max_output_tokens_per_request == DEEPSEEK_V4_FLASH_MAX_OUTPUT_TOKENS
    assert definition.transport_wall_time_seconds == 14_400
    assert all(plan.network_budget.task_wall_time_seconds == 480 for plan in plans)


def test_h0_and_h1_differ_only_in_registered_tool_surface():
    definition = load_campaign_definition(_definition_path())
    config = build_campaign_run_config(workspace_identity="integration-tree")
    first, second = build_episode_plans(definition, config)[:2]

    assert first.case.case_sha256 == second.case.case_sha256
    assert first.pairing_sha256 == second.pairing_sha256
    assert first.messages == second.messages
    assert first.network_budget == second.network_budget
    assert first.hypothesis.tool_schema_sha256 != second.hypothesis.tool_schema_sha256


def test_public_record_removes_private_reasoning_and_credentials():
    raw = {
        "message": {
            "role": "assistant",
            "content": "Visible English answer.",
            "reasoning_content": "private trace",
            "tool_calls": [{"id": "call-1"}],
        },
        "authorization": "Bearer private-value",
        "usage": {"input_tokens": 3, "reasoning_tokens": 7},
    }

    public = sanitize_public_record(raw)

    assert public["message"]["content"] == "Visible English answer."
    assert public["message"]["tool_calls"] == [{"id": "call-1"}]
    assert "reasoning_content" not in public["message"]
    assert "authorization" not in public
    assert public["usage"] == {"input_tokens": 3, "reasoning_tokens": 7}
