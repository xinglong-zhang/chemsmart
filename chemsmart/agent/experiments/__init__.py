"""Reproducible, provider-backed agent experiments."""

from chemsmart.agent.experiments.deepseek_program_management import (
    CAMPAIGN_SYSTEM_PROMPT_V1,
    CampaignArm,
    CampaignDefinitionV1,
    CampaignEpisodePlanV1,
    CampaignEpisodeReceiptV1,
    CampaignHostInputsFactory,
    CampaignRunConfigV1,
    CampaignRunReceiptV1,
    DeepSeekProgramManagementCampaignRunner,
    DeterministicOracleReceiptV1,
    build_campaign_run_config,
    build_episode_plans,
    evaluate_episode_observations,
    load_campaign_definition,
    sanitize_public_record,
)
from chemsmart.agent.experiments.program_management_context import (
    CampaignArtifactSlotV1,
    CampaignHostFixtureV1,
    CampaignPublicArtifactV1,
    CampaignPublicContextV1,
    CampaignPublicReceiptRefV1,
    CampaignToolInputV1,
)
from chemsmart.agent.experiments.program_management_fixtures import (
    ProgramManagementHostFixtureFactoryV1,
)

__all__ = [
    "CAMPAIGN_SYSTEM_PROMPT_V1",
    "CampaignArm",
    "CampaignArtifactSlotV1",
    "CampaignDefinitionV1",
    "CampaignEpisodePlanV1",
    "CampaignEpisodeReceiptV1",
    "CampaignHostInputsFactory",
    "CampaignHostFixtureV1",
    "CampaignPublicArtifactV1",
    "CampaignPublicContextV1",
    "CampaignPublicReceiptRefV1",
    "CampaignRunConfigV1",
    "CampaignRunReceiptV1",
    "CampaignToolInputV1",
    "DeepSeekProgramManagementCampaignRunner",
    "DeterministicOracleReceiptV1",
    "ProgramManagementHostFixtureFactoryV1",
    "build_campaign_run_config",
    "build_episode_plans",
    "evaluate_episode_observations",
    "load_campaign_definition",
    "sanitize_public_record",
]
