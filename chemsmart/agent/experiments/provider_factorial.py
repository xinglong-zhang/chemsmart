"""Provider-bound D x F x C blocks over the existing PySCF experiment IR.

The historical Qwen/PySCF schemas remain the replay format for arms, cases,
and provider-free preparations.  This additive block binds those immutable
records to one exact provider profile and selection without renaming or
migrating previously recorded Qwen evidence.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable

from chemsmart.agent._contracts import (
    ContractError,
    canonical_sha256,
    require_identifier,
    require_sha256,
)
from chemsmart.agent.experiments.qwen_pyscf_dfc import (
    QwenDfcArmV1,
    QwenExperimentPreparationV1,
    QwenPyscfCaseSpecV1,
)
from chemsmart.agent.provider_config import (
    ALIBABA_TOKEN_PLAN_MODEL,
    ALIBABA_TOKEN_PLAN_PROVIDER,
    AgentProviderSelectionV1,
)


PRIMARY_DFC_ARM_ORDER_V1 = (
    "d0-ffull-c0",
    "d1-fcausal-c1",
    "d0-fcausal-c1",
    "d1-ffull-c0",
    "d0-fcausal-c0",
    "d1-ffull-c1",
    "d0-ffull-c1",
    "d1-fcausal-c0",
)


@dataclass(frozen=True)
class ProviderFactorialBlockV1:
    """One complete, ordered 2^3 block for one provider and one case."""

    schema_version: str
    block_id: str
    provider_profile_name: str
    provider_id: str
    model_id: str
    reasoning_effort: str
    provider_profile_sha256: str
    provider_selection_sha256: str
    case_id: str
    case_sha256: str
    repeat_index: int
    deterministic_oracle_id: str
    source_sha256s: tuple[str, ...]
    arm_bindings: tuple[tuple[str, str, str], ...]
    prompt_sha256: str
    tool_schema_sha256: str
    task_order_sha256: str
    token_budget: int
    tool_call_budget: int
    wall_time_seconds: int
    max_concurrency: int
    block_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.provider-factorial-block.v1":
            raise ContractError("unsupported provider factorial block")
        if require_identifier(self.block_id, "block_id") != self.block_id:
            raise ContractError("factorial block ID is not canonical")
        for name, value in (
            ("provider_profile_name", self.provider_profile_name),
            ("provider_id", self.provider_id),
            ("model_id", self.model_id),
            ("reasoning_effort", self.reasoning_effort),
            ("case_id", self.case_id),
            ("deterministic_oracle_id", self.deterministic_oracle_id),
        ):
            if not str(value).strip():
                raise ContractError(f"{name} is required")
        if self.repeat_index < 0:
            raise ContractError("factorial repeat index must be non-negative")
        for name in (
            "provider_profile_sha256",
            "provider_selection_sha256",
            "case_sha256",
            "prompt_sha256",
            "tool_schema_sha256",
            "task_order_sha256",
        ):
            require_sha256(getattr(self, name), name)
        if not self.source_sha256s or self.source_sha256s != tuple(
            sorted(set(self.source_sha256s))
        ):
            raise ContractError("factorial sources must be non-empty and canonical")
        for digest in self.source_sha256s:
            require_sha256(digest, "source_sha256")
        if tuple(item[0] for item in self.arm_bindings) != (
            PRIMARY_DFC_ARM_ORDER_V1
        ):
            raise ContractError("factorial block does not use the frozen arm order")
        if len(self.arm_bindings) != len(set(self.arm_bindings)):
            raise ContractError("factorial block repeats an arm binding")
        for arm_id, arm_sha256, config_sha256 in self.arm_bindings:
            if not arm_id:
                raise ContractError("factorial arm identity is required")
            require_sha256(arm_sha256, "arm_sha256")
            require_sha256(config_sha256, "config_sha256")
        if min(
            self.token_budget,
            self.tool_call_budget,
            self.wall_time_seconds,
        ) < 1:
            raise ContractError("factorial budgets must be positive")
        if not 1 <= self.max_concurrency <= 4:
            raise ContractError("factorial concurrency must be within [1, 4]")
        body = _without_field(self, "block_sha256")
        if self.block_sha256 != canonical_sha256(body):
            raise ContractError("provider factorial block digest mismatch")

    def config_sha256_for_arm(self, arm_id: str) -> str:
        """Return the frozen configuration digest for one canonical arm."""

        matches = tuple(
            config_sha256
            for observed_arm_id, _, config_sha256 in self.arm_bindings
            if observed_arm_id == arm_id
        )
        if len(matches) != 1:
            raise ContractError("arm is outside the provider factorial block")
        return matches[0]

    def public_record(self) -> dict[str, Any]:
        return {
            **_without_field(self, "block_sha256"),
            "block_sha256": self.block_sha256,
        }


def build_provider_factorial_block(
    *,
    block_id: str,
    selection: AgentProviderSelectionV1,
    case: QwenPyscfCaseSpecV1,
    arms: Iterable[QwenDfcArmV1],
    preparations: Iterable[QwenExperimentPreparationV1],
    source_sha256s: Iterable[str] = (),
) -> ProviderFactorialBlockV1:
    """Freeze eight observed preparations under one exact provider profile."""

    if not isinstance(selection, AgentProviderSelectionV1):
        raise ContractError("factorial provider selection is not typed")
    if selection.fallback_profiles:
        raise ContractError("factorial blocks forbid provider fallbacks")
    profile = selection.active_profile
    runtime_config = profile.runtime_config()
    if (
        runtime_config.provider != profile.provider
        or runtime_config.model != profile.model
        or runtime_config.reasoning_effort != profile.reasoning_effort
    ):
        raise ContractError("provider profile differs from its runtime adapter")
    if not isinstance(case, QwenPyscfCaseSpecV1):
        raise ContractError("factorial case is not typed")

    arm_rows = tuple(arms)
    if tuple(item.arm_id for item in arm_rows) != PRIMARY_DFC_ARM_ORDER_V1:
        raise ContractError("factorial arms must use the frozen primary order")
    if len({item.arm_sha256 for item in arm_rows}) != 8:
        raise ContractError("factorial block requires eight distinct arms")

    preparation_rows = tuple(preparations)
    if len(preparation_rows) != 8:
        raise ContractError("factorial block requires eight preparations")
    by_arm_id = {item.arm_id: item for item in preparation_rows}
    if len(by_arm_id) != 8 or tuple(by_arm_id) != PRIMARY_DFC_ARM_ORDER_V1:
        # Mapping insertion order is part of the recorded dispatch order.
        raise ContractError("factorial preparations differ from the frozen order")

    repeats = {item.repeat_index for item in preparation_rows}
    if len(repeats) != 1:
        raise ContractError("factorial preparations mix repeat indexes")
    prompt_sha256s: set[str] = set()
    tool_schema_sha256s: set[str] = set()
    task_order_sha256s: set[str] = set()
    token_budgets: set[int] = set()
    tool_call_budgets: set[int] = set()
    wall_time_budgets: set[int] = set()
    concurrency_values: set[int] = set()
    bindings: list[tuple[str, str, str]] = []
    observed_sources = set(source_sha256s)
    observed_sources.update(case.source_sha256s)

    for arm in arm_rows:
        preparation = by_arm_id[arm.arm_id]
        config = preparation.experiment_config
        legacy_episode_id = (
            f"{case.case_id}.{arm.arm_id}.r{preparation.repeat_index}"
        )
        expected_episode_id = (
            legacy_episode_id
            if profile.provider == ALIBABA_TOKEN_PLAN_PROVIDER
            and profile.model == ALIBABA_TOKEN_PLAN_MODEL
            else f"{profile.provider}.{legacy_episode_id}"
        )
        if preparation.episode_id != expected_episode_id:
            raise ContractError(
                "factorial preparation lacks its provider episode namespace"
            )
        if (
            preparation.case_sha256 != case.case_sha256
            or preparation.arm_sha256 != arm.arm_sha256
            or preparation.provider_profile_sha256 != profile.profile_sha256
        ):
            raise ContractError("factorial preparation binding differs")
        if (
            config.decomposition != arm.decomposition
            or config.feedback_projection != arm.feedback_projection
            or config.critic != arm.critic
            or config.max_concurrency != arm.max_concurrency
        ):
            raise ContractError("factorial configuration differs from its arm")
        if (
            config.provider_id != profile.provider
            or config.model_id != profile.model
            or config.reasoning_effort != profile.reasoning_effort
        ):
            raise ContractError("factorial configuration uses another provider")
        prompt_sha256s.add(config.prompt_sha256)
        tool_schema_sha256s.add(config.tool_schema_sha256)
        task_order_sha256s.add(config.task_order_sha256)
        token_budgets.add(config.token_budget)
        tool_call_budgets.add(config.tool_call_budget)
        wall_time_budgets.add(config.wall_time_seconds)
        concurrency_values.add(config.max_concurrency)
        observed_sources.update(preparation.artifact_sha256s)
        bindings.append(
            (arm.arm_id, arm.arm_sha256, config.config_sha256)
        )

    common_values = (
        prompt_sha256s,
        tool_schema_sha256s,
        task_order_sha256s,
        token_budgets,
        tool_call_budgets,
        wall_time_budgets,
        concurrency_values,
    )
    if any(len(values) != 1 for values in common_values):
        raise ContractError("factorial preparations do not share frozen inputs")
    body = {
        "schema_version": "chemsmart.provider-factorial-block.v1",
        "block_id": str(block_id).strip().lower(),
        "provider_profile_name": profile.profile_name,
        "provider_id": profile.provider,
        "model_id": profile.model,
        "reasoning_effort": profile.reasoning_effort,
        "provider_profile_sha256": profile.profile_sha256,
        "provider_selection_sha256": selection.selection_sha256,
        "case_id": case.case_id,
        "case_sha256": case.case_sha256,
        "repeat_index": next(iter(repeats)),
        "deterministic_oracle_id": case.deterministic_oracle_id,
        "source_sha256s": tuple(sorted(observed_sources)),
        "arm_bindings": tuple(bindings),
        "prompt_sha256": next(iter(prompt_sha256s)),
        "tool_schema_sha256": next(iter(tool_schema_sha256s)),
        "task_order_sha256": next(iter(task_order_sha256s)),
        "token_budget": next(iter(token_budgets)),
        "tool_call_budget": next(iter(tool_call_budgets)),
        "wall_time_seconds": next(iter(wall_time_budgets)),
        "max_concurrency": next(iter(concurrency_values)),
    }
    return ProviderFactorialBlockV1(
        **body, block_sha256=canonical_sha256(body)
    )


def validate_provider_factorial_binding(
    *,
    block: ProviderFactorialBlockV1,
    selection: AgentProviderSelectionV1,
    case: QwenPyscfCaseSpecV1,
    arm: QwenDfcArmV1,
    config_sha256: str,
) -> None:
    """Reject cross-provider, cross-case, or cross-arm episode reuse."""

    if selection.fallback_profiles:
        raise ContractError("factorial blocks forbid provider fallbacks")
    profile = selection.active_profile
    if (
        block.provider_profile_sha256 != profile.profile_sha256
        or block.provider_selection_sha256 != selection.selection_sha256
        or block.provider_id != profile.provider
        or block.model_id != profile.model
        or block.reasoning_effort != profile.reasoning_effort
    ):
        raise ContractError("factorial block uses another provider selection")
    if block.case_sha256 != case.case_sha256 or block.case_id != case.case_id:
        raise ContractError("factorial block uses another case")
    expected = block.config_sha256_for_arm(arm.arm_id)
    if expected != require_sha256(config_sha256, "config_sha256"):
        raise ContractError("factorial block uses another arm configuration")


def _without_field(value: Any, field_name: str) -> dict[str, Any]:
    return {
        key: item
        for key, item in value.__dict__.items()
        if key != field_name
    }


__all__ = [
    "PRIMARY_DFC_ARM_ORDER_V1",
    "ProviderFactorialBlockV1",
    "build_provider_factorial_block",
    "validate_provider_factorial_binding",
]
