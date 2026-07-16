"""Bounded rule-driven repair policy for observable runtime failures."""

from __future__ import annotations

from enum import Enum
from typing import Iterable

from pydantic import BaseModel, ConfigDict


class RepairAction(str, Enum):
    DETERMINISTIC_REPAIR = "deterministic_repair"
    ASK_USER = "ask_user"
    TERMINATE = "terminate"
    REVIEW = "review"


class RepairDecision(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    action: RepairAction
    rule_ids: tuple[str, ...]
    missing_slots: tuple[str, ...] = ()
    reason: str


_ASK_RULE_SLOTS = {
    "cmd.runtime.input_not_found": ("input_path",),
    "cmd.runtime.project_not_found": ("project",),
    "cmd.runtime.server_invalid": ("server",),
    "workflow.project.not_found": ("project",),
}
_TERMINATE_RULES = {
    "cmd.runtime.dependency_missing",
    "cmd.runtime.runner_unavailable",
    "cmd.semantic.timeout",
}
_DETERMINISTIC_RULES = {
    "cmd.runtime.cli_value_error",
    "cmd.semantic.option_order",
    "cmd.semantic.orca_aux_basis_short_flag",
    "cmd.semantic.strict_parser",
    "yaml.update.stringified_json_invalid",
}


def decide_repair(
    rule_ids: Iterable[str],
    *,
    repeated_count: int = 1,
) -> RepairDecision:
    rules = tuple(dict.fromkeys(str(item) for item in rule_ids if str(item)))
    if repeated_count >= 2:
        return RepairDecision(
            action=RepairAction.TERMINATE,
            rule_ids=rules,
            reason="the same failure repeated without a state change",
        )
    missing_slots = tuple(
        dict.fromkeys(
            slot for rule in rules for slot in _ASK_RULE_SLOTS.get(rule, ())
        )
    )
    if missing_slots:
        return RepairDecision(
            action=RepairAction.ASK_USER,
            rule_ids=rules,
            missing_slots=missing_slots,
            reason="a required user-owned value is missing or unresolved",
        )
    if any(rule in _TERMINATE_RULES for rule in rules):
        return RepairDecision(
            action=RepairAction.TERMINATE,
            rule_ids=rules,
            reason="the environment cannot be repaired by changing chemistry intent",
        )
    if rules and all(rule in _DETERMINISTIC_RULES for rule in rules):
        return RepairDecision(
            action=RepairAction.DETERMINISTIC_REPAIR,
            rule_ids=rules,
            reason="all failures have one schema-preserving repair",
        )
    return RepairDecision(
        action=RepairAction.REVIEW,
        rule_ids=rules,
        reason="repair may change scientific intent and requires explicit review",
    )


__all__ = ["RepairAction", "RepairDecision", "decide_repair"]
