"""Provider-independent evaluation primitives for agentic harness trials."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any, Iterable

from chemsmart.agent.harness.command_semantics import (
    CommandSemanticResult,
    evaluate_command_semantics,
)
from chemsmart.agent.harness.intent import IntentResult, IntentSpec, evaluate_intent


class OutcomeClass(str, Enum):
    DIRECT_PASS = "direct_pass"
    REPAIR_PASS = "repair_pass"
    VALID_ASK = "valid_ask"
    SPURIOUS_ASK = "spurious_ask"
    FORMAT_SCHEMA_FAILURE = "format_schema_failure"
    CLI_RUNTIME_FAILURE = "cli_runtime_failure"
    GENERATED_INPUT_FAILURE = "generated_input_invariant_failure"
    INTENT_DRIFT = "intent_drift"
    YAML_STATE_FAILURE = "yaml_state_failure"
    TERMINAL_ENVIRONMENT_FAILURE = "terminal_environment_failure"


_GENERATED_INPUT_RULE_IDS = frozenset(
    {
        "cmd.semantic.generated_input_missing",
        "cmd.semantic.generated_route_missing",
        "cmd.semantic.submit_generated_input_not_observed",
    }
)
_GENERATED_INPUT_RULE_PREFIXES = (
    "input.",
    "gaussian.ts.",
    "gaussian.freq.",
    "gaussian.irc.",
    "orca.ts.",
    "orca.freq.",
    "orca.irc.",
)
_FORMAT_RULE_IDS = frozenset(
    {
        "cmd.runtime.cli_value_error",
        "cmd.semantic.option_order",
        "cmd.semantic.orca_aux_basis_short_flag",
        "cmd.semantic.strict_parser",
    }
)
_TERMINAL_RULE_IDS = frozenset(
    {
        "cmd.runtime.dependency_missing",
        "cmd.runtime.input_not_found",
        "cmd.runtime.runner_unavailable",
        "cmd.runtime.server_invalid",
        "cmd.semantic.timeout",
    }
)


@dataclass(frozen=True)
class HarnessCase:
    case_id: str
    family: str
    turns: tuple[str, ...]
    expected_outcome: str
    intent: IntentSpec
    fixture: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        value = asdict(self)
        value["intent"] = self.intent.to_dict()
        return value


@dataclass(frozen=True)
class CaseEvaluation:
    case_id: str
    outcome: OutcomeClass
    passed: bool
    semantic: CommandSemanticResult | None = None
    intent: IntentResult | None = None
    failed_rule_ids: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        return {
            "case_id": self.case_id,
            "outcome": self.outcome.value,
            "passed": self.passed,
            "semantic": self.semantic.to_dict() if self.semantic else None,
            "intent": self.intent.to_dict() if self.intent else None,
            "failed_rule_ids": list(self.failed_rule_ids),
        }


def load_case_matrix(path: str | Path) -> list[HarnessCase]:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    cases: list[HarnessCase] = []
    for family in payload.get("families", []):
        base_intent = dict(family.get("intent") or {})
        for variant in family.get("variants", []):
            intent = _deep_merge(base_intent, variant.get("intent") or {})
            cases.append(
                HarnessCase(
                    case_id=f"{family['id']}::{variant['id']}",
                    family=str(family["id"]),
                    turns=tuple(str(item) for item in variant.get("turns") or ()),
                    expected_outcome=str(variant.get("expected_outcome") or "direct_pass"),
                    intent=IntentSpec.from_dict(intent),
                    fixture=dict(family.get("fixture") or {}),
                )
            )
    return cases


def evaluate_case_command(
    case: HarnessCase,
    command: str,
    *,
    cwd: str | Path | None = None,
    repaired: bool = False,
) -> CaseEvaluation:
    semantic = evaluate_command_semantics(command, cwd=cwd)
    intent = evaluate_intent(command, case.intent, cwd=str(cwd) if cwd else None)
    semantic_rules = tuple(semantic.failed_rule_ids)
    intent_rules = tuple(intent.failed_rule_ids)
    if semantic.verdict == "reject":
        outcome = classify_semantic_failure(semantic_rules)
        return CaseEvaluation(case.case_id, outcome, False, semantic, intent, semantic_rules)
    if intent.verdict == "reject":
        return CaseEvaluation(case.case_id, OutcomeClass.INTENT_DRIFT, False, semantic, intent, intent_rules)
    outcome = OutcomeClass.REPAIR_PASS if repaired else OutcomeClass.DIRECT_PASS
    return CaseEvaluation(case.case_id, outcome, True, semantic, intent)


def classify_agent_result(
    *,
    status: str,
    expected_outcome: str,
    semantic_rule_ids: Iterable[str] = (),
    intent_rule_ids: Iterable[str] = (),
    repaired: bool = False,
) -> OutcomeClass:
    semantic = tuple(semantic_rule_ids)
    intent = tuple(intent_rule_ids)
    if status == "needs_clarification":
        return OutcomeClass.VALID_ASK if expected_outcome == "valid_ask" else OutcomeClass.SPURIOUS_ASK
    if intent:
        return OutcomeClass.INTENT_DRIFT
    if semantic:
        return classify_semantic_failure(semantic)
    if status not in {"ready", "ok"}:
        return OutcomeClass.FORMAT_SCHEMA_FAILURE
    return OutcomeClass.REPAIR_PASS if repaired else OutcomeClass.DIRECT_PASS


def reliability_metrics(rows: Iterable[dict[str, Any]], *, k: int = 3) -> dict[str, Any]:
    by_case: dict[str, list[bool]] = {}
    for row in rows:
        by_case.setdefault(str(row["case_id"]), []).append(bool(row["passed"]))
    trials = [value[:k] for value in by_case.values() if value]
    all_attempts = [item for values in trials for item in values]
    return {
        "case_count": len(trials),
        "trial_count": len(all_attempts),
        "pass_at_1": (sum(values[0] for values in trials) / len(trials)) if trials else 0.0,
        f"pass_at_{k}": (sum(any(values) for values in trials) / len(trials)) if trials else 0.0,
        f"pass_power_{k}": (sum(len(values) == k and all(values) for values in trials) / len(trials)) if trials else 0.0,
    }


def classify_semantic_failure(rule_ids: Iterable[str]) -> OutcomeClass:
    """Classify stable harness rules without relying on provider error prose."""

    rules = tuple(str(rule_id) for rule_id in rule_ids)
    if any(
        rule in _GENERATED_INPUT_RULE_IDS
        or rule.startswith(_GENERATED_INPUT_RULE_PREFIXES)
        for rule in rules
    ):
        return OutcomeClass.GENERATED_INPUT_FAILURE
    if any("project" in rule or "yaml" in rule for rule in rules):
        return OutcomeClass.YAML_STATE_FAILURE
    if any(rule in _TERMINAL_RULE_IDS or rule.startswith("terminal.") for rule in rules):
        return OutcomeClass.TERMINAL_ENVIRONMENT_FAILURE
    if any(rule in _FORMAT_RULE_IDS or rule.startswith("cmd.contract.") for rule in rules):
        return OutcomeClass.FORMAT_SCHEMA_FAILURE
    return OutcomeClass.CLI_RUNTIME_FAILURE


def _deep_merge(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    merged = dict(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = _deep_merge(merged[key], value)
        else:
            merged[key] = value
    return merged


__all__ = [
    "CaseEvaluation",
    "HarnessCase",
    "OutcomeClass",
    "classify_agent_result",
    "classify_semantic_failure",
    "evaluate_case_command",
    "load_case_matrix",
    "reliability_metrics",
]
