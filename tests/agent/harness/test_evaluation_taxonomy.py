from __future__ import annotations

import pytest

from chemsmart.agent.harness.evaluation import (
    OutcomeClass,
    classify_semantic_failure,
    reliability_metrics,
)


@pytest.mark.parametrize(
    ("rule_id", "expected"),
    [
        ("input.gaussian.scan.directive", OutcomeClass.GENERATED_INPUT_FAILURE),
        ("gaussian.ts.route", OutcomeClass.GENERATED_INPUT_FAILURE),
        ("orca.freq.route", OutcomeClass.GENERATED_INPUT_FAILURE),
        ("cmd.semantic.generated_input_missing", OutcomeClass.GENERATED_INPUT_FAILURE),
        ("cmd.runtime.project_not_found", OutcomeClass.YAML_STATE_FAILURE),
        ("project_yaml.runtime_loader", OutcomeClass.YAML_STATE_FAILURE),
        ("cmd.runtime.server_invalid", OutcomeClass.TERMINAL_ENVIRONMENT_FAILURE),
        ("cmd.runtime.dependency_missing", OutcomeClass.TERMINAL_ENVIRONMENT_FAILURE),
        ("cmd.semantic.strict_parser", OutcomeClass.FORMAT_SCHEMA_FAILURE),
        ("cmd.contract.job_subcommand_required", OutcomeClass.FORMAT_SCHEMA_FAILURE),
        ("cmd.runtime.builder_error", OutcomeClass.CLI_RUNTIME_FAILURE),
    ],
)
def test_semantic_failure_taxonomy_uses_stable_rule_ids(rule_id, expected):
    assert classify_semantic_failure([rule_id]) is expected


def test_semantic_failure_taxonomy_prioritizes_generated_input_evidence():
    assert classify_semantic_failure(
        ["cmd.runtime.builder_error", "input.orca.neb.nimages"]
    ) is OutcomeClass.GENERATED_INPUT_FAILURE


def test_reliability_pass_at_one_uses_only_first_attempt_per_case():
    rows = [
        {"case_id": "a", "passed": True},
        {"case_id": "a", "passed": False},
        {"case_id": "a", "passed": False},
        {"case_id": "b", "passed": False},
        {"case_id": "b", "passed": True},
        {"case_id": "b", "passed": True},
    ]

    metrics = reliability_metrics(rows, k=3)

    assert metrics["pass_at_1"] == 0.5
    assert metrics["pass_at_3"] == 1.0
    assert metrics["pass_power_3"] == 0.0
