from __future__ import annotations

from dataclasses import replace

import pytest

from chemsmart.agent._contracts import ContractError, canonical_sha256
from chemsmart.agent.experiments.qwen_pyscf_dfc import (
    ComplexityGateInputV1,
    all_dfc_arms,
    bind_harness_experiment_config,
    build_episode_plans_from_preparations,
    build_episode_plans,
    build_qwen_experiment_preparation,
    build_qwen_dfc_arm,
    evaluate_complexity_gate,
)
from chemsmart.agent.experiments.qwen_pyscf_fixtures import (
    qwen_pyscf_cases_v1,
)
from chemsmart.agent.feedback import project_tool_feedback


def test_feedback_full_is_unchanged_and_causal_uses_action_signature():
    result = {
        "schema_version": "chemsmart.tool-result.v1",
        "tool": "preview_command",
        "status": "ok",
        "result": {
            "program": "pyscf",
            "jobtype": "sp",
            "method": "b3lyp",
            "basis": "def2-svp",
            "charge": 0,
            "multiplicity": 1,
            "findings": [],
            "receipt_sha256": "a" * 64,
            "script_text": "host-generated-view\n" * 1000,
            "inventory": [f"irrelevant-{index}" for index in range(100)],
        },
    }

    full = project_tool_feedback(
        tool="preview_command", result=result, mode="full-v1"
    )
    causal = project_tool_feedback(
        tool="preview_command", result=result, mode="causal-v1"
    )

    assert full.content == result
    assert causal.content != result
    assert "canonical_result" not in causal.content
    assert causal.content["result"]["method"] == "b3lyp"
    assert "script_text" not in causal.content["result"]
    assert "inventory" not in causal.content["result"]
    assert causal.action_signature is not None
    assert causal.receipt.causal_action_signature_sha256 == (
        causal.action_signature.signature_sha256
    )
    assert causal.receipt.verdict == "causal_action_signature_preserved"
    assert causal.receipt.provider_feedback_bytes < (
        causal.receipt.canonical_bytes
    )
    assert "$.result.script_text" in causal.receipt.omitted_paths
    assert causal.receipt.bytes_reduced > 0


def test_feedback_rejects_an_unregistered_projection():
    with pytest.raises(ContractError, match="unsupported"):
        project_tool_feedback(tool="x", result={}, mode="summary-v0")

    with pytest.raises(ContractError, match="registered tool projection"):
        project_tool_feedback(
            tool="unregistered_tool",
            result={
                "schema_version": "chemsmart.tool-result.v1",
                "tool": "unregistered_tool",
                "status": "ok",
                "result": {},
            },
            mode="causal-v1",
        )


def test_causal_feedback_surfaces_nested_scientific_failure():
    result = {
        "schema_version": "chemsmart.tool-result.v1",
        "tool": "preview_command",
        "status": "ok",
        "result": {
            "safe_preview": {
                "status": "failed",
                "rule_ids": ["preview.click_invocation_failed"],
            },
            "validator": {"status": "invalid"},
        },
    }

    causal = project_tool_feedback(
        tool="preview_command", result=result, mode="causal-v1"
    )

    assert causal.content["status"] == "failed"
    assert causal.content["transport_status"] == "ok"
    assert causal.action_signature is not None
    assert causal.action_signature.rule_ids == (
        "preview.click_invocation_failed",
    )


def test_all_dfc_arms_are_orthogonal_and_qwen_only():
    rows = all_dfc_arms(max_concurrency=4)

    assert len(rows) == 8
    assert len({item.arm_sha256 for item in rows}) == 8
    assert {item.max_concurrency for item in rows} == {4}
    assert {item.feedback_projection for item in rows} == {
        "full-v1",
        "causal-v1",
    }


def test_dfc_config_rejects_a_silent_model_or_role_change():
    arm = build_qwen_dfc_arm(
        decomposition=True,
        feedback_projection="full-v1",
        critic=True,
    )
    config = bind_harness_experiment_config(
        arm=arm,
        experiment_id="qp-dev-001.d1-ffull-c1.r0",
        prompt_sha256=canonical_sha256("prompt"),
        tool_schema_sha256=canonical_sha256("tools"),
        task_order_sha256=canonical_sha256("order"),
        token_budget=1_000_000,
        tool_call_budget=256,
        wall_time_seconds=5400,
    )

    assert config.provider_id == "alibaba-token-plan"
    assert config.model_id == "qwen3.8-max"
    assert config.reasoning_effort == "xhigh"
    with pytest.raises(ContractError, match="digest"):
        replace(config, model_id="qwen3.8-max-preview")
    with pytest.raises(ContractError, match="roles"):
        replace(arm, specialist_roles=())


def test_complexity_gate_only_dispatches_bounded_roles_when_eligible():
    enabled = build_qwen_dfc_arm(
        decomposition=True,
        feedback_projection="causal-v1",
        critic=False,
    )
    disabled = build_qwen_dfc_arm(
        decomposition=False,
        feedback_projection="causal-v1",
        critic=False,
    )
    signals = ComplexityGateInputV1(
        multi_stage=True, producer_artifact_edges=True
    )

    active = evaluate_complexity_gate(enabled, signals)
    inactive = evaluate_complexity_gate(disabled, signals)

    assert active.activated is True
    assert active.reason_ids == ("multi_stage", "producer_artifact_edges")
    assert active.requested_roles == (
        "scientific-specialist",
        "pyscf-specialist",
        "dag-specialist",
    )
    assert inactive.activated is False
    assert inactive.requested_roles == ()


def test_case_registry_covers_required_families_without_engine_authority():
    cases = qwen_pyscf_cases_v1()
    families = {item.family for item in cases}

    assert {item.split for item in cases} == {"development", "transfer"}
    assert len(cases) == len({item.case_id for item in cases})
    assert {
        "specified-dft-sp",
        "hf-branch",
        "missing-method",
        "missing-state",
        "complete-solvent",
        "incomplete-solvent",
        "workflow-edges",
        "closed-shell-tda",
        "closed-shell-tddft",
        "gpu-unavailable",
        "unsupported-ts-irc",
        "unsupported-methods",
        "future-artifact-boundary",
        "paraphrase-roundtrip",
    }.issubset(families)
    assert all("execute" not in item.expected_observation.lower() for item in cases)


def test_episode_plan_is_paired_counterbalanced_and_count_unbounded():
    cases = qwen_pyscf_cases_v1()[:2]
    arms = all_dfc_arms()
    plans = build_episode_plans(
        cases=cases,
        arms=arms,
        repeats=3,
        prompt_sha256=canonical_sha256("prompt-v1"),
        tool_schema_sha256=canonical_sha256("tools-v1"),
        task_order_sha256=canonical_sha256("order-v1"),
        token_budget=1_000_000,
        tool_call_budget=256,
        wall_time_seconds=5400,
    )

    assert len(plans) == 2 * 8 * 3
    assert tuple(item.order_index for item in plans) == tuple(range(len(plans)))
    assert all(item.engine_calls == 0 for item in plans)
    assert len({item.hypothesis.hypothesis_id for item in plans}) == len(plans)
    assert all(item.hypothesis.distinct_from_prior for item in plans)
    paired = {}
    for item in plans:
        paired.setdefault((item.case_sha256, item.repeat_index), set()).add(
            item.pairing_key
        )
    assert all(len(values) == 1 for values in paired.values())
    repeats_in_order = tuple(item.repeat_index for item in plans)
    assert repeats_in_order == (
        (0,) * (len(cases) * len(arms))
        + (1,) * (len(cases) * len(arms))
        + (2,) * (len(cases) * len(arms))
    )


def test_observed_preparation_is_the_plan_configuration_authority():
    case = qwen_pyscf_cases_v1()[0]
    arm = build_qwen_dfc_arm(
        decomposition=False,
        feedback_projection="causal-v1",
        critic=True,
    )
    preparation = build_qwen_experiment_preparation(
        case=case,
        arm=arm,
        repeat_index=2,
        task_spec_sha256=canonical_sha256("task"),
        artifact_sha256s=(canonical_sha256("geometry"),),
        provider_profile_sha256=canonical_sha256("profile"),
        prompt_sha256=canonical_sha256("observed-base-messages"),
        tool_schema_sha256=canonical_sha256("observed-tools"),
        task_order_sha256=canonical_sha256("observed-order"),
        token_budget=1_000_000,
        tool_call_budget=256,
        wall_time_seconds=5400,
    )

    plans = build_episode_plans_from_preparations(
        cases=(case,), arms=(arm,), preparations=(preparation,)
    )

    assert len(plans) == 1
    assert plans[0].experiment_config is preparation.experiment_config
    assert plans[0].hypothesis.prompt_sha256 == (
        preparation.experiment_config.prompt_sha256
    )
    assert preparation.preparation_sha256 in (
        plans[0].hypothesis.source_sha256s
    )
