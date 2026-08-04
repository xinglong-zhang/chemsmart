from __future__ import annotations

from chemsmart.agent._contracts import canonical_sha256
from chemsmart.agent.experiments.qwen_pyscf_analysis import (
    analyze_qwen_pyscf_campaign,
)
from chemsmart.agent.experiments.qwen_pyscf_campaign import (
    QwenPyscfEpisodeLedgerV1,
    QwenPyscfSplitLedgerV1,
)


def _episode(
    arm: str,
    *,
    verdict: str = "pass",
    safety=(),
    factor_status: str = "valid",
) -> QwenPyscfEpisodeLedgerV1:
    factor_findings = (
        ("experiment.factor.specialist_set_incomplete",)
        if factor_status == "invalid"
        else ()
    )
    body = {
        "schema_version": "chemsmart.qwen-episode-ledger.v1",
        "episode_id": f"QP-DEV-006.{arm}.r0",
        "split": "development",
        "case_sha256": canonical_sha256("case-QP-DEV-006"),
        "plan_sha256": canonical_sha256({"plan": arm}),
        "hypothesis_sha256": canonical_sha256({"hypothesis": arm}),
        "experiment_config_sha256": canonical_sha256({"config": arm}),
        "source_sha256": canonical_sha256("source"),
        "workspace_binding_sha256": canonical_sha256({"workspace": arm}),
        "result_sha256": canonical_sha256({"result": arm}),
        "grade_sha256": canonical_sha256({"grade": arm}),
        "session_terminal_state": "complete",
        "scientific_state": "previewed",
        "verdict": verdict,
        "failure_class": (
            "experiment_factor_invalid" if factor_status == "invalid" else ""
        ),
        "factor_realization_status": factor_status,
        "factor_realization_findings": factor_findings,
        "safety_violations": tuple(sorted(safety)),
    }
    return QwenPyscfEpisodeLedgerV1(
        **body, ledger_sha256=canonical_sha256(body)
    )


def _split(*episodes: QwenPyscfEpisodeLedgerV1) -> QwenPyscfSplitLedgerV1:
    body = {
        "schema_version": "chemsmart.qwen-split-ledger.v1",
        "campaign_window_sha256": canonical_sha256("window"),
        "freeze_manifest_sha256": canonical_sha256("freeze"),
        "split": "development",
        "started_utc": "2026-08-04T00:00:00Z",
        "observed_at_utc": "2026-08-04T01:00:00Z",
        "termination_reason": "split_complete",
        "episode_ledgers": episodes,
        "last_hypothesis_sha256": episodes[-1].hypothesis_sha256,
    }
    return QwenPyscfSplitLedgerV1(
        **body, ledger_sha256=canonical_sha256(body)
    )


def _v2_receipt(
    *, tool: str, mode: str, canonical_bytes: int, provider_bytes: int
):
    body = {
        "schema_version": "chemsmart.tool-feedback-projection-receipt.v2",
        "tool": tool,
        "mode": mode,
        "projection_contract": (
            "chemsmart.full-tool-result.v1"
            if mode == "full-v1"
            else "chemsmart.causal-action-projection.v1"
        ),
        "canonical_result_sha256": canonical_sha256(
            {"tool": tool, "mode": mode, "kind": "canonical"}
        ),
        "provider_feedback_sha256": canonical_sha256(
            {"tool": tool, "mode": mode, "kind": "provider"}
        ),
        "causal_feedback_sha256": canonical_sha256(
            {"tool": tool, "mode": mode, "kind": "causal"}
        ),
        "causal_action_signature_sha256": canonical_sha256(
            {"tool": tool, "mode": mode, "kind": "action"}
        ),
        "omitted_paths": (),
        "canonical_bytes": canonical_bytes,
        "provider_feedback_bytes": provider_bytes,
        "bytes_reduced": canonical_bytes - provider_bytes,
        "reduction_basis_points": round(
            (canonical_bytes - provider_bytes) * 10_000 / canonical_bytes
        ),
        "verdict": (
            "full_result_preserved"
            if mode == "full-v1"
            else "causal_action_signature_preserved"
        ),
    }
    return {**body, "receipt_sha256": canonical_sha256(body)}


def _v1_receipt(
    *, tool: str, mode: str, canonical_bytes: int, provider_bytes: int
):
    body = {
        "schema_version": "chemsmart.tool-feedback-equivalence.v1",
        "tool": tool,
        "mode": mode,
        "canonical_result_sha256": canonical_sha256(
            {"old": tool, "mode": mode, "kind": "canonical"}
        ),
        "provider_feedback_sha256": canonical_sha256(
            {"old": tool, "mode": mode, "kind": "provider"}
        ),
        "causal_feedback_sha256": canonical_sha256(
            {"old": tool, "mode": mode, "kind": "causal"}
        ),
        "semantic_signature_sha256": canonical_sha256(
            {"old": tool, "mode": mode, "kind": "semantic"}
        ),
        "omitted_paths": (),
        "canonical_bytes": canonical_bytes,
        "provider_feedback_bytes": provider_bytes,
        "verdict": "causal_semantics_preserved",
    }
    return {**body, "receipt_sha256": canonical_sha256(body)}


def _observation(
    episode: QwenPyscfEpisodeLedgerV1,
    *receipts,
    input_tokens: int = 1000,
    attempts: int = 2,
    activated: bool = False,
    critic: bool = False,
    specialists=(),
):
    critic_record = {
        "status": "complete" if critic else "not_enabled",
        "findings": (
            ({"severity": "critical", "rule_id": "seeded.critical"},)
            if critic
            else ()
        ),
    }
    body = {
        "schema_version": "chemsmart.live-harness-experiment-observations.v1",
        "experiment_config_sha256": episode.experiment_config_sha256,
        "gate": {"activated": activated},
        "specialists": tuple(specialists),
        "merge": {"status": "merged" if activated else "not_dispatched"},
        "critic": critic_record,
        "usage": {
            "provider_attempts": attempts,
            "input_tokens": input_tokens,
            "output_tokens": 100,
            "reasoning_tokens": 50,
            "wall_time_millis": 250,
            "successful_tool_calls": len(receipts),
            "failed_tool_calls": 0,
        },
    }
    experiment = {**body, "record_sha256": canonical_sha256(body)}
    return {
        "experiment_observations": experiment,
        "feedback_receipts": tuple(receipts),
    }


def _pair(result, factor: str):
    return next(item for item in result.pair_analyses if item.factor == factor)


def test_analyzer_reports_complete_v2_full_vs_causal_pair_and_usage():
    full = _episode("d0-ffull-c0")
    causal = _episode("d0-fcausal-c0")
    result = analyze_qwen_pyscf_campaign(
        split_ledgers=(_split(full, causal),),
        live_observations={
            full.episode_id: _observation(
                full,
                {
                    "feedback_equivalence_receipt": _v2_receipt(
                        tool="preview_command",
                        mode="full-v1",
                        canonical_bytes=1000,
                        provider_bytes=1000,
                    )
                },
                input_tokens=1000,
                attempts=3,
            ),
            causal.episode_id: _observation(
                causal,
                _v2_receipt(
                    tool="preview_command",
                    mode="causal-v1",
                    canonical_bytes=1000,
                    provider_bytes=300,
                ),
                input_tokens=700,
                attempts=2,
            ),
        },
    )

    pair = _pair(result, "F")
    assert pair.status == "complete"
    assert pair.effect_scope == "feedback-projection-effect"
    assert pair.scientific_comparison == "tied_pass"
    assert pair.feedback_contract_family == "causal-action-v2"
    assert dict(pair.metric_deltas)["input_tokens"] == -300
    assert dict(pair.metric_deltas)["model_attempts"] == -1
    assert dict(pair.provider_feedback_byte_deltas_by_tool) == {
        "preview_command": -700
    }
    causal_row = next(
        item
        for item in result.episode_analyses
        if item.feedback_projection == "causal-v1"
    )
    assert causal_row.reasoning_tokens == 50
    assert causal_row.feedback_bytes_by_tool[0].signature_kind == (
        "causal_action_signature_sha256"
    )


def test_analyzer_accepts_receipts_embedded_in_public_experiment_record():
    full = _episode("d0-ffull-c0")
    observation = _observation(
        full,
        _v2_receipt(
            tool="preview_command",
            mode="full-v1",
            canonical_bytes=100,
            provider_bytes=100,
        ),
    )
    nested = dict(observation["experiment_observations"])
    nested["feedback_receipts"] = observation["feedback_receipts"]
    nested.pop("record_sha256")
    nested["record_sha256"] = canonical_sha256(nested)

    result = analyze_qwen_pyscf_campaign(
        split_ledgers=(_split(full),),
        live_observations={
            full.episode_id: {"experiment_observations": nested}
        },
    )

    row = result.episode_analyses[0]
    assert row.observation_status == "valid"
    assert row.feedback_coverage_status == "complete"
    assert row.feedback_contract_families == ("causal-action-v2",)


def test_critic_pair_is_detection_only_not_a_scientific_improvement():
    without_critic = _episode("d0-ffull-c0")
    with_critic = _episode("d0-ffull-c1")
    result = analyze_qwen_pyscf_campaign(
        split_ledgers=(_split(without_critic, with_critic),),
        live_observations={
            without_critic.episode_id: _observation(
                without_critic,
                _v2_receipt(
                    tool="preview_command",
                    mode="full-v1",
                    canonical_bytes=100,
                    provider_bytes=100,
                ),
            ),
            with_critic.episode_id: _observation(
                with_critic,
                _v2_receipt(
                    tool="preview_command",
                    mode="full-v1",
                    canonical_bytes=100,
                    provider_bytes=100,
                ),
            ),
        },
    )

    pair = _pair(result, "C")
    assert pair.status == "complete"
    assert pair.effect_scope == "critic-detection-only"
    assert pair.scientific_comparison == "inconclusive"
    assert dict(pair.metric_deltas)["model_attempts"] == 0


def test_incomplete_pairs_are_explicitly_inconclusive():
    full = _episode("d0-ffull-c0")
    result = analyze_qwen_pyscf_campaign(
        split_ledgers=(_split(full),),
        live_observations={
            full.episode_id: _observation(
                full,
                _v2_receipt(
                    tool="preview_command",
                    mode="full-v1",
                    canonical_bytes=100,
                    provider_bytes=100,
                ),
            )
        },
    )

    pair = _pair(result, "F")
    assert pair.status == "inconclusive"
    assert pair.scientific_comparison == "inconclusive"
    assert "analysis.pair.candidate_missing" in pair.findings
    assert pair.metric_deltas == ()


def test_old_semantic_and_new_action_receipts_are_never_pooled():
    full = _episode("d0-ffull-c0")
    causal = _episode("d0-fcausal-c0")
    result = analyze_qwen_pyscf_campaign(
        split_ledgers=(_split(full, causal),),
        live_observations={
            full.episode_id: _observation(
                full,
                _v1_receipt(
                    tool="preview_command",
                    mode="full-v1",
                    canonical_bytes=900,
                    provider_bytes=900,
                ),
            ),
            causal.episode_id: _observation(
                causal,
                _v2_receipt(
                    tool="preview_command",
                    mode="causal-v1",
                    canonical_bytes=900,
                    provider_bytes=250,
                ),
            ),
        },
    )

    pair = _pair(result, "F")
    assert pair.status == "inconclusive"
    assert pair.feedback_contract_family == ""
    assert "analysis.pair.feedback_contract_mismatch" in pair.findings
    families = {
        item.feedback_contract_families[0]
        for item in result.episode_analyses
    }
    assert families == {"historical-semantic-v1", "causal-action-v2"}


def test_mixed_feedback_versions_in_one_episode_invalidate_observation():
    full = _episode("d0-ffull-c0")
    causal = _episode("d0-fcausal-c0")
    result = analyze_qwen_pyscf_campaign(
        split_ledgers=(_split(full, causal),),
        live_observations={
            full.episode_id: _observation(
                full,
                _v2_receipt(
                    tool="preview_command",
                    mode="full-v1",
                    canonical_bytes=100,
                    provider_bytes=100,
                ),
            ),
            causal.episode_id: _observation(
                causal,
                _v1_receipt(
                    tool="validate_project_yaml",
                    mode="causal-v1",
                    canonical_bytes=500,
                    provider_bytes=450,
                ),
                _v2_receipt(
                    tool="preview_command",
                    mode="causal-v1",
                    canonical_bytes=500,
                    provider_bytes=100,
                ),
            ),
        },
    )

    candidate = next(
        item
        for item in result.episode_analyses
        if item.feedback_projection == "causal-v1"
    )
    assert candidate.observation_status == "invalid"
    assert "analysis.feedback_contract_mixed_versions" in (
        candidate.observation_findings
    )
    assert len(candidate.feedback_bytes_by_tool) == 2
    assert _pair(result, "F").status == "inconclusive"


def test_strict_safety_factor_and_specialist_critic_outcomes_are_separate():
    red = _episode(
        "d0-ffull-c0",
        safety=("safety.false_ready",),
    )
    invalid = _episode(
        "d1-ffull-c1",
        verdict="pass",
        factor_status="invalid",
    )
    result = analyze_qwen_pyscf_campaign(
        split_ledgers=(_split(red, invalid),),
        live_observations={
            red.episode_id: _observation(
                red,
                _v2_receipt(
                    tool="preview_command",
                    mode="full-v1",
                    canonical_bytes=100,
                    provider_bytes=100,
                ),
            ),
            invalid.episode_id: _observation(
                invalid,
                _v2_receipt(
                    tool="preview_command",
                    mode="full-v1",
                    canonical_bytes=100,
                    provider_bytes=100,
                ),
                activated=True,
                critic=True,
                specialists=(
                    {"role": "dag-specialist", "status": "complete"},
                ),
            ),
        },
    )

    red_row = next(item for item in result.episode_analyses if not item.critic)
    invalid_row = next(item for item in result.episode_analyses if item.critic)
    assert red_row.strict_scientific_verdict == "fail"
    assert red_row.safety_status == "red"
    assert invalid_row.strict_scientific_verdict == "inconclusive"
    assert invalid_row.factor_realization_status == "invalid"
    assert invalid_row.specialist_outcomes == (("dag-specialist", "complete"),)
    assert invalid_row.critic_status == "complete"
    assert invalid_row.critic_findings == 1
    assert invalid_row.critic_critical_findings == 1
    assert result.safety_red_episodes == 1
    assert result.factor_invalid_episodes == 1
