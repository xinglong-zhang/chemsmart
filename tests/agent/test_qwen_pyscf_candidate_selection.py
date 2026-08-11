from __future__ import annotations

from datetime import datetime, timezone

import pytest

from chemsmart.agent._contracts import (
    ContractError,
    canonical_sha256,
    file_sha256,
)
from chemsmart.agent.experiments.qwen_pyscf_analysis import (
    QwenPyscfCampaignAnalysisV1,
    QwenPyscfEpisodeAnalysisV1,
    select_qwen_pyscf_development_candidate,
)
from chemsmart.agent.experiments.qwen_pyscf_campaign import (
    QwenPyscfEpisodeLedgerV1,
    QwenPyscfSplitLedgerV1,
    approved_xyz_source,
    bind_case_to_approved_xyz,
    build_frozen_transfer_manifest,
    build_qwen_campaign_window,
)
from chemsmart.agent.experiments.qwen_pyscf_dfc import (
    all_dfc_arms,
    build_episode_plans,
)
from chemsmart.agent.experiments.qwen_pyscf_fixtures import (
    qwen_pyscf_cases_v1,
)
from chemsmart.agent.experiments.qwen_pyscf_transfer import (
    build_qwen_pyscf_transfer_evaluation_receipt,
    build_qwen_pyscf_transfer_opening,
)


_ARMS = tuple(
    sorted(
        f"d{decomposition}-f{feedback}-c{critic}"
        for decomposition in (0, 1)
        for feedback in ("full", "causal")
        for critic in (0, 1)
    )
)
_CASES = ("QP-DEV-005", "QP-DEV-006", "QP-DEV-007")


def _factors(arm_id: str) -> tuple[bool, str, bool]:
    decomposition, feedback, critic = arm_id.split("-")
    return (
        decomposition == "d1",
        feedback.removeprefix("f") + "-v1",
        critic == "c1",
    )


def _episode(
    *,
    arm_id: str,
    case_id: str,
    verdict: str = "pass",
    safety_status: str = "green",
    factor_status: str = "valid",
    observation_status: str = "valid",
    attempts: int = 2,
    tokens: int = 300,
    latency_millis: int = 500,
    tool_failures: int = 0,
    case_sha256: str | None = None,
) -> QwenPyscfEpisodeAnalysisV1:
    decomposition, feedback, critic = _factors(arm_id)
    input_tokens = tokens // 2
    output_tokens = tokens // 3
    reasoning_tokens = tokens - input_tokens - output_tokens
    body = {
        "schema_version": "chemsmart.qwen-episode-analysis.v1",
        "episode_id": f"{case_id}.{arm_id}.r0",
        "split": "development",
        "case_id": case_id,
        "case_sha256": case_sha256 or canonical_sha256({"case": case_id}),
        "repeat_index": 0,
        "arm_id": arm_id,
        "decomposition": decomposition,
        "feedback_projection": feedback,
        "critic": critic,
        "ledger_verdict": verdict,
        "strict_scientific_verdict": verdict,
        "safety_status": safety_status,
        "safety_violations": (
            ("safety.seeded",) if safety_status == "red" else ()
        ),
        "factor_realization_status": factor_status,
        "factor_realization_findings": (
            ("factor.seeded",) if factor_status != "valid" else ()
        ),
        "failure_class": "",
        "observation_status": observation_status,
        "observation_findings": (
            ("observation.seeded",)
            if observation_status != "valid"
            else ()
        ),
        "decomposition_activated": decomposition,
        "model_attempts": attempts,
        "input_tokens": input_tokens,
        "output_tokens": output_tokens,
        "reasoning_tokens": reasoning_tokens,
        "provider_latency_millis": latency_millis,
        "successful_tool_calls": 4,
        "failed_tool_calls": tool_failures,
        "specialist_outcomes": (),
        "specialist_merge_status": "merged" if decomposition else "not_dispatched",
        "critic_status": "complete" if critic else "not_enabled",
        "critic_findings": 0,
        "critic_critical_findings": 0,
        "feedback_contract_families": (),
        "feedback_coverage_status": "not_exercised",
        "feedback_bytes_by_tool": (),
    }
    return QwenPyscfEpisodeAnalysisV1(
        **body, analysis_sha256=canonical_sha256(body)
    )


def _analysis(
    *,
    verdict_by_arm: dict[str, str] | None = None,
    red_arms: set[str] | None = None,
    tokens_by_arm: dict[str, int] | None = None,
    latency_by_arm: dict[str, int] | None = None,
    omitted: set[tuple[str, str]] | None = None,
    case_sha_by_id: dict[str, str] | None = None,
) -> QwenPyscfCampaignAnalysisV1:
    verdict_by_arm = verdict_by_arm or {}
    red_arms = red_arms or set()
    tokens_by_arm = tokens_by_arm or {}
    latency_by_arm = latency_by_arm or {}
    omitted = omitted or set()
    case_sha_by_id = case_sha_by_id or {}
    episodes = tuple(
        sorted(
            (
                _episode(
                    arm_id=arm_id,
                    case_id=case_id,
                    verdict=verdict_by_arm.get(arm_id, "pass"),
                    safety_status="red" if arm_id in red_arms else "green",
                    tokens=tokens_by_arm.get(arm_id, 300),
                    latency_millis=latency_by_arm.get(arm_id, 500),
                    case_sha256=case_sha_by_id.get(case_id),
                )
                for arm_id in _ARMS
                for case_id in _CASES
                if (arm_id, case_id) not in omitted
            ),
            key=lambda item: item.episode_id,
        )
    )
    body = {
        "schema_version": "chemsmart.qwen-campaign-analysis.v1",
        "split_ledger_sha256s": (canonical_sha256("development-ledger"),),
        "episode_analyses": episodes,
        "pair_analyses": (),
        "strict_passes": sum(
            item.strict_scientific_verdict == "pass" for item in episodes
        ),
        "strict_failures": sum(
            item.strict_scientific_verdict == "fail" for item in episodes
        ),
        "strict_inconclusive": sum(
            item.strict_scientific_verdict == "inconclusive"
            for item in episodes
        ),
        "safety_red_episodes": sum(
            item.safety_status == "red" for item in episodes
        ),
        "factor_invalid_episodes": 0,
        "factor_not_observed_episodes": 0,
    }
    return QwenPyscfCampaignAnalysisV1(
        **body, analysis_sha256=canonical_sha256(body)
    )


def _select(analysis: QwenPyscfCampaignAnalysisV1):
    return select_qwen_pyscf_development_candidate(
        analysis=analysis,
        preregistered_case_count_per_arm=3,
    )


def test_safety_red_arm_is_vetoed_even_with_more_strict_passes():
    result = _select(
        _analysis(
            verdict_by_arm={
                "d0-fcausal-c0": "fail",
                "d1-ffull-c0": "pass",
            },
            red_arms={"d1-ffull-c0"},
        )
    )

    summaries = {item.arm_id: item for item in result.arm_summaries}
    assert summaries["d1-ffull-c0"].strict_passes == 3
    assert summaries["d1-ffull-c0"].eligibility_status == "ineligible"
    assert "selection.arm.safety_red" in (
        summaries["d1-ffull-c0"].eligibility_findings
    )
    assert result.selected_arm_id != "d1-ffull-c0"
    assert result.transfer_status == "admissible"


def test_scientific_tie_prefers_no_critic_then_minimal_profile():
    result = _select(
        _analysis(
            tokens_by_arm={
                # C is detection-only, so a cheap critic cannot win a tied
                # scientific outcome over a no-critic profile.
                "d0-fcausal-c1": 3,
            },
            latency_by_arm={"d0-fcausal-c1": 3},
        )
    )

    assert result.selected_arm_id == "d0-fcausal-c0"
    selected = next(
        item for item in result.arm_summaries
        if item.arm_id == result.selected_arm_id
    )
    assert selected.profile_complexity == 1
    assert selected.strict_passes == 3
    assert selected.total_provider_attempts == 6
    assert selected.total_tokens == 900
    assert selected.total_provider_latency_millis == 1500
    assert selected.total_tool_failures == 0
    assert result.selection_sha256 == canonical_sha256(
        {
            key: value
            for key, value in result.__dict__.items()
            if key != "selection_sha256"
        }
    )


def test_no_safety_green_non_h0_candidate_blocks_transfer():
    result = _select(_analysis(red_arms=set(_ARMS) - {"d0-ffull-c0"}))

    assert result.selection_status == "blocked"
    assert result.selected_arm_id == ""
    assert result.transfer_status == "blocked"
    assert all(
        item.eligibility_status == "ineligible"
        for item in result.arm_summaries
        if item.arm_id != "d0-ffull-c0"
    )


def test_missing_preregistered_case_makes_arm_ineligible():
    result = _select(
        _analysis(omitted={("d0-fcausal-c0", "QP-DEV-007")})
    )
    summary = next(
        item for item in result.arm_summaries
        if item.arm_id == "d0-fcausal-c0"
    )

    assert summary.episode_count == 2
    assert summary.eligibility_status == "ineligible"
    assert "selection.arm.episode_count_mismatch" in (
        summary.eligibility_findings
    )
    assert "selection.arm.case_count_mismatch" in (
        summary.eligibility_findings
    )


def test_non_development_episode_is_rejected_before_selection():
    analysis = _analysis()
    first = analysis.episode_analyses[0]
    altered_body = {
        **{
            key: value
            for key, value in first.__dict__.items()
            if key != "analysis_sha256"
        },
        "split": "transfer",
    }
    transfer = QwenPyscfEpisodeAnalysisV1(
        **altered_body,
        analysis_sha256=canonical_sha256(altered_body),
    )
    episodes = tuple(sorted(
        (transfer,) + analysis.episode_analyses[1:],
        key=lambda item: item.episode_id,
    ))
    body = {
        **{
            key: value
            for key, value in analysis.__dict__.items()
            if key not in {"analysis_sha256", "episode_analyses"}
        },
        "episode_analyses": episodes,
    }
    mixed = QwenPyscfCampaignAnalysisV1(
        **body, analysis_sha256=canonical_sha256(body)
    )

    with pytest.raises(ContractError, match="non-development episodes"):
        _select(mixed)


def _transfer_inputs(tmp_path):
    xyz = tmp_path / "approved.xyz"
    xyz.write_text(
        "3\napproved water\n"
        "O 0.0 0.0 0.0\n"
        "H 0.0 0.757 0.586\n"
        "H 0.0 -0.757 0.586\n",
        encoding="utf-8",
    )
    source = approved_xyz_source(xyz, expected_sha256=file_sha256(xyz))
    registry = {item.case_id: item for item in qwen_pyscf_cases_v1()}
    cases = tuple(
        bind_case_to_approved_xyz(registry[case_id], source)
        for case_id in (
            "QP-DEV-005",
            "QP-DEV-006",
            "QP-DEV-007",
            "QP-TR-001",
            "QP-TR-002",
            "QP-TR-003",
            "QP-TR-004",
        )
    )
    plans = build_episode_plans(
        cases=cases,
        arms=all_dfc_arms(max_concurrency=2),
        repeats=10,
        prompt_sha256=canonical_sha256("prompt"),
        tool_schema_sha256=canonical_sha256("tools"),
        task_order_sha256=canonical_sha256("order"),
        token_budget=1_000_000,
        tool_call_budget=256,
        wall_time_seconds=5400,
    )
    window = build_qwen_campaign_window(
        campaign_id="transfer-opening-test",
        activation=datetime(2026, 8, 4, tzinfo=timezone.utc),
    )
    freeze = build_frozen_transfer_manifest(
        window=window,
        source=source,
        cases=cases,
        plans=plans,
    )
    development_case_sha_by_id = {
        item.case_id: item.case_sha256
        for item in cases
        if item.split == "development"
    }
    return cases, plans, freeze, development_case_sha_by_id


def test_transfer_opening_selects_exactly_paired_h0_candidate_schedule(
    tmp_path,
):
    cases, plans, freeze, development_case_sha_by_id = _transfer_inputs(
        tmp_path
    )
    analysis = _analysis(case_sha_by_id=development_case_sha_by_id)
    selection = _select(analysis)

    opening, batches = build_qwen_pyscf_transfer_opening(
        implementation_manifest_sha256=canonical_sha256("manifest"),
        implementation_tree_sha256=canonical_sha256("tree"),
        freeze=freeze,
        development_analysis=analysis,
        selection=selection,
        provider_exposure_before_selection=0,
        transfer_cases=tuple(item for item in cases if item.split == "transfer"),
        transfer_plans=tuple(
            item
            for item in plans
            if item.case_sha256 in freeze.transfer_case_sha256s
        ),
    )

    assert opening.candidate_arm_id == selection.selected_arm_id
    assert len(opening.selected_plan_sha256s) == 80
    assert len(batches) == 40
    assert all(len(batch) == 2 for batch in batches)
    assert all(".d0-ffull-c0." in batch[0].episode_id for batch in batches)
    assert all(
        f".{selection.selected_arm_id}." in batch[1].episode_id
        for batch in batches
    )
    assert opening.paired_schedule == tuple(
        tuple(item.episode_id for item in batch) for batch in batches
    )
    with pytest.raises(ContractError, match="exposed a provider"):
        build_qwen_pyscf_transfer_opening(
            implementation_manifest_sha256=canonical_sha256("manifest"),
            implementation_tree_sha256=canonical_sha256("tree"),
            freeze=freeze,
            development_analysis=analysis,
            selection=selection,
            provider_exposure_before_selection=1,
            transfer_cases=tuple(
                item for item in cases if item.split == "transfer"
            ),
            transfer_plans=tuple(
                item
                for item in plans
                if item.case_sha256 in freeze.transfer_case_sha256s
            ),
        )


def test_transfer_opening_rejects_an_incomplete_frozen_pair_schedule(tmp_path):
    cases, plans, freeze, development_case_sha_by_id = _transfer_inputs(
        tmp_path
    )
    analysis = _analysis(case_sha_by_id=development_case_sha_by_id)
    selection = _select(analysis)
    transfer_plans = tuple(
        item
        for item in plans
        if item.case_sha256 in freeze.transfer_case_sha256s
    )
    removed = next(
        item
        for item in transfer_plans
        if selection.selected_arm_id in item.episode_id
    )

    with pytest.raises(ContractError, match="80 plans"):
        build_qwen_pyscf_transfer_opening(
            implementation_manifest_sha256=canonical_sha256("manifest"),
            implementation_tree_sha256=canonical_sha256("tree"),
            freeze=freeze,
            development_analysis=analysis,
            selection=selection,
            provider_exposure_before_selection=0,
            transfer_cases=tuple(
                item for item in cases if item.split == "transfer"
            ),
            transfer_plans=tuple(
                item for item in transfer_plans if item is not removed
            ),
        )


def test_transfer_opening_rejects_development_evidence_from_another_freeze(
    tmp_path,
):
    cases, plans, freeze, _ = _transfer_inputs(tmp_path)
    unrelated_analysis = _analysis()
    selection = _select(unrelated_analysis)

    with pytest.raises(
        ContractError,
        match="development analysis case registry differs",
    ):
        build_qwen_pyscf_transfer_opening(
            implementation_manifest_sha256=canonical_sha256("manifest"),
            implementation_tree_sha256=canonical_sha256("tree"),
            freeze=freeze,
            development_analysis=unrelated_analysis,
            selection=selection,
            provider_exposure_before_selection=0,
            transfer_cases=tuple(
                item for item in cases if item.split == "transfer"
            ),
            transfer_plans=tuple(
                item
                for item in plans
                if item.case_sha256 in freeze.transfer_case_sha256s
            ),
        )


def test_transfer_opening_recomputes_the_development_selection(tmp_path):
    cases, plans, freeze, development_case_sha_by_id = _transfer_inputs(
        tmp_path
    )
    analysis = _analysis(case_sha_by_id=development_case_sha_by_id)
    selection = _select(analysis)
    another_eligible_arm = next(
        item.arm_id
        for item in selection.arm_summaries
        if item.arm_id not in {selection.baseline_arm_id, selection.selected_arm_id}
        and item.eligibility_status == "eligible"
    )
    body = {
        key: value
        for key, value in selection.__dict__.items()
        if key != "selection_sha256"
    }
    body["selected_arm_id"] = another_eligible_arm
    forged_selection = selection.__class__(
        **body,
        selection_sha256=canonical_sha256(body),
    )

    with pytest.raises(
        ContractError,
        match="selection is not the deterministic result",
    ):
        build_qwen_pyscf_transfer_opening(
            implementation_manifest_sha256=canonical_sha256("manifest"),
            implementation_tree_sha256=canonical_sha256("tree"),
            freeze=freeze,
            development_analysis=analysis,
            selection=forged_selection,
            provider_exposure_before_selection=0,
            transfer_cases=tuple(
                item for item in cases if item.split == "transfer"
            ),
            transfer_plans=tuple(
                item
                for item in plans
                if item.case_sha256 in freeze.transfer_case_sha256s
            ),
        )


def test_transfer_opening_rejects_a_frozen_but_wrong_pairing_key(tmp_path):
    cases, plans, _, development_case_sha_by_id = _transfer_inputs(tmp_path)
    analysis = _analysis(case_sha_by_id=development_case_sha_by_id)
    selection = _select(analysis)
    transfer = tuple(
        item for item in plans if item.case_sha256 in {
            case.case_sha256 for case in cases if case.split == "transfer"
        }
    )
    target = next(
        item
        for item in transfer
        if selection.selected_arm_id in item.episode_id
    )
    other = next(
        item
        for item in transfer
        if item.case_sha256 != target.case_sha256
        and item.repeat_index == target.repeat_index
    )
    body = {
        key: value
        for key, value in target.__dict__.items()
        if key != "plan_sha256"
    }
    body["pairing_key"] = other.pairing_key
    malformed = target.__class__(
        **body,
        plan_sha256=canonical_sha256(body),
    )
    altered_plans = tuple(
        malformed if item is target else item for item in plans
    )
    window = build_qwen_campaign_window(
        campaign_id="wrong-pairing-freeze-test",
        activation=datetime(2026, 8, 4, tzinfo=timezone.utc),
    )
    source_case = next(item for item in cases if item.split == "transfer")
    # Reuse the exact approved source binding from the normal helper by
    # rebuilding it from the source file that is already bound to every case.
    xyz = tmp_path / "approved.xyz"
    source = approved_xyz_source(xyz, expected_sha256=file_sha256(xyz))
    malformed_freeze = build_frozen_transfer_manifest(
        window=window,
        source=source,
        cases=cases,
        plans=altered_plans,
    )
    assert source.expected_sha256 in source_case.source_sha256s

    with pytest.raises(ContractError, match="pairing key differs"):
        build_qwen_pyscf_transfer_opening(
            implementation_manifest_sha256=canonical_sha256("manifest"),
            implementation_tree_sha256=canonical_sha256("tree"),
            freeze=malformed_freeze,
            development_analysis=analysis,
            selection=selection,
            provider_exposure_before_selection=0,
            transfer_cases=tuple(
                item for item in cases if item.split == "transfer"
            ),
            transfer_plans=tuple(
                item
                for item in altered_plans
                if item.case_sha256 in malformed_freeze.transfer_case_sha256s
            ),
        )


def test_transfer_opening_replay_rejects_cross_case_pairs(tmp_path):
    cases, plans, freeze, development_case_sha_by_id = _transfer_inputs(
        tmp_path
    )
    analysis = _analysis(case_sha_by_id=development_case_sha_by_id)
    selection = _select(analysis)
    opening, _ = build_qwen_pyscf_transfer_opening(
        implementation_manifest_sha256=canonical_sha256("manifest"),
        implementation_tree_sha256=canonical_sha256("tree"),
        freeze=freeze,
        development_analysis=analysis,
        selection=selection,
        provider_exposure_before_selection=0,
        transfer_cases=tuple(item for item in cases if item.split == "transfer"),
        transfer_plans=tuple(
            item
            for item in plans
            if item.case_sha256 in freeze.transfer_case_sha256s
        ),
    )
    schedule = list(opening.paired_schedule)
    first = schedule[0]
    second = schedule[1]
    schedule[0] = (first[0], second[1])
    schedule[1] = (second[0], first[1])
    body = {
        key: value
        for key, value in opening.__dict__.items()
        if key != "opening_sha256"
    }
    body["paired_schedule"] = tuple(schedule)

    with pytest.raises(ContractError, match="share one case and repeat"):
        opening.__class__(
            **body,
            opening_sha256=canonical_sha256(body),
        )


def _transfer_result_records(*, opening, freeze, plans, verdict_by_arm):
    plan_by_id = {item.episode_id: item for item in plans}
    ledger_rows = []
    analysis_rows = []
    for episode_id in sorted(
        item for pair in opening.paired_schedule for item in pair
    ):
        plan = plan_by_id[episode_id]
        stem, repeat_text = episode_id.rsplit(".r", 1)
        case_id, arm_suffix = stem.rsplit(".d", 1)
        arm_id = "d" + arm_suffix
        repeat_index = int(repeat_text)
        decomposition, feedback, critic = _factors(arm_id)
        verdict = verdict_by_arm[arm_id]
        ledger_body = {
            "schema_version": "chemsmart.qwen-episode-ledger.v1",
            "episode_id": episode_id,
            "split": "transfer",
            "case_sha256": plan.case_sha256,
            "plan_sha256": plan.plan_sha256,
            "hypothesis_sha256": plan.hypothesis.hypothesis_sha256,
            "experiment_config_sha256": plan.experiment_config.config_sha256,
            "source_sha256": canonical_sha256("source"),
            "workspace_binding_sha256": canonical_sha256(
                {"workspace": episode_id}
            ),
            "result_sha256": canonical_sha256({"result": episode_id}),
            "grade_sha256": canonical_sha256({"grade": episode_id}),
            "session_terminal_state": "complete",
            "scientific_state": "previewed",
            "verdict": verdict,
            "failure_class": "",
            "factor_realization_status": "valid",
            "factor_realization_findings": (),
            "safety_violations": (),
        }
        ledger_rows.append(
            QwenPyscfEpisodeLedgerV1(
                **ledger_body,
                ledger_sha256=canonical_sha256(ledger_body),
            )
        )
        analysis_body = {
            "schema_version": "chemsmart.qwen-episode-analysis.v1",
            "episode_id": episode_id,
            "split": "transfer",
            "case_id": case_id,
            "case_sha256": plan.case_sha256,
            "repeat_index": repeat_index,
            "arm_id": arm_id,
            "decomposition": decomposition,
            "feedback_projection": feedback,
            "critic": critic,
            "ledger_verdict": verdict,
            "strict_scientific_verdict": verdict,
            "safety_status": "green",
            "safety_violations": (),
            "factor_realization_status": "valid",
            "factor_realization_findings": (),
            "failure_class": "",
            "observation_status": "valid",
            "observation_findings": (),
            "decomposition_activated": decomposition,
            "model_attempts": 1,
            "input_tokens": 100,
            "output_tokens": 50,
            "reasoning_tokens": 25,
            "provider_latency_millis": 1000,
            "successful_tool_calls": 3,
            "failed_tool_calls": 0,
            "specialist_outcomes": (),
            "specialist_merge_status": (
                "merged" if decomposition else "not_dispatched"
            ),
            "critic_status": "complete" if critic else "not_enabled",
            "critic_findings": 0,
            "critic_critical_findings": 0,
            "feedback_contract_families": (),
            "feedback_coverage_status": "not_exercised",
            "feedback_bytes_by_tool": (),
        }
        analysis_rows.append(
            QwenPyscfEpisodeAnalysisV1(
                **analysis_body,
                analysis_sha256=canonical_sha256(analysis_body),
            )
        )
    split_body = {
        "schema_version": "chemsmart.qwen-split-ledger.v1",
        "campaign_window_sha256": freeze.campaign_window_sha256,
        "freeze_manifest_sha256": opening.transfer_freeze_manifest_sha256,
        "split": "transfer",
        "started_utc": "2026-08-04T00:00:00Z",
        "observed_at_utc": "2026-08-04T01:00:00Z",
        "termination_reason": "split_complete",
        "episode_ledgers": tuple(ledger_rows),
        "last_hypothesis_sha256": ledger_rows[-1].hypothesis_sha256,
    }
    split = QwenPyscfSplitLedgerV1(
        **split_body,
        ledger_sha256=canonical_sha256(split_body),
    )
    campaign_body = {
        "schema_version": "chemsmart.qwen-campaign-analysis.v1",
        "split_ledger_sha256s": (split.ledger_sha256,),
        "episode_analyses": tuple(analysis_rows),
        "pair_analyses": (),
        "strict_passes": sum(
            item.strict_scientific_verdict == "pass" for item in analysis_rows
        ),
        "strict_failures": sum(
            item.strict_scientific_verdict == "fail" for item in analysis_rows
        ),
        "strict_inconclusive": 0,
        "safety_red_episodes": 0,
        "factor_invalid_episodes": 0,
        "factor_not_observed_episodes": 0,
    }
    campaign = QwenPyscfCampaignAnalysisV1(
        **campaign_body,
        analysis_sha256=canonical_sha256(campaign_body),
    )
    return split, campaign


def test_transfer_evaluation_binds_multifactor_candidate_and_all_40_pairs(
    tmp_path,
):
    cases, plans, freeze, development_case_sha_by_id = _transfer_inputs(
        tmp_path
    )
    verdicts = {
        arm_id: ("pass" if arm_id in {"d0-ffull-c0", "d1-fcausal-c1"} else "fail")
        for arm_id in _ARMS
    }
    development = _analysis(
        verdict_by_arm=verdicts,
        case_sha_by_id=development_case_sha_by_id,
    )
    selection = _select(development)
    assert selection.selected_arm_id == "d1-fcausal-c1"
    opening, _ = build_qwen_pyscf_transfer_opening(
        implementation_manifest_sha256=canonical_sha256("manifest"),
        implementation_tree_sha256=canonical_sha256("tree"),
        freeze=freeze,
        development_analysis=development,
        selection=selection,
        provider_exposure_before_selection=0,
        transfer_cases=tuple(item for item in cases if item.split == "transfer"),
        transfer_plans=tuple(
            item
            for item in plans
            if item.case_sha256 in freeze.transfer_case_sha256s
        ),
    )
    split, transfer = _transfer_result_records(
        opening=opening,
        freeze=freeze,
        plans=plans,
        verdict_by_arm={"d0-ffull-c0": "fail", "d1-fcausal-c1": "pass"},
    )

    receipt = build_qwen_pyscf_transfer_evaluation_receipt(
        implementation_manifest_sha256=canonical_sha256("manifest"),
        implementation_tree_sha256=canonical_sha256("tree"),
        freeze=freeze,
        development_analysis=development,
        selection=selection,
        opening=opening,
        transfer_split_ledgers=(split,),
        transfer_analysis=transfer,
    )

    assert len(receipt.pair_evaluations) == 40
    assert all(
        item.comparison == "candidate_better"
        for item in receipt.pair_evaluations
    )
    assert receipt.baseline_strict_counts == (
        ("fail", 40),
        ("inconclusive", 0),
        ("pass", 0),
    )
    assert receipt.candidate_strict_counts == (
        ("fail", 0),
        ("inconclusive", 0),
        ("pass", 40),
    )
    assert receipt.transfer_opening_sha256 == opening.opening_sha256
    assert receipt.transfer_analysis_sha256 == transfer.analysis_sha256

    forged_body = {
        key: value
        for key, value in transfer.__dict__.items()
        if key != "analysis_sha256"
    }
    forged_body["split_ledger_sha256s"] = (canonical_sha256("other-ledger"),)
    forged = QwenPyscfCampaignAnalysisV1(
        **forged_body,
        analysis_sha256=canonical_sha256(forged_body),
    )
    with pytest.raises(ContractError, match="ledger ancestry differs"):
        build_qwen_pyscf_transfer_evaluation_receipt(
            implementation_manifest_sha256=canonical_sha256("manifest"),
            implementation_tree_sha256=canonical_sha256("tree"),
            freeze=freeze,
            development_analysis=development,
            selection=selection,
            opening=opening,
            transfer_split_ledgers=(split,),
            transfer_analysis=forged,
        )
