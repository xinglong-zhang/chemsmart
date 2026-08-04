from __future__ import annotations

from datetime import datetime, timedelta, timezone
from pathlib import Path
import threading
import time
from types import SimpleNamespace

import pytest

from chemsmart.agent._contracts import ContractError, canonical_sha256, file_sha256
from chemsmart.agent.experiments.qwen_pyscf_campaign import (
    QwenPyscfCampaignControllerV1,
    approved_xyz_source,
    bind_case_to_approved_xyz,
    build_frozen_transfer_manifest,
    build_qwen_campaign_window,
    prepare_qwen_campaign_plans,
)
from chemsmart.agent.experiments.qwen_pyscf_dfc import (
    build_episode_plans_from_preparations,
    build_qwen_dfc_arm,
    build_qwen_experiment_preparation,
)
from chemsmart.agent.experiments.qwen_pyscf_fixtures import qwen_pyscf_cases_v1


_ACTIVATION = datetime(2026, 8, 4, 0, 0, tzinfo=timezone.utc)


class _MutableClock:
    def __init__(self, value: datetime) -> None:
        self.value = value

    def __call__(self) -> datetime:
        return self.value

    def sleep(self, seconds: float) -> None:
        self.value += timedelta(seconds=seconds)


class _ProviderFailure(RuntimeError):
    def __init__(self, error_class: str, retry_after_seconds=None) -> None:
        super().__init__("private provider detail must not enter the ledger")
        self.error_class = error_class
        self.retry_after_seconds = retry_after_seconds


def _write_provider_config(path: Path, *, fallback: bool = False) -> None:
    fallback_block = "\nfallback: [qwen-secondary]" if fallback else ""
    secondary = (
        "\n  qwen-secondary:\n"
        "    type: openai-chat-completions\n"
        "    api_key_env: ALIBABA_TOKEN_PLAN_KEY\n"
        "    model: qwen3.8-max\n"
        "    base_url: https://token-plan.ap-southeast-1.maas.aliyuncs.com/"
        "compatible-mode/v1\n"
        "    reasoning_effort: xhigh\n"
        "    preserve_thinking: true\n"
        if fallback
        else ""
    )
    path.write_text(
        "active: alibaba-token-plan"
        + fallback_block
        + "\nproviders:\n"
        "  alibaba-token-plan:\n"
        "    type: openai-chat-completions\n"
        "    api_key_env: ALIBABA_TOKEN_PLAN_KEY\n"
        "    model: qwen3.8-max\n"
        "    base_url: https://token-plan.ap-southeast-1.maas.aliyuncs.com/"
        "compatible-mode/v1\n"
        "    reasoning_effort: xhigh\n"
        "    preserve_thinking: true\n"
        + secondary,
        encoding="utf-8",
    )


def _campaign_inputs(
    tmp_path: Path,
    *,
    repeats: int = 1,
    concurrency: int = 1,
    decomposition: bool = False,
    critic: bool = False,
):
    xyz = tmp_path / "source.xyz"
    xyz.write_text(
        "3\napproved water\nO 0.0 0.0 0.0\nH 0.0 0.757 0.586\nH 0.0 -0.757 0.586\n",
        encoding="utf-8",
    )
    source = approved_xyz_source(xyz, expected_sha256=file_sha256(xyz))
    registry = {item.case_id: item for item in qwen_pyscf_cases_v1()}
    cases = tuple(
        bind_case_to_approved_xyz(registry[case_id], source)
        for case_id in ("QP-DEV-003", "QP-TR-001")
    )
    arm = build_qwen_dfc_arm(
        decomposition=decomposition,
        feedback_projection="full-v1",
        critic=critic,
        max_concurrency=concurrency,
    )
    preparations = tuple(
        build_qwen_experiment_preparation(
            case=case,
            arm=arm,
            repeat_index=repeat_index,
            task_spec_sha256=canonical_sha256(
                {"case": case.case_id, "source": source.expected_sha256}
            ),
            artifact_sha256s=(source.expected_sha256,),
            provider_profile_sha256=canonical_sha256("profile-v1"),
            prompt_sha256=canonical_sha256(
                {"base-prompt": case.case_id, "arm": arm.arm_id}
            ),
            tool_schema_sha256=canonical_sha256("tools-v1"),
            task_order_sha256=canonical_sha256("fixed-order-v1"),
            token_budget=1_000_000,
            tool_call_budget=256,
            wall_time_seconds=5400,
        )
        for repeat_index in range(repeats)
        for case in cases
    )
    plans = build_episode_plans_from_preparations(
        cases=cases,
        arms=(arm,),
        preparations=preparations,
    )
    window = build_qwen_campaign_window(
        campaign_id="qwen-pyscf-24h-test", activation=_ACTIVATION
    )
    freeze = build_frozen_transfer_manifest(
        window=window, source=source, cases=cases, plans=plans
    )
    provider_config = tmp_path / "agent.yaml"
    _write_provider_config(provider_config)
    secret_file = tmp_path / "api.env"
    secret_file.write_text("ALIBABA_TOKEN_PLAN_KEY=unused-by-test\n", encoding="utf-8")
    return {
        "source": source,
        "cases": cases,
        "plans": plans,
        "window": window,
        "freeze": freeze,
        "provider_config": provider_config,
        "secret_file": secret_file,
        "workspace_root": tmp_path / "workspaces",
    }


def _controller(
    inputs,
    *,
    runner,
    clock,
    concurrency=1,
    resume_hook=None,
    outcome_observer=None,
    authorize_transfer=False,
):
    kwargs = {
        "window": inputs["window"],
        "freeze": inputs["freeze"],
        "source": inputs["source"],
        "cases": inputs["cases"],
        "provider_config_file": inputs["provider_config"],
        "secret_file": inputs["secret_file"],
        "workspace_root": inputs["workspace_root"],
        "max_concurrency": concurrency,
        "runner": runner,
        "clock": clock,
        "outcome_observer": outcome_observer,
        "authorized_transfer_plan_sha256s": (
            inputs["freeze"].transfer_plan_sha256s
            if authorize_transfer
            else ()
        ),
    }
    if resume_hook is not None:
        kwargs["resume_hook"] = resume_hook
    if hasattr(clock, "sleep"):
        kwargs["sleeper"] = clock.sleep
    return QwenPyscfCampaignControllerV1(**kwargs)


def _planned_result(**kwargs):
    workspace = Path(kwargs["workspace"])
    private = workspace / ".chemsmart-agent" / "runs" / "fake"
    private.mkdir(parents=True, exist_ok=True)
    (private / "events.jsonl").write_text("private raw event\n", encoding="utf-8")
    config = kwargs["experiment_config"]
    activated = bool(config.decomposition)
    requested_roles = (
        ("dag-specialist", "pyscf-specialist", "scientific-specialist")
        if activated
        else ()
    )
    specialists = tuple(
        {
            "role": role,
            "status": "complete",
            "result_sha256": canonical_sha256({"result": role}),
            "outcome_sha256": canonical_sha256({"outcome": role}),
        }
        for role in requested_roles
    )
    def output_envelope(role):
        body = {
            "schema_version": "chemsmart.specialist-output-envelope.v1",
            "mode": "strict_json",
            "raw_text_sha256": canonical_sha256({"raw": role}),
            "normalized_object_sha256": canonical_sha256(
                {"normalized": role}
            ),
            "ignored_prefix_bytes": 0,
            "ignored_suffix_bytes": 0,
        }
        return {**body, "receipt_sha256": canonical_sha256(body)}

    provider_sessions = tuple(
        {
            "role": role,
            "session_id": f"fresh-{role}",
            "feedback_projection": config.feedback_projection,
            "output_envelope_receipt": output_envelope(role),
        }
        for role in requested_roles
    )
    critic = {
        "status": "not_enabled",
        "review_sha256": "",
        "outcome_sha256": "",
        "nonsecret_error_class": "",
    }
    if config.critic:
        critic = {
            "status": "complete",
            "review_sha256": canonical_sha256("critic-review"),
            "outcome_sha256": canonical_sha256("critic-outcome"),
            "nonsecret_error_class": "",
        }
        provider_sessions = (
            *provider_sessions,
            {
                "role": "critic",
                "session_id": "fresh-critic",
                "feedback_projection": config.feedback_projection,
                "output_envelope_receipt": output_envelope("critic"),
            },
        )
    body = {
        # This synthetic controller fixture predates Runtime V2 host-oracle
        # bundles.  Production live results may not set this marker.
        "legacy_transcript_fixture": True,
        "terminal_state": "planned",
        "final_text": "A consequential value is missing.",
        "public_transcript": (
            {
                "role": "tool",
                "content": {
                    "schema_version": "chemsmart.tool-result.v1",
                    "tool": "record_scientific_decision",
                    "status": "ok",
                    "result": {
                        "record_sha256": "a" * 64,
                        "uncertainties": (
                            "The method is unspecified and the electronic "
                            "charge and multiplicity may also be unknown.",
                        ),
                        "diagnostics": (
                            "Clarification is required before materialization.",
                        ),
                    },
                },
            },
        ),
        "experiment_observations": {
            "experiment_config_sha256": (
                kwargs["experiment_config"].config_sha256
            ),
            "feedback_projection": config.feedback_projection,
            "observed_experiment_config_sha256": (
                kwargs["experiment_config"].config_sha256
            ),
            "preparation_sha256": canonical_sha256(
                {
                    "case": kwargs["experiment_case"].case_id,
                    "repeat": kwargs["experiment_repeat_index"],
                    "config": kwargs["experiment_config"].config_sha256,
                }
            ),
            "gate": {
                "activated": activated,
                "requested_roles": requested_roles,
            },
            "specialists": specialists,
            "merge": {
                "status": "merged" if activated else "not_dispatched",
            },
            "critic": critic,
            "provider_sessions": provider_sessions,
            "nonsecret_specialist_error_class": "",
        },
        "successful_tool_calls": 1,
        "failed_tool_calls": 0,
    }
    body["result_sha256"] = canonical_sha256(body)
    return body


def _rewrite_observations(result, **changes):
    result["experiment_observations"] = {
        **result["experiment_observations"],
        **changes,
    }
    result["result_sha256"] = canonical_sha256(
        {key: value for key, value in result.items() if key != "result_sha256"}
    )
    return result


def test_window_is_exactly_24_hours_and_never_dispatches_outside_it(tmp_path):
    inputs = _campaign_inputs(tmp_path)
    calls = []
    clock = _MutableClock(_ACTIVATION - timedelta(seconds=1))
    controller = _controller(
        inputs, runner=lambda **kwargs: calls.append(kwargs), clock=clock
    )
    plans = tuple(
        item
        for item in inputs["plans"]
        if inputs["cases"][0].case_sha256 == item.case_sha256
    )

    before = controller.run_split(split="development", plans=plans)
    clock.value = _ACTIVATION + timedelta(hours=24)
    at_deadline = controller.run_split(split="development", plans=plans)

    assert before.termination_reason == "not_activated"
    assert at_deadline.termination_reason == "deadline_reached"
    assert calls == []
    assert (
        datetime.fromisoformat(inputs["window"].deadline_utc.replace("Z", "+00:00"))
        - datetime.fromisoformat(inputs["window"].activation_utc.replace("Z", "+00:00"))
    ) == timedelta(hours=24)


def test_episode_is_not_dispatched_when_its_frozen_budget_crosses_deadline(
    tmp_path,
):
    inputs = _campaign_inputs(tmp_path)
    calls = []
    clock = _MutableClock(
        _ACTIVATION + timedelta(hours=24, seconds=-5399)
    )
    controller = _controller(
        inputs, runner=lambda **kwargs: calls.append(kwargs), clock=clock
    )
    plans = tuple(
        item
        for item in inputs["plans"]
        if inputs["cases"][0].case_sha256 == item.case_sha256
    )

    ledger = controller.run_split(split="development", plans=plans)

    assert ledger.termination_reason == "deadline_reached"
    assert ledger.episode_ledgers == ()
    assert calls == []


def test_episode_is_not_dispatched_at_final_hour_reserve_boundary(tmp_path):
    inputs = _campaign_inputs(tmp_path)
    calls = []
    # The frozen episode allowance is 5,400 s; exactly 3,600 s must remain
    # after it for the Goal's final local validation hour.
    clock = _MutableClock(
        _ACTIVATION + timedelta(hours=24, seconds=-(5400 + 3600))
    )
    controller = _controller(
        inputs, runner=lambda **kwargs: calls.append(kwargs), clock=clock
    )
    plans = tuple(
        item
        for item in inputs["plans"]
        if inputs["cases"][0].case_sha256 == item.case_sha256
    )

    ledger = controller.run_split(split="development", plans=plans)

    assert ledger.termination_reason == "deadline_reached"
    assert ledger.episode_ledgers == ()
    assert calls == []


def test_development_run_is_qwen_only_preview_only_and_path_free(tmp_path):
    inputs = _campaign_inputs(tmp_path)
    observed = []

    def runner(**kwargs):
        observed.append(kwargs)
        return _planned_result(**kwargs)

    controller = _controller(
        inputs,
        runner=runner,
        clock=_MutableClock(_ACTIVATION + timedelta(seconds=1)),
    )
    plans = tuple(
        item
        for item in inputs["plans"]
        if item.case_sha256 == inputs["cases"][0].case_sha256
    )
    ledger = controller.run_split(split="development", plans=plans)

    assert ledger.termination_reason == "split_complete"
    assert len(ledger.episode_ledgers) == 1
    assert ledger.episode_ledgers[0].verdict == "pass"
    assert observed[0]["provider"] == "alibaba-token-plan"
    assert observed[0]["execution_enabled"] is False
    assert observed[0]["approval_file"] is None
    workspace = Path(observed[0]["workspace"])
    assert file_sha256(workspace / "approved.xyz") == inputs["source"].expected_sha256
    assert (workspace / ".chemsmart-agent" / "runs" / "fake" / "events.jsonl").is_file()
    public = ledger.public_summary_json()
    assert str(tmp_path) not in public
    assert "private raw event" not in public
    assert ledger.episode_ledgers[0].result_sha256
    assert ledger.episode_ledgers[0].grade_sha256


def test_activated_decomposition_requires_every_requested_specialist_and_merge(
    tmp_path,
):
    inputs = _campaign_inputs(tmp_path, decomposition=True)

    def runner(**kwargs):
        result = _planned_result(**kwargs)
        return _rewrite_observations(
            result,
            specialists=(),
            merge={"status": "failed"},
            provider_sessions=(
                {"role": "dag-specialist", "session_id": "failed-dag"},
            ),
            nonsecret_specialist_error_class="ContractError",
        )

    controller = _controller(
        inputs,
        runner=runner,
        clock=_MutableClock(_ACTIVATION + timedelta(seconds=1)),
    )
    development = tuple(
        plan
        for plan in inputs["plans"]
        if plan.case_sha256 == inputs["cases"][0].case_sha256
    )

    ledger = controller.run_split(split="development", plans=development)
    episode = ledger.episode_ledgers[0]

    assert ledger.termination_reason == "split_complete"
    assert episode.failure_class == "experiment_factor_invalid"
    assert episode.factor_realization_status == "invalid"
    assert episode.verdict == "inconclusive"
    assert episode.safety_violations == ()
    assert "experiment.factor.specialist_set_incomplete" in (
        episode.factor_realization_findings
    )
    assert "experiment.factor.worker_sessions_incomplete" in (
        episode.factor_realization_findings
    )
    assert "experiment.factor.merge_not_realized" in (
        episode.factor_realization_findings
    )


def test_activated_decomposition_accepts_exact_worker_set_and_merge(tmp_path):
    inputs = _campaign_inputs(tmp_path, decomposition=True)
    controller = _controller(
        inputs,
        runner=_planned_result,
        clock=_MutableClock(_ACTIVATION + timedelta(seconds=1)),
    )
    development = tuple(
        plan
        for plan in inputs["plans"]
        if plan.case_sha256 == inputs["cases"][0].case_sha256
    )

    ledger = controller.run_split(split="development", plans=development)
    episode = ledger.episode_ledgers[0]

    assert episode.factor_realization_status == "valid"
    assert episode.factor_realization_findings == ()
    assert episode.failure_class == ""
    assert episode.verdict == "pass"


def test_normalized_worker_envelope_is_valid_but_not_raw_schema_success(
    tmp_path,
):
    inputs = _campaign_inputs(tmp_path, decomposition=True)
    observed_modes = []

    def runner(**kwargs):
        result = _planned_result(**kwargs)
        sessions = list(result["experiment_observations"]["provider_sessions"])
        first = dict(sessions[0])
        receipt = dict(first["output_envelope_receipt"])
        receipt["mode"] = "single_json_object_extracted"
        receipt["ignored_prefix_bytes"] = 18
        receipt["receipt_sha256"] = canonical_sha256(
            {key: value for key, value in receipt.items() if key != "receipt_sha256"}
        )
        first["output_envelope_receipt"] = receipt
        sessions[0] = first
        observed_modes.extend(
            item["output_envelope_receipt"]["mode"] for item in sessions
        )
        return _rewrite_observations(result, provider_sessions=tuple(sessions))

    controller = _controller(
        inputs,
        runner=runner,
        clock=_MutableClock(_ACTIVATION + timedelta(seconds=1)),
    )
    development = tuple(
        plan
        for plan in inputs["plans"]
        if plan.case_sha256 == inputs["cases"][0].case_sha256
    )
    episode = controller.run_split(
        split="development", plans=development
    ).episode_ledgers[0]

    assert episode.factor_realization_status == "valid"
    assert observed_modes.count("single_json_object_extracted") == 1
    assert observed_modes.count("strict_json") == 2


def test_tampered_worker_envelope_invalidates_factor_realization(tmp_path):
    inputs = _campaign_inputs(tmp_path, decomposition=True)

    def runner(**kwargs):
        result = _planned_result(**kwargs)
        sessions = list(result["experiment_observations"]["provider_sessions"])
        first = dict(sessions[0])
        receipt = dict(first["output_envelope_receipt"])
        receipt["receipt_sha256"] = "0" * 64
        first["output_envelope_receipt"] = receipt
        sessions[0] = first
        return _rewrite_observations(result, provider_sessions=tuple(sessions))

    controller = _controller(
        inputs,
        runner=runner,
        clock=_MutableClock(_ACTIVATION + timedelta(seconds=1)),
    )
    development = tuple(
        plan
        for plan in inputs["plans"]
        if plan.case_sha256 == inputs["cases"][0].case_sha256
    )
    episode = controller.run_split(
        split="development", plans=development
    ).episode_ledgers[0]

    assert episode.factor_realization_status == "invalid"
    assert "experiment.factor.output_envelope_receipt_invalid" in (
        episode.factor_realization_findings
    )
    assert episode.verdict == "inconclusive"


def test_critic_on_requires_one_fresh_review_session(tmp_path):
    inputs = _campaign_inputs(tmp_path, critic=True)

    def runner(**kwargs):
        result = _planned_result(**kwargs)
        return _rewrite_observations(
            result,
            critic={
                "status": "failed",
                "review_sha256": "",
                "outcome_sha256": "",
                "nonsecret_error_class": "ContractError",
            },
            provider_sessions=(),
        )

    controller = _controller(
        inputs,
        runner=runner,
        clock=_MutableClock(_ACTIVATION + timedelta(seconds=1)),
    )
    development = tuple(
        plan
        for plan in inputs["plans"]
        if plan.case_sha256 == inputs["cases"][0].case_sha256
    )

    episode = controller.run_split(
        split="development", plans=development
    ).episode_ledgers[0]

    assert episode.factor_realization_status == "invalid"
    assert episode.failure_class == "experiment_factor_invalid"
    assert "experiment.factor.critic_not_realized" in (
        episode.factor_realization_findings
    )
    assert "experiment.factor.critic_session_count" in (
        episode.factor_realization_findings
    )


def test_critic_on_accepts_one_fresh_review_session(tmp_path):
    inputs = _campaign_inputs(tmp_path, critic=True)
    controller = _controller(
        inputs,
        runner=_planned_result,
        clock=_MutableClock(_ACTIVATION + timedelta(seconds=1)),
    )
    development = tuple(
        plan
        for plan in inputs["plans"]
        if plan.case_sha256 == inputs["cases"][0].case_sha256
    )

    episode = controller.run_split(
        split="development", plans=development
    ).episode_ledgers[0]

    assert episode.factor_realization_status == "valid"
    assert episode.factor_realization_findings == ()
    assert episode.failure_class == ""


def test_critic_off_rejects_an_unconfigured_critic_session(tmp_path):
    inputs = _campaign_inputs(tmp_path, critic=False)

    def runner(**kwargs):
        result = _planned_result(**kwargs)
        return _rewrite_observations(
            result,
            critic={
                "status": "complete",
                "review_sha256": canonical_sha256("unexpected-review"),
                "outcome_sha256": canonical_sha256("unexpected-outcome"),
                "nonsecret_error_class": "",
            },
            provider_sessions=(
                {"role": "critic", "session_id": "unexpected-critic"},
            ),
        )

    controller = _controller(
        inputs,
        runner=runner,
        clock=_MutableClock(_ACTIVATION + timedelta(seconds=1)),
    )
    development = tuple(
        plan
        for plan in inputs["plans"]
        if plan.case_sha256 == inputs["cases"][0].case_sha256
    )

    episode = controller.run_split(
        split="development", plans=development
    ).episode_ledgers[0]

    assert episode.factor_realization_status == "invalid"
    assert episode.failure_class == "experiment_factor_invalid"
    assert "experiment.factor.unexpected_critic" in (
        episode.factor_realization_findings
    )
    assert "experiment.factor.unexpected_critic_session" in (
        episode.factor_realization_findings
    )


def test_transfer_is_closed_until_an_exact_frozen_plan_set_is_authorized(
    tmp_path,
):
    inputs = _campaign_inputs(tmp_path)
    controller = _controller(
        inputs,
        runner=_planned_result,
        clock=_MutableClock(_ACTIVATION + timedelta(seconds=1)),
    )
    transfer = tuple(
        item
        for item in inputs["plans"]
        if item.case_sha256 == inputs["cases"][1].case_sha256
    )

    with pytest.raises(ContractError, match="has not been opened"):
        controller.run_split(split="transfer", plans=transfer)


def test_transfer_accepts_only_exact_frozen_plans(tmp_path):
    inputs = _campaign_inputs(tmp_path, repeats=2)
    controller = _controller(
        inputs,
        runner=_planned_result,
        clock=_MutableClock(_ACTIVATION + timedelta(seconds=1)),
        authorize_transfer=True,
    )
    transfer = tuple(
        item
        for item in inputs["plans"]
        if item.case_sha256 == inputs["cases"][1].case_sha256
    )

    ledger = controller.run_split(split="transfer", plans=transfer)

    assert len(ledger.episode_ledgers) == 2
    with pytest.raises(ContractError, match="crosses"):
        controller.run_split(split="transfer", plans=(next(
            item
            for item in inputs["plans"]
            if item.case_sha256 == inputs["cases"][0].case_sha256
        ),))

    unfrozen_arm = build_qwen_dfc_arm(
            decomposition=False,
            feedback_projection="full-v1",
            critic=False,
            max_concurrency=1,
    )
    unfrozen_preparation = build_qwen_experiment_preparation(
        case=inputs["cases"][1],
        arm=unfrozen_arm,
        repeat_index=2,
        task_spec_sha256=canonical_sha256("unfrozen-task"),
        artifact_sha256s=(inputs["source"].expected_sha256,),
        provider_profile_sha256=canonical_sha256("profile-v1"),
        prompt_sha256=canonical_sha256(
            {
                "base-prompt": inputs["cases"][1].case_id,
                "arm": unfrozen_arm.arm_id,
            }
        ),
        tool_schema_sha256=canonical_sha256("tools-v1"),
        task_order_sha256=canonical_sha256("fixed-order-v1"),
        token_budget=1_000_000,
        tool_call_budget=256,
        wall_time_seconds=5400,
    )
    unfrozen = build_episode_plans_from_preparations(
        cases=inputs["cases"],
        arms=(unfrozen_arm,),
        preparations=(unfrozen_preparation,),
    )
    new_transfer_repeat = next(
        item
        for item in unfrozen
        if item.case_sha256 == inputs["cases"][1].case_sha256
        and item.repeat_index == 2
    )
    with pytest.raises(ContractError, match="not frozen"):
        controller.run_split(
            split="transfer", plans=(new_transfer_repeat,)
        )


def test_rate_limit_uses_resume_hook_then_continues_with_next_hypothesis(tmp_path):
    inputs = _campaign_inputs(tmp_path, repeats=2)
    clock = _MutableClock(_ACTIVATION + timedelta(seconds=1))
    calls = 0

    def runner(**kwargs):
        nonlocal calls
        calls += 1
        if calls == 1:
            raise _ProviderFailure("rate_limited", retry_after_seconds=7)
        return _planned_result(**kwargs)

    controller = _controller(inputs, runner=runner, clock=clock)
    development = tuple(
        item
        for item in inputs["plans"]
        if item.case_sha256 == inputs["cases"][0].case_sha256
    )
    ledger = controller.run_split(split="development", plans=development)

    assert ledger.termination_reason == "split_complete"
    assert tuple(item.failure_class for item in ledger.episode_ledgers) == (
        "rate_limited",
        "",
    )
    assert len({item.hypothesis_sha256 for item in ledger.episode_ledgers}) == 2
    assert clock.value == _ACTIVATION + timedelta(seconds=8)
    assert "private provider detail" not in ledger.public_summary_json()


def test_quota_window_with_retry_after_resumes_before_deadline(tmp_path):
    inputs = _campaign_inputs(tmp_path, repeats=2)
    clock = _MutableClock(_ACTIVATION + timedelta(seconds=1))
    calls = 0

    def runner(**kwargs):
        nonlocal calls
        calls += 1
        if calls == 1:
            raise _ProviderFailure("quota_exhausted", retry_after_seconds=11)
        return _planned_result(**kwargs)

    controller = _controller(inputs, runner=runner, clock=clock)
    development = tuple(
        item
        for item in inputs["plans"]
        if item.case_sha256 == inputs["cases"][0].case_sha256
    )

    ledger = controller.run_split(split="development", plans=development)

    assert ledger.termination_reason == "split_complete"
    assert tuple(item.failure_class for item in ledger.episode_ledgers) == (
        "quota_exhausted",
        "",
    )
    assert clock.value == _ACTIVATION + timedelta(seconds=12)


def test_concurrency_uses_distinct_mutable_workspaces(tmp_path):
    inputs = _campaign_inputs(tmp_path, repeats=2, concurrency=2)
    clock = _MutableClock(_ACTIVATION + timedelta(seconds=1))
    lock = threading.Lock()
    active = 0
    max_active = 0
    workspaces = []

    def runner(**kwargs):
        nonlocal active, max_active
        with lock:
            active += 1
            max_active = max(max_active, active)
            workspaces.append(Path(kwargs["workspace"]))
        time.sleep(0.02)
        with lock:
            active -= 1
        return _planned_result(**kwargs)

    controller = _controller(
        inputs, runner=runner, clock=clock, concurrency=2
    )
    development = tuple(
        item
        for item in inputs["plans"]
        if item.case_sha256 == inputs["cases"][0].case_sha256
    )
    ledger = controller.run_split(split="development", plans=development)

    assert len(ledger.episode_ledgers) == 2
    assert max_active == 2
    assert len(set(workspaces)) == 2
    assert all(
        file_sha256(path / "approved.xyz") == inputs["source"].expected_sha256
        for path in workspaces
    )


def test_concurrent_terminal_failure_does_not_drop_a_completed_peer(tmp_path):
    inputs = _campaign_inputs(tmp_path, repeats=2, concurrency=2)

    def runner(**kwargs):
        if kwargs["experiment_repeat_index"] == 0:
            raise _ProviderFailure("quota_exhausted")
        return _planned_result(**kwargs)

    controller = _controller(
        inputs,
        runner=runner,
        clock=_MutableClock(_ACTIVATION + timedelta(seconds=1)),
        concurrency=2,
    )
    development = tuple(
        item
        for item in inputs["plans"]
        if item.case_sha256 == inputs["cases"][0].case_sha256
    )

    ledger = controller.run_split(split="development", plans=development)

    assert ledger.termination_reason == "quota_exhausted"
    assert len(ledger.episode_ledgers) == 2
    assert tuple(item.failure_class for item in ledger.episode_ledgers) == (
        "quota_exhausted",
        "",
    )


def test_outcome_observer_receives_the_controller_bound_result_and_grade(
    tmp_path,
):
    inputs = _campaign_inputs(tmp_path)
    observed = []
    controller = _controller(
        inputs,
        runner=_planned_result,
        clock=_MutableClock(_ACTIVATION + timedelta(seconds=1)),
        outcome_observer=lambda **row: observed.append(row),
    )
    development = tuple(
        item
        for item in inputs["plans"]
        if item.case_sha256 == inputs["cases"][0].case_sha256
    )

    split = controller.run_split(split="development", plans=development)

    assert len(observed) == 1
    row = observed[0]
    assert row["result"]["result_sha256"] == row["ledger"].result_sha256
    assert row["grade"].grade_sha256 == row["ledger"].grade_sha256
    assert row["ledger"] == split.episode_ledgers[0]


def test_controller_rejects_a_configured_provider_fallback(tmp_path):
    inputs = _campaign_inputs(tmp_path)
    _write_provider_config(inputs["provider_config"], fallback=True)

    with pytest.raises(ContractError, match="no fallback"):
        _controller(
            inputs,
            runner=_planned_result,
            clock=_MutableClock(_ACTIVATION + timedelta(seconds=1)),
        )


def test_campaign_planner_freezes_provider_free_probe_contracts(
    tmp_path, monkeypatch
):
    import chemsmart.agent.live_session as live

    inputs = _campaign_inputs(tmp_path)

    snapshot = SimpleNamespace(snapshot_sha256=canonical_sha256("host-snapshot"))
    snapshot_builds = []
    observed_snapshots = []

    def build_snapshot(**kwargs):
        snapshot_builds.append(kwargs)
        return snapshot

    def probe(**kwargs):
        observed_snapshots.append(kwargs["campaign_preparation_snapshot"])
        coordinate = Path(kwargs["workspace"]) / "approved.xyz"
        return build_qwen_experiment_preparation(
            case=kwargs["experiment_case"],
            arm=kwargs["experiment_arm"],
            repeat_index=kwargs["experiment_repeat_index"],
            task_spec_sha256=canonical_sha256(kwargs["task"]),
            artifact_sha256s=(file_sha256(coordinate),),
            provider_profile_sha256=canonical_sha256("profile"),
            prompt_sha256=canonical_sha256(
                {"case": kwargs["experiment_case"].case_id}
            ),
            tool_schema_sha256=canonical_sha256("live-tools"),
            task_order_sha256=canonical_sha256("live-order"),
            token_budget=1_000_000,
            tool_call_budget=256,
            wall_time_seconds=5400,
            host_snapshot_sha256=(
                kwargs["campaign_preparation_snapshot"].snapshot_sha256
            ),
        )

    monkeypatch.setattr(
        live, "build_campaign_preparation_host_snapshot", build_snapshot
    )
    monkeypatch.setattr(live, "probe_live_experiment_preparation", probe)
    arm = build_qwen_dfc_arm(
        decomposition=False,
        feedback_projection="causal-v1",
        critic=False,
    )
    plans, preparations, observed_snapshot = prepare_qwen_campaign_plans(
        source=inputs["source"],
        cases=inputs["cases"],
        arms=(arm,),
        repeats=1,
        provider_config_file=inputs["provider_config"],
        preparation_workspace_root=tmp_path / "preparation-workspaces",
    )

    assert len(plans) == len(preparations) == 2
    assert len(snapshot_builds) == 1
    assert observed_snapshot is snapshot
    assert observed_snapshots == [snapshot, snapshot]
    by_episode = {item.episode_id: item for item in preparations}
    assert all(
        plan.experiment_config.config_sha256
        == by_episode[plan.episode_id].experiment_config.config_sha256
        for plan in plans
    )
    assert all(
        snapshot.snapshot_sha256 in plan.hypothesis.source_sha256s
        for plan in plans
    )


def test_controller_treats_observed_config_drift_as_a_safety_red_line(tmp_path):
    inputs = _campaign_inputs(tmp_path)

    def runner(**kwargs):
        result = _planned_result(**kwargs)
        result["experiment_observations"] = {
            **result["experiment_observations"],
            "observed_experiment_config_sha256": "f" * 64,
        }
        result["result_sha256"] = canonical_sha256(
            {
                key: value
                for key, value in result.items()
                if key != "result_sha256"
            }
        )
        return result

    controller = _controller(
        inputs,
        runner=runner,
        clock=_MutableClock(_ACTIVATION + timedelta(seconds=1)),
    )
    development = tuple(
        item
        for item in inputs["plans"]
        if item.case_sha256 == inputs["cases"][0].case_sha256
    )

    ledger = controller.run_split(split="development", plans=development)

    assert ledger.termination_reason == "safety_red_line"
    assert ledger.episode_ledgers[0].failure_class == (
        "experiment_integrity_violation"
    )
    assert ledger.episode_ledgers[0].safety_violations == (
        "safety.experiment_config_mismatch",
    )


@pytest.mark.parametrize(
    ("provider_error", "termination", "public_failure"),
    (
        ("credential_invalid", "credential_revoked", "credential_revoked"),
        ("quota_exhausted", "quota_exhausted", "quota_exhausted"),
        ("provider_5xx", "provider_paused", "transient_transport"),
    ),
)
def test_provider_failures_use_a_closed_nonsecret_taxonomy(
    tmp_path, provider_error, termination, public_failure
):
    root = tmp_path / provider_error
    root.mkdir()
    inputs = _campaign_inputs(root)

    def runner(**kwargs):
        del kwargs
        raise _ProviderFailure(provider_error)

    controller = _controller(
        inputs,
        runner=runner,
        clock=_MutableClock(_ACTIVATION + timedelta(seconds=1)),
    )
    development = tuple(
        item
        for item in inputs["plans"]
        if item.case_sha256 == inputs["cases"][0].case_sha256
    )
    ledger = controller.run_split(split="development", plans=development)

    assert ledger.termination_reason == termination
    assert ledger.episode_ledgers[0].failure_class == public_failure
    assert "private provider detail" not in ledger.public_summary_json()
