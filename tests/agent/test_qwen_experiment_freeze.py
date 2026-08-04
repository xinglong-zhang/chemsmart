from __future__ import annotations

from types import SimpleNamespace

import pytest

from chemsmart.agent._contracts import ContractError, canonical_sha256
from chemsmart.agent.experiments.host_oracle import HostOracleInputBundleV1
from chemsmart.agent.experiments.qwen_pyscf_dfc import (
    bind_harness_experiment_config,
    build_qwen_dfc_arm,
)
from chemsmart.agent.experiments.qwen_pyscf_fixtures import (
    qwen_pyscf_cases_v1,
)
from chemsmart.agent.provider_config import (
    ALIBABA_TOKEN_PLAN_ENDPOINT,
    AgentProviderProfileV1,
)
from chemsmart.agent.tool_specs import AgentToolSurfaceV1
from chemsmart.agent.live_session import (
    build_campaign_preparation_host_snapshot,
    probe_live_experiment_preparation,
    run_live_agent_session,
)


def _empty_host_oracle_bundle(**kwargs) -> HostOracleInputBundleV1:
    body = {
        "schema_version": "chemsmart.host-oracle-input-bundle.v1",
        "session_id": kwargs["session_id"],
        "event_stream_head_sha256": kwargs["event_stream_head_sha256"],
        "observations": (),
        "successful_tool_calls": kwargs["successful_tool_calls"],
        "failed_tool_calls": kwargs["failed_tool_calls"],
        "tool_counts": (),
        "tool_actions_sha256": canonical_sha256(()),
    }
    return HostOracleInputBundleV1(
        **body, bundle_sha256=canonical_sha256(body)
    )


def _profile() -> AgentProviderProfileV1:
    body = {
        "schema_version": "chemsmart.agent-provider-profile.v1",
        "profile_name": "alibaba-token-plan",
        "provider": "alibaba-token-plan",
        "wire_protocol": "openai-chat-completions",
        "api_key_env": "ALIBABA_TOKEN_PLAN_KEY",
        "model": "qwen3.8-max",
        "endpoint": ALIBABA_TOKEN_PLAN_ENDPOINT,
        "reasoning_effort": "xhigh",
        "preserve_thinking": True,
        "context_tokens": 1_000_000,
        "max_output_tokens": 262_144,
    }
    return AgentProviderProfileV1(
        **body, profile_sha256=canonical_sha256(body)
    )


def _surface() -> AgentToolSurfaceV1:
    tools = ()
    return AgentToolSurfaceV1(
        schema_version="chemsmart.agent-tool-surface.v1",
        profile="command_compiled_frontier",
        tool_definitions=tools,
        tool_schema_sha256=canonical_sha256(tools),
    )


def _patch_local_preparation(monkeypatch, live):
    profile = _profile()
    selection = SimpleNamespace(
        active_profile=profile, fallback_profiles=(), selection_sha256="1" * 64
    )
    monkeypatch.setattr(
        live, "load_agent_provider_selection", lambda *args, **kwargs: selection
    )
    monkeypatch.setattr(
        live,
        "load_program_capabilities",
        lambda: SimpleNamespace(registry_sha256="2" * 64),
    )
    monkeypatch.setattr(
        live,
        "build_live_click_schema",
        lambda: SimpleNamespace(schema_sha256="3" * 64),
    )
    monkeypatch.setattr(live, "_bootstrap_conformance", lambda **kwargs: ((), ()))
    monkeypatch.setattr(live, "_observe_environments", lambda: ((), (), ()))
    monkeypatch.setattr(
        live, "build_command_compiled_tool_surface", lambda registry: _surface()
    )
    return profile


def _local_snapshot(monkeypatch, live, workspace):
    profile = _patch_local_preparation(monkeypatch, live)
    xyz = workspace / "approved.xyz"
    xyz.write_text(
        "3\nwater\nO 0 0 0\nH 0 0.7 0.5\nH 0 -0.7 0.5\n",
        encoding="utf-8",
    )
    snapshot = build_campaign_preparation_host_snapshot(
        provider="alibaba-token-plan",
        provider_config_file=workspace / "unused.yaml",
        workspace=workspace,
    )
    return profile, xyz, snapshot


def test_probe_observes_exact_config_without_secret_or_provider(tmp_path, monkeypatch):
    import chemsmart.agent.live_session as live

    _patch_local_preparation(monkeypatch, live)
    monkeypatch.setattr(
        live,
        "load_secret_lease",
        lambda **kwargs: pytest.fail("preparation probe loaded a credential"),
    )
    xyz = tmp_path / "approved.xyz"
    xyz.write_text(
        "3\nwater\nO 0 0 0\nH 0 0.7 0.5\nH 0 -0.7 0.5\n",
        encoding="utf-8",
    )
    case = qwen_pyscf_cases_v1()[0]
    arm = build_qwen_dfc_arm(
        decomposition=False,
        feedback_projection="causal-v1",
        critic=False,
    )

    preparation = probe_live_experiment_preparation(
        task=case.task,
        provider="alibaba-token-plan",
        provider_config_file=tmp_path / "unused.yaml",
        workspace=tmp_path,
        experiment_arm=arm,
        experiment_case=case,
        experiment_repeat_index=0,
    )

    assert preparation.provider_calls == 0
    assert preparation.engine_calls == 0
    assert preparation.approval_files == 0
    assert preparation.experiment_config.prompt_sha256
    assert preparation.experiment_config.tool_schema_sha256 == (
        _surface().tool_schema_sha256
    )


def test_campaign_snapshot_rejects_another_artifact(tmp_path, monkeypatch):
    import chemsmart.agent.live_session as live

    _profile_value, xyz, snapshot = _local_snapshot(monkeypatch, live, tmp_path)
    xyz.write_text(
        "3\nchanged\nO 0 0 0\nH 0 0.8 0.5\nH 0 -0.8 0.5\n",
        encoding="utf-8",
    )
    case = qwen_pyscf_cases_v1()[0]
    arm = build_qwen_dfc_arm(
        decomposition=False,
        feedback_projection="full-v1",
        critic=False,
    )

    with pytest.raises(ContractError, match="artifact mismatch"):
        probe_live_experiment_preparation(
            task=case.task,
            provider="alibaba-token-plan",
            provider_config_file=tmp_path / "unused.yaml",
            workspace=tmp_path,
            experiment_arm=arm,
            experiment_case=case,
            experiment_repeat_index=0,
            campaign_preparation_snapshot=snapshot,
        )


def test_campaign_snapshot_rejects_provider_profile_drift(tmp_path, monkeypatch):
    import chemsmart.agent.live_session as live

    profile, _xyz, snapshot = _local_snapshot(monkeypatch, live, tmp_path)
    profile_body = {
        key: value
        for key, value in profile.__dict__.items()
        if key != "profile_sha256"
    }
    profile_body["context_tokens"] -= 1
    drifted = AgentProviderProfileV1(
        **profile_body,
        profile_sha256=canonical_sha256(profile_body),
    )
    monkeypatch.setattr(
        live,
        "load_agent_provider_selection",
        lambda *args, **kwargs: SimpleNamespace(
            active_profile=drifted,
            fallback_profiles=(),
            selection_sha256="5" * 64,
        ),
    )
    case = qwen_pyscf_cases_v1()[0]
    arm = build_qwen_dfc_arm(
        decomposition=False,
        feedback_projection="full-v1",
        critic=False,
    )

    with pytest.raises(ContractError, match="provider profile mismatch"):
        probe_live_experiment_preparation(
            task=case.task,
            provider="alibaba-token-plan",
            provider_config_file=tmp_path / "unused.yaml",
            workspace=tmp_path,
            experiment_arm=arm,
            experiment_case=case,
            experiment_repeat_index=0,
            campaign_preparation_snapshot=snapshot,
        )


def test_campaign_snapshot_rejects_schema_object_drift(tmp_path, monkeypatch):
    import chemsmart.agent.live_session as live

    _profile_value, _xyz, snapshot = _local_snapshot(
        monkeypatch, live, tmp_path
    )
    snapshot.live_schema.schema_sha256 = "6" * 64
    case = qwen_pyscf_cases_v1()[0]
    arm = build_qwen_dfc_arm(
        decomposition=False,
        feedback_projection="full-v1",
        critic=False,
    )

    with pytest.raises(ContractError, match="live schema mismatch"):
        probe_live_experiment_preparation(
            task=case.task,
            provider="alibaba-token-plan",
            provider_config_file=tmp_path / "unused.yaml",
            workspace=tmp_path,
            experiment_arm=arm,
            experiment_case=case,
            experiment_repeat_index=0,
            campaign_preparation_snapshot=snapshot,
        )


def test_initial_coordinator_prompt_is_factor_and_repeat_blind(
    tmp_path, monkeypatch
):
    import chemsmart.agent.live_session as live

    _patch_local_preparation(monkeypatch, live)
    xyz = tmp_path / "approved.xyz"
    xyz.write_text(
        "3\nwater\nO 0 0 0\nH 0 0.7 0.5\nH 0 -0.7 0.5\n",
        encoding="utf-8",
    )
    case = qwen_pyscf_cases_v1()[0]
    preparations = tuple(
        probe_live_experiment_preparation(
            task=case.task,
            provider="alibaba-token-plan",
            provider_config_file=tmp_path / "unused.yaml",
            workspace=tmp_path,
            experiment_arm=build_qwen_dfc_arm(
                decomposition=decomposition,
                feedback_projection=feedback,
                critic=critic,
            ),
            experiment_case=case,
            experiment_repeat_index=repeat,
        )
        for decomposition, feedback, critic, repeat in (
            (False, "full-v1", False, 0),
            (True, "causal-v1", True, 2),
        )
    )

    assert len(
        {item.experiment_config.prompt_sha256 for item in preparations}
    ) == 1


def test_live_session_rejects_prompt_drift_before_specialist_or_provider(
    tmp_path, monkeypatch
):
    import chemsmart.agent.live_session as live

    _patch_local_preparation(monkeypatch, live)
    xyz = tmp_path / "approved.xyz"
    xyz.write_text(
        "3\nwater\nO 0 0 0\nH 0 0.7 0.5\nH 0 -0.7 0.5\n",
        encoding="utf-8",
    )
    case = qwen_pyscf_cases_v1()[0]
    arm = build_qwen_dfc_arm(
        decomposition=False,
        feedback_projection="full-v1",
        critic=False,
    )
    observed = probe_live_experiment_preparation(
        task=case.task,
        provider="alibaba-token-plan",
        provider_config_file=tmp_path / "unused.yaml",
        workspace=tmp_path,
        experiment_arm=arm,
        experiment_case=case,
        experiment_repeat_index=0,
    )
    drifted = bind_harness_experiment_config(
        arm=arm,
        experiment_id=observed.episode_id,
        prompt_sha256=canonical_sha256("another prompt"),
        tool_schema_sha256=(
            observed.experiment_config.tool_schema_sha256
        ),
        task_order_sha256=(
            observed.experiment_config.task_order_sha256
        ),
        token_budget=observed.experiment_config.token_budget,
        tool_call_budget=observed.experiment_config.tool_call_budget,
        wall_time_seconds=observed.experiment_config.wall_time_seconds,
    )
    calls = {"lease": 0, "specialist": 0}

    class _Host:
        def __init__(self, **kwargs):
            self.surface = kwargs["tool_surface"]
            self.scientific_decisions = {}

    class _Campaign:
        @classmethod
        def start(cls, **kwargs):
            calls["specialist"] += 1
            return cls()

    monkeypatch.setattr(live, "CommandCompiledToolHostV1", _Host)
    monkeypatch.setattr(live, "LiveSpecialistCampaignV1", _Campaign)
    monkeypatch.setattr(
        live,
        "load_secret_lease",
        lambda **kwargs: calls.__setitem__("lease", calls["lease"] + 1),
    )

    with pytest.raises(ContractError, match="prompt_sha256"):
        run_live_agent_session(
            task=case.task,
            provider="alibaba-token-plan",
            provider_config_file=tmp_path / "unused.yaml",
            secret_file=tmp_path / "unused.env",
            workspace=tmp_path,
            execution_enabled=False,
            approval_file=None,
            experiment_arm=arm,
            experiment_case=case,
            experiment_repeat_index=0,
            experiment_config=drifted,
        )

    assert calls == {"lease": 0, "specialist": 0}


def test_exact_frozen_config_is_returned_in_public_observations(
    tmp_path, monkeypatch
):
    import chemsmart.agent.live_session as live

    _patch_local_preparation(monkeypatch, live)
    xyz = tmp_path / "approved.xyz"
    xyz.write_text(
        "3\nwater\nO 0 0 0\nH 0 0.7 0.5\nH 0 -0.7 0.5\n",
        encoding="utf-8",
    )
    local_probe_calls = {"conformance": 0, "environment": 0}

    def conformance(**kwargs):
        del kwargs
        local_probe_calls["conformance"] += 1
        return (), ()

    def environments():
        local_probe_calls["environment"] += 1
        return (), (), ()

    monkeypatch.setattr(live, "_bootstrap_conformance", conformance)
    monkeypatch.setattr(live, "_observe_environments", environments)
    snapshot = build_campaign_preparation_host_snapshot(
        provider="alibaba-token-plan",
        provider_config_file=tmp_path / "unused.yaml",
        workspace=tmp_path,
    )
    assert local_probe_calls == {"conformance": 1, "environment": 1}

    def repeated_probe(**kwargs):
        del kwargs
        pytest.fail("episode repeated campaign conformance")

    monkeypatch.setattr(live, "_bootstrap_conformance", repeated_probe)
    monkeypatch.setattr(
        live,
        "_observe_environments",
        lambda: pytest.fail("episode repeated campaign environment probe"),
    )
    case = qwen_pyscf_cases_v1()[0]
    arm = build_qwen_dfc_arm(
        decomposition=False,
        feedback_projection="full-v1",
        critic=False,
    )
    preparation = probe_live_experiment_preparation(
        task=case.task,
        provider="alibaba-token-plan",
        provider_config_file=tmp_path / "unused.yaml",
        workspace=tmp_path,
        experiment_arm=arm,
        experiment_case=case,
        experiment_repeat_index=0,
        campaign_preparation_snapshot=snapshot,
    )

    class _Host:
        def __init__(self, **kwargs):
            self.surface = kwargs["tool_surface"]
            self.scientific_decisions = {}

    class _Campaign:
        def __init__(self, config):
            self.config = config

        @classmethod
        def start(cls, **kwargs):
            return cls(kwargs["experiment_config"])

        def coordinator_advisory_record(self):
            return {}

        def run_critic(self, **kwargs):
            del kwargs

        def public_observation_record(self, *, coordinator_usage):
            del coordinator_usage
            body = {
                "schema_version": (
                    "chemsmart.live-harness-experiment-observations.v1"
                ),
                "experiment_config_sha256": self.config.config_sha256,
                "specialists": (),
                "critic": {"status": "not_enabled"},
            }
            return {**body, "record_sha256": canonical_sha256(body)}

    class _Runner:
        def __init__(self, **kwargs):
            del kwargs

        def run(self, **kwargs):
            del kwargs
            return SimpleNamespace(
                terminal_state="planned",
                final_text="Planning stopped before execution.",
                public_transcript=(),
                public_transcript_sha256=canonical_sha256(()),
                successful_tool_calls=0,
                failed_tool_calls=0,
                event_stream_head_sha256="4" * 64,
            )

    monkeypatch.setattr(live, "CommandCompiledToolHostV1", _Host)
    monkeypatch.setattr(live, "LiveSpecialistCampaignV1", _Campaign)
    monkeypatch.setattr(live, "UnifiedSessionRunner", _Runner)
    monkeypatch.setattr(
        live, "build_host_oracle_input_bundle", _empty_host_oracle_bundle
    )
    monkeypatch.setattr(
        live,
        "load_secret_lease",
        lambda **kwargs: SimpleNamespace(provider=kwargs["provider"]),
    )

    result = run_live_agent_session(
        task=case.task,
        provider="alibaba-token-plan",
        provider_config_file=tmp_path / "unused.yaml",
        secret_file=tmp_path / "unused.env",
        workspace=tmp_path,
        execution_enabled=False,
        approval_file=None,
        experiment_arm=arm,
        experiment_case=case,
        experiment_repeat_index=0,
        experiment_config=preparation.experiment_config,
        campaign_preparation_snapshot=snapshot,
    )

    assert result.experiment_observations[
        "experiment_config_sha256"
    ] == preparation.experiment_config.config_sha256
    assert result.experiment_observations[
        "observed_experiment_config_sha256"
    ] == preparation.experiment_config.config_sha256
    assert result.experiment_observations["preparation_sha256"] == (
        preparation.preparation_sha256
    )
    assert preparation.host_snapshot_sha256 == snapshot.snapshot_sha256
    assert result.experiment_observations["host_snapshot_sha256"] == (
        snapshot.snapshot_sha256
    )
    assert local_probe_calls == {"conformance": 1, "environment": 1}
