from __future__ import annotations

from dataclasses import asdict
import json
import threading
import time
from types import SimpleNamespace

import pytest

from chemsmart.agent._contracts import ContractError
from chemsmart.agent._contracts import canonical_json, canonical_sha256
from chemsmart.agent.experiments.host_oracle import HostOracleInputBundleV1
from chemsmart.agent.experiments.qwen_pyscf_dfc import (
    bind_harness_experiment_config,
    build_qwen_dfc_arm,
)
from chemsmart.agent.experiments.qwen_pyscf_fixtures import (
    qwen_pyscf_cases_v1,
)
from chemsmart.agent.execution import build_scientific_decision_record
from chemsmart.agent.feedback import project_tool_feedback
from chemsmart.agent.live_specialists import (
    LiveSpecialistCampaignV1,
    _bounded_public_json_object,
    build_experiment_seed_plan,
    build_f_invariant_critic_candidate,
    derive_specialist_budget,
)
from chemsmart.agent.live_session import (
    _bind_feedback_receipts,
    run_live_agent_session,
)
from chemsmart.agent.provider_config import (
    ALIBABA_TOKEN_PLAN_CONTEXT_TOKENS,
    ALIBABA_TOKEN_PLAN_ENDPOINT,
    ALIBABA_TOKEN_PLAN_MAX_OUTPUT_TOKENS,
    AgentProviderProfileV1,
)
from chemsmart.agent.tool_specs import AgentToolSurfaceV1


def test_specialist_output_envelope_distinguishes_raw_and_normalized_json():
    strict, strict_receipt = _bounded_public_json_object('{"status":"complete"}')
    assert strict == {"status": "complete"}
    assert strict_receipt.mode == "strict_json"
    assert strict_receipt.ignored_prefix_bytes == 0
    assert strict_receipt.ignored_suffix_bytes == 0

    normalized, normalized_receipt = _bounded_public_json_object(
        'Evidence gathered.\n{"status":"complete"}\n'
    )
    assert normalized == strict
    assert normalized_receipt.mode == "single_json_object_extracted"
    assert normalized_receipt.ignored_prefix_bytes > 0
    assert normalized_receipt.ignored_suffix_bytes > 0
    assert normalized_receipt.raw_text_sha256 != strict_receipt.raw_text_sha256
    assert (
        normalized_receipt.normalized_object_sha256
        == strict_receipt.normalized_object_sha256
    )


@pytest.mark.parametrize(
    "text",
    (
        '```json\n{"status":"complete"}\n```',
        'first {"status":"complete"} second {"status":"blocked"}',
        'Evidence gathered.\n{"status":"complete"}\ntrailing prose',
        'python exploit.py\n{"status":"complete"}',
        'rm -rf /tmp/example\n{"status":"complete"}',
        'export TOKEN=value\n{"status":"complete"}',
        '! B3LYP def2-SVP\n{"status":"complete"}',
        'python3.10 exploit.py.\n{"status":"complete"}',
        '/usr/bin/python exploit.py.\n{"status":"complete"}',
        'env python exploit.py.\n{"status":"complete"}',
        'g16 job.com.\n{"status":"complete"}',
        'orca job.inp.\n{"status":"complete"}',
        'xtb molecule.xyz.\n{"status":"complete"}',
        '{"status":"failed","status":"complete"}',
        '{"status":"complete","nested":{"value":1,"value":2}}',
        'no object here',
    ),
)
def test_specialist_output_envelope_rejects_ambiguous_wrappers(text):
    with pytest.raises(ContractError):
        _bounded_public_json_object(text)


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


def _critic_candidate(
    *,
    candidate_id="candidate-water",
    task_spec_sha256="d" * 64,
):
    bundle = _empty_host_oracle_bundle(
        session_id="coordinator-session",
        event_stream_head_sha256="f" * 64,
        successful_tool_calls=0,
        failed_tool_calls=0,
    )
    return build_f_invariant_critic_candidate(
        candidate_id=candidate_id,
        task_spec_sha256=task_spec_sha256,
        host_oracle_input_bundle=bundle,
        coordinator_public_decisions=(),
    )


def _tool(name: str) -> dict:
    return {
        "type": "function",
        "function": {
            "name": name,
            "description": name,
            "parameters": {
                "type": "object",
                "properties": {},
                "additionalProperties": False,
            },
        },
    }


def _surface() -> AgentToolSurfaceV1:
    tools = tuple(
        _tool(name)
        for name in (
            "inspect_program_capability",
            "inspect_program_environment",
            "read_project_yaml",
            "validate_project_yaml",
            "inspect_calculation_artifact",
            "render_project_yaml",
            "repair_command",
            "execute_approved_program_node",
        )
    )
    return AgentToolSurfaceV1(
        schema_version="chemsmart.agent-tool-surface.v1",
        profile="test-coordinator",
        tool_definitions=tools,
        tool_schema_sha256=canonical_sha256(tools),
    )


def _profile() -> AgentProviderProfileV1:
    body = {
        "schema_version": "chemsmart.agent-provider-profile.v1",
        "profile_name": "qwen-token-plan",
        "provider": "alibaba-token-plan",
        "wire_protocol": "openai-chat-completions",
        "api_key_env": "ALIBABA_TOKEN_PLAN_KEY",
        "model": "qwen3.8-max",
        "endpoint": ALIBABA_TOKEN_PLAN_ENDPOINT,
        "reasoning_effort": "xhigh",
        "preserve_thinking": True,
        "context_tokens": ALIBABA_TOKEN_PLAN_CONTEXT_TOKENS,
        "max_output_tokens": ALIBABA_TOKEN_PLAN_MAX_OUTPUT_TOKENS,
    }
    return AgentProviderProfileV1(
        **body, profile_sha256=canonical_sha256(body)
    )


def _config(arm):
    return bind_harness_experiment_config(
        arm=arm,
        experiment_id="qp-dev-006.d1-ffull-c1.r0",
        prompt_sha256="a" * 64,
        tool_schema_sha256=_surface().tool_schema_sha256,
        task_order_sha256="b" * 64,
        token_budget=1_000_000,
        tool_call_budget=250,
        wall_time_seconds=5_000,
    )


class _BaseHost:
    def __init__(self, surface):
        self.surface = surface

    def record_seeded_evidence(self, turn_id):
        return None

    def dispatch(self, *, turn_id, tool_name, arguments):
        return {"tool": tool_name, "status": "ok"}


class _FakeRunner:
    created = []
    scientific_field_path = "scientific.suitability_assessment"
    scientific_value = (
        "Closed-shell B3LYP is suitable for this state; numerical accuracy "
        "remains unverified."
    )
    dag_field_path = "workflow.edges"
    dag_value = []

    def __init__(self, *, host, event_store, credential_lease, provider_config):
        self.host = host
        self.event_store = event_store
        self.provider_config = provider_config
        self.__class__.created.append(self)

    def run(self, *, messages, envelope, hypothesis, network_budget,
            feedback_projection):
        system = messages[0]["content"]
        self.feedback_projection = feedback_projection
        assert self.provider_config.model == "qwen3.8-max"
        assert envelope.session_id == self.event_store.session_id
        assert network_budget.engine_calls == 0
        assert feedback_projection in {"full-v1", "causal-v1"}
        exposed = {
            item["function"]["name"] for item in self.host.surface.tool_definitions
        }
        assert "render_project_yaml" not in exposed
        assert "repair_command" not in exposed
        assert "execute_approved_program_node" not in exposed
        if "read-only computational-chemistry critic" in system:
            critic_context = json.loads(messages[1]["content"])
            assert "expected_observation" not in canonical_json(critic_context)
            output = {
                "status": "complete",
                "public_summary": "The candidate remains advisory.",
                "findings": [
                    {
                        "rule_id": "authority-separation",
                        "severity": "info",
                        "expected": "The host retains readiness authority.",
                        "observed": "The host retains readiness authority.",
                        "disposition": "open",
                    }
                ],
            }
        elif "scientific-specialist" in system:
            output = {
                "status": "complete",
                "public_summary": "Preserve the requested DFT method.",
                "field_proposals": [
                    {
                        "field_path": self.scientific_field_path,
                        "value": self.scientific_value,
                    }
                ],
                "unresolved_fields": [],
            }
        elif "pyscf-specialist" in system:
            output = {
                "status": "complete",
                "public_summary": "Use the requested basis literal.",
                "field_proposals": [
                    {"field_path": "project.sp.basis", "value": "def2-svp"}
                ],
                "unresolved_fields": [],
            }
        else:
            output = {
                "status": "complete",
                "public_summary": "OPT produces the downstream geometry.",
                "field_proposals": [
                    {
                        "field_path": self.dag_field_path,
                        "value": self.dag_value,
                    }
                ],
                "unresolved_fields": [],
            }
        text = canonical_json(output)
        transcript = (
            messages[0],
            messages[1],
            {"role": "assistant", "content": text},
        )
        return SimpleNamespace(
            terminal_state="blocked",
            final_text=text,
            public_transcript=transcript,
            public_transcript_sha256=canonical_sha256(transcript),
            successful_tool_calls=0,
            failed_tool_calls=0,
            event_stream_head_sha256="c" * 64,
        )


def _lease_loader(**kwargs):
    return SimpleNamespace(provider=kwargs["provider"])


def test_seed_plan_keeps_initial_sp_and_opt_as_siblings():
    case = next(
        item for item in qwen_pyscf_cases_v1()
        if item.case_id == "QP-DEV-006"
    )
    plan = build_experiment_seed_plan(
        case=case,
        task_spec_sha256="d" * 64,
        artifact_sha256s=("e" * 64,),
    )

    assert tuple(item.node_id for item in plan.nodes) == (
        "candidate-sp",
        "candidate-opt",
        "candidate-hess",
    )
    assert tuple(
        (item.source_node_id, item.target_node_id, item.edge_kind)
        for item in plan.edges
    ) == (("candidate-opt", "candidate-hess", "data"),)

    gpu_case = next(
        item for item in qwen_pyscf_cases_v1()
        if item.case_id == "QP-DEV-008"
    )
    gpu_plan = build_experiment_seed_plan(
        case=gpu_case,
        task_spec_sha256="d" * 64,
        artifact_sha256s=("e" * 64,),
    )
    assert tuple((item.stage, item.engine) for item in gpu_plan.nodes) == (
        ("sp", "gpu"),
    )


def test_budget_is_derived_from_all_active_participant_slots():
    arm = build_qwen_dfc_arm(
        decomposition=True, feedback_projection="full-v1", critic=True
    )
    budget = derive_specialist_budget(
        experiment_config=_config(arm), participant_slots=5
    )

    assert budget.token_budget == 200_000
    assert budget.tool_call_budget == 50
    assert budget.wall_time_seconds == 1_000


def test_critic_toggle_does_not_redistribute_specialist_budget(tmp_path):
    case = next(
        item for item in qwen_pyscf_cases_v1()
        if item.case_id == "QP-DEV-006"
    )
    plan = build_experiment_seed_plan(
        case=case,
        task_spec_sha256="d" * 64,
        artifact_sha256s=("e" * 64,),
    )
    budgets = []
    for critic in (False, True):
        arm = build_qwen_dfc_arm(
            decomposition=True,
            feedback_projection="full-v1",
            critic=critic,
        )
        campaign = LiveSpecialistCampaignV1.start(
            arm=arm,
            experiment_config=_config(arm),
            plan=plan,
            coordinator_session_id=f"coordinator-{critic}",
            public_context={"task": case.task},
            source_sha256s=case.source_sha256s,
            artifact_sha256s=("e" * 64,),
            base_tool_surface=_surface(),
            provider_profile=_profile(),
            secret_file=tmp_path / "unused.env",
            run_directory=tmp_path / str(critic),
            host_builder=lambda event_store, request: _BaseHost(_surface()),
            runner_factory=_FakeRunner,
            lease_loader=_lease_loader,
        )
        budgets.append(campaign.budget)

    assert budgets[0] == budgets[1]
    assert budgets[0].token_budget == 250_000
    assert budgets[0].tool_call_budget == 62
    assert budgets[0].wall_time_seconds == 1_250


def test_live_campaign_uses_fresh_qwen_runtime_sessions_and_public_receipts(
    tmp_path,
):
    _FakeRunner.created.clear()
    case = next(
        item for item in qwen_pyscf_cases_v1()
        if item.case_id == "QP-DEV-006"
    )
    plan = build_experiment_seed_plan(
        case=case,
        task_spec_sha256="d" * 64,
        artifact_sha256s=("e" * 64,),
    )
    arm = build_qwen_dfc_arm(
        decomposition=True, feedback_projection="full-v1", critic=True
    )
    surface = _surface()
    campaign = LiveSpecialistCampaignV1.start(
        arm=arm,
        experiment_config=_config(arm),
        plan=plan,
        coordinator_session_id="coordinator-session",
        public_context={"task": case.task},
        source_sha256s=(),
        artifact_sha256s=("e" * 64,),
        base_tool_surface=surface,
        provider_profile=_profile(),
        secret_file=tmp_path / "unused.env",
        run_directory=tmp_path,
        host_builder=lambda event_store, request: _BaseHost(surface),
        runner_factory=_FakeRunner,
        lease_loader=_lease_loader,
    )

    assert campaign.gate.activated is True
    assert len(campaign.outcomes) == 3
    assert campaign.merge is not None
    assert campaign.merge.status == "merged"
    assert len(_FakeRunner.created) == 3
    assert tuple(
        item.host.surface.profile.removesuffix("-read-only")
        for item in _FakeRunner.created
    ) == campaign.gate.requested_roles
    assert {
        item.feedback_projection for item in _FakeRunner.created
    } == {"full-v1"}
    assert len(
        {item.event_store.session_id for item in _FakeRunner.created}
    ) == 3
    advisory = campaign.coordinator_advisory_record()
    assert advisory["authority"].startswith("advisory_only")
    assert advisory["merge"]["status"] == "merged"
    merged = {
        item["field_path"]: item["value"]
        for item in advisory["merge"]["merged_fields"]
    }
    assert ";" in merged["scientific.suitability_assessment"]
    assert advisory["merge"]["validation_findings"] == ()

    campaign.run_critic(
        coordinator_session_id="coordinator-session",
        candidate=_critic_candidate(),
        public_context={"task": case.task},
        source_sha256s=(),
        artifact_sha256s=("e" * 64,),
    )
    record = campaign.public_observation_record()

    assert len(_FakeRunner.created) == 4
    assert record["critic"]["status"] == "complete"
    assert record["critic"]["findings"][0]["disposition"] == "open"
    assert record["critic"]["candidate_sha256"]
    assert record["critic"]["review"]["schema_version"] == (
        "chemsmart.critic-review-record.v1"
    )
    assert record["usage"]["provider_sessions"] == 4
    scientific_row = next(
        item for item in record["specialists"]
        if item["role"] == "scientific-specialist"
    )
    assert scientific_row["result_packet"]["public_summary"] == (
        "Preserve the requested DFT method."
    )
    assert scientific_row["candidate"]["proposals"][0][
        "field_path"
    ] == "scientific.suitability_assessment"
    assert (
        record["specialist_dispatch"]["provider_concurrency_limit"] == 1
    )
    assert (
        record["specialist_dispatch"]["provider_concurrency_observed"] == 1
    )
    assert record["usage"]["worker"]["provider_concurrency_observed"] == 1
    assert "provider_wall_time_sum_millis" in record["usage"]["worker"]
    assert "observed_wall_time_millis" in record["usage"]["worker"]
    assert all(
        "path" not in key.lower()
        for row in record["provider_sessions"]
        for key in row
    )
    assert record["record_sha256"] == canonical_sha256(
        {key: value for key, value in record.items() if key != "record_sha256"}
    )


def test_live_campaign_preserves_role_local_specialist_validation_finding(
    tmp_path,
):
    class _MalformedFinalRoleRunner(_FakeRunner):
        created = []
        dag_field_path = "workflow.dependency_assessment"
        dag_value = "python run.py; echo finished"

    case = next(
        item for item in qwen_pyscf_cases_v1()
        if item.case_id == "QP-DEV-006"
    )
    plan = build_experiment_seed_plan(
        case=case,
        task_spec_sha256="d" * 64,
        artifact_sha256s=("e" * 64,),
    )
    arm = build_qwen_dfc_arm(
        decomposition=True,
        feedback_projection="full-v1",
        critic=False,
        max_concurrency=3,
    )
    campaign = LiveSpecialistCampaignV1.start(
        arm=arm,
        experiment_config=_config(arm),
        plan=plan,
        coordinator_session_id="coordinator-session",
        public_context={"task": case.task},
        source_sha256s=(),
        artifact_sha256s=("e" * 64,),
        base_tool_surface=_surface(),
        provider_profile=_profile(),
        secret_file=tmp_path / "unused.env",
        run_directory=tmp_path,
        host_builder=lambda event_store, request: _BaseHost(_surface()),
        runner_factory=_MalformedFinalRoleRunner,
        lease_loader=_lease_loader,
    )

    advisory = campaign.coordinator_advisory_record()
    record = campaign.public_observation_record()
    finding = advisory["merge"]["validation_findings"][0]

    assert tuple(
        outcome.result_packet.role for outcome in campaign.outcomes
    ) == ("pyscf-specialist", "scientific-specialist")
    assert campaign.merge is not None
    assert advisory["merge"]["status"] == "merged"
    assert advisory["merge"]["nonsecret_error_class"]
    assert len(advisory["specialists"]) == 2
    assert len(advisory["gate"]["requested_roles"]) == 3
    assert finding["role"] == "dag-specialist"
    assert finding["rule_id"] == "specialist.advisory.shell-authority.v1"
    assert record["specialist_validation_findings"] == (finding,)
    assert record["nonsecret_specialist_error_class"]
    assert len(_MalformedFinalRoleRunner.created) == 3
    assert (
        record["specialist_dispatch"]["provider_concurrency_limit"] == 3
    )


def test_parallel_specialists_overlap_but_merge_and_observations_are_role_sorted(
    tmp_path,
):
    class _OverlappingRunner(_FakeRunner):
        created = []
        barrier = threading.Barrier(3)
        state_lock = threading.Lock()
        completion_order = []

        def run(self, **kwargs):
            role = self.host.surface.profile.removesuffix("-read-only")
            self.__class__.barrier.wait(timeout=5)
            time.sleep(
                {
                    "dag-specialist": 0.0,
                    "pyscf-specialist": 0.04,
                    "scientific-specialist": 0.08,
                }[role]
            )
            result = super().run(**kwargs)
            with self.__class__.state_lock:
                self.__class__.completion_order.append(role)
            return result

    lease_lock = threading.Lock()
    lease_objects = []

    def lease_loader(**kwargs):
        lease = SimpleNamespace(
            provider=kwargs["provider"], lease_id=object()
        )
        with lease_lock:
            lease_objects.append(lease)
        return lease

    case = next(
        item for item in qwen_pyscf_cases_v1()
        if item.case_id == "QP-DEV-006"
    )
    plan = build_experiment_seed_plan(
        case=case,
        task_spec_sha256="d" * 64,
        artifact_sha256s=("e" * 64,),
    )
    arm = build_qwen_dfc_arm(
        decomposition=True,
        feedback_projection="full-v1",
        critic=False,
        max_concurrency=3,
    )
    campaign = LiveSpecialistCampaignV1.start(
        arm=arm,
        experiment_config=_config(arm),
        plan=plan,
        coordinator_session_id="coordinator-session",
        public_context={"task": case.task},
        source_sha256s=(),
        artifact_sha256s=("e" * 64,),
        base_tool_surface=_surface(),
        provider_profile=_profile(),
        secret_file=tmp_path / "unused.env",
        run_directory=tmp_path,
        host_builder=lambda event_store, request: _BaseHost(_surface()),
        runner_factory=_OverlappingRunner,
        lease_loader=lease_loader,
    )
    record = campaign.public_observation_record()

    assert _OverlappingRunner.completion_order == [
        "dag-specialist",
        "pyscf-specialist",
        "scientific-specialist",
    ]
    assert tuple(
        item.result_packet.role for item in campaign.outcomes
    ) == (
        "dag-specialist",
        "pyscf-specialist",
        "scientific-specialist",
    )
    assert tuple(
        item["role"] for item in record["provider_sessions"]
    ) == (
        "dag-specialist",
        "pyscf-specialist",
        "scientific-specialist",
    )
    assert record["specialist_output_envelopes"]["strict_json_count"] == 3
    assert record["specialist_output_envelopes"]["normalized_count"] == 0
    assert all(
        item["output_envelope_receipt"]["mode"] == "strict_json"
        for item in record["provider_sessions"]
    )
    assert campaign.merge is not None
    assert campaign.merge.status == "merged"
    assert record["specialist_dispatch"]["provider_concurrency_limit"] == 3
    assert record["specialist_dispatch"]["provider_concurrency_observed"] == 3
    assert record["specialist_dispatch"]["observed_wall_time_millis"] > 0
    assert (
        record["specialist_dispatch"][
            "provider_session_wall_time_sum_millis"
        ]
        > record["specialist_dispatch"]["observed_wall_time_millis"]
    )
    assert len(lease_objects) == 3
    assert len({id(item) for item in lease_objects}) == 3
    assert len(
        {item.event_store.session_id for item in _OverlappingRunner.created}
    ) == 3
    worker_directories = tuple(
        (tmp_path / "specialist-sessions").iterdir()
    )
    assert len(worker_directories) == 3


def test_live_experiment_wires_advisory_before_coordinator_and_critic_after(
    monkeypatch, tmp_path,
):
    import chemsmart.agent.live_session as live

    case = next(
        item for item in qwen_pyscf_cases_v1()
        if item.case_id == "QP-DEV-006"
    )
    arm = build_qwen_dfc_arm(
        decomposition=True, feedback_projection="causal-v1", critic=True
    )
    profile = _profile()
    selection = SimpleNamespace(
        active_profile=profile,
        fallback_profiles=(),
    )
    surface = _surface()
    registry = SimpleNamespace(registry_sha256="1" * 64)
    live_schema = SimpleNamespace(schema_sha256="2" * 64)

    class _CoordinatorHost:
        def __init__(self, **kwargs):
            self.surface = kwargs["tool_surface"]
            self.artifacts = dict(kwargs["artifacts"])
            self.scientific_identities = {}
            self.settings_objects = {}
            self.run_receipts = {}
            self.scientific_claim_evidence = {}
            self.functional_equivalence_receipts = {}
            self.substitution_approvals = {}
            self.capabilities = {}
            self.environments = {}
            self.program_bindings = {}
            self.engine_bindings = {}
            self.project_validations = {}
            self.scientific_decisions = {}

    class _Campaign:
        critic_context = None
        config_sha256 = ""

        @classmethod
        def start(cls, **kwargs):
            assert kwargs["provider_profile"].model == "qwen3.8-max"
            assert kwargs["experiment_config"].decomposition is True
            cls.config_sha256 = kwargs["experiment_config"].config_sha256
            return cls()

        def coordinator_advisory_record(self):
            return {
                "schema_version": "chemsmart.coordinator-specialist-advisory.v1",
                "authority": "advisory_only",
                "gate": {"activated": True},
                "specialists": (),
                "merge": {"status": "merged"},
            }

        def run_critic(self, **kwargs):
            self.__class__.critic_context = kwargs
            assert "expected_observation" not in canonical_json(
                kwargs["public_context"]
            )

        def public_observation_record(self, *, coordinator_usage):
            assert coordinator_usage["successful_tool_calls"] == 0
            body = {
                "schema_version": (
                    "chemsmart.live-harness-experiment-observations.v1"
                ),
                "experiment_config_sha256": self.config_sha256,
                "gate": {"activated": True},
                "specialists": (),
                "critic": {"status": "complete"},
                "provider_sessions": (),
                "usage": {"provider_sessions": 4},
            }
            return {**body, "record_sha256": canonical_sha256(body)}

    class _CoordinatorRunner:
        observed_messages = None

        def __init__(self, **kwargs):
            pass

        def run(self, *, messages, **kwargs):
            self.__class__.observed_messages = messages
            context = json.loads(messages[1]["content"])
            assert context["specialist_advisory"]["authority"] == (
                "advisory_only"
            )
            transcript = tuple(messages) + (
                {"role": "assistant", "content": "Advisory reviewed."},
            )
            return SimpleNamespace(
                terminal_state="planned",
                final_text="Advisory reviewed.",
                public_transcript=transcript,
                public_transcript_sha256=canonical_sha256(transcript),
                successful_tool_calls=0,
                failed_tool_calls=0,
                event_stream_head_sha256="3" * 64,
            )

    xyz = tmp_path / "water.xyz"
    xyz.write_text(
        "3\nwater\nO 0.0 0.0 0.0\nH 0.0 0.0 1.0\nH 1.0 0.0 0.0\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(
        live, "load_agent_provider_selection", lambda *args, **kwargs: selection
    )
    monkeypatch.setattr(live, "load_program_capabilities", lambda: registry)
    monkeypatch.setattr(live, "build_live_click_schema", lambda: live_schema)
    monkeypatch.setattr(live, "_bootstrap_conformance", lambda **kwargs: ((), ()))
    monkeypatch.setattr(live, "_observe_environments", lambda: ((), (), ()))
    monkeypatch.setattr(
        live, "build_command_compiled_tool_surface", lambda registry: surface
    )
    monkeypatch.setattr(live, "CommandCompiledToolHostV1", _CoordinatorHost)
    monkeypatch.setattr(live, "LiveSpecialistCampaignV1", _Campaign)
    monkeypatch.setattr(live, "UnifiedSessionRunner", _CoordinatorRunner)
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
        provider="qwen-token-plan",
        provider_config_file=tmp_path / "unused.yaml",
        secret_file=tmp_path / "unused.env",
        workspace=tmp_path,
        execution_enabled=False,
        approval_file=None,
        experiment_arm=arm,
        experiment_case=case,
        experiment_repeat_index=0,
    )

    assert result.terminal_state == "planned"
    assert result.experiment_observations["critic"]["status"] == "complete"
    assert result.experiment_observations["feedback_receipts"] == ()
    assert _Campaign.critic_context is not None
    candidate = _Campaign.critic_context["candidate"]
    assert candidate.task_spec_sha256 == result.task_spec_sha256
    serialized = canonical_json(candidate.review_record())
    assert "public_transcript" not in serialized
    assert "final_text" not in serialized
    assert "feedback_projection" not in serialized


def test_feedback_receipt_binding_persists_a_real_tool_receipt():
    projected = project_tool_feedback(
        tool="inspect_program_capability",
        result={
            "schema_version": "chemsmart.agent-tool-result.v1",
            "tool": "inspect_program_capability",
            "status": "ok",
            "result": {"declared": True},
        },
        mode="full-v1",
    )
    event = SimpleNamespace(
        kind="tool_succeeded",
        payload={
            "feedback_equivalence_receipt": asdict(projected.receipt),
        },
    )

    bound = _bind_feedback_receipts(
        observations={
            "schema_version": (
                "chemsmart.live-harness-experiment-observations.v1"
            ),
            "record_sha256": "0" * 64,
        },
        events=(event,),
    )

    assert len(bound["feedback_receipts"]) == 1
    assert bound["feedback_receipts"][0]["receipt_sha256"] == (
        projected.receipt.receipt_sha256
    )
    assert bound["record_sha256"] == canonical_sha256(
        {
            key: value
            for key, value in bound.items()
            if key != "record_sha256"
        }
    )


def test_feedback_receipt_binding_preserves_workers_before_coordinator():
    worker = project_tool_feedback(
        tool="inspect_program_environment",
        result={
            "schema_version": "chemsmart.agent-tool-result.v1",
            "tool": "inspect_program_environment",
            "status": "ok",
            "result": {"observed_status": "available"},
        },
        mode="causal-v1",
    )
    coordinator = project_tool_feedback(
        tool="inspect_program_capability",
        result={
            "schema_version": "chemsmart.agent-tool-result.v1",
            "tool": "inspect_program_capability",
            "status": "ok",
            "result": {"declared": True},
        },
        mode="causal-v1",
    )
    event = SimpleNamespace(
        kind="tool_succeeded",
        payload={
            "feedback_equivalence_receipt": asdict(coordinator.receipt),
        },
    )

    bound = _bind_feedback_receipts(
        observations={
            "schema_version": (
                "chemsmart.live-harness-experiment-observations.v1"
            ),
            "feedback_receipts": (asdict(worker.receipt),),
            "record_sha256": "0" * 64,
        },
        events=(event,),
    )

    assert tuple(
        item["receipt_sha256"] for item in bound["feedback_receipts"]
    ) == (
        worker.receipt.receipt_sha256,
        coordinator.receipt.receipt_sha256,
    )


def test_critic_candidate_digest_is_invariant_to_feedback_projection():
    decision = build_scientific_decision_record(
        decision_id="water.plan",
        task_spec_sha256="d" * 64,
        stage_order=("sp",),
        assumptions=("The supplied geometry is fixed.",),
        method_rationale="Use the explicit host-resolved method.",
    )
    bundle = _empty_host_oracle_bundle(
        session_id="coordinator-session",
        event_stream_head_sha256="f" * 64,
        successful_tool_calls=0,
        failed_tool_calls=0,
    )
    candidates = {
        projection: build_f_invariant_critic_candidate(
            candidate_id="candidate-water",
            task_spec_sha256="d" * 64,
            host_oracle_input_bundle=bundle,
            coordinator_public_decisions=(decision,),
        )
        for projection in ("full-v1", "causal-v1")
    }

    assert candidates["full-v1"].candidate_sha256 == (
        candidates["causal-v1"].candidate_sha256
    )
    record = candidates["full-v1"].review_record()
    assert record["coordinator_public_decisions"][0]["record_sha256"] == (
        decision.record_sha256
    )
    assert "feedback_projection" not in canonical_json(record)
