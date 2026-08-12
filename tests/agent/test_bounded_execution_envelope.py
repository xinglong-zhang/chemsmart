from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

from click.testing import CliRunner
import pytest
import yaml

from chemsmart.agent._contracts import (
    ContractError,
    TrustedArtifactRefV1,
    file_sha256,
)
from chemsmart.agent.capabilities import (
    SupportLevel,
    build_program_component_conformance_receipt,
)
from chemsmart.agent.execution import (
    build_producer_edge_rule,
    handoff_optimized_native_geometry,
    handoff_optimized_xtb_geometry,
)
from chemsmart.agent.execution_envelope import (
    BoundedExecutionEnvelopeV1,
    load_bounded_execution_envelope,
)
from chemsmart.agent.runtime.event_store import RuntimeEventStore
from chemsmart.agent.tool_runtime import (
    CommandCompiledToolHostV1,
    _CommandContext,
    _program_process_environment,
)
from chemsmart.agent.tool_specs import build_approved_execution_tool_surface
from chemsmart.agent.workflows import (
    MaterializedNodeV1,
    ScientificWorkflowEdgeV2,
    ScientificWorkflowNodeV2,
    build_materialized_workflow,
    build_scientific_workflow_plan,
)
from chemsmart.cli.agent import agent


def _payload(tmp_path: Path) -> dict:
    return {
        "schema_version": "chemsmart.bounded-execution-envelope.v1",
        "mode": "continuous-local",
        "allowed_program_engines": {"pyscf": ["cpu"], "xtb": ["cpu"]},
        "resources": {
            "execution_target": "run",
            "cores": 8,
            "memory_gb": 28,
            "gpu_count": 0,
            "scratch_policy": "server",
            "node_timeout_seconds": 600,
        },
        "episode_wall_time_seconds": 3600,
        "postprocess_reserve_seconds": 300,
        "max_engine_calls": 4,
        "scratch_root": str(tmp_path / "scratch"),
    }


def _write_envelope(tmp_path: Path, payload: dict | None = None) -> Path:
    path = tmp_path / "envelope.yaml"
    path.write_text(
        yaml.safe_dump(payload or _payload(tmp_path), sort_keys=True),
        encoding="utf-8",
    )
    return path


def _artifact(
    path: Path, *, artifact_id: str, kind: str
) -> TrustedArtifactRefV1:
    return TrustedArtifactRefV1(
        artifact_id=artifact_id,
        kind=kind,
        sha256=file_sha256(path),
        size_bytes=path.stat().st_size,
        path=str(path.resolve()),
        cli_value=str(path.resolve()),
    )


def test_loader_accepts_science_free_resource_and_program_bounds(tmp_path):
    observed = load_bounded_execution_envelope(_write_envelope(tmp_path))

    assert observed.allows("pyscf", "cpu") is True
    assert observed.allows("orca", "cpu") is False
    assert observed.resources.cores == 8
    assert observed.resources.memory_gb == 28
    assert observed.max_engine_calls == 4


@pytest.mark.parametrize(
    ("mutate", "message"),
    (
        (
            lambda value: value.update(
                {"workflow": {"nodes": ["answer-specific"]}}
            ),
            "unsupported fields: workflow",
        ),
        (
            lambda value: value["resources"].update(
                {"execution_target": "sub"}
            ),
            "requires target run",
        ),
        (
            lambda value: value["resources"].update(
                {"scratch_policy": "none"}
            ),
            "requires the server scratch policy",
        ),
        (
            lambda value: value.update({"scratch_root": "/"}),
            "narrow absolute path",
        ),
        (
            lambda value: value["resources"].update({"gpu_count": 1}),
            "GPU resources require",
        ),
        (
            lambda value: value["allowed_program_engines"].update(
                {"pyscf": ["gpu"]}
            ),
            "non-CPU engine requires",
        ),
        (
            lambda value: value.update({"max_engine_calls": 0}),
            "max_engine_calls must be positive",
        ),
    ),
)
def test_loader_rejects_malformed_or_unsafe_envelopes(
    tmp_path, mutate, message
):
    payload = _payload(tmp_path)
    mutate(payload)

    with pytest.raises(ContractError, match=message):
        load_bounded_execution_envelope(_write_envelope(tmp_path, payload))


class _Result:
    def public_summary_json(self):
        return "{}"


def _cli_args(tmp_path: Path) -> list[str]:
    workspace = tmp_path / "workspace"
    workspace.mkdir()
    secret = tmp_path / "secret.env"
    secret.write_text("KEY=value\n", encoding="utf-8")
    return [
        "run",
        "--task",
        "preview this calculation",
        "--secret-file",
        str(secret),
        "--workspace",
        str(workspace),
    ]


def test_agent_run_without_authorization_remains_preview_only(
    tmp_path, monkeypatch
):
    captured = {}

    def fake_run(**kwargs):
        captured.update(kwargs)
        return _Result()

    import chemsmart.agent.live_session as live_session

    monkeypatch.setattr(live_session, "run_live_agent_session", fake_run)
    result = CliRunner().invoke(agent, _cli_args(tmp_path))

    assert result.exit_code == 0, result.output
    assert captured["execution_enabled"] is False
    assert captured["approval_file"] is None
    assert captured["execution_envelope_file"] is None


def test_agent_run_rejects_approval_and_envelope_together(tmp_path):
    approval = tmp_path / "approval.json"
    approval.write_text("{}", encoding="utf-8")
    args = _cli_args(tmp_path) + [
        "--approval-file",
        str(approval),
        "--execution-envelope",
        str(_write_envelope(tmp_path)),
    ]

    result = CliRunner().invoke(agent, args)

    assert result.exit_code == 2
    assert "mutually exclusive" in result.output


def test_existing_approval_path_is_forwarded_unchanged(tmp_path, monkeypatch):
    captured = {}
    approval = tmp_path / "approval.json"
    approval.write_text("{}", encoding="utf-8")

    def fake_run(**kwargs):
        captured.update(kwargs)
        return _Result()

    import chemsmart.agent.live_session as live_session

    monkeypatch.setattr(live_session, "run_live_agent_session", fake_run)
    result = CliRunner().invoke(
        agent,
        _cli_args(tmp_path) + ["--approval-file", str(approval)],
    )

    assert result.exit_code == 0, result.output
    assert captured["execution_enabled"] is True
    assert captured["approval_file"] == approval
    assert captured["execution_envelope_file"] is None


def test_envelope_path_is_forwarded_as_bounded_execution(
    tmp_path, monkeypatch
):
    captured = {}
    envelope = _write_envelope(tmp_path)

    def fake_run(**kwargs):
        captured.update(kwargs)
        return _Result()

    import chemsmart.agent.live_session as live_session

    monkeypatch.setattr(live_session, "run_live_agent_session", fake_run)
    result = CliRunner().invoke(
        agent,
        _cli_args(tmp_path) + ["--execution-envelope", str(envelope)],
    )

    assert result.exit_code == 0, result.output
    assert captured["execution_enabled"] is True
    assert captured["approval_file"] is None
    assert captured["execution_envelope_file"] == envelope


def test_deferred_admission_builds_existing_approval_contracts(tmp_path):
    envelope = load_bounded_execution_envelope(_write_envelope(tmp_path))
    plan = build_scientific_workflow_plan(
        workflow_id="water-sp",
        task_spec_sha256="a" * 64,
        scientific_identity_sha256="b" * 64,
        nodes=(
            ScientificWorkflowNodeV2(
                node_id="sp-initial",
                stage="sp",
                requested_program="pyscf",
                program="pyscf",
                engine="cpu",
                project_role="project.water",
                unresolved_fields=(),
            ),
        ),
        edges=(),
    )
    materialized = build_materialized_workflow(
        plan=plan,
        live_cli_schema_sha256="c" * 64,
        resource_sha256=envelope.resources.resource_sha256,
        nodes=(
            MaterializedNodeV1(
                node_id="sp-initial",
                input_artifact_sha256="d" * 64,
                project_artifact_sha256="e" * 64,
                project_validation_receipt_sha256="f" * 64,
                environment_receipt_sha256="1" * 64,
                invocation_sha256="2" * 64,
                preflight_receipt_sha256="3" * 64,
                invocation_identity_sha256="4" * 64,
                state="previewed",
            ),
        ),
        unresolved_node_ids=(),
        status="previewed",
    )
    project = SimpleNamespace(artifact_id="project.water", sha256="e" * 64)
    geometry = SimpleNamespace(artifact_id="geometry.water", sha256="d" * 64)
    identity = SimpleNamespace(
        binding_sha256="b" * 64, charge=0, multiplicity=1
    )
    context = _CommandContext(
        proposal=SimpleNamespace(program="pyscf", jobtype="sp"),
        capability=SimpleNamespace(),
        program_binding=SimpleNamespace(),
        engine_binding=SimpleNamespace(
            engine="cpu", environment_receipt_sha256="1" * 64
        ),
        project_artifact=project,
        project_validation=SimpleNamespace(
            status="valid", settings_sha256="5" * 64
        ),
        input_artifact=geometry,
        scientific_identity=identity,
    )
    host = object.__new__(CommandCompiledToolHostV1)
    host.bounded_execution_envelope = envelope
    host.execution_resources = envelope.resources
    host.execution_receipts = {}
    host._bounded_execution_started_at = 0.0
    host.scientific_workflow_plans = {plan.plan_sha256: plan}
    host.task_spec_sha256s = frozenset({plan.task_spec_sha256})
    host.materialized_workflows = {
        materialized.materialized_sha256: materialized
    }
    host.approved_workspace = tmp_path.resolve()
    host.stationary_point_policy = None
    host.workflow_drafts = {}
    host.workflow_execution_approval = None
    host.frozen_workflow_approval = None
    host._latest_invocation_for_node = lambda _node_id: (
        SimpleNamespace(auxiliary_input_bindings=()),
        context,
    )
    host._environment_identity_for = lambda _receipt: "6" * 64
    host._require_bounded_launch_budget = lambda: None

    host._admit_bounded_workflow(node_id="sp-initial")

    assert host.workflow_execution_approval.workflow_sha256 == plan.plan_sha256
    assert host.workflow_execution_approval.node(
        "sp-initial"
    ).settings_sha256 == ("5" * 64)
    assert host.frozen_workflow_approval.plan_sha256 == plan.plan_sha256
    assert host.frozen_workflow_approval.execution_admissible is True


def test_deferred_admission_enforces_program_allowlist(tmp_path):
    envelope = load_bounded_execution_envelope(_write_envelope(tmp_path))
    assert isinstance(envelope, BoundedExecutionEnvelopeV1)
    assert envelope.allows("gaussian", "cpu") is False


def test_bounded_overlay_keeps_allowed_red_program_reference_only(
    tmp_path,
    fake_capability_registry,
    fake_click_schema,
):
    payload = _payload(tmp_path)
    payload["allowed_program_engines"] = {
        "demo": ["cpu"],
        "optional": ["cpu"],
    }
    envelope = load_bounded_execution_envelope(
        _write_envelope(tmp_path, payload)
    )
    conformance = build_program_component_conformance_receipt(
        program="demo",
        registry_sha256=fake_capability_registry.registry_sha256,
        live_cli_schema_sha256=fake_click_schema.schema_sha256,
        fixture_bundle_sha256="1" * 64,
        covered_jobtypes=("sp",),
        covered_engines=("cpu",),
        covered_engine_job_pairs=(("cpu", "sp"),),
        compiler_receipt_sha256="2" * 64,
        preview_receipt_sha256="3" * 64,
        preflight_receipt_sha256="4" * 64,
        verifier_receipt_sha256="5" * 64,
        compiler_status="passed",
        preview_status="passed",
        preflight_status="passed",
        verifier_status="passed",
    )

    host = CommandCompiledToolHostV1(
        event_store=RuntimeEventStore(
            tmp_path / "events.jsonl", session_id="bounded-overlay"
        ),
        component_conformance_receipts=(conformance,),
        tool_surface=build_approved_execution_tool_surface(
            fake_capability_registry
        ),
        registry=fake_capability_registry,
        live_schema=fake_click_schema,
        approved_workspace=tmp_path,
        execution_resources=envelope.resources,
        bounded_execution_envelope=envelope,
    )

    assert host.overlay.get("demo").support_level is SupportLevel.AVAILABLE
    assert (
        host.overlay.get("optional").support_level
        is SupportLevel.REFERENCE_ONLY
    )


def test_planning_time_reduces_effective_node_timeout(tmp_path, monkeypatch):
    envelope = load_bounded_execution_envelope(_write_envelope(tmp_path))
    host = object.__new__(CommandCompiledToolHostV1)
    host.bounded_execution_envelope = envelope
    host.execution_resources = envelope.resources
    host.execution_receipts = {}
    host._bounded_execution_started_at = 1000.0
    monkeypatch.setattr(
        "chemsmart.agent.tool_runtime.time.monotonic", lambda: 3800.0
    )

    observed = host._require_bounded_launch_budget()

    # 3600 episode - 2800 planning - 300 analysis reserve.
    assert observed == pytest.approx(500.0)


def test_completed_node_can_be_followed_by_shorter_bounded_launch(
    tmp_path, monkeypatch
):
    envelope = load_bounded_execution_envelope(_write_envelope(tmp_path))
    host = object.__new__(CommandCompiledToolHostV1)
    host.bounded_execution_envelope = envelope
    host.execution_resources = envelope.resources
    host.execution_receipts = {"opt-initial": SimpleNamespace()}
    host._bounded_execution_started_at = 1000.0
    now = iter((3900.0, 4299.0))
    monkeypatch.setattr(
        "chemsmart.agent.tool_runtime.time.monotonic", lambda: next(now)
    )

    # One completed engine call does not force the next node to request the
    # original 600 seconds: 3600 - 2900 - 300 leaves a 400-second window.
    assert host._require_bounded_launch_budget() == pytest.approx(400.0)
    # Later, the same episode can still launch a node with the 1-second usable
    # window while retaining the full 300-second analysis reserve.
    assert host._require_bounded_launch_budget() == pytest.approx(1.0)


def test_bounded_launch_refuses_to_consume_postprocessing_reserve(
    tmp_path, monkeypatch
):
    envelope = load_bounded_execution_envelope(_write_envelope(tmp_path))
    host = object.__new__(CommandCompiledToolHostV1)
    host.bounded_execution_envelope = envelope
    host.execution_resources = envelope.resources
    host.execution_receipts = {"opt-initial": SimpleNamespace()}
    host._bounded_execution_started_at = 1000.0
    monkeypatch.setattr(
        "chemsmart.agent.tool_runtime.time.monotonic", lambda: 4300.5
    )

    with pytest.raises(ContractError, match="postprocessing reserve"):
        host._require_bounded_launch_budget()


def test_session_rejects_allowlist_outside_live_registry(tmp_path):
    from chemsmart.agent.live_session import (
        _validate_bounded_envelope_against_registry,
    )

    payload = _payload(tmp_path)
    payload["allowed_program_engines"] = {"invented": ["cpu"]}
    envelope = load_bounded_execution_envelope(
        _write_envelope(tmp_path, payload)
    )
    registry = SimpleNamespace(get=lambda _program: None)

    with pytest.raises(ContractError, match="unknown program 'invented'"):
        _validate_bounded_envelope_against_registry(
            envelope, registry=registry
        )


def test_session_rejects_program_without_continuous_result_validator(tmp_path):
    from chemsmart.agent.live_session import (
        _validate_bounded_envelope_against_registry,
    )

    payload = _payload(tmp_path)
    payload["allowed_program_engines"] = {"nciplot": ["cpu"]}
    envelope = load_bounded_execution_envelope(
        _write_envelope(tmp_path, payload)
    )
    capability = SimpleNamespace(
        engines=("cpu",), execution_engine_job_pairs=()
    )
    registry = SimpleNamespace(get=lambda _program: capability)

    with pytest.raises(ContractError, match="lacks a result validator"):
        _validate_bounded_envelope_against_registry(
            envelope, registry=registry
        )


def test_program_environment_excludes_provider_secret(monkeypatch):
    monkeypatch.setenv("LAB_DEEPSEEK_API_KEY", "must-not-reach-engine")
    monkeypatch.setenv("GAUSSIAN_LICENSE_FILE", "/licenses/g16")

    observed = _program_process_environment(
        overrides={"OMP_NUM_THREADS": "8"},
        remove=("LAB_DEEPSEEK_API_KEY",),
    )

    assert "LAB_DEEPSEEK_API_KEY" not in observed
    assert observed["GAUSSIAN_LICENSE_FILE"] == "/licenses/g16"
    assert observed["OMP_NUM_THREADS"] == "8"


@pytest.mark.parametrize(
    ("program", "consumer_stage"),
    (("xtb", "hess"), ("gaussian", "sp")),
)
def test_future_node_context_joins_same_engine_stages_by_capability(
    program, consumer_stage
):
    identity = SimpleNamespace(
        binding_sha256="a" * 64, charge=-1, multiplicity=2
    )
    geometry = SimpleNamespace(artifact_id="geometry.initial", sha256="b" * 64)
    producer_context = _CommandContext(
        proposal=SimpleNamespace(program="pyscf", jobtype="opt"),
        capability=SimpleNamespace(),
        program_binding=SimpleNamespace(),
        engine_binding=SimpleNamespace(
            engine="cpu", environment_receipt_sha256="c" * 64
        ),
        project_artifact=SimpleNamespace(
            artifact_id="project.opt", sha256="d" * 64
        ),
        project_validation=SimpleNamespace(
            status="valid", settings_sha256="e" * 64
        ),
        input_artifact=geometry,
        scientific_identity=identity,
    )
    plan = build_scientific_workflow_plan(
        workflow_id=f"radical-opt-{consumer_stage}",
        task_spec_sha256="f" * 64,
        scientific_identity_sha256=identity.binding_sha256,
        nodes=(
            ScientificWorkflowNodeV2(
                node_id="opt-initial",
                stage="opt",
                requested_program=program,
                program=program,
                engine="cpu",
                project_role="project.opt",
                unresolved_fields=(),
            ),
            ScientificWorkflowNodeV2(
                node_id=f"{consumer_stage}-optimized",
                stage=consumer_stage,
                requested_program=program,
                program=program,
                engine="cpu",
                project_role=f"project.{consumer_stage}",
                unresolved_fields=(),
                support_state="unresolved_future",
            ),
        ),
        edges=(
            ScientificWorkflowEdgeV2(
                edge_id=f"opt-to-{consumer_stage}",
                source_node_id="opt-initial",
                target_node_id=f"{consumer_stage}-optimized",
                edge_kind="data",
                artifact_class="geometry_xyz",
                producer_output_id="optimized-geometry",
                consumer_input_id="geometry",
            ),
        ),
    )
    project = SimpleNamespace(
        artifact_id=f"project.{consumer_stage}", sha256="1" * 64
    )
    producer_capability = SimpleNamespace(
        receipt_sha256="2" * 64,
        status=SimpleNamespace(value="supported"),
        query=SimpleNamespace(program=program, jobtype="opt", engine="cpu"),
    )
    consumer_capability = SimpleNamespace(
        receipt_sha256="3" * 64,
        status=SimpleNamespace(value="supported"),
        query=SimpleNamespace(
            program=program, jobtype=consumer_stage, engine="cpu"
        ),
    )
    alternate_capability = SimpleNamespace(
        receipt_sha256="4" * 64,
        status=SimpleNamespace(value="supported"),
        query=SimpleNamespace(
            program=program, jobtype=consumer_stage, engine="gpu"
        ),
    )
    validation = SimpleNamespace(
        project_artifact_id=project.artifact_id,
        project_sha256=project.sha256,
        capability_receipt_sha256=consumer_capability.receipt_sha256,
        program=program,
        jobtype=consumer_stage,
        status="valid",
    )
    alternate_validation = SimpleNamespace(
        project_artifact_id=project.artifact_id,
        project_sha256=project.sha256,
        capability_receipt_sha256=alternate_capability.receipt_sha256,
        program=program,
        jobtype=consumer_stage,
        status="valid",
    )
    producer_engine = SimpleNamespace(
        program=program,
        selected_engine="cpu",
        execution_ready=True,
        capability_receipt_sha256=producer_capability.receipt_sha256,
        program_binding_sha256="5" * 64,
        environment_receipt_sha256="6" * 64,
    )
    consumer_engine = SimpleNamespace(
        program=program,
        selected_engine="cpu",
        execution_ready=True,
        capability_receipt_sha256=consumer_capability.receipt_sha256,
        program_binding_sha256="7" * 64,
        environment_receipt_sha256="8" * 64,
    )
    alternate_engine = SimpleNamespace(
        program=program,
        selected_engine="gpu",
        execution_ready=True,
        capability_receipt_sha256=alternate_capability.receipt_sha256,
        program_binding_sha256="9" * 64,
        environment_receipt_sha256="0" * 64,
    )
    host = object.__new__(CommandCompiledToolHostV1)
    host.workflow_drafts = {
        "draft": SimpleNamespace(
            workflow_id=plan.workflow_id,
            nodes=(
                SimpleNamespace(
                    node_id=f"{consumer_stage}-optimized",
                    project_role=project.artifact_id,
                    inputs=(SimpleNamespace(producer_node_id="opt-initial"),),
                ),
            ),
        )
    }
    host.artifacts = {project.artifact_id: project}
    host.project_validations = {
        "validation": validation,
        "alternate-validation": alternate_validation,
    }
    host.engine_bindings = {
        "producer": producer_engine,
        "consumer": consumer_engine,
        "alternate": alternate_engine,
    }
    host.capabilities = {
        capability.receipt_sha256: capability
        for capability in (
            producer_capability,
            consumer_capability,
            alternate_capability,
        )
    }
    host.program_bindings = {
        binding.program_binding_sha256: SimpleNamespace()
        for binding in (producer_engine, consumer_engine, alternate_engine)
    }
    host._latest_invocation_for_node = lambda _node_id: (
        SimpleNamespace(),
        producer_context,
    )

    context = host._bounded_node_context(
        plan=plan,
        planned_node=plan.nodes[1],
        data_target_ids={f"{consumer_stage}-optimized"},
    )

    assert context.scientific_identity is identity
    assert context.proposal.charge == -1
    assert context.proposal.multiplicity == 2
    assert (
        context.proposal.scientific_identity_sha256 == identity.binding_sha256
    )
    assert context.capability is consumer_capability
    assert context.engine_binding is consumer_engine
    assert context.project_validation is validation


def test_bounded_readiness_defers_exact_optimized_geometry_consumer(tmp_path):
    envelope = load_bounded_execution_envelope(_write_envelope(tmp_path))
    plan = build_scientific_workflow_plan(
        workflow_id="water-opt-hess",
        task_spec_sha256="a" * 64,
        scientific_identity_sha256="b" * 64,
        nodes=(
            ScientificWorkflowNodeV2(
                node_id="opt-initial",
                stage="opt",
                requested_program="xtb",
                program="xtb",
                engine="cpu",
                project_role="project.opt",
                unresolved_fields=(),
            ),
            ScientificWorkflowNodeV2(
                node_id="hess-optimized",
                stage="hess",
                requested_program="xtb",
                program="xtb",
                engine="cpu",
                project_role="project.hess",
                unresolved_fields=(),
                support_state="unresolved_future",
            ),
        ),
        edges=(
            ScientificWorkflowEdgeV2(
                edge_id="opt-to-hess",
                source_node_id="opt-initial",
                target_node_id="hess-optimized",
                edge_kind="data",
                artifact_class="geometry_xyz",
                producer_output_id="optimized-geometry",
                consumer_input_id="geometry",
            ),
        ),
    )
    host = object.__new__(CommandCompiledToolHostV1)
    host.bounded_execution_envelope = envelope
    host._preflight_by_node = {
        "opt-initial": SimpleNamespace(
            plan_state="previewed", critical_finding_sha256s=()
        )
    }

    readiness = host._approval_readiness(plan)

    assert readiness["approvable"] is True
    assert readiness["blocking_node_ids"] == ()
    assert readiness["deferred_node_ids"] == ("hess-optimized",)
    states = {node["node_id"]: node for node in readiness["nodes"]}
    assert states["opt-initial"]["approval_state"] == "previewed"
    assert states["hess-optimized"]["approval_state"] == "deferred_admissible"
    assert states["hess-optimized"]["blocks_approval"] is False

    host.bounded_execution_envelope = None
    exact_readiness = host._approval_readiness(plan)

    assert exact_readiness["approvable"] is False
    assert exact_readiness["blocking_node_ids"] == ("hess-optimized",)
    assert exact_readiness["deferred_node_ids"] == ()


def test_bounded_server_profile_preserves_program_environment_and_timeout(
    tmp_path, monkeypatch
):
    import chemsmart.agent.live_session as live_session

    monkeypatch.setattr(
        live_session,
        "_active_server_program_blocks",
        lambda: {
            "ORCA": {
                "EXEFOLDER": "/opt/orca",
                "LOCAL_RUN": False,
                "ENVARS": (
                    "export PATH=/opt/openmpi/bin:$PATH\n"
                    "export LD_LIBRARY_PATH=/opt/openmpi/lib:$LD_LIBRARY_PATH\n"
                ),
            },
            "GAUSSIAN": {
                "EXEFOLDER": "/opt/g16",
                "LOCAL_RUN": True,
                "ENVARS": "export g16root=/opt/g16\n",
            },
        },
    )
    payload = _payload(tmp_path)
    payload["resources"]["node_timeout_seconds"] = 7201
    payload["episode_wall_time_seconds"] = 12000
    envelope = load_bounded_execution_envelope(
        _write_envelope(tmp_path, payload)
    )
    run_directory = tmp_path / "run"
    run_directory.mkdir()

    profile = live_session._write_execution_server_profile(
        run_directory,
        envelope.resources,
        scratch_root=tmp_path / "scratch",
    )
    observed = yaml.safe_load(profile.read_text(encoding="utf-8"))

    assert observed["SERVER"]["NUM_HOURS"] == 3
    assert "openmpi/bin" in observed["ORCA"]["ENVARS"]
    assert "LD_LIBRARY_PATH" in observed["ORCA"]["ENVARS"]
    assert "export g16root=/opt/g16" in observed["GAUSSIAN"]["ENVARS"]
    assert (
        f"export GAUSS_SCRDIR={tmp_path / 'scratch'}"
        in (observed["GAUSSIAN"]["ENVARS"])
    )


def test_xtb_optimized_geometry_handoff_preserves_identity_and_state(tmp_path):
    input_path = tmp_path / "input.xyz"
    input_path.write_text(
        "3\nwater radical cation\n"
        "O 0.0 0.0 0.1\nH -0.7 0.0 -0.4\nH 0.7 0.0 -0.4\n",
        encoding="utf-8",
    )
    result_path = tmp_path / "xtbopt.xyz"
    result_path.write_text(
        "3\noptimized\nO 0.0 0.0 0.12\nH -0.76 0.0 -0.47\nH 0.76 0.0 -0.47\n",
        encoding="utf-8",
    )
    input_artifact = _artifact(
        input_path, artifact_id="geometry.initial", kind="geometry_xyz"
    )
    result_artifact = _artifact(
        result_path, artifact_id="result.opt.geometry", kind="geometry_xyz"
    )
    edge = build_producer_edge_rule(
        producer_node_id="opt-initial",
        consumer_node_id="hess-optimized",
        artifact_kind="geometry_xyz",
        selection_rule="validated_optimized_geometry",
    )
    receipt = SimpleNamespace(
        validated=True,
        node_id="opt-initial",
        receipt_sha256="a" * 64,
        output_artifacts=(result_artifact,),
    )

    geometry, handoff = handoff_optimized_xtb_geometry(
        producer_receipt=receipt,
        result_artifact=result_artifact,
        input_artifact=input_artifact,
        producer_edge=edge,
        approved_workspace=tmp_path,
        geometry_artifact_id="geometry.xtb.optimized",
        expected_charge=1,
        expected_multiplicity=2,
    )

    assert handoff.symbols == ("O", "H", "H")
    assert (handoff.charge, handoff.multiplicity) == (1, 2)
    assert geometry.sha256 == file_sha256(geometry.path)
    assert "charge=1; multiplicity=2" in Path(geometry.path).read_text(
        encoding="utf-8"
    )


def test_xtb_optimized_geometry_handoff_rejects_atom_reordering(tmp_path):
    input_path = tmp_path / "input.xyz"
    input_path.write_text(
        "3\nwater\nO 0 0 0.1\nH -0.7 0 -0.4\nH 0.7 0 -0.4\n",
        encoding="utf-8",
    )
    result_path = tmp_path / "xtbopt.xyz"
    result_path.write_text(
        "3\nreordered\nH -0.7 0 -0.4\nO 0 0 0.1\nH 0.7 0 -0.4\n",
        encoding="utf-8",
    )
    input_artifact = _artifact(
        input_path, artifact_id="geometry.initial", kind="geometry_xyz"
    )
    result_artifact = _artifact(
        result_path, artifact_id="result.opt.geometry", kind="geometry_xyz"
    )
    edge = build_producer_edge_rule(
        producer_node_id="opt-initial",
        consumer_node_id="hess-optimized",
        artifact_kind="geometry_xyz",
        selection_rule="validated_optimized_geometry",
    )
    receipt = SimpleNamespace(
        validated=True,
        node_id="opt-initial",
        receipt_sha256="a" * 64,
        output_artifacts=(result_artifact,),
    )

    with pytest.raises(ContractError, match="atom identity or atom order"):
        handoff_optimized_xtb_geometry(
            producer_receipt=receipt,
            result_artifact=result_artifact,
            input_artifact=input_artifact,
            producer_edge=edge,
            approved_workspace=tmp_path,
            geometry_artifact_id="geometry.xtb.optimized",
            expected_charge=0,
            expected_multiplicity=1,
        )


@pytest.mark.parametrize(
    (
        "program",
        "relative_output",
        "expected_symbols",
        "charge",
        "multiplicity",
    ),
    (
        (
            "gaussian",
            "tests/data/GaussianTests/outputs/acetone.log",
            ("O", "C", "C", "C", "H", "H", "H", "H", "H", "H"),
            0,
            1,
        ),
        (
            "orca",
            "tests/data/ORCATests/outputs/sn2_ts.out",
            ("C", "Cl", "H", "H", "H", "F"),
            -1,
            1,
        ),
    ),
)
def test_native_optimized_geometry_handoff_uses_parser_stationary_geometry(
    tmp_path,
    program,
    relative_output,
    expected_symbols,
    charge,
    multiplicity,
):
    output_path = Path(relative_output).resolve()
    if program == "gaussian":
        from chemsmart.io.gaussian.output import Gaussian16Output

        molecule = Gaussian16Output(str(output_path)).optimized_structure
    else:
        from chemsmart.io.orca.output import ORCAOutput

        molecule = ORCAOutput(str(output_path)).molecule
    input_path = tmp_path / f"{program}-input.xyz"
    lines = [str(len(expected_symbols)), "same molecule before optimization"]
    for symbol, position in zip(expected_symbols, molecule.positions):
        x, y, z = (float(value) + 0.01 for value in position)
        lines.append(f"{symbol} {x} {y} {z}")
    input_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    input_artifact = _artifact(
        input_path,
        artifact_id=f"geometry.{program}.initial",
        kind="geometry_xyz",
    )
    result_artifact = _artifact(
        output_path,
        artifact_id=f"result.{program}.optimized",
        kind=f"{program}_output",
    )
    edge = build_producer_edge_rule(
        producer_node_id="opt-initial",
        consumer_node_id="sp-optimized",
        artifact_kind="geometry_xyz",
        selection_rule="validated_optimized_geometry",
    )
    receipt = SimpleNamespace(
        validated=True,
        node_id="opt-initial",
        receipt_sha256="a" * 64,
        output_artifacts=(result_artifact,),
    )

    geometry, handoff = handoff_optimized_native_geometry(
        program=program,
        producer_receipt=receipt,
        result_artifact=result_artifact,
        input_artifact=input_artifact,
        producer_edge=edge,
        approved_workspace=tmp_path,
        geometry_artifact_id=f"geometry.{program}.optimized",
        expected_charge=charge,
        expected_multiplicity=multiplicity,
    )

    assert handoff.symbols == expected_symbols
    assert (handoff.charge, handoff.multiplicity) == (charge, multiplicity)
    assert geometry.sha256 == file_sha256(geometry.path)


def test_gaussian_execution_validator_accepts_real_optimum():
    opt_path = Path("tests/data/GaussianTests/outputs/acetone.log").resolve()
    opt = _artifact(
        opt_path, artifact_id="result.gaussian.opt", kind="gaussian_output"
    )

    evaluation = CommandCompiledToolHostV1._evaluate_execution_outputs(
        program="gaussian",
        jobtype="opt",
        charge=0,
        multiplicity=1,
        expected_settings={
            "functional": "M06-2X",
            "basis": "def2-SVP",
            "freq": False,
        },
        output_artifacts=(opt,),
        exit_status=0,
    )

    assert evaluation.validated is True
    assert evaluation.validator_id == "gaussian-result-validator"
    row = evaluation.observations["gaussian"]["outputs"][0]
    assert row["normal_termination"] is True
    assert row["optimization_converged"] is True
    assert row["energy_hartree"] == pytest.approx(-192.919415566)


def test_gaussian_execution_validator_requires_requested_td_roots():
    td_path = Path(
        "tests/data/GaussianTests/tddft/tddft_r1s50_gas_radical_anion.log"
    ).resolve()
    td = _artifact(
        td_path, artifact_id="result.gaussian.td", kind="gaussian_output"
    )

    evaluation = CommandCompiledToolHostV1._evaluate_execution_outputs(
        program="gaussian",
        jobtype="td",
        charge=-1,
        multiplicity=2,
        expected_settings={
            "functional": "CAM-B3LYP",
            "basis": "gen",
            "nstates": 50,
        },
        output_artifacts=(td,),
        exit_status=0,
    )
    too_many_roots = CommandCompiledToolHostV1._evaluate_execution_outputs(
        program="gaussian",
        jobtype="td",
        charge=-1,
        multiplicity=2,
        expected_settings={
            "functional": "CAM-B3LYP",
            "basis": "gen",
            "nstates": 51,
        },
        output_artifacts=(td,),
        exit_status=0,
    )

    assert evaluation.validated is True
    assert (
        evaluation.observations["gaussian"]["outputs"][0]["transition_count"]
        == 50
    )
    assert "gaussian.result.transitions_missing" in too_many_roots.findings


def test_gaussian_link_validator_requires_stable_link_and_target_result():
    link_path = Path(
        "tests/data/GaussianTests/outputs/link/"
        "dppeFeCl2_opt_quintet_link_opt_link.log"
    ).resolve()
    artifact = _artifact(
        link_path,
        artifact_id="result.gaussian.stable-link",
        kind="gaussian_output",
    )

    evaluation = CommandCompiledToolHostV1._evaluate_execution_outputs(
        program="gaussian",
        jobtype="link",
        charge=0,
        multiplicity=5,
        expected_settings={
            "functional": "MN15",
            "basis": "def2-SVP",
            "jobtype": "opt",
            "freq": True,
            "stable": "opt",
        },
        output_artifacts=(artifact,),
        exit_status=0,
    )

    assert evaluation.validated is True
    row = evaluation.observations["gaussian"]["outputs"][0]
    assert row["optimization_converged"] is True
    assert row["wavefunction_stability_history"][-1] == (
        "stable_under_considered_perturbations"
    )
