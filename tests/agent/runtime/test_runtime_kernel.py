from __future__ import annotations

import json
from dataclasses import dataclass

import pytest

from chemsmart.agent.harness.workflow_state import (
    current_workflow_state,
    hydrate_workflow_state,
    reset_workflow_state,
    workflow_state_scope,
)
from chemsmart.agent.runtime.contracts import (
    AgentAction,
    AgentDecision,
    ProviderRole,
    RuntimeV2Mode,
    TaskPhase,
)
from chemsmart.agent.runtime.event_store import (
    EventStoreCorruptionError,
    RuntimeEventStore,
)
from chemsmart.agent.runtime.events import EventKind
from chemsmart.agent.runtime.lifecycle import ToolExposureViolation
from chemsmart.agent.runtime.orchestrator import (
    RuntimeController,
    provider_role,
    route_initial_phase,
)
from chemsmart.agent.runtime.reducer import reduce_events
from chemsmart.agent.runtime.repair_policy import RepairAction, decide_repair
from chemsmart.agent.runtime.tool_catalog import ToolCatalog


@dataclass(frozen=True)
class _Tool:
    name: str

    def openai_tool_def(self):
        return {
            "type": "function",
            "function": {"name": self.name, "parameters": {"type": "object"}},
        }


class _Registry:
    def __init__(self):
        names = {
            "extract_project_protocol",
            "render_project_yaml",
            "validate_project_yaml",
            "critic_project_yaml",
            "write_project_yaml",
            "read_project_yaml",
            "update_project_yaml",
            "search_basis_sets",
            "synthesize_command",
            "repair_command",
            "execute_chemsmart_command",
            "build_job",
            "dry_run_input",
        }
        self._tools = {name: _Tool(name) for name in names}

    def list_tools(self):
        return list(self._tools.values())

    def get_tool(self, name):
        return self._tools.get(name)

    def tool_defs_for_provider(self, provider_name, tools):
        del provider_name
        return [tool.openai_tool_def() for tool in tools]


def test_event_store_is_idempotent_and_replayable(tmp_path):
    store = RuntimeEventStore(tmp_path / "events.jsonl")
    first = store.append(
        session_id="s1",
        turn_id="bootstrap",
        kind=EventKind.SESSION_STARTED,
        payload={"cwd": str(tmp_path)},
        idempotency_key="session-start",
    )
    duplicate = store.append(
        session_id="s1",
        turn_id="bootstrap",
        kind=EventKind.SESSION_STARTED,
        payload={"cwd": "ignored"},
        idempotency_key="session-start",
    )
    store.append(
        session_id="s1",
        turn_id="turn_0001",
        kind=EventKind.TURN_STARTED,
        payload={"request": "optimize", "phase": "synthesis"},
    )

    events = store.load()
    state = reduce_events(events)

    assert duplicate.event_id == first.event_id
    assert len(events) == 2
    assert state.session_id == "s1"
    assert state.phase is TaskPhase.SYNTHESIS


def test_event_store_detects_payload_tampering(tmp_path):
    path = tmp_path / "events.jsonl"
    store = RuntimeEventStore(path)
    store.append(
        session_id="s1",
        turn_id="bootstrap",
        kind=EventKind.SESSION_STARTED,
        payload={"cwd": str(tmp_path)},
    )
    row = json.loads(path.read_text())
    row["payload"]["cwd"] = "/tampered"
    path.write_text(json.dumps(row) + "\n")

    with pytest.raises(EventStoreCorruptionError, match="invalid hash"):
        store.load()


def test_replay_preserves_inflight_tool_after_crash(tmp_path):
    store = RuntimeEventStore(tmp_path / "events.jsonl")
    store.append(
        session_id="s1",
        turn_id="bootstrap",
        kind=EventKind.SESSION_STARTED,
        payload={"cwd": str(tmp_path)},
    )
    store.append(
        session_id="s1",
        turn_id="turn_0001",
        kind=EventKind.TURN_STARTED,
        payload={"request": "optimize", "phase": "synthesis"},
    )
    store.append(
        session_id="s1",
        turn_id="turn_0001",
        kind=EventKind.TOOL_STARTED,
        payload={"request_id": "call-1", "tool": "synthesize_command"},
    )

    state = reduce_events(RuntimeEventStore(store.path).load())

    assert state.active_tool_calls == {"call-1": "synthesize_command"}


def test_durable_project_state_hydrates_session_scoped_compatibility_store(
    tmp_path,
):
    reset_workflow_state()
    with workflow_state_scope("session-1", cwd=tmp_path):
        hydrate_workflow_state(
            {
                "cwd": str(tmp_path),
                "project": {
                    "name": "co2",
                    "program": "gaussian",
                    "path": str(tmp_path / "co2.yaml"),
                    "sha256": "abc123",
                },
                "previous_command": "chemsmart run gaussian -p co2 ...",
            },
            cwd=tmp_path,
        )
        restored = current_workflow_state(tmp_path)

    assert restored.project is not None
    assert restored.project.name == "co2"
    assert restored.previous_command.startswith("chemsmart run gaussian")
    reset_workflow_state()


def test_durable_project_state_overrides_inherited_workspace_default(tmp_path):
    reset_workflow_state()
    with workflow_state_scope("default", cwd=tmp_path):
        hydrate_workflow_state(
            {
                "cwd": str(tmp_path),
                "project": {
                    "name": "workspace-default",
                    "program": "gaussian",
                    "path": str(tmp_path / "workspace-default.yaml"),
                },
            },
            cwd=tmp_path,
            overwrite=True,
        )

    with workflow_state_scope("resumed-session", cwd=tmp_path):
        restored = hydrate_workflow_state(
            {
                "cwd": str(tmp_path),
                "project": {
                    "name": "durable-project",
                    "program": "orca",
                    "path": str(tmp_path / "durable-project.yaml"),
                    "sha256": "durable-hash",
                },
            },
            cwd=tmp_path,
            overwrite=True,
        )

    assert restored.project is not None
    assert restored.project.name == "durable-project"
    assert restored.project.program == "orca"
    assert restored.project.sha256 == "durable-hash"
    reset_workflow_state()


def test_phase_catalog_limits_controller_and_local_tool_surfaces():
    catalog = ToolCatalog(_Registry())
    controller = catalog.select(
        phase=TaskPhase.PROJECT,
        provider_role=ProviderRole.CONTROLLER,
    )
    specialist = catalog.select(
        phase=TaskPhase.PROJECT,
        provider_role=ProviderRole.SYNTHESIS_SPECIALIST,
    )

    assert len(controller.direct) == 5
    assert "write_project_yaml" not in controller.direct
    assert "read_project_yaml" in controller.direct
    assert "critic_project_yaml" in controller.deferred
    assert specialist.direct == ("synthesize_command", "repair_command")
    assert "build_job" in controller.hidden


def test_router_preserves_local_specialist_and_write_boundary():
    assert provider_role("local-mlx") is ProviderRole.SYNTHESIS_SPECIALIST
    assert (
        route_initial_phase(
            "Write the validated project YAML now.",
            role=ProviderRole.CONTROLLER,
        )
        is TaskPhase.PROJECT_WRITE
    )
    assert (
        route_initial_phase(
            "Write a new Gaussian project YAML with B3LYP/def2-SVP.",
            role=ProviderRole.CONTROLLER,
        )
        is TaskPhase.PROJECT
    )
    assert (
        route_initial_phase(
            "Create a project YAML from this method section.",
            role=ProviderRole.CONTROLLER,
        )
        is TaskPhase.PROJECT
    )
    assert (
        route_initial_phase(
            "Read project YAML co2 and explain the settings.",
            role=ProviderRole.CONTROLLER,
        )
        is TaskPhase.PROJECT_READ
    )
    assert (
        route_initial_phase(
            "Read the workspace project YAML named co2 and validate it.",
            role=ProviderRole.CONTROLLER,
        )
        is TaskPhase.PROJECT_READ
    )
    assert (
        route_initial_phase(
            "Create a project YAML.",
            role=ProviderRole.SYNTHESIS_SPECIALIST,
        )
        is TaskPhase.SYNTHESIS
    )


def test_active_lifecycle_rejects_unexposed_internal_tool(tmp_path):
    controller = RuntimeController(
        session_dir=tmp_path,
        session_id="s1",
        registry=_Registry(),
        mode=RuntimeV2Mode.ACTIVE,
    )
    controller.start_turn(
        request="Optimize water with Gaussian.",
        turn_index=1,
        provider_name="openai",
        cwd=str(tmp_path),
    )

    with pytest.raises(ToolExposureViolation, match="not exposed"):
        controller.lifecycle().before_tool(
            request_id="call-1",
            tool_name="build_job",
            arguments={},
        )


def test_shadow_lifecycle_records_violation_without_blocking(tmp_path):
    controller = RuntimeController(
        session_dir=tmp_path,
        session_id="s1",
        registry=_Registry(),
        mode=RuntimeV2Mode.SHADOW,
    )
    controller.start_turn(
        request="Optimize water with Gaussian.",
        turn_index=1,
        provider_name="openai",
        cwd=str(tmp_path),
    )

    controller.lifecycle().before_tool(
        request_id="call-1",
        tool_name="build_job",
        arguments={},
    )

    assert controller.state.shadow_violations == ["runtime.tool.not_exposed"]
    assert controller.state.active_tool_calls == {"call-1": "build_job"}


def test_lifecycle_records_command_project_and_artifact_receipts(tmp_path):
    project = tmp_path / "demo.yaml"
    project.write_text("gas:\n  functional: b3lyp\n")
    generated = tmp_path / "water.com"
    generated.write_text("# b3lyp/6-31g(d) opt\n")
    controller = RuntimeController(
        session_dir=tmp_path / "session",
        session_id="s1",
        registry=_Registry(),
        mode=RuntimeV2Mode.ACTIVE,
    )
    controller.start_turn(
        request="Prepare a Gaussian optimization.",
        turn_index=1,
        provider_name="openai",
        cwd=str(tmp_path),
    )
    lifecycle = controller.lifecycle()
    lifecycle.before_tool(
        request_id="call-1",
        tool_name="synthesize_command",
        arguments={"request": "optimize"},
    )
    lifecycle.after_tool(
        request_id="call-1",
        tool_name="synthesize_command",
        result={
            "command": "chemsmart run gaussian -p demo -f water.xyz -c 0 -m 1 opt",
            "semantic": {"verdict": "ok"},
            "generated_input": str(generated),
            "state_delta": {
                "project": {
                    "selected": True,
                    "project": "demo",
                    "program": "gaussian",
                    "path": str(project),
                    "sha256": "abc123",
                }
            },
        },
    )
    controller.complete()

    replayed = RuntimeController(
        session_dir=tmp_path / "session",
        session_id="s1",
        registry=_Registry(),
        mode=RuntimeV2Mode.ACTIVE,
    )

    assert replayed.state.previous_command.startswith("chemsmart run gaussian")
    assert replayed.state.active_project is not None
    assert replayed.state.active_project.name == "demo"
    assert replayed.state.artifacts[0].sha256
    assert replayed.state.phase is TaskPhase.COMPLETE


def test_lifecycle_records_cli_grounded_direct_dry_run_command(tmp_path):
    command = (
        "chemsmart run gaussian -p water_demo -f h2o.xyz -c 0 -m 1 "
        "scan --coordinates '[[1,2]]' --num-steps 10 --step-size 0.05"
    )
    controller = RuntimeController(
        session_dir=tmp_path / "session",
        session_id="s1",
        registry=_Registry(),
        mode=RuntimeV2Mode.SHADOW,
    )
    controller.start_turn(
        request="Prepare a Gaussian scan.",
        turn_index=1,
        provider_name="openai",
        cwd=str(tmp_path),
    )
    lifecycle = controller.lifecycle()
    lifecycle.before_tool(
        request_id="call-1",
        tool_name="dry_run_input",
        arguments={"job": "job_1"},
    )
    lifecycle.after_tool(
        request_id="call-1",
        tool_name="dry_run_input",
        result={"command": command, "cli_grounded": True},
    )

    assert controller.state.previous_command == command


@pytest.mark.parametrize(
    ("rules", "repeat", "action"),
    [
        (["cmd.runtime.project_not_found"], 1, RepairAction.ASK_USER),
        (["cmd.runtime.dependency_missing"], 1, RepairAction.TERMINATE),
        (["cmd.semantic.option_order"], 1, RepairAction.DETERMINISTIC_REPAIR),
        (["input.gaussian.tddft.root"], 1, RepairAction.REVIEW),
        (["cmd.semantic.option_order"], 2, RepairAction.TERMINATE),
    ],
)
def test_repair_policy_respects_scientific_ownership(rules, repeat, action):
    assert decide_repair(rules, repeated_count=repeat).action is action


def test_controller_rejects_invalid_phase_transition(tmp_path):
    controller = RuntimeController(
        session_dir=tmp_path,
        session_id="s1",
        registry=_Registry(),
        mode=RuntimeV2Mode.ACTIVE,
    )
    controller.start_turn(
        request="Optimize water.",
        turn_index=1,
        provider_name="openai",
        cwd=str(tmp_path),
    )
    decision = AgentDecision(
        action=AgentAction.BUILD_PROJECT,
        phase=TaskPhase.PROJECT,
        summary="Switch to project authoring.",
    )

    with pytest.raises(ValueError, match="invalid runtime transition"):
        controller.validate_decision(decision)


def test_project_authoring_completion_requires_validated_render(tmp_path):
    controller = RuntimeController(
        session_dir=tmp_path,
        session_id="s1",
        registry=_Registry(),
        mode=RuntimeV2Mode.ACTIVE,
    )
    controller.start_turn(
        request="Create a Gaussian project YAML for water.",
        turn_index=1,
        provider_name="openai",
        cwd=str(tmp_path),
    )

    assert controller.complete() is False
    assert controller.state.phase is TaskPhase.BLOCKED
    assert controller.state.blocked_reason == "runtime.project.render_required"


def test_project_authoring_completion_accepts_validated_render(tmp_path):
    controller = RuntimeController(
        session_dir=tmp_path,
        session_id="s1",
        registry=_Registry(),
        mode=RuntimeV2Mode.ACTIVE,
    )
    controller.start_turn(
        request="Create a Gaussian project YAML for water.",
        turn_index=1,
        provider_name="openai",
        cwd=str(tmp_path),
    )
    lifecycle = controller.lifecycle()
    lifecycle.before_tool(
        request_id="call-1",
        tool_name="render_project_yaml",
        arguments={"protocol": {}},
    )
    lifecycle.after_tool(
        request_id="call-1",
        tool_name="render_project_yaml",
        result={
            "ok": True,
            "yaml_text": "gas: {}",
            "validation": {"verdict": "ok"},
        },
    )

    assert controller.complete() is True
    assert controller.state.phase is TaskPhase.COMPLETE
    assert "no workspace file was written" in controller.completion_notice()
    assert "/write-project" in controller.completion_notice()
