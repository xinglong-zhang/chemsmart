from __future__ import annotations

import subprocess
from dataclasses import dataclass
from pathlib import Path

import pytest

import chemsmart.agent.tools_command as tools_command
from chemsmart.agent.harness.command_semantics import CommandSemanticResult
from chemsmart.agent.harness.intent import IntentSpec
from chemsmart.agent.harness.workflow_state import (
    select_workspace_project,
    workflow_state_scope,
)


@pytest.fixture(autouse=True)
def _reset_command_tools_state():
    tools_command.reset_command_tools_state()
    yield
    tools_command.reset_command_tools_state()


@dataclass
class Completed:
    returncode: int = 0
    stdout: str = ""
    stderr: str = ""


def test_execute_command_rejects_before_subprocess(monkeypatch):
    calls = {"run": 0}

    monkeypatch.setattr(
        tools_command,
        "_gate_command",
        lambda command: CommandSemanticResult(
            verdict="reject",
            command=command,
        ),
    )

    def fake_run(*args, **kwargs):
        calls["run"] += 1
        return Completed()

    monkeypatch.setattr(tools_command.subprocess, "run", fake_run)

    result = tools_command.execute_chemsmart_command(
        "chemsmart run gaussian -f examples/h2o.xyz -c 0 -m 1 sp"
    )

    assert result["ok"] is False
    assert result["status"] == "rejected"
    assert calls["run"] == 0


def test_public_synthesis_trace_excludes_provider_private_reasoning():
    trace = tools_command._public_synthesis_trace(
        {
            "status": "ready",
            "command": (
                "chemsmart run gaussian -p demo -f h2o.xyz -c 0 -m 1 opt"
            ),
            "decision_trace": {
                "decision_summary": "The user requested a geometry optimization.",
                "evidence": ["A structure file and spin state were supplied."],
            },
        },
        CommandSemanticResult(
            verdict="ok",
            command="chemsmart run gaussian -p demo -f h2o.xyz -c 0 -m 1 opt",
        ),
    )

    assert "geometry optimization" in trace
    assert "Selected gaussian opt command." in trace
    assert "Deterministic semantic gate verdict: ok." in trace
    assert "private" not in trace.lower()


def test_execute_command_test_mode_adds_fake_runner_flags(monkeypatch):
    seen = {}

    monkeypatch.setattr(
        tools_command,
        "_gate_command",
        lambda command: CommandSemanticResult(
            verdict="ok",
            command=command,
        ),
    )

    def fake_run(argv, **kwargs):
        seen["argv"] = argv
        return Completed(returncode=0, stdout="ok\n")

    monkeypatch.setattr(tools_command.subprocess, "run", fake_run)

    result = tools_command.execute_chemsmart_command(
        "chemsmart run gaussian -f examples/h2o.xyz -c 0 -m 1 sp",
        test=True,
    )

    assert result["ok"] is True
    assert "--no-verbose" in seen["argv"]
    assert "--fake" in seen["argv"]
    assert "--no-scratch" in seen["argv"]


def test_execute_command_test_mode_adds_submit_test_flags(monkeypatch):
    seen = {}

    monkeypatch.setattr(
        tools_command,
        "_gate_command",
        lambda command: CommandSemanticResult(
            verdict="ok",
            command=command,
        ),
    )

    def fake_run(argv, **kwargs):
        seen["argv"] = argv
        return Completed(returncode=0)

    monkeypatch.setattr(tools_command.subprocess, "run", fake_run)

    tools_command.execute_chemsmart_command(
        "chemsmart sub gaussian -f examples/h2o.xyz -c 0 -m 1 sp",
        test=True,
    )

    assert "--test" in seen["argv"]
    assert "--fake" in seen["argv"]


def test_submit_terminal_state_records_server_yaml_and_scheduler(
    monkeypatch, tmp_path
):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("HOME", str(tmp_path))
    server_dir = tmp_path / ".chemsmart" / "server"
    server_dir.mkdir(parents=True)
    (server_dir / "mock-pbs.yaml").write_text(
        "SERVER:\n  SCHEDULER: PBS\n", encoding="utf-8"
    )
    project_dir = tmp_path / ".chemsmart" / "gaussian"
    project_dir.mkdir(parents=True)
    (project_dir / "demo.yaml").write_text(
        "gas:\n  functional: b3lyp\n  basis: def2svp\n"
        "solv:\n  functional: b3lyp\n  basis: def2svp\n",
        encoding="utf-8",
    )
    select_workspace_project("demo", "gaussian")

    monkeypatch.setattr(
        tools_command,
        "_gate_command",
        lambda command: CommandSemanticResult(verdict="ok", command=command),
    )
    def fake_run(*args, **kwargs):
        (tmp_path / "chemsmart_sub_h2o_opt.sh").write_text(
            "#!/bin/sh\n#PBS -N h2o\n", encoding="utf-8"
        )
        return Completed(returncode=0, stdout="job-1\n")

    monkeypatch.setattr(tools_command.subprocess, "run", fake_run)

    command = (
        "chemsmart sub -s mock-pbs gaussian -p demo "
        "-f h2o.xyz -c 0 -m 1 opt"
    )
    tools_command.register_command_intent(
        command,
        IntentSpec.from_dict(
            {
                "action": "sub",
                "program": "gaussian",
                "kind": "gaussian.opt",
                "project": "demo",
                "server": "mock-pbs",
                "input_path": "h2o.xyz",
                "charge": 0,
                "multiplicity": 1,
            }
        ),
    )
    result = tools_command.execute_chemsmart_command(command)

    state = result["terminal_state"]
    assert state["server"] == "mock-pbs"
    assert state["scheduler"] == "PBS"
    assert state["status"] == "passed"
    assert state["all_passed"] is True


def test_submit_terminal_state_rejects_stale_submit_script(
    monkeypatch, tmp_path
):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("HOME", str(tmp_path))
    server_dir = tmp_path / ".chemsmart" / "server"
    server_dir.mkdir(parents=True)
    (server_dir / "mock-pbs.yaml").write_text(
        "SERVER:\n  SCHEDULER: PBS\n", encoding="utf-8"
    )
    project_dir = tmp_path / ".chemsmart" / "gaussian"
    project_dir.mkdir(parents=True)
    (project_dir / "demo.yaml").write_text(
        "gas:\n  functional: b3lyp\n  basis: def2svp\n"
        "solv:\n  functional: b3lyp\n  basis: def2svp\n",
        encoding="utf-8",
    )
    (tmp_path / "chemsmart_sub_stale.sh").write_text(
        "#!/bin/sh\n#PBS -N stale\n", encoding="utf-8"
    )
    select_workspace_project("demo", "gaussian")
    command = (
        "chemsmart sub -s mock-pbs gaussian -p demo "
        "-f h2o.xyz -c 0 -m 1 opt"
    )
    tools_command.register_command_intent(
        command,
        IntentSpec.from_dict(
            {
                "action": "sub",
                "program": "gaussian",
                "kind": "gaussian.opt",
                "project": "demo",
                "server": "mock-pbs",
                "input_path": "h2o.xyz",
                "charge": 0,
                "multiplicity": 1,
            }
        ),
    )
    monkeypatch.setattr(
        tools_command,
        "_gate_command",
        lambda value: CommandSemanticResult(verdict="ok", command=value),
    )
    monkeypatch.setattr(
        tools_command.subprocess,
        "run",
        lambda *args, **kwargs: Completed(returncode=0),
    )

    result = tools_command.execute_chemsmart_command(command)
    assertions = {
        row["id"]: row for row in result["terminal_state"]["assertions"]
    }

    assert result["ok"] is True
    assert result["terminal_state"]["status"] == "failed"
    assert assertions["sub.submit_script_exists"]["status"] == "fail"


def test_submit_terminal_state_requires_independent_intent_provenance(
    monkeypatch, tmp_path
):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.setattr(
        tools_command,
        "_gate_command",
        lambda command: CommandSemanticResult(verdict="ok", command=command),
    )
    monkeypatch.setattr(
        tools_command.subprocess,
        "run",
        lambda *args, **kwargs: Completed(returncode=0),
    )

    command = (
        "chemsmart sub -s mock-pbs gaussian -p demo "
        "-f changed.xyz -c 0 -m 1 opt"
    )
    result = tools_command.execute_chemsmart_command(command)
    assertions = {
        row["id"]: row for row in result["terminal_state"]["assertions"]
    }

    assert result["ok"] is True
    assert result["terminal_state"]["status"] == "failed"
    assert assertions["sub.intent_provenance_present"]["status"] == "fail"


def test_submit_terminal_state_rejects_cached_input_path_drift(
    monkeypatch, tmp_path
):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.setattr(
        tools_command,
        "_gate_command",
        lambda command: CommandSemanticResult(verdict="ok", command=command),
    )
    monkeypatch.setattr(
        tools_command.subprocess,
        "run",
        lambda *args, **kwargs: Completed(returncode=0),
    )

    command = (
        "chemsmart sub -s mock-pbs gaussian -p demo "
        "-f changed.xyz -c 0 -m 1 opt"
    )
    with workflow_state_scope("drift-test"):
        tools_command.register_command_intent(
            command,
            IntentSpec.from_dict(
                {
                    "action": "sub",
                    "program": "gaussian",
                    "kind": "gaussian.opt",
                    "project": "demo",
                    "server": "mock-pbs",
                    "input_path": "requested.xyz",
                    "charge": 0,
                    "multiplicity": 1,
                }
            ),
        )
        result = tools_command.execute_chemsmart_command(command)

    failed = {
        row["id"]
        for row in result["terminal_state"]["assertions"]
        if row["status"] == "fail"
    }
    assert "sub.intent_filename" in failed
    assert result["terminal_state"]["status"] == "failed"


def test_execute_command_timeout_returns_structured_result(monkeypatch):
    monkeypatch.setattr(
        tools_command,
        "_gate_command",
        lambda command: CommandSemanticResult(
            verdict="ok",
            command=command,
        ),
    )

    def fake_run(argv, **kwargs):
        raise subprocess.TimeoutExpired(
            cmd=argv, timeout=5, output="partial out", stderr="partial err"
        )

    monkeypatch.setattr(tools_command.subprocess, "run", fake_run)

    result = tools_command.execute_chemsmart_command(
        "chemsmart run gaussian -f examples/h2o.xyz -c 0 -m 1 sp",
        timeout_s=5,
    )

    assert result["ok"] is False
    assert result["status"] == "timeout"
    assert result["returncode"] is None
    assert "partial out" in result["stdout_tail"]
    assert "partial err" in result["stderr_tail"]


def test_gate_cache_invalidated_by_workspace_yaml_change(
    monkeypatch, tmp_path
):
    # A cached ok verdict must NOT survive a workspace project YAML edit —
    # execution after update_project_yaml has to re-run the semantic gate.
    monkeypatch.chdir(tmp_path)
    calls = {"gate": 0}

    def fake_gate(command, **kwargs):
        calls["gate"] += 1
        return CommandSemanticResult(verdict="ok", command=command)

    monkeypatch.setattr(tools_command, "evaluate_command_semantics", fake_gate)

    command = "chemsmart run gaussian -p x -f a.xyz -c 0 -m 1 opt"
    ws = Path(".chemsmart/gaussian")
    ws.mkdir(parents=True)
    yaml_path = ws / "x.yaml"
    yaml_path.write_text(
        "gas:\n  functional: b3lyp\n  basis: def2svp\n"
        "solv:\n  functional: b3lyp\n  basis: def2svp\n"
    )

    first = tools_command._gate_command(command)
    second = tools_command._gate_command(command)
    assert first.verdict == second.verdict == "ok"
    assert calls["gate"] == 1  # second call hit the cache

    # Simulate update_project_yaml: change the file content (mtime/size).
    yaml_path.write_text(
        "gas:\n  functional: m062x\n  basis: def2tzvp\n"
        "solv:\n  functional: m062x\n  basis: def2tzvp\n"
    )

    tools_command._gate_command(command)
    assert calls["gate"] == 2  # fingerprint change forced a re-gate


def test_session_refreshes_default_project_each_call(monkeypatch, tmp_path):
    # The shared SynthesisSession must pick up a workspace project written
    # AFTER the session was created (the build-yaml -> synthesize chain).
    monkeypatch.chdir(tmp_path)

    class DummyProvider:
        name = "openai"

        def chat(self, messages, **kwargs):
            raise AssertionError("no chat expected")

    monkeypatch.setattr(
        tools_command,
        "build_chemsmart_cli_schema",
        lambda: {"subcommands": {}},
    )
    import chemsmart.agent.providers as providers_mod

    monkeypatch.setattr(
        providers_mod, "get_provider", lambda *a, **k: DummyProvider()
    )

    session_before = tools_command._session_for_cwd()
    assert session_before.default_project == ""

    ws = Path(".chemsmart/gaussian")
    ws.mkdir(parents=True)
    (ws / "demo.yaml").write_text(
        "gas:\n  functional: b3lyp\n  basis: def2svp\n"
        "solv:\n  functional: b3lyp\n  basis: def2svp\n"
    )

    session_after = tools_command._session_for_cwd()
    assert session_after is session_before  # same cached session object
    assert session_after.default_project == "demo"


def _mini_full_schema():
    def node(name, **subcommands):
        return {
            "name": name,
            "description": None,
            "options": [],
            "subcommands": dict(subcommands),
        }

    gaussian = node(
        "gaussian",
        opt=node("opt", qmmm=node("qmmm")),
        sp=node("sp"),
        ts=node("ts"),
    )
    orca = node("orca", opt=node("opt"), sp=node("sp"), ts=node("ts"))
    return node(
        "chemsmart",
        run=node("run", gaussian=gaussian, orca=orca),
        sub=node("sub", gaussian=gaussian, orca=orca),
        agent=node("agent"),
    )


@pytest.fixture
def _pruning_session(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)

    class DummyProvider:
        name = "openai"

        def chat(self, messages, **kwargs):
            raise AssertionError("no chat expected")

    monkeypatch.setattr(
        tools_command, "build_chemsmart_cli_schema", _mini_full_schema
    )
    import chemsmart.agent.providers as providers_mod

    monkeypatch.setattr(
        providers_mod, "get_provider", lambda *a, **k: DummyProvider()
    )
    return tools_command._session_for_cwd()


def test_synthesize_command_uses_pruned_schema_and_restores(
    _pruning_session,
):
    session = _pruning_session
    full_schema = session.schema
    seen = {}

    def fake_prepare(request):
        seen["schema"] = session.schema
        return {
            "status": "ready",
            "command": "chemsmart run gaussian -f a.xyz -c 0 -m 1 opt",
            "explanation": "ok",
            "confidence": "high",
        }

    session.prepare_command = fake_prepare

    payload = tools_command.synthesize_command(
        "optimize a.xyz with gaussian, charge 0 multiplicity 1"
    )

    assert payload["schema_variant"] == "run/gaussian[opt,sp]"
    pruned = seen["schema"]
    assert pruned is not full_schema
    assert set(pruned["subcommands"]) == {"run"}
    programs = pruned["subcommands"]["run"]["subcommands"]
    assert set(programs) == {"gaussian"}
    assert set(programs["gaussian"]["subcommands"]) == {"opt", "sp"}
    # The shared session must leave with the full schema for other callers.
    assert session.schema is full_schema


def test_synthesize_command_falls_back_to_full_schema_on_infeasible(
    _pruning_session,
):
    session = _pruning_session
    full_schema = session.schema
    calls = []

    def fake_prepare(request):
        calls.append(session.schema)
        if len(calls) == 1:
            return {"status": "infeasible", "command": "", "explanation": ""}
        return {
            "status": "ready",
            "command": "chemsmart run gaussian -f a.xyz -c 0 -m 1 ts",
            "explanation": "ok",
            "confidence": "high",
        }

    session.prepare_command = fake_prepare

    payload = tools_command.synthesize_command(
        "transition state of a.xyz with gaussian"
    )

    assert len(calls) == 2
    assert calls[0] is not full_schema  # first attempt pruned
    assert calls[1] is full_schema  # retry with the full schema
    assert payload["schema_variant"] == "full-fallback"
    assert payload["status"] == "ready"
    assert session.schema is full_schema


def test_repair_command_runs_against_full_schema(
    _pruning_session, monkeypatch
):
    session = _pruning_session
    full_schema = session.schema
    monkeypatch.setattr(
        tools_command,
        "_gate_command",
        lambda command: CommandSemanticResult(
            verdict="reject", command=command
        ),
    )
    seen = {}

    def fake_repair(request, result):
        seen["schema"] = session.schema
        return {
            "status": "ready",
            "command": "chemsmart run gaussian -f a.xyz -c 0 -m 1 opt",
        }

    session._repair_ready_result = fake_repair
    # Simulate a stale pruned schema left behind by a crashed caller.
    session.schema = {"subcommands": {}}

    tools_command.repair_command(
        "chemsmart run gaussian -f a.xyz opt", failure="missing charge"
    )

    assert seen["schema"] is full_schema
