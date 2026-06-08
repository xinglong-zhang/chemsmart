from __future__ import annotations

from pathlib import Path

from chemsmart.agent.core import AgentSession
from chemsmart.agent.registry import ToolRegistry

from ._agent_session_helpers import (
    FakeProvider,
    critic_ok,
    critic_reject,
    planner_plan,
)


def test_warn_with_override_executes_but_warn_without_override_refuses(
    monkeypatch,
    single_molecule_xyz_file,
    tmp_path: Path,
):
    import chemsmart.agent.tools as agent_tools

    calls = {"count": 0}

    def fake_validate_runtime(job, server=None):
        return {
            "ok": "ok",
            "local_ok": True,
            "local_issues": [],
            "remote_unknown": [],
        }

    def fake_dry_run_input(job):
        return {
            "inputfile": str(Path(job.folder) / "override.com"),
            "content": "%chk=override.chk\nmissing route line\n",
        }

    def fake_run_local(job):
        calls["count"] += 1
        return {
            "ok": True,
            "returncode": 0,
            "stdout_path": str(Path(job.folder) / "override.stdout"),
            "stderr_path": str(Path(job.folder) / "override.stderr"),
            "output_summary": {},
        }

    monkeypatch.setattr(agent_tools, "validate_runtime", fake_validate_runtime)
    monkeypatch.setattr(agent_tools, "dry_run_input", fake_dry_run_input)
    monkeypatch.setattr(agent_tools, "run_local", fake_run_local)

    provider = FakeProvider(
        [
            planner_plan(
                single_molecule_xyz_file, "override_case", "run_local"
            ),
            critic_ok(),
        ]
    )
    session = AgentSession(
        provider=provider,
        registry=ToolRegistry.default(),
        session_root=tmp_path / "no_override",
    )
    blocked = session.run("warn without override", dry_submit=True)

    assert blocked["blocked"] is True
    assert calls["count"] == 0

    provider = FakeProvider(
        [
            planner_plan(
                single_molecule_xyz_file, "override_case_2", "run_local"
            ),
            critic_ok(),
        ]
    )
    session = AgentSession(
        provider=provider,
        registry=ToolRegistry.default(),
        session_root=tmp_path / "with_override",
    )
    allowed = session.run(
        "warn with override",
        dry_submit=True,
        allow_critic_override=True,
    )

    assert allowed["blocked"] is False
    assert calls["count"] == 1


def test_reject_never_executes_even_with_override(
    monkeypatch,
    single_molecule_xyz_file,
    tmp_path: Path,
):
    import chemsmart.agent.tools as agent_tools

    calls = {"count": 0}

    def fake_validate_runtime(job, server=None):
        return {
            "ok": "fail",
            "local_ok": False,
            "local_issues": ["server.executable_path missing"],
            "remote_unknown": [],
        }

    def fake_run_local(job):
        calls["count"] += 1
        return {
            "ok": True,
            "returncode": 0,
            "stdout_path": str(Path(job.folder) / "reject.stdout"),
            "stderr_path": str(Path(job.folder) / "reject.stderr"),
            "output_summary": {},
        }

    monkeypatch.setattr(agent_tools, "validate_runtime", fake_validate_runtime)
    monkeypatch.setattr(agent_tools, "run_local", fake_run_local)

    provider = FakeProvider(
        [
            planner_plan(single_molecule_xyz_file, "reject_case", "run_local"),
            critic_reject("bad plan"),
        ]
    )
    session = AgentSession(
        provider=provider,
        registry=ToolRegistry.default(),
        session_root=tmp_path,
    )
    result = session.run(
        "reject regardless of override",
        dry_submit=True,
        allow_remote_unknown=True,
        allow_critic_override=True,
    )

    assert result["blocked"] is True
    assert result["critic_verdict"].verdict == "reject"
    assert calls["count"] == 0
