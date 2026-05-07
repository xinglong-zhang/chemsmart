from __future__ import annotations

import os
from pathlib import Path

from chemsmart.agent.core import AgentSession
from chemsmart.agent.registry import ToolRegistry

from ._agent_session_helpers import FakeProvider, critic_ok, planner_plan


def test_session_resume_uses_original_absolute_cwd(
    monkeypatch,
    single_molecule_xyz_file,
    tmp_path: Path,
):
    import chemsmart.agent.tools as agent_tools

    dir_a = tmp_path / "dir_a"
    dir_b = tmp_path / "dir_b"
    dir_a.mkdir()
    dir_b.mkdir()

    provider = FakeProvider(
        [
            planner_plan(single_molecule_xyz_file, "resume_case", "run_local"),
            critic_ok(),
        ]
    )
    recorded = {}

    def fake_run_local(job):
        recorded["cwd"] = os.getcwd()
        recorded["job_folder"] = job.folder
        return {
            "ok": True,
            "returncode": 0,
            "stdout_path": str(Path(job.folder) / "resume.stdout"),
            "stderr_path": str(Path(job.folder) / "resume.stderr"),
            "output_summary": {},
        }

    monkeypatch.setattr(agent_tools, "run_local", fake_run_local)

    old_cwd = os.getcwd()
    os.chdir(dir_a)
    try:
        session = AgentSession(
            provider=provider,
            registry=ToolRegistry.default(),
            session_root=tmp_path / "sessions",
        )
        first = session.run(
            "prepare resume case",
            dry_submit=True,
            allow_remote_unknown=False,
        )
        assert first["blocked"] is True
        session_id = first["session_id"]
    finally:
        os.chdir(old_cwd)

    os.chdir(dir_b)
    try:
        resumed = AgentSession.resume(
            session_id,
            provider=provider,
            registry=ToolRegistry.default(),
            session_root=tmp_path / "sessions",
            dry_submit=True,
            allow_remote_unknown=True,
        )
    finally:
        os.chdir(old_cwd)

    assert resumed["blocked"] is False
    assert Path(recorded["cwd"]).resolve() == dir_a.resolve()
    assert Path(recorded["job_folder"]).resolve() == dir_a.resolve()
