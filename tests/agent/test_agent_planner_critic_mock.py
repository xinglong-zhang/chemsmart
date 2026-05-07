from __future__ import annotations

from pathlib import Path

import pytest

from chemsmart.agent.core import AgentSession, Plan
from chemsmart.agent.registry import ToolRegistry

from ._agent_session_helpers import FakeProvider, critic_ok, planner_plan


@pytest.mark.parametrize(
    ("user_request", "label"),
    [
        ("optimize examples/h2o.xyz", "plan_h2o"),
        ("optimize then frequency examples/co2.xyz", "plan_co2"),
        ("optimize methane with Gaussian", "plan_ch4"),
    ],
)
def test_planner_emits_valid_plan_for_canonical_requests(
    single_molecule_xyz_file,
    tmp_path: Path,
    user_request: str,
    label: str,
):
    provider = FakeProvider([planner_plan(single_molecule_xyz_file, label)])
    session = AgentSession(
        provider=provider,
        registry=ToolRegistry.default(),
        session_root=tmp_path,
    )

    plan = session._planner_call(user_request)

    assert isinstance(plan, Plan)
    assert plan.steps[0].tool == "build_molecule"
    assert plan.steps[-1].tool == "validate_runtime"


def test_critic_blocks_malformed_input(
    monkeypatch,
    single_molecule_xyz_file,
    tmp_path: Path,
):
    import chemsmart.agent.tools as agent_tools

    provider = FakeProvider(
        [
            planner_plan(
                single_molecule_xyz_file, "malformed_case", "run_local"
            ),
            critic_ok(),
        ]
    )

    def fake_validate_runtime(job, server=None):
        return {
            "ok": "ok",
            "local_ok": True,
            "local_issues": [],
            "remote_unknown": [],
        }

    def fake_dry_run_input(job):
        return {
            "inputfile": str(Path(job.folder) / "malformed.com"),
            "content": "%chk=bad.chk\nmissing route line\n",
        }

    monkeypatch.setattr(agent_tools, "validate_runtime", fake_validate_runtime)
    monkeypatch.setattr(agent_tools, "dry_run_input", fake_dry_run_input)

    session = AgentSession(
        provider=provider,
        registry=ToolRegistry.default(),
        session_root=tmp_path,
    )
    result = session.run("optimize malformed input", dry_submit=True)

    assert result["blocked"] is True
    assert result["critic_verdict"].verdict == "warn"
    assert (
        "Gaussian route line missing or malformed"
        in result["critic_verdict"].issues
    )


def test_critic_refuses_partial_without_allow_remote_unknown(
    monkeypatch,
    single_molecule_xyz_file,
    tmp_path: Path,
):
    import chemsmart.agent.tools as agent_tools

    provider = FakeProvider(
        [
            planner_plan(
                single_molecule_xyz_file, "partial_case", "run_local"
            ),
            critic_ok(),
        ]
    )
    run_local_calls = {"count": 0}

    def fake_run_local(job):
        run_local_calls["count"] += 1
        return {
            "ok": True,
            "returncode": 0,
            "stdout_path": str(Path(job.folder) / "run.stdout"),
            "stderr_path": str(Path(job.folder) / "run.stderr"),
            "output_summary": {},
        }

    monkeypatch.setattr(agent_tools, "run_local", fake_run_local)

    session = AgentSession(
        provider=provider,
        registry=ToolRegistry.default(),
        session_root=tmp_path,
    )
    result = session.run("optimize with partial runtime", dry_submit=True)

    assert result["blocked"] is True
    assert result["critic_verdict"].verdict == "warn"
    assert run_local_calls["count"] == 0
    assert "server.queue required" in result["critic_verdict"].issues
