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


def _planner_ts_plan(filepath: str, label: str):
    return {
        "steps": [
            {
                "tool": "build_molecule",
                "args": {"filepath": filepath},
                "rationale": "Load the structure.",
            },
            {
                "tool": "build_gaussian_settings",
                "args": {"functional": "b3lyp", "basis": "6-31g*"},
                "rationale": "Use a simple DFT level.",
            },
            {
                "tool": "build_job",
                "args": {
                    "kind": "gaussian.ts",
                    "molecule": "$step1",
                    "settings": "$step2",
                    "label": label,
                },
                "rationale": "Assemble the Gaussian TS job.",
            },
            {
                "tool": "dry_run_input",
                "args": {"job": "$step3"},
                "rationale": "Write the input file for inspection.",
            },
            {
                "tool": "validate_runtime",
                "args": {"job": "$step3", "server": None},
                "rationale": "Check local/runtime prerequisites.",
            },
        ],
        "rationale": "Prepare a Gaussian transition-state workflow.",
        "estimated_cost": "low",
    }


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


def test_harness_rejects_malformed_ts_route_even_when_critic_ok(
    monkeypatch,
    single_molecule_xyz_file,
    tmp_path: Path,
):
    import json

    import chemsmart.agent.tools as agent_tools

    def fake_validate_runtime(job, server=None):
        return {
            "ok": "ok",
            "local_ok": True,
            "local_issues": [],
            "remote_unknown": [],
        }

    def fake_dry_run_input(job):
        return {
            "inputfile": str(Path(job.folder) / "bad_ts.com"),
            "content": (
                "%chk=bad_ts.chk\n"
                "# b3lyp/6-31g* opt=(ts,calcfc,noeigentest,['ts'])\n\n"
            ),
        }

    monkeypatch.setattr(agent_tools, "validate_runtime", fake_validate_runtime)
    monkeypatch.setattr(agent_tools, "dry_run_input", fake_dry_run_input)

    provider = FakeProvider(
        [
            _planner_ts_plan(single_molecule_xyz_file, "bad_ts"),
            critic_ok(),
        ]
    )
    session = AgentSession(
        provider=provider,
        registry=ToolRegistry.default(),
        session_root=tmp_path,
    )

    result = session.run(
        "find a transition state",
        dry_submit=True,
        allow_critic_override=True,
    )

    assert result["blocked"] is True
    assert result["critic_verdict"].verdict == "reject"
    assert any(
        issue.startswith("gaussian.ts.route:")
        for issue in result["critic_verdict"].issues
    )

    session_dir = Path(result["session_dir"])
    harness_result = json.loads(
        (session_dir / "harness_result.json").read_text()
    )
    metadata = json.loads((session_dir / "session_metadata.json").read_text())
    assert harness_result["verdict"] == "reject"
    assert harness_result["failed_rule_ids"] == ["gaussian.ts.route"]
    assert metadata["harness_verdict"] == "reject"
    assert metadata["harness_issue_count"] >= 1
