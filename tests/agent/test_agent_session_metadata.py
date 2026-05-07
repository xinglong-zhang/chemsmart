from __future__ import annotations

import json
from pathlib import Path

import pytest

from chemsmart.agent.core import AgentSession
from chemsmart.agent.registry import ToolRegistry

from ._agent_session_helpers import FakeProvider, planner_plan


def _read_decision_log(session_dir: Path) -> list[dict]:
    return [
        json.loads(line)
        for line in (session_dir / "decision_log.jsonl")
        .read_text()
        .splitlines()
    ]


def _patch_runtime_ok(monkeypatch):
    import chemsmart.agent.tools as agent_tools

    def fake_validate_runtime(job, server=None):
        return {
            "ok": "ok",
            "local_ok": True,
            "local_issues": [],
            "remote_unknown": [],
        }

    monkeypatch.setattr(agent_tools, "validate_runtime", fake_validate_runtime)


def _patch_run_local(monkeypatch):
    import chemsmart.agent.tools as agent_tools

    def fake_run_local(job):
        return {
            "ok": True,
            "returncode": 0,
            "stdout_path": str(Path(job.folder) / "session.stdout"),
            "stderr_path": str(Path(job.folder) / "session.stderr"),
            "output_summary": {},
        }

    monkeypatch.setattr(agent_tools, "run_local", fake_run_local)


def test_session_metadata_json_is_written_with_required_keys(
    monkeypatch,
    single_molecule_xyz_file,
    tmp_path: Path,
):
    _patch_runtime_ok(monkeypatch)

    provider = FakeProvider(
        [
            planner_plan(single_molecule_xyz_file, "metadata_case"),
            {
                "verdict": "ok",
                "confidence": 0.83,
                "issues": [],
                "rationale": "Looks good.",
            },
        ]
    )
    session = AgentSession(
        provider=provider,
        registry=ToolRegistry.default(),
        session_root=tmp_path / "sessions",
    )

    result = session.run(
        f"optimize {single_molecule_xyz_file}", dry_submit=True
    )

    assert result["blocked"] is False
    metadata = json.loads(
        (Path(result["session_dir"]) / "session_metadata.json").read_text()
    )
    assert {
        "session_id",
        "request",
        "intent",
        "plan_steps",
        "executed_steps",
        "input_file",
        "critic_verdict",
        "critic_confidence",
        "blocked",
        "wall_time_ms",
        "started_at",
        "ended_at",
    } <= metadata.keys()
    assert metadata["intent"] == "opt"
    assert metadata["plan_steps"] == 5
    assert metadata["critic_confidence"] == 0.83
    assert metadata["input_file"].endswith("metadata_case.com")


def test_tool_error_is_logged_before_runtime_error(
    single_molecule_xyz_file,
    tmp_path: Path,
):
    provider = FakeProvider(
        [
            {
                "steps": [
                    {
                        "tool": "build_molecule",
                        "args": {"filepath": single_molecule_xyz_file},
                        "rationale": "Load structure.",
                    },
                    {
                        "tool": "build_gaussian_settings",
                        "args": {
                            "functional": "b3lyp",
                            "basis": "6-31g*",
                        },
                        "rationale": "Make settings.",
                    },
                    {
                        "tool": "build_job",
                        "args": {
                            "kind": "gaussian.frequency",
                            "molecule": "$step1",
                            "settings": "$step2",
                            "label": "error_case",
                        },
                        "rationale": "Trigger a tool error.",
                    },
                ],
                "rationale": "Bad plan for testing.",
                "estimated_cost": "low",
            }
        ]
    )
    session = AgentSession(
        provider=provider,
        registry=ToolRegistry.default(),
        session_root=tmp_path / "sessions",
    )

    with pytest.raises(RuntimeError, match="validation error"):
        session.run("frequency test failure", dry_submit=True)

    assert session.session_dir is not None
    entries = _read_decision_log(session.session_dir)
    kinds = [entry["kind"] for entry in entries]
    assert "tool_error" in kinds
    assert kinds.index("tool_error") < kinds.index("session_summary")
    tool_error = next(
        entry for entry in entries if entry["kind"] == "tool_error"
    )
    assert tool_error["payload"]["tool"] == "build_job"
    assert tool_error["payload"]["error_type"] == "ValidationError"


def test_session_summary_is_last_entry_with_tools_called(
    monkeypatch,
    single_molecule_xyz_file,
    tmp_path: Path,
):
    _patch_runtime_ok(monkeypatch)
    _patch_run_local(monkeypatch)

    provider = FakeProvider(
        [
            planner_plan(
                single_molecule_xyz_file, "summary_case", "run_local"
            ),
            {
                "verdict": "ok",
                "confidence": 0.72,
                "issues": [],
                "rationale": "Looks good.",
            },
        ]
    )
    session = AgentSession(
        provider=provider,
        registry=ToolRegistry.default(),
        session_root=tmp_path / "sessions",
    )

    result = session.run("optimize then run locally", dry_submit=True)

    assert result["blocked"] is False
    entries = _read_decision_log(Path(result["session_dir"]))
    assert entries[-1]["kind"] == "session_summary"
    assert entries[-1]["payload"]["tools_called"] == [
        "build_molecule",
        "build_gaussian_settings",
        "build_job",
        "dry_run_input",
        "validate_runtime",
        "run_local",
    ]


def test_critic_confidence_is_logged_in_critic_verdict_entry(
    monkeypatch,
    single_molecule_xyz_file,
    tmp_path: Path,
):
    _patch_runtime_ok(monkeypatch)

    provider = FakeProvider(
        [
            planner_plan(single_molecule_xyz_file, "confidence_case"),
            {
                "verdict": "ok",
                "confidence": 0.91,
                "issues": [],
                "rationale": "High confidence.",
            },
        ]
    )
    session = AgentSession(
        provider=provider,
        registry=ToolRegistry.default(),
        session_root=tmp_path / "sessions",
    )

    result = session.run("optimize with confidence logging", dry_submit=True)

    verdict_entry = next(
        entry
        for entry in _read_decision_log(Path(result["session_dir"]))
        if entry["kind"] == "critic_verdict"
    )
    assert verdict_entry["payload"]["confidence"] == 0.91
