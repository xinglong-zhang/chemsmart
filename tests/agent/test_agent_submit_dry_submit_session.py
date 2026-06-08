from __future__ import annotations

import json
from pathlib import Path

from chemsmart.agent.core import AgentSession
from chemsmart.agent.registry import ToolRegistry

from ._agent_session_helpers import FakeProvider, critic_ok


def _orca_submit_plan(filepath: str) -> dict:
    return {
        "steps": [
            {
                "tool": "build_molecule",
                "args": {"filepath": filepath},
                "rationale": "Load the requested structure.",
            },
            {
                "tool": "build_orca_settings",
                "args": {
                    "ab_initio": "DLPNO-CCSD(T)",
                    "functional": None,
                    "basis": "cc-pVTZ",
                },
                "rationale": "Use the requested correlated ORCA level.",
            },
            {
                "tool": "build_job",
                "args": {
                    "kind": "orca.sp",
                    "molecule": "$step1",
                    "settings": "$step2",
                    "label": "orca_sp_case",
                },
                "rationale": "Assemble the ORCA single-point job.",
            },
            {
                "tool": "dry_run_input",
                "args": {"job": "$step3"},
                "rationale": "Render the ORCA input before execution.",
            },
            {
                "tool": "validate_runtime",
                "args": {"job": "$step3", "server": None},
                "rationale": "Confirm runtime prerequisites before submit.",
            },
            {
                "tool": "submit_hpc",
                "args": {"job": "$step3"},
                "rationale": "Submit the prepared ORCA job after review.",
            },
        ],
        "rationale": "Prepare an ORCA single-point submission workflow.",
        "estimated_cost": "medium",
    }


def _read_decision_log(session_dir: Path) -> list[dict]:
    return [
        json.loads(line)
        for line in (session_dir / "decision_log.jsonl")
        .read_text()
        .splitlines()
    ]


def _patch_runtime_and_submit(monkeypatch, submit_calls: list[dict]) -> None:
    import chemsmart.agent.tools as agent_tools

    def fake_validate_runtime(job, server=None):
        return {
            "ok": "ok",
            "local_ok": True,
            "local_issues": [],
            "remote_unknown": [],
        }

    def fake_submit_hpc(job, server=None, transport=None, execute=False):
        submit_calls.append(
            {
                "label": job.label,
                "folder": job.folder,
                "server": server,
                "transport": transport,
                "execute": execute,
            }
        )
        return {
            "transport": "FakeTransport",
            "script_path": str(Path(job.folder) / "submit.sh"),
            "script_bytes": None,
            "command_executed": None,
            "job_id": "job-123" if execute else None,
            "duplicate_check": {"duplicate": False, "message": None},
        }

    monkeypatch.setattr(agent_tools, "validate_runtime", fake_validate_runtime)
    monkeypatch.setattr(agent_tools, "submit_hpc", fake_submit_hpc)


def test_orca_dry_submit_skips_submit_hpc(
    monkeypatch,
    single_molecule_xyz_file,
    tmp_path: Path,
):
    monkeypatch.chdir(tmp_path)
    submit_calls: list[dict] = []
    _patch_runtime_and_submit(monkeypatch, submit_calls)
    provider = FakeProvider(
        [_orca_submit_plan(single_molecule_xyz_file), critic_ok()]
    )
    session = AgentSession(
        provider=provider,
        registry=ToolRegistry.default(),
        session_root=tmp_path / "sessions",
    )

    result = session.run("ORCA dry-submit regression", dry_submit=True)

    assert result["blocked"] is False
    assert submit_calls == []
    assert result["completed_steps"] == 5
    assert Path(result["dry_run_result"]["inputfile"]).suffix == ".inp"
    entries = _read_decision_log(Path(result["session_dir"]))
    assert not any(
        entry["kind"] in {"tool_call", "tool_preview"}
        and entry["payload"]["tool"] == "submit_hpc"
        for entry in entries
    )
    assert any(
        entry["kind"] == "tool_skipped"
        and entry["payload"]["tool"] == "submit_hpc"
        for entry in entries
    )


def test_execute_actually_calls_submit_hpc(
    monkeypatch,
    single_molecule_xyz_file,
    tmp_path: Path,
):
    monkeypatch.chdir(tmp_path)
    submit_calls: list[dict] = []
    _patch_runtime_and_submit(monkeypatch, submit_calls)
    provider = FakeProvider(
        [_orca_submit_plan(single_molecule_xyz_file), critic_ok()]
    )
    session = AgentSession(
        provider=provider,
        registry=ToolRegistry.default(),
        session_root=tmp_path / "sessions",
    )

    result = session.run("ORCA execute regression", dry_submit=False)

    assert result["blocked"] is False
    assert any(call["execute"] for call in submit_calls)
    assert result["completed_steps"] == 6


def test_dry_submit_resume_execute_calls_submit_hpc(
    monkeypatch,
    single_molecule_xyz_file,
    tmp_path: Path,
):
    monkeypatch.chdir(tmp_path)
    submit_calls: list[dict] = []
    _patch_runtime_and_submit(monkeypatch, submit_calls)
    provider = FakeProvider(
        [_orca_submit_plan(single_molecule_xyz_file), critic_ok()]
    )
    session_root = tmp_path / "sessions"
    session = AgentSession(
        provider=provider,
        registry=ToolRegistry.default(),
        session_root=session_root,
    )

    first = session.run("ORCA resume regression", dry_submit=True)
    resumed = AgentSession.resume(
        first["session_id"],
        provider=provider,
        registry=ToolRegistry.default(),
        session_root=session_root,
        dry_submit=False,
    )

    assert first["blocked"] is False
    assert resumed["blocked"] is False
    assert any(call["execute"] for call in submit_calls)
    assert resumed["completed_steps"] == 6
