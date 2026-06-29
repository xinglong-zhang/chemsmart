from __future__ import annotations

from pathlib import Path

import pytest

from chemsmart.agent import dry_run_input, submit_hpc

from ._submit_helpers import (
    build_gaussian_submit_job,
    build_submit_server,
    invoke_manual_submit_cli,
    read_submit_script_bytes,
)


@pytest.mark.parametrize(
    ("case_name", "kind", "jobtype", "freq"),
    [
        ("gaussian-opt", "gaussian.opt", "opt", False),
        ("gaussian-opt-freq", "gaussian.opt", "opt", True),
        ("gaussian-ts", "gaussian.ts", "ts", False),
    ],
)
def test_submit_hpc_matches_manual_sub_cli_submit_script_bytes(
    monkeypatch,
    tmp_path: Path,
    single_molecule_xyz_file,
    case_name: str,
    kind: str,
    jobtype: str,
    freq: bool,
):
    server = build_submit_server(monkeypatch, tmp_path, server_name="wave5")

    agent_job = build_gaussian_submit_job(
        single_molecule_xyz_file,
        tmp_path / f"{case_name}-agent",
        kind=kind,
        label=case_name,
        freq=freq,
    )
    manual_job = build_gaussian_submit_job(
        single_molecule_xyz_file,
        tmp_path / f"{case_name}-manual",
        kind=kind,
        label=case_name,
        freq=freq,
    )

    agent_result = submit_hpc(agent_job, server=server, execute=False)
    manual_input = dry_run_input(manual_job)
    invoke_manual_submit_cli(
        server_name="wave5",
        job_folder=Path(manual_job.folder),
        inputfile=manual_input["inputfile"],
        label=manual_job.label,
        jobtype=jobtype,
    )

    manual_bytes = read_submit_script_bytes(manual_job, server)

    assert agent_result["transport"] == "LocalDryRunTransport"
    assert agent_result["duplicate_check"] == {
        "duplicate": False,
        "message": None,
    }
    assert isinstance(agent_result["script_bytes"], bytes)
    assert agent_result["script_bytes"] == manual_bytes
