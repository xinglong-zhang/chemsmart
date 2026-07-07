from __future__ import annotations

from pathlib import Path

from chemsmart.agent import MockTransport, submit_hpc
from chemsmart.settings.server import Server

from ._submit_helpers import build_gaussian_submit_job, build_submit_server


def test_submit_hpc_returns_duplicate_check_and_skips_transport(
    monkeypatch,
    tmp_path: Path,
    single_molecule_xyz_file,
):
    server = build_submit_server(
        monkeypatch, tmp_path, server_name="dup-server"
    )
    job = build_gaussian_submit_job(
        single_molecule_xyz_file,
        tmp_path / "duplicate-job",
        kind="gaussian.opt",
        label="duplicate_job",
    )
    transport = MockTransport()

    def fake_duplicate_check(submit_job):
        raise SystemExit(f"Duplicate job NOT submitted: {submit_job.label}")

    monkeypatch.setattr(Server, "_check_running_jobs", fake_duplicate_check)

    result = submit_hpc(
        job,
        server=server,
        transport=transport,
        execute=True,
    )

    assert result["script_path"] is None
    assert result["script_bytes"] is None
    assert result["command_executed"] is None
    assert result["job_id"] is None
    assert result["duplicate_check"]["duplicate"] is True
    assert (
        "Duplicate job NOT submitted" in result["duplicate_check"]["message"]
    )
    assert transport.calls == []
