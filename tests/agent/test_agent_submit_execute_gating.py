from __future__ import annotations

from pathlib import Path

from chemsmart.agent import MockTransport, submit_hpc

from ._submit_helpers import build_gaussian_submit_job, build_submit_server


def test_submit_hpc_execute_false_forces_local_dry_run_transport(
    monkeypatch,
    tmp_path: Path,
    single_molecule_xyz_file,
):
    server = build_submit_server(
        monkeypatch, tmp_path, server_name="gate-server"
    )
    job = build_gaussian_submit_job(
        single_molecule_xyz_file,
        tmp_path / "gate-false",
        kind="gaussian.opt",
        label="gate_false",
    )
    transport = MockTransport()

    result = submit_hpc(
        job,
        server=server,
        transport=transport,
        execute=False,
    )

    assert result["transport"] == "LocalDryRunTransport"
    assert result["job_id"] is None
    assert transport.calls == []


def test_submit_hpc_execute_true_uses_provided_transport(
    monkeypatch,
    tmp_path: Path,
    single_molecule_xyz_file,
):
    server = build_submit_server(
        monkeypatch, tmp_path, server_name="gate-server"
    )
    job = build_gaussian_submit_job(
        single_molecule_xyz_file,
        tmp_path / "gate-true",
        kind="gaussian.opt",
        label="gate_true",
    )
    transport = MockTransport()

    result = submit_hpc(
        job,
        server=server,
        transport=transport,
        execute=True,
    )

    assert result["transport"] == "MockTransport"
    assert result["job_id"] == "mock-job-0001"
    assert len(transport.calls) == 1
    assert transport.calls[0]["working_dir"] == str(Path(job.folder).resolve())


def test_submit_hpc_uses_sole_configured_server_when_omitted(
    monkeypatch,
    tmp_path: Path,
    single_molecule_xyz_file,
):
    build_submit_server(monkeypatch, tmp_path, server_name="default-server")
    job = build_gaussian_submit_job(
        single_molecule_xyz_file,
        tmp_path / "default-server-job",
        kind="gaussian.opt",
        label="default_server_job",
    )

    result = submit_hpc(
        job,
        execute=False,
    )

    assert result["transport"] == "LocalDryRunTransport"
    assert result["script_path"] is not None
    assert Path(result["script_path"]).exists()
