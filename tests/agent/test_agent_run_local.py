from __future__ import annotations

from pathlib import Path

from chemsmart.agent.tools import (
    build_gaussian_settings,
    build_job,
    build_molecule,
    run_local,
)


def _build_job(single_molecule_xyz_file, tmp_path: Path):
    molecule = build_molecule(single_molecule_xyz_file)
    settings = build_gaussian_settings("B3LYP", "6-31G*")
    job = build_job(
        "gaussian.opt",
        molecule=molecule,
        settings=settings,
        label="agent_run_local",
    )
    job.set_folder(str(tmp_path))
    return job


def test_run_local_parses_gaussian_output_summary(
    monkeypatch,
    tmp_path: Path,
    single_molecule_xyz_file,
    gaussian_co2_opt_outfile,
):
    job = _build_job(single_molecule_xyz_file, tmp_path)

    def fake_run():
        Path(job.outputfile).write_text(
            Path(gaussian_co2_opt_outfile).read_text()
        )
        print("job finished")

    monkeypatch.setattr(job, "run", fake_run)

    result = run_local(job)

    assert result["ok"] is True
    assert result["returncode"] == 0
    assert Path(result["stdout_path"]).exists()
    assert Path(result["stderr_path"]).exists()
    assert result["output_summary"]["energy"] is not None
    assert result["output_summary"]["converged"] is True
    assert result["output_summary"]["imag_freqs"] == []
    assert result["output_summary"]["optimized_geometry_count"] > 0


def test_run_local_returns_failure_when_job_run_raises(
    monkeypatch,
    tmp_path: Path,
    single_molecule_xyz_file,
):
    job = _build_job(single_molecule_xyz_file, tmp_path)

    def fake_run():
        raise RuntimeError("boom")

    monkeypatch.setattr(job, "run", fake_run)

    result = run_local(job)

    assert result["ok"] is False
    assert result["returncode"] != 0
    assert Path(result["stderr_path"]).exists()
    assert "boom" in Path(result["stderr_path"]).read_text()


def test_run_local_returns_empty_summary_for_malformed_output(
    monkeypatch,
    tmp_path: Path,
    single_molecule_xyz_file,
):
    job = _build_job(single_molecule_xyz_file, tmp_path)

    def fake_run():
        Path(job.outputfile).write_text("not a valid gaussian output\n")

    monkeypatch.setattr(job, "run", fake_run)

    result = run_local(job)

    assert result["ok"] is True
    assert result["returncode"] == 0
    assert result["output_summary"] == {}
