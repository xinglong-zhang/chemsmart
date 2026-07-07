from __future__ import annotations

import os
from pathlib import Path

import yaml

from chemsmart.agent.tools import (
    build_gaussian_settings,
    build_job,
    validate_runtime,
)
from chemsmart.settings.server import Server
from chemsmart.settings.user import ChemsmartUserSettings


def _write_server_yaml(server_file: Path, config: dict) -> None:
    with open(server_file, "w") as file:
        yaml.safe_dump(config, file, sort_keys=False)


def _build_gaussian_job(single_molecule_xyz_file, tmp_path: Path):
    from chemsmart.agent.tools import build_molecule

    molecule = build_molecule(single_molecule_xyz_file)
    settings = build_gaussian_settings("B3LYP", "6-31G*")
    job = build_job(
        "gaussian.opt",
        molecule=molecule,
        settings=settings,
        label="agent_validate_runtime",
    )
    job.set_folder(str(tmp_path))
    return job


def _patch_user_config(monkeypatch, config_dir: Path) -> None:
    import chemsmart.jobs.runner as runner_module
    import chemsmart.settings.executable as executable_module
    import chemsmart.settings.server as server_module

    monkeypatch.setattr(
        ChemsmartUserSettings,
        "USER_CONFIG_DIR",
        str(config_dir),
    )
    runner_module.user_settings = ChemsmartUserSettings()
    executable_module.user_settings = ChemsmartUserSettings()
    server_module.user_settings = ChemsmartUserSettings()


def test_validate_runtime_fails_when_executable_path_missing(
    monkeypatch,
    tmp_path: Path,
    single_molecule_xyz_file,
):
    config_dir = tmp_path / ".chemsmart"
    server_dir = config_dir / "server"
    scratch_dir = tmp_path / "scratch"
    server_dir.mkdir(parents=True)
    scratch_dir.mkdir()
    _patch_user_config(monkeypatch, config_dir)

    server_file = server_dir / "missing-exe.yaml"
    _write_server_yaml(
        server_file,
        {
            "SERVER": {
                "SCHEDULER": "PBS",
                "QUEUE_NAME": "normal",
                "PROJECT": "proj123",
                "SCRATCH_DIR": str(scratch_dir),
            },
            "GAUSSIAN": {
                "EXEFOLDER": str(tmp_path / "does-not-exist"),
                "LOCAL_RUN": False,
            },
        },
    )

    job = _build_gaussian_job(single_molecule_xyz_file, tmp_path)
    result = validate_runtime(job, server=Server.from_yaml(str(server_file)))

    assert result["ok"] == "fail"
    assert result["local_ok"] is False
    assert "server.executable_path missing" in result["local_issues"]


def test_validate_runtime_modules_server_is_partial(
    monkeypatch,
    tmp_path: Path,
    single_molecule_xyz_file,
):
    config_dir = tmp_path / ".chemsmart"
    server_dir = config_dir / "server"
    scratch_dir = tmp_path / "scratch"
    server_dir.mkdir(parents=True)
    scratch_dir.mkdir()
    _patch_user_config(monkeypatch, config_dir)

    server_file = server_dir / "modules-only.yaml"
    _write_server_yaml(
        server_file,
        {
            "SERVER": {
                "SCHEDULER": "PBS",
                "QUEUE_NAME": "normal",
                "PROJECT": "proj123",
                "SCRATCH_DIR": str(scratch_dir),
            },
            "GAUSSIAN": {
                "EXEFOLDER": None,
                "LOCAL_RUN": False,
                "MODULES": "module load gaussian/g16\n",
            },
        },
    )

    job = _build_gaussian_job(single_molecule_xyz_file, tmp_path)
    result = validate_runtime(job, server=Server.from_yaml(str(server_file)))

    assert result["ok"] == "partial"
    assert result["local_ok"] is True
    assert result["local_issues"] == []
    assert "module load succeeds on HPC" in result["remote_unknown"]
    assert "scratch_dir writable on HPC" in result["remote_unknown"]


def test_validate_runtime_fails_when_scratch_dir_is_unresolved(
    monkeypatch,
    tmp_path: Path,
    single_molecule_xyz_file,
):
    config_dir = tmp_path / ".chemsmart"
    server_dir = config_dir / "server"
    bin_dir = tmp_path / "bin"
    server_dir.mkdir(parents=True)
    bin_dir.mkdir()
    _patch_user_config(monkeypatch, config_dir)

    fake_g16 = bin_dir / "g16"
    fake_g16.write_text("#!/bin/sh\nexit 0\n")
    fake_g16.chmod(0o755)
    monkeypatch.setenv(
        "PATH",
        f"{bin_dir}{os.pathsep}{os.environ.get('PATH', '')}",
    )

    server_file = server_dir / "unresolved-scratch.yaml"
    _write_server_yaml(
        server_file,
        {
            "SERVER": {
                "SCHEDULER": "PBS",
                "QUEUE_NAME": "normal",
                "PROJECT": "proj123",
                "SCRATCH_DIR": "$CHEMSMART_SCRATCH/runtime",
            },
            "GAUSSIAN": {
                "EXEFOLDER": str(bin_dir),
                "LOCAL_RUN": False,
            },
        },
    )

    job = _build_gaussian_job(single_molecule_xyz_file, tmp_path)
    result = validate_runtime(job, server=Server.from_yaml(str(server_file)))

    assert result["ok"] == "fail"
    assert result["local_ok"] is False
    assert "server.scratch_dir unresolved" in result["local_issues"]


def test_validate_runtime_without_server_is_partial(
    tmp_path: Path,
    single_molecule_xyz_file,
):
    job = _build_gaussian_job(single_molecule_xyz_file, tmp_path)

    result = validate_runtime(job, server=None)

    assert result["ok"] == "partial"
    assert result["local_ok"] is True
    assert result["local_issues"] == []
    assert "server.queue required" in result["remote_unknown"]
    assert "server.account required" in result["remote_unknown"]
    assert "server.scratch_dir required" in result["remote_unknown"]
    assert (
        "server.modules_or_executable_path required"
        in result["remote_unknown"]
    )
