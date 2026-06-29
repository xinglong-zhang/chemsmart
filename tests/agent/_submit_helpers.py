from __future__ import annotations

import contextlib
import os
from pathlib import Path

import yaml
from click.testing import CliRunner

from chemsmart.agent.tools import (
    build_gaussian_settings,
    build_job,
    build_molecule,
)
from chemsmart.cli.sub import sub as sub_cli
from chemsmart.settings.server import Server
from chemsmart.settings.user import ChemsmartUserSettings


def patch_user_config(monkeypatch, config_dir: Path) -> None:
    import chemsmart.jobs.runner as runner_module
    import chemsmart.settings.executable as executable_module
    import chemsmart.settings.server as server_module
    import chemsmart.settings.submitters as submitters_module

    monkeypatch.setattr(
        ChemsmartUserSettings,
        "USER_CONFIG_DIR",
        str(config_dir),
    )
    runner_module.user_settings = ChemsmartUserSettings()
    executable_module.user_settings = ChemsmartUserSettings()
    server_module.user_settings = ChemsmartUserSettings()
    submitters_module.user_settings = ChemsmartUserSettings()


def write_server_yaml(
    config_dir: Path,
    server_name: str,
    executable_dir: Path,
) -> Path:
    server_dir = config_dir / "server"
    server_dir.mkdir(parents=True, exist_ok=True)
    server_file = server_dir / f"{server_name}.yaml"
    with open(server_file, "w") as file:
        yaml.safe_dump(
            {
                "SERVER": {
                    "SCHEDULER": "PBS",
                    "QUEUE_NAME": "normal",
                    "NUM_HOURS": 24,
                    "MEM_GB": 32,
                    "NUM_CORES": 8,
                    "NUM_GPUS": 0,
                    "NUM_THREADS": 8,
                    "SUBMIT_COMMAND": "qsub",
                    "PROJECT": "proj123",
                    "SCRATCH_DIR": str(config_dir / "scratch"),
                },
                "GAUSSIAN": {
                    "EXEFOLDER": str(executable_dir),
                    "LOCAL_RUN": False,
                },
            },
            file,
            sort_keys=False,
        )
    return server_file


def write_gaussian_project_yaml(
    config_dir: Path, project_name="agent_contract"
) -> Path:
    gaussian_dir = config_dir / "gaussian"
    gaussian_dir.mkdir(parents=True, exist_ok=True)
    defaults_file = gaussian_dir / "defaults.yaml"
    defaults_file.write_text(
        "\n".join(
            [
                "chk: True",
                "freq: True",
                "jobtype: Null",
                "title: 'Gaussian job'",
                "ab_initio: Null",
                "functional: Null",
                "basis: Null",
                "semiempirical: Null",
                "charge: Null",
                "multiplicity: Null",
                "route_to_be_written: Null",
                "solvent_model: Null",
                "solvent_id: Null",
                "dieze_tag: Null",
                "additional_opt_options_in_route: Null",
                "additional_route_parameters: Null",
                "modred: Null",
                "gen_genecp: Null",
                "heavy_elements: Null",
                "heavy_elements_basis: Null",
                "light_elements_basis: Null",
                "custom_solvent: Null",
                "append_additional_info: Null",
                "forces: False",
            ]
        )
        + "\n"
    )
    project_file = gaussian_dir / f"{project_name}.yaml"
    with open(project_file, "w") as file:
        yaml.safe_dump(
            {
                "gas": {
                    "functional": "B3LYP",
                    "basis": "6-31G(d)",
                },
                "solv": {
                    "functional": "B3LYP",
                    "basis": "6-31G(d)",
                },
            },
            file,
            sort_keys=False,
        )
    return project_file


def build_submit_server(
    monkeypatch, tmp_path: Path, server_name="wave5-server"
):
    config_dir = tmp_path / ".chemsmart"
    executable_dir = tmp_path / "bin"
    scratch_dir = config_dir / "scratch"
    executable_dir.mkdir(parents=True)
    scratch_dir.mkdir(parents=True)
    (executable_dir / "g16").write_text("#!/bin/sh\nexit 0\n")
    (executable_dir / "g16").chmod(0o755)
    patch_user_config(monkeypatch, config_dir)
    monkeypatch.setattr(
        Server, "_check_running_jobs", staticmethod(lambda job: None)
    )
    write_server_yaml(config_dir, server_name, executable_dir)
    write_gaussian_project_yaml(config_dir)
    return Server.from_servername(server_name)


def build_gaussian_submit_job(
    single_molecule_xyz_file,
    tmp_path: Path,
    *,
    kind: str,
    label: str,
    functional="B3LYP",
    basis="6-31G(d)",
    freq: bool = False,
):
    molecule = build_molecule(single_molecule_xyz_file)
    settings = build_gaussian_settings(functional, basis)
    settings.freq = freq
    job = build_job(
        kind,
        molecule=molecule,
        settings=settings,
        label=label,
    )
    job.set_folder(str(tmp_path))
    return job


def read_submit_script_bytes(job, server: Server) -> bytes:
    submitter = server.get_submitter(job)
    script_path = Path(job.folder) / submitter.submit_script
    return script_path.read_bytes()


def invoke_manual_submit_cli(
    server_name: str,
    job_folder: Path,
    inputfile: str,
    label: str,
    jobtype: str,
):
    runner = CliRunner()
    args = [
        "-s",
        server_name,
        "--test",
        "gaussian",
        "-p",
        "agent_contract",
        "-f",
        inputfile,
        "-l",
        label,
        jobtype,
    ]
    with pushd(job_folder):
        result = runner.invoke(sub_cli, args, catch_exceptions=False)
    assert result.exit_code == 0, result.output
    return result


@contextlib.contextmanager
def pushd(path: Path):
    cwd = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(cwd)
