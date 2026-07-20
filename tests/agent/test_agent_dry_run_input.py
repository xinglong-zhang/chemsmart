from __future__ import annotations

import shlex
from pathlib import Path

from chemsmart.agent.harness.workflow_state import (
    select_workspace_project,
    workflow_state_scope,
)
from chemsmart.agent.tools import build_job, build_molecule, dry_run_input
from chemsmart.settings.gaussian import GaussianProjectSettings


class TestDryRunInput:
    def test_gaussian_dry_run_matches_golden_output(
        self,
        tmp_path: Path,
        single_molecule_xyz_file,
        gaussian_yaml_settings_gas_solv_project_name,
        gaussian_jobrunner_no_scratch,
        gaussian_written_opt_file,
    ):
        gaussian_jobrunner_no_scratch.scratch_dir = str(tmp_path)
        settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_gas_solv_project_name
        ).opt_settings()
        settings.charge = 0
        settings.multiplicity = 1
        molecule = build_molecule(single_molecule_xyz_file)
        job = build_job(
            "gaussian.opt",
            molecule=molecule,
            settings=settings,
            label="gaussian_opt",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        job.set_folder(str(tmp_path))

        result = dry_run_input(job)

        assert result["content"] == Path(gaussian_written_opt_file).read_text()

    def test_dry_run_is_idempotent_and_returns_absolute_path(
        self,
        tmp_path: Path,
        single_molecule_xyz_file,
        gaussian_yaml_settings_gas_solv_project_name,
        gaussian_jobrunner_no_scratch,
    ):
        gaussian_jobrunner_no_scratch.scratch_dir = str(tmp_path)
        settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_gas_solv_project_name
        ).opt_settings()
        settings.charge = 0
        settings.multiplicity = 1
        molecule = build_molecule(single_molecule_xyz_file)
        job = build_job(
            "gaussian.opt",
            molecule=molecule,
            settings=settings,
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        job.set_folder(str(tmp_path))

        first_result = dry_run_input(job)
        second_result = dry_run_input(job)

        assert first_result["content"] == second_result["content"]
        assert Path(first_result["inputfile"]).is_absolute()
        assert first_result["cli_grounded"] is True
        assert first_result["command"].startswith("chemsmart run gaussian ")

    def test_dry_run_only_touches_tmp_path(
        self,
        tmp_path: Path,
        single_molecule_xyz_file,
        gaussian_yaml_settings_gas_solv_project_name,
        gaussian_jobrunner_no_scratch,
        monkeypatch,
    ):
        monkeypatch.chdir(tmp_path)
        gaussian_jobrunner_no_scratch.scratch_dir = str(tmp_path)
        settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_gas_solv_project_name
        ).opt_settings()
        settings.charge = 0
        settings.multiplicity = 1
        molecule = build_molecule(single_molecule_xyz_file)
        job = build_job(
            "gaussian.opt",
            molecule=molecule,
            settings=settings,
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        job.set_folder(str(tmp_path))

        result = dry_run_input(job)

        created_paths = {
            path.relative_to(tmp_path).as_posix()
            for path in tmp_path.rglob("*")
            if path.is_file()
        }
        assert Path(result["inputfile"]).parent == tmp_path
        assert created_paths == {Path(result["inputfile"]).name}

    def test_gaussian_scan_dry_run_includes_cli_grounding_command(
        self,
        tmp_path: Path,
        single_molecule_xyz_file,
    ):
        from chemsmart.agent.tools import build_gaussian_settings

        molecule = build_molecule(single_molecule_xyz_file)
        settings = build_gaussian_settings(
            "B3LYP",
            "6-31G(d)",
            scan_definition="B 1 2 S 10 0.05\nB 1 3 F",
        )
        job = build_job(
            "gaussian.scan",
            molecule=molecule,
            settings=settings,
            label="h2o_scan",
        )
        job.set_folder(str(tmp_path))

        result = dry_run_input(job)

        assert result["cli_grounded"] is True
        assert result["command"] == (
            "chemsmart run gaussian -c 0 -m 1 -x B3LYP -b '6-31G(d)' "
            # The builder shell-quotes each argument, which is visible only
            # when the fixture path contains characters a POSIX shell would
            # treat specially - as a Windows path's backslashes do.
            f"-f {shlex.quote(str(single_molecule_xyz_file))} "
            "-l h2o_scan scan --coordinates "
            "'[[1,2]]' --num-steps 10 --step-size 0.05 "
            "--constrained-coordinates '[[1,3]]'"
        )

    def test_gaussian_scan_repairs_route_constraint_and_preserves_project(
        self,
        tmp_path: Path,
        single_molecule_xyz_file,
        monkeypatch,
    ):
        from chemsmart.agent.tools import build_gaussian_settings

        monkeypatch.chdir(tmp_path)
        project_path = tmp_path / ".chemsmart" / "gaussian" / "water_demo.yaml"
        project_path.parent.mkdir(parents=True)
        project_path.write_text(
            "gas:\n  functional: b3lyp empiricaldispersion=gd3bj\n"
            "  basis: def2svp\n  freq: true\n",
            encoding="utf-8",
        )
        with workflow_state_scope("project-scan", cwd=tmp_path):
            selected = select_workspace_project(
                "water_demo", "gaussian", cwd=tmp_path
            )
            assert selected["selected"] is True
            molecule = build_molecule(single_molecule_xyz_file)
            settings = build_gaussian_settings(
                "b3lyp empiricaldispersion=gd3bj",
                "def2svp",
                scan_definition="B 1 2 S 10 0.05",
                additional_opt_options_in_route="B 1 3 F",
            )
            job = build_job(
                "gaussian.scan",
                molecule=molecule,
                settings=settings,
                label="h2o_bond_scan",
            )
            job.set_folder(str(tmp_path))
            result = dry_run_input(job)

        route = next(
            line.strip()
            for line in result["content"].splitlines()
            if line.lstrip().startswith("#")
        )
        assert "B,1,3,F" not in route
        assert "B 1 2 S 10 0.05" in result["content"]
        assert "B 1 3 F" in result["content"]
        assert result["command"] == (
            "chemsmart run gaussian -p water_demo -c 0 -m 1 "
            f"-f {shlex.quote(str(single_molecule_xyz_file))} "
            "-l h2o_bond_scan scan "
            "--coordinates '[[1,2]]' --num-steps 10 --step-size 0.05 "
            "--constrained-coordinates '[[1,3]]'"
        )
