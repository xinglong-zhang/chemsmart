from __future__ import annotations

from pathlib import Path

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
