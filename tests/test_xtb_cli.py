import os
from pathlib import Path
from unittest.mock import MagicMock, patch

from click.testing import CliRunner

from chemsmart.cli.run import run
from chemsmart.cli.sub import sub
from chemsmart.cli.xtb.xtb import xtb
from chemsmart.jobs.xtb.hess import XTBHessJob
from chemsmart.jobs.xtb.opt import XTBOptJob
from chemsmart.jobs.xtb.singlepoint import XTBSinglePointJob


def _invoke_xtb_and_capture_settings(job_path, args):
    runner = CliRunner()
    mock_job = MagicMock()
    with patch(job_path, return_value=mock_job) as mock_job_cls:
        result = runner.invoke(
            xtb,
            args,
            obj={"jobrunner": MagicMock()},
            catch_exceptions=False,
        )
        settings = mock_job_cls.call_args[1]["settings"]
    return result, settings, mock_job_cls


class TestXTBHelp:
    def test_run_help_lists_xtb(self):
        result = CliRunner().invoke(run, ["--help"])
        assert result.exit_code == 0, result.output
        assert "xtb" in result.output

    def test_sub_help_lists_xtb(self):
        result = CliRunner().invoke(sub, ["--help"])
        assert result.exit_code == 0, result.output
        assert "xtb" in result.output

    def test_xtb_help_lists_jobs(self):
        result = CliRunner().invoke(run, ["xtb", "--help"])
        assert result.exit_code == 0, result.output
        assert "opt" in result.output
        assert "sp" in result.output
        assert "hess" in result.output

    def test_xtb_job_help(self):
        for jobtype in ("opt", "sp", "hess"):
            result = CliRunner().invoke(run, ["xtb", jobtype, "--help"])
            assert result.exit_code == 0, result.output
            assert "--skip-completed" in result.output


class TestXTBCLISettings:
    def test_sp_settings_merge(self, single_molecule_xyz_file):
        result, settings, mock_job_cls = _invoke_xtb_and_capture_settings(
            "chemsmart.jobs.xtb.singlepoint.XTBSinglePointJob",
            [
                "-p",
                "test",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "-1",
                "-m",
                "2",
                "-g",
                "gfn1",
                "-sm",
                "alpb",
                "-si",
                "water",
                "--grad",
                "sp",
            ],
        )
        assert result.exit_code == 0, result.output
        assert mock_job_cls.call_args[0] == ()
        assert settings.jobtype == "sp"
        assert settings.charge == -1
        assert settings.multiplicity == 2
        assert settings.gfn_version == "gfn1"
        assert settings.solvent_model == "alpb"
        assert settings.solvent_id == "water"
        assert settings.grad is True

    def test_opt_instantiates_xtb_opt_job(self, single_molecule_xyz_file):
        result, settings, _ = _invoke_xtb_and_capture_settings(
            "chemsmart.jobs.xtb.opt.XTBOptJob",
            [
                "-p",
                "test",
                "-f",
                single_molecule_xyz_file,
                "opt",
                "--optimization-level",
                "loose",
            ],
        )
        assert result.exit_code == 0, result.output
        assert settings.jobtype == "opt"
        assert settings.optimization_level == "loose"

    def test_hess_instantiates_xtb_hess_job(self, single_molecule_xyz_file):
        result, settings, _ = _invoke_xtb_and_capture_settings(
            "chemsmart.jobs.xtb.hess.XTBHessJob",
            ["-p", "test", "-f", single_molecule_xyz_file, "hess"],
        )
        assert result.exit_code == 0, result.output
        assert settings.jobtype == "hess"


class TestXTBJobClasses:
    def test_xtb_job_types_are_distinct(self):
        assert XTBOptJob.TYPE == "xtbopt"
        assert XTBSinglePointJob.TYPE == "xtbsp"
        assert XTBHessJob.TYPE == "xtbhess"


class TestXTBProjectOptional:
    """xTB needs no project YAML: defaults come from XTBJobSettings."""

    def test_from_project_none_returns_gfn2_defaults(self):
        from chemsmart.settings.xtb import XTBProjectSettings

        settings = XTBProjectSettings.from_project(None).opt_settings()

        assert settings.gfn_version == "gfn2"
        assert settings.jobtype == "opt"

    def test_from_project_missing_name_still_raises(self):
        import pytest

        from chemsmart.settings.xtb import XTBProjectSettings

        with pytest.raises(FileNotFoundError):
            XTBProjectSettings.from_project("does-not-exist")

    def test_run_without_project_uses_defaults(self, single_molecule_xyz_file):
        # No -p on the command line: the opt leaf must still build a job from
        # GFN2 defaults rather than failing to resolve a project.
        result, settings, _cls = _invoke_xtb_and_capture_settings(
            "chemsmart.jobs.xtb.opt.XTBOptJob",
            ["-f", str(single_molecule_xyz_file), "opt"],
        )

        assert result.exit_code == 0, result.output
        assert settings.gfn_version == "gfn2"
        assert settings.jobtype == "opt"


class TestXTBSubmission:
    def test_sub_test_writes_xtb_scripts(
        self, tmp_path, server_yaml_file, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        xyz_file = (
            Path(__file__).resolve().parents[1]
            / "examples"
            / "xtb"
            / "water.xyz"
        )
        result = CliRunner().invoke(
            sub,
            [
                "-s",
                server_yaml_file,
                "--test",
                "xtb",
                "-p",
                "test",
                "-f",
                str(xyz_file),
                "sp",
            ],
            catch_exceptions=False,
        )
        assert result.exit_code == 0, result.output
        run_script = tmp_path / "chemsmart_run_water_sp.py"
        submit_script = tmp_path / "chemsmart_sub_water_sp.sh"
        assert run_script.exists()
        assert submit_script.exists()
        assert "'--test'" not in run_script.read_text()
        assert "'xtb'" in run_script.read_text()
        assert (
            "conda activate ~/anaconda3/envs/chemsmart"
            in submit_script.read_text()
        )


class TestXTBExamples:
    def test_example_xyz_files_exist(self):
        examples_dir = Path(__file__).resolve().parents[1] / "examples" / "xtb"
        for filename in ("water.xyz", "methane.xyz"):
            path = examples_dir / filename
            assert path.exists()
            assert os.path.getsize(path) > 0
