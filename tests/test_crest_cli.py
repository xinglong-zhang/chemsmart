from pathlib import Path
from unittest.mock import MagicMock, patch

from click.testing import CliRunner

from chemsmart.cli.crest.crest import crest
from chemsmart.cli.run import run
from chemsmart.cli.sub import sub
from chemsmart.jobs.crest import CRESTConformerSearchJob


def _invoke_crest_and_capture_settings(job_path, args):
    runner = CliRunner()
    mock_job = MagicMock()
    with patch(job_path, return_value=mock_job) as mock_job_cls:
        result = runner.invoke(
            crest,
            args,
            obj={"jobrunner": MagicMock()},
            catch_exceptions=False,
        )
        settings = mock_job_cls.call_args[1]["settings"]
    return result, settings, mock_job_cls


class TestCRESTHelp:
    def test_run_help_lists_crest(self):
        result = CliRunner().invoke(run, ["--help"])
        assert result.exit_code == 0, result.output
        assert "crest" in result.output

    def test_sub_help_lists_crest(self):
        result = CliRunner().invoke(sub, ["--help"])
        assert result.exit_code == 0, result.output
        assert "crest" in result.output

    def test_crest_help_lists_subcommands(self):
        result = CliRunner().invoke(run, ["crest", "--help"])
        assert result.exit_code == 0, result.output
        assert "conformers" in result.output


class TestCRESTCLISettings:
    def test_conformers_settings_merge(self, single_molecule_xyz_file):
        result, settings, mock_job_cls = _invoke_crest_and_capture_settings(
            "chemsmart.jobs.crest.conformers.CRESTConformerSearchJob",
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
                "-O",
                "loose",
                "-sm",
                "alpb",
                "-si",
                "water",
                "--ewin",
                "8.0",
                "conformers",
            ],
        )
        assert result.exit_code == 0, result.output
        assert settings.jobtype == "conformers"
        assert settings.charge == -1
        assert settings.multiplicity == 2
        assert settings.gfn_version == "gfn1"
        assert settings.optimization_level == "loose"
        assert settings.solvent_model == "alpb"
        assert settings.solvent_id == "water"
        assert settings.energy_window == 8.0

    def test_conformers_default_settings(self, single_molecule_xyz_file):
        result, settings, _ = _invoke_crest_and_capture_settings(
            "chemsmart.jobs.crest.conformers.CRESTConformerSearchJob",
            [
                "-p",
                "test",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "conformers",
            ],
        )
        assert result.exit_code == 0, result.output
        assert settings.jobtype == "conformers"
        assert settings.gfn_version == "gfn2"
        assert settings.constraints is None

    def test_constrained_settings(self, single_molecule_xyz_file):
        result, settings, _ = _invoke_crest_and_capture_settings(
            "chemsmart.jobs.crest.conformers.CRESTConformerSearchJob",
            [
                "-p",
                "test",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "conformers",
                "-c",
                "[[1,2],[1,3]]",
                "--force-constant",
                "0.5",
            ],
        )
        assert result.exit_code == 0, result.output
        assert settings.constraints == [[1, 2], [1, 3]]
        assert settings.force_constant == 0.5


class TestCRESTJobClasses:
    def test_crest_job_types_are_distinct(self):
        assert CRESTConformerSearchJob.TYPE == "crestconformers"


class TestCRESTSubmission:
    def test_sub_test_writes_crest_scripts(
        self, tmp_path, server_yaml_file, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        xyz_file = (
            Path(__file__).resolve().parent
            / "data"
            / "StructuresTests"
            / "xyz"
            / "1-mer.xyz"
        )
        result = CliRunner().invoke(
            sub,
            [
                "-s",
                server_yaml_file,
                "--test",
                "crest",
                "-p",
                "test",
                "-f",
                str(xyz_file),
                "-c",
                "0",
                "-m",
                "1",
                "conformers",
            ],
            catch_exceptions=False,
        )
        assert result.exit_code == 0, result.output
        run_script = tmp_path / "chemsmart_run_1-mer_conformers.py"
        submit_script = tmp_path / "chemsmart_sub_1-mer_conformers.sh"
        assert run_script.exists()
        assert submit_script.exists()
        assert "'--test'" not in run_script.read_text()
        assert "'crest'" in run_script.read_text()
        assert (
            "conda activate ~/anaconda3/envs/chemsmart"
            in submit_script.read_text()
        )
