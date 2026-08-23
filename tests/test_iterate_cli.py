"""
Direct tests for the ``chemsmart run iterate`` CLI group in
``chemsmart.cli.iterate.iterate``.

Mocks out ``IterateJob``/``IterateJobRunner``/``generate_template``/
``validate_config`` so these tests exercise only the CLI argument
parsing, validation, and result-callback wiring -- not the underlying
job/optimization pipeline (covered separately by the ``jobs.iterate``
unit tests).
"""

import importlib
from unittest.mock import MagicMock

import pytest
from click.testing import CliRunner

from chemsmart.cli.iterate.iterate import iterate

# The `chemsmart.cli.iterate` package's __init__.py does
# ``from .iterate import iterate``, which rebinds the package attribute
# ``iterate`` to the Click group object rather than the submodule.
# String-based ``mocker.patch("chemsmart.cli.iterate.iterate.X")`` would
# therefore resolve to the Click group, not the submodule, and raise
# AttributeError. Look the submodule up via ``sys.modules`` (through
# ``importlib.import_module``) and patch attributes on that object.
iterate_module = importlib.import_module("chemsmart.cli.iterate.iterate")


@pytest.fixture()
def valid_config_toml(tmp_path):
    path = tmp_path / "config.toml"
    path.write_text("[[skeletons]]\nfile = 'a.xyz'\n")
    return str(path)


def _patch_job(mocker, run_return="out.xyz"):
    mock_job = MagicMock()
    mock_job.run.return_value = run_return
    mocker.patch.object(iterate_module, "IterateJob", return_value=mock_job)
    mocker.patch.object(iterate_module, "IterateJobRunner")
    mocker.patch.object(iterate_module, "IterateJobSettings")
    mocker.patch.object(
        iterate_module,
        "validate_config",
        return_value={"skeletons": [1], "substituents": [1]},
    )
    return mock_job


class TestGenerateTemplateOption:
    def test_generate_template_writes_file_and_exits(self, mocker, tmp_path):
        out_path = str(tmp_path / "template.toml")
        mock_generate = mocker.patch.object(
            iterate_module, "generate_template", return_value=out_path
        )
        runner = CliRunner()
        result = runner.invoke(iterate, ["-g", out_path])
        assert result.exit_code == 0
        assert f"Generated template: {out_path}" in result.output
        mock_generate.assert_called_once_with(out_path, overwrite=False)


class TestOutputModeValidation:
    def test_outputfile_forbidden_with_separate_outputs(
        self, mocker, valid_config_toml
    ):
        _patch_job(mocker)
        runner = CliRunner()
        result = runner.invoke(
            iterate,
            [
                "-f",
                valid_config_toml,
                "--separate-outputs",
                "-o",
                "custom_out",
            ],
        )
        assert result.exit_code != 0
        assert "not allowed when '--separate-outputs'" in result.output

    def test_directory_defaults_to_cwd_with_separate_outputs(
        self, mocker, valid_config_toml
    ):
        mock_job = _patch_job(mocker)
        runner = CliRunner()
        result = runner.invoke(
            iterate,
            ["-f", valid_config_toml, "--separate-outputs"],
        )
        assert result.exit_code == 0
        mock_job.run.assert_called_once()

    def test_directory_explicitly_given_with_separate_outputs(
        self, mocker, valid_config_toml, tmp_path
    ):
        mock_job = _patch_job(mocker)
        runner = CliRunner()
        result = runner.invoke(
            iterate,
            [
                "-f",
                valid_config_toml,
                "--separate-outputs",
                "-d",
                str(tmp_path),
            ],
        )
        assert result.exit_code == 0
        mock_job.run.assert_called_once()

    def test_directory_forbidden_without_separate_outputs(
        self, mocker, valid_config_toml, tmp_path
    ):
        _patch_job(mocker)
        runner = CliRunner()
        result = runner.invoke(
            iterate,
            ["-f", valid_config_toml, "-d", str(tmp_path)],
        )
        assert result.exit_code != 0
        assert "not allowed when '--no-separate-outputs'" in result.output


class TestFilenameValidation:
    def test_missing_filename_raises(self, mocker):
        _patch_job(mocker)
        runner = CliRunner()
        result = runner.invoke(iterate, [])
        assert result.exit_code != 0
        assert "configuration file is required" in result.output

    def test_nonexistent_file_raises(self, mocker, tmp_path):
        _patch_job(mocker)
        runner = CliRunner()
        missing = str(tmp_path / "missing.toml")
        result = runner.invoke(iterate, ["-f", missing])
        assert result.exit_code != 0
        assert "does not exist" in result.output

    def test_non_toml_extension_raises(self, mocker, tmp_path):
        _patch_job(mocker)
        bad_file = tmp_path / "config.txt"
        bad_file.write_text("not toml")
        runner = CliRunner()
        result = runner.invoke(iterate, ["-f", str(bad_file)])
        assert result.exit_code != 0
        assert "must be a configuration file" in result.output


class TestEmptyConfigFile:
    def test_empty_toml_file_treated_as_empty_dict(self, mocker, tmp_path):
        _patch_job(mocker)
        mock_validate = mocker.patch.object(
            iterate_module,
            "validate_config",
            return_value={"skeletons": [], "substituents": []},
        )
        empty_file = tmp_path / "empty.toml"
        empty_file.write_text("")
        runner = CliRunner()
        result = runner.invoke(iterate, ["-f", str(empty_file)])
        assert result.exit_code == 0
        mock_validate.assert_called_once()
        called_raw_config = mock_validate.call_args[0][0]
        assert called_raw_config == {}

    def test_none_unwrap_result_becomes_empty_dict(
        self, mocker, valid_config_toml
    ):
        # tomlkit's TOMLDocument.unwrap() always returns a dict, never
        # None, even for an empty file -- but the CLI defensively
        # normalizes a None result to {} in case that ever changes (or
        # a differently-behaved TOML backend is swapped in). Directly
        # mock tomlkit.load(...).unwrap() to return None to exercise
        # that normalization branch.
        _patch_job(mocker)
        mock_validate = mocker.patch.object(
            iterate_module,
            "validate_config",
            return_value={"skeletons": [], "substituents": []},
        )
        mock_doc = MagicMock()
        mock_doc.unwrap.return_value = None
        mocker.patch.object(
            iterate_module.tomlkit, "load", return_value=mock_doc
        )
        runner = CliRunner()
        result = runner.invoke(iterate, ["-f", valid_config_toml])
        assert result.exit_code == 0
        called_raw_config = mock_validate.call_args[0][0]
        assert called_raw_config == {}


class TestResultCallback:
    def test_successful_run_reports_output_path(
        self, mocker, valid_config_toml, tmp_path
    ):
        out_path = tmp_path / "result.xyz"
        out_path.write_text("dummy")
        _patch_job(mocker, run_return=str(out_path))
        runner = CliRunner()
        result = runner.invoke(iterate, ["-f", valid_config_toml])
        assert result.exit_code == 0
        assert "Output saved to:" in result.output

    def test_no_structures_generated_warns(self, mocker, valid_config_toml):
        _patch_job(mocker, run_return=None)
        runner = CliRunner()
        result = runner.invoke(iterate, ["-f", valid_config_toml])
        assert result.exit_code == 0

    def test_job_run_exception_becomes_click_exception(
        self, mocker, valid_config_toml
    ):
        mock_job = _patch_job(mocker)
        mock_job.run.side_effect = RuntimeError("boom")
        runner = CliRunner()
        result = runner.invoke(iterate, ["-f", valid_config_toml])
        assert result.exit_code != 0
        assert "boom" in str(result.exception) or "boom" in result.output
