"""
CLI tests for the ``thermochemistry`` group and its ``boltzmann``
subcommand (chemsmart/cli/thermochemistry/thermochemistry.py and
chemsmart/cli/thermochemistry/boltzmann.py).

Each test uses :class:`click.testing.CliRunner` to invoke the
``thermochemistry`` group and :mod:`unittest.mock` to intercept
``BoltzmannAverageThermochemistryJob`` (always appending the harmless
``boltzmann`` subcommand token) so the group's own directory/filename
parsing logic can be exercised without triggering the result callback's
real computation, and without writing any real .dat output files as a
side effect.
"""

import os
import shutil
from unittest.mock import MagicMock, patch

from click.testing import CliRunner

from chemsmart.cli.thermochemistry.thermochemistry import (
    resolve_entropy_cutoff,
    thermochemistry,
)
from chemsmart.jobs.thermochemistry.job import ThermochemistryJob

BOLTZMANN_JOB = "chemsmart.jobs.thermochemistry.boltzmann.BoltzmannAverageThermochemistryJob"


def run_boltzmann_and_capture_kwargs(cli_args):
    """Run the ``thermochemistry ... boltzmann`` CLI with a patched
    job class and capture the constructor call kwargs, plus whether
    compute_boltzmann_averages/show_results were invoked."""
    runner = CliRunner()
    with patch(BOLTZMANN_JOB) as mock_job_cls:
        mock_instance = MagicMock()
        mock_job_cls.return_value = mock_instance
        result = runner.invoke(
            thermochemistry,
            cli_args,
            obj={},
            catch_exceptions=True,
        )
        call_kwargs = None
        if mock_job_cls.call_args is not None:
            call_kwargs = mock_job_cls.call_args.kwargs
    return result, call_kwargs, mock_instance


class TestThermochemistryBoltzmannCommand:
    """See BUGS_FOUND.md #1: the CLI now passes a real `filename` to
    BoltzmannAverageThermochemistryJob (previously omitted entirely,
    which crashed on every invocation)."""

    def test_single_file_passes_it_as_filename(self, gaussian_co2_opt_outfile):
        result, kwargs, instance = run_boltzmann_and_capture_kwargs(
            ["-f", gaussian_co2_opt_outfile, "-T", "298.15", "boltzmann"]
        )
        assert result.exit_code == 0, result.output
        assert kwargs["files"] == (gaussian_co2_opt_outfile,)
        assert kwargs["filename"] == gaussian_co2_opt_outfile
        assert kwargs["energy_type"] == "gibbs"
        instance.compute_boltzmann_averages.assert_called_once()
        instance.show_results.assert_called_once()

    def test_multiple_files_uses_first_as_filename(
        self, gaussian_co2_opt_outfile
    ):
        result, kwargs, _ = run_boltzmann_and_capture_kwargs(
            [
                "-f",
                gaussian_co2_opt_outfile,
                "-f",
                gaussian_co2_opt_outfile,
                "-T",
                "298.15",
                "boltzmann",
            ]
        )
        assert result.exit_code == 0, result.output
        assert kwargs["files"] == (
            gaussian_co2_opt_outfile,
            gaussian_co2_opt_outfile,
        )
        assert kwargs["filename"] == gaussian_co2_opt_outfile

    def test_energy_type_electronic_option(self, gaussian_co2_opt_outfile):
        result, kwargs, _ = run_boltzmann_and_capture_kwargs(
            [
                "-f",
                gaussian_co2_opt_outfile,
                "-T",
                "298.15",
                "boltzmann",
                "-w",
                "electronic",
            ]
        )
        assert result.exit_code == 0, result.output
        assert kwargs["energy_type"] == "electronic"

    def test_rerun_completed_flag(self, gaussian_co2_opt_outfile):
        result, kwargs, _ = run_boltzmann_and_capture_kwargs(
            ["-f", gaussian_co2_opt_outfile, "-T", "298.15", "boltzmann", "-R"]
        )
        assert result.exit_code == 0, result.output
        assert kwargs["skip_completed"] is False

    def test_no_files_leaves_filename_none(self):
        """No -f given at all: files is None, so the fixed
        `filename=files[0] if files else None` guard falls back to
        None rather than raising IndexError."""
        result, kwargs, _ = run_boltzmann_and_capture_kwargs(
            ["-T", "298.15", "boltzmann"]
        )
        assert result.exit_code == 0, result.output
        assert kwargs["files"] is None
        assert kwargs["filename"] is None

    def test_no_files_real_job_class_raises(self):
        """Without the patched job class, a no-files invocation still
        surfaces the underlying ThermochemistryJob filename validation
        (since Boltzmann needs at least one real file to average)."""
        runner = CliRunner()
        result = runner.invoke(
            thermochemistry,
            ["-T", "298.15", "boltzmann"],
            obj={},
            catch_exceptions=True,
        )
        assert result.exit_code != 0
        assert isinstance(result.exception, ValueError)
        assert "'filename' must be provided" in str(result.exception)

    def test_outputfile_parameter_is_always_none(
        self, gaussian_co2_opt_outfile
    ):
        """Regression test for BUGS_FOUND.md #39: -o/--outputfile is
        only defined at the thermochemistry group scope and is never
        forwarded into boltzmann()'s own (option-less) `outputfile`
        parameter."""
        result, kwargs, _ = run_boltzmann_and_capture_kwargs(
            [
                "-f",
                gaussian_co2_opt_outfile,
                "-T",
                "298.15",
                "-o",
                "custom_out.dat",
                "boltzmann",
            ]
        )
        assert result.exit_code == 0, result.output
        assert kwargs["outputfile"] is None

    def test_outputfile_option_not_recognized_on_subcommand(
        self, gaussian_co2_opt_outfile
    ):
        runner = CliRunner()
        result = runner.invoke(
            thermochemistry,
            [
                "-f",
                gaussian_co2_opt_outfile,
                "-T",
                "298.15",
                "boltzmann",
                "-o",
                "custom_out.dat",
            ],
            obj={},
            catch_exceptions=True,
        )
        assert result.exit_code != 0
        assert "No such option" in result.output


class TestResolveEntropyCutoff:
    """Direct unit tests for the pure resolve_entropy_cutoff helper."""

    def test_both_given_raises(self):
        import pytest

        with pytest.raises(ValueError, match="Cannot specify both"):
            resolve_entropy_cutoff(100.0, 100.0)

    def test_truhlar_only(self):
        assert resolve_entropy_cutoff(None, 100.0) == (100.0, "truhlar")

    def test_grimme_only(self):
        assert resolve_entropy_cutoff(100.0, None) == (100.0, "grimme")

    def test_neither_given(self):
        assert resolve_entropy_cutoff(None, None) == (None, None)


def run_thermochemistry_directory_mode(cli_args, obj=None):
    """Invoke `thermochemistry ... boltzmann` with a patched Boltzmann
    job class so the group's own directory/filename-resolution logic
    runs (populating ctx.obj) without ever reaching the result
    callback's real per-job computation."""
    runner = CliRunner()
    obj = obj if obj is not None else {}
    with patch(BOLTZMANN_JOB) as mock_job_cls:
        mock_job_cls.return_value = MagicMock()
        result = runner.invoke(
            thermochemistry,
            [*cli_args, "boltzmann"],
            obj=obj,
            catch_exceptions=True,
        )
    return result, obj


class TestThermochemistryGroupDirectoryMode:
    """Covers the `thermochemistry` group's own -d/--directory
    handling: xTB per-subdirectory discovery (valid/empty/invalid
    calculation dirs), program- and filetype-based file discovery, and
    the invalid-file skip-and-continue behavior."""

    def test_xtb_directory_discovers_valid_calculation_subdirs_only(
        self, tmp_path, xtb_water_outfolder, xtb_test_directory
    ):
        shutil.copytree(xtb_water_outfolder, tmp_path / "valid")
        (tmp_path / "empty").mkdir()
        (tmp_path / "multi").mkdir()
        shutil.copy(
            os.path.join(xtb_water_outfolder, "water_ohess.out"),
            tmp_path / "multi" / "water_ohess.out",
        )
        co2_out = os.path.join(
            xtb_test_directory, "outputs", "co2_ohess", "co2_ohess.out"
        )
        shutil.copy(co2_out, tmp_path / "multi" / "co2_ohess.out")

        result, obj = run_thermochemistry_directory_mode(
            ["-d", str(tmp_path), "-p", "xtb", "-T", "298.15"]
        )
        assert result.exit_code == 0, result.output
        assert obj["filenames"] == [
            str(tmp_path / "valid" / "water_ohess.out")
        ]

    def test_xtb_directory_with_no_valid_calculations_raises(self, tmp_path):
        (tmp_path / "empty").mkdir()
        result, _ = run_thermochemistry_directory_mode(
            ["-d", str(tmp_path), "-p", "xtb", "-T", "298.15"]
        )
        assert result.exit_code != 0
        assert "No xTB output files found" in str(result.exception)

    def test_directory_and_program_only(self, tmp_path):
        shutil.copy(
            "tests/data/GaussianTests/outputs/co2.log", tmp_path / "co2.log"
        )
        shutil.copy(
            "tests/data/GaussianTests/outputs/he.log", tmp_path / "he.log"
        )
        result, obj = run_thermochemistry_directory_mode(
            ["-d", str(tmp_path), "-p", "gaussian", "-T", "298.15"]
        )
        assert result.exit_code == 0, result.output
        assert len(obj["filenames"]) == 2

    def test_directory_and_filetype_only(self, tmp_path):
        shutil.copy(
            "tests/data/GaussianTests/outputs/co2.log", tmp_path / "co2.log"
        )
        result, obj = run_thermochemistry_directory_mode(
            ["-d", str(tmp_path), "-t", "log", "-T", "298.15"]
        )
        assert result.exit_code == 0, result.output
        assert len(obj["filenames"]) == 1

    def test_directory_and_program_and_filetype(self, tmp_path):
        shutil.copy(
            "tests/data/GaussianTests/outputs/co2.log", tmp_path / "co2.log"
        )
        result, obj = run_thermochemistry_directory_mode(
            [
                "-d",
                str(tmp_path),
                "-p",
                "gaussian",
                "-t",
                "log",
                "-T",
                "298.15",
            ]
        )
        assert result.exit_code == 0, result.output
        assert len(obj["filenames"]) == 1

    def test_directory_without_program_or_filetype_raises(self, tmp_path):
        result, _ = run_thermochemistry_directory_mode(
            ["-d", str(tmp_path), "-T", "298.15"]
        )
        assert result.exit_code != 0
        assert "Must specify --program or --filetype" in str(result.exception)

    def test_directory_filetype_no_matches_raises(self, tmp_path):
        shutil.copy(
            "tests/data/GaussianTests/outputs/co2.log", tmp_path / "co2.log"
        )
        result, _ = run_thermochemistry_directory_mode(
            ["-d", str(tmp_path), "-t", "xyz", "-T", "298.15"]
        )
        assert result.exit_code != 0
        assert "No output files matching" in str(result.exception)

    def test_unsupported_filetype_in_filenames_mode_raises(
        self, single_molecule_xyz_file
    ):
        result, _ = run_thermochemistry_directory_mode(
            ["-f", single_molecule_xyz_file, "-T", "298.15"]
        )
        assert result.exit_code != 0
        assert "Unsupported output file type" in str(result.exception)

    def test_invalid_file_in_directory_is_skipped_not_fatal(self, tmp_path):
        """One file in the directory fails ThermochemistryJob
        construction; it should be logged and skipped rather than
        aborting the whole directory scan."""
        shutil.copy(
            "tests/data/GaussianTests/outputs/co2.log", tmp_path / "co2.log"
        )
        shutil.copy(
            "tests/data/GaussianTests/outputs/he.log", tmp_path / "he.log"
        )
        orig_from_filename = ThermochemistryJob.from_filename.__func__

        def fake_from_filename(cls, filename, **kwargs):
            if "he" in filename:
                raise ValueError("simulated parse failure")
            return orig_from_filename(cls, filename, **kwargs)

        obj = {}
        with patch.object(
            ThermochemistryJob,
            "from_filename",
            classmethod(fake_from_filename),
        ):
            result, obj = run_thermochemistry_directory_mode(
                ["-d", str(tmp_path), "-t", "log", "-T", "298.15"], obj=obj
            )
        assert result.exit_code == 0, result.output
        assert len(obj["jobs"]) == 1
        assert obj["jobs"][0].label == "co2"

    def test_directory_and_filenames_both_given_raises(self, tmp_path):
        result, _ = run_thermochemistry_directory_mode(
            ["-d", str(tmp_path), "-f", "a.log", "-T", "298.15"]
        )
        assert result.exit_code != 0
        assert "Cannot specify both --directory and --filenames" in str(
            result.exception
        )


class TestThermochemistryResultCallback:
    """Covers thermochemistry_process_pipeline's no-subcommand path:
    running every created job's compute_thermochemistry(), then either
    each job's own show_results() (no --outputfile) or only the first
    job's show_results() (combined --outputfile)."""

    def _make_directory(self, tmp_path):
        shutil.copy(
            "tests/data/GaussianTests/outputs/co2.log", tmp_path / "co2.log"
        )
        shutil.copy(
            "tests/data/GaussianTests/outputs/he.log", tmp_path / "he.log"
        )

    def test_no_outputfile_calls_show_results_on_every_job(self, tmp_path):
        self._make_directory(tmp_path)
        runner = CliRunner()
        with (
            patch.object(
                ThermochemistryJob, "compute_thermochemistry"
            ) as mock_compute,
            patch.object(ThermochemistryJob, "show_results") as mock_show,
        ):
            result = runner.invoke(
                thermochemistry,
                ["-d", str(tmp_path), "-t", "log", "-T", "298.15"],
                obj={},
                catch_exceptions=True,
            )
        assert result.exit_code == 0, result.output
        assert mock_compute.call_count == 2
        assert mock_show.call_count == 2

    def test_outputfile_calls_show_results_once_on_first_job(self, tmp_path):
        self._make_directory(tmp_path)
        runner = CliRunner()
        with (
            patch.object(
                ThermochemistryJob, "compute_thermochemistry"
            ) as mock_compute,
            patch.object(ThermochemistryJob, "show_results") as mock_show,
        ):
            result = runner.invoke(
                thermochemistry,
                [
                    "-d",
                    str(tmp_path),
                    "-t",
                    "log",
                    "-T",
                    "298.15",
                    "-o",
                    "combined.dat",
                ],
                obj={},
                catch_exceptions=True,
            )
        assert result.exit_code == 0, result.output
        assert mock_compute.call_count == 2
        assert mock_show.call_count == 1

    def test_compute_thermochemistry_failure_is_logged_not_fatal(
        self, tmp_path
    ):
        self._make_directory(tmp_path)
        runner = CliRunner()
        with (
            patch.object(
                ThermochemistryJob,
                "compute_thermochemistry",
                side_effect=RuntimeError("boom"),
            ),
            patch.object(ThermochemistryJob, "show_results") as mock_show,
        ):
            result = runner.invoke(
                thermochemistry,
                ["-d", str(tmp_path), "-t", "log", "-T", "298.15"],
                obj={},
                catch_exceptions=True,
            )
        assert result.exit_code == 0, result.output
        assert mock_show.call_count == 2
