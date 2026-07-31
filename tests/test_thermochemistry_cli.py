"""
CLI tests for the ``thermochemistry`` group's ``boltzmann`` subcommand
(chemsmart/cli/thermochemistry/boltzmann.py).

Each test uses :class:`click.testing.CliRunner` to invoke the
``thermochemistry`` group and :mod:`unittest.mock` to intercept
``BoltzmannAverageThermochemistryJob`` so the merged call kwargs can be
inspected without running an actual Boltzmann-averaging computation.
"""

from unittest.mock import MagicMock, patch

from click.testing import CliRunner

from chemsmart.cli.thermochemistry.thermochemistry import thermochemistry

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
