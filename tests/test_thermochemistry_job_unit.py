"""
Direct unit tests for :class:`ThermochemistryJob`.

The CLI-level propagation of options is already covered by
``TestThermochemistryCLI`` in ``test_thermochemistry.py``. These tests
exercise the job class directly: constructor validation, path
properties, completion checks, and the ``compute_thermochemistry``/
``show_results`` methods (with the underlying ``Thermochemistry``
analysis mocked out).
"""

import os
from unittest.mock import MagicMock, patch

import pytest

from chemsmart.jobs.thermochemistry.job import ThermochemistryJob
from chemsmart.jobs.thermochemistry.settings import (
    ThermochemistryJobSettings,
)


class TestThermochemistryJobConstruction:
    def test_rejects_unsupported_file_extension(self):
        with pytest.raises(ValueError, match="Unsupported file extension"):
            ThermochemistryJob(filename="molecule.txt")

    def test_rejects_invalid_settings_type(self, gaussian_co2_opt_outfile):
        with pytest.raises(ValueError, match="Settings must be instance"):
            ThermochemistryJob(
                filename=gaussian_co2_opt_outfile, settings="not-a-settings"
            )

    def test_rejects_invalid_molecule_type(self, gaussian_co2_opt_outfile):
        with pytest.raises(ValueError, match="Molecule must be instance"):
            ThermochemistryJob(
                filename=gaussian_co2_opt_outfile, molecule="not-a-molecule"
            )

    def test_label_defaults_to_filename_stem(self, gaussian_co2_opt_outfile):
        job = ThermochemistryJob(filename=gaussian_co2_opt_outfile)
        assert job.label == "co2"

    def test_default_settings_used_when_none_given(
        self, gaussian_co2_opt_outfile
    ):
        job = ThermochemistryJob(filename=gaussian_co2_opt_outfile)
        assert isinstance(job.settings, ThermochemistryJobSettings)

    def test_settings_class_returns_expected_type(self):
        assert (
            ThermochemistryJob.settings_class() is ThermochemistryJobSettings
        )


class TestThermochemistryJobPaths:
    def test_inputfile_is_absolute_path(self, gaussian_co2_opt_outfile):
        job = ThermochemistryJob(filename=gaussian_co2_opt_outfile)
        assert job.inputfile == os.path.abspath(gaussian_co2_opt_outfile)

    def test_inputfile_none_when_no_filename(self):
        job = ThermochemistryJob.__new__(ThermochemistryJob)
        job.filename = None
        assert job.inputfile is None

    def test_outputfile_and_errfile_use_label_and_folder(
        self, gaussian_co2_opt_outfile
    ):
        job = ThermochemistryJob(
            filename=gaussian_co2_opt_outfile, label="my_label"
        )
        assert job.outputfile == os.path.join(job.folder, "my_label.dat")
        assert job.errfile == os.path.join(job.folder, "my_label.err")


class TestThermochemistryJobCompletion:
    def test_job_is_complete_false_when_missing(
        self, gaussian_co2_opt_outfile, tmp_path
    ):
        job = ThermochemistryJob(filename=gaussian_co2_opt_outfile)
        job.settings.outputfile = str(tmp_path / "does_not_exist.dat")
        assert job._job_is_complete() is False

    def test_job_is_complete_true_when_present(
        self, gaussian_co2_opt_outfile, tmp_path
    ):
        job = ThermochemistryJob(filename=gaussian_co2_opt_outfile)
        outfile = tmp_path / "exists.dat"
        outfile.write_text("done")
        job.settings.outputfile = str(outfile)
        assert job._job_is_complete() is True

    def test_output_returns_none_when_missing(
        self, gaussian_co2_opt_outfile, tmp_path
    ):
        job = ThermochemistryJob(filename=gaussian_co2_opt_outfile)
        job.settings.outputfile = str(tmp_path / "missing.dat")
        assert job._output() is None

    def test_output_returns_abspath_when_present(
        self, gaussian_co2_opt_outfile, tmp_path
    ):
        job = ThermochemistryJob(filename=gaussian_co2_opt_outfile)
        outfile = tmp_path / "present.dat"
        outfile.write_text("done")
        job.settings.outputfile = str(outfile)
        assert job._output() == os.path.abspath(str(outfile))


class TestThermochemistryJobRun:
    def test_run_delegates_to_jobrunner(self, gaussian_co2_opt_outfile):
        mock_runner = MagicMock()
        job = ThermochemistryJob(
            filename=gaussian_co2_opt_outfile, jobrunner=mock_runner
        )
        job._run()
        mock_runner.run.assert_called_once_with(job)


class TestThermochemistryJobComputeAndShow:
    def test_compute_thermochemistry_requires_filename(self):
        job = ThermochemistryJob.__new__(ThermochemistryJob)
        job.filename = None
        job.settings = ThermochemistryJobSettings()
        with pytest.raises(ValueError, match="No input file provided"):
            job.compute_thermochemistry()

    def test_compute_thermochemistry_calls_analysis_and_logs_results(
        self, gaussian_co2_opt_outfile, tmp_path
    ):
        job = ThermochemistryJob(filename=gaussian_co2_opt_outfile)
        job.settings.outputfile = str(tmp_path / "out.dat")

        mock_thermo = MagicMock()
        mock_thermo.compute_thermochemistry.return_value = (
            "structure",
            -100.0,
            0.01,
            -99.9,
            -99.9,
            -0.02,
            -0.02,
            -99.92,
            -99.92,
        )
        with patch(
            "chemsmart.jobs.thermochemistry.job.Thermochemistry",
            return_value=mock_thermo,
        ) as mock_cls:
            job.compute_thermochemistry()

        mock_cls.assert_called_once()
        mock_thermo.compute_thermochemistry.assert_called_once()
        mock_thermo.log_results_to_file.assert_called_once()

    def test_compute_thermochemistry_reraises_on_failure(
        self, gaussian_co2_opt_outfile, tmp_path
    ):
        job = ThermochemistryJob(filename=gaussian_co2_opt_outfile)
        job.settings.outputfile = str(tmp_path / "out.dat")

        with patch(
            "chemsmart.jobs.thermochemistry.job.Thermochemistry",
            side_effect=RuntimeError("parse failed"),
        ):
            with pytest.raises(RuntimeError, match="parse failed"):
                job.compute_thermochemistry()

    def test_show_results_prints_output_file_contents(
        self, gaussian_co2_opt_outfile, tmp_path, capsys
    ):
        job = ThermochemistryJob(filename=gaussian_co2_opt_outfile)
        outfile = tmp_path / "results.dat"
        outfile.write_text("some thermochemistry results\n")
        job.settings.outputfile = str(outfile)

        job.show_results()

        captured = capsys.readouterr()
        assert "some thermochemistry results" in captured.out

    def test_show_results_no_op_when_outputfile_missing(
        self, gaussian_co2_opt_outfile, tmp_path, capsys
    ):
        job = ThermochemistryJob(filename=gaussian_co2_opt_outfile)
        job.settings.outputfile = str(tmp_path / "missing.dat")

        job.show_results()

        captured = capsys.readouterr()
        assert captured.out == ""


class TestThermochemistryJobFromFilename:
    def test_from_filename_reads_molecule_and_builds_job(
        self, gaussian_co2_opt_outfile
    ):
        mock_jobrunner = MagicMock()
        job = ThermochemistryJob.from_filename(
            filename=gaussian_co2_opt_outfile, jobrunner=mock_jobrunner
        )
        assert isinstance(job, ThermochemistryJob)
        assert job.filename == gaussian_co2_opt_outfile
        assert job.jobrunner is mock_jobrunner
        assert job.label == "co2"

    def test_from_filename_creates_jobrunner_when_not_given(
        self, gaussian_co2_opt_outfile
    ):
        mock_runner_instance = MagicMock()
        with patch(
            "chemsmart.jobs.thermochemistry.job.JobRunner.from_job",
            return_value=mock_runner_instance,
        ) as mock_from_job:
            job = ThermochemistryJob.from_filename(
                filename=gaussian_co2_opt_outfile
            )

        mock_from_job.assert_called_once()
        assert job.jobrunner is mock_runner_instance
