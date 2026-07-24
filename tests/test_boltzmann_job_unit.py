"""
Direct unit tests for :class:`BoltzmannAverageThermochemistryJob`.

Note: the class's parent ``ThermochemistryJob.__init__`` unconditionally
calls ``filename.endswith(...)`` to validate the file extension, so a
``filename`` must always be supplied even though ``BoltzmannAverage
ThermochemistryJob`` itself operates on ``files`` (plural). The
``chemsmart sub thermochemistry boltzmann`` CLI command
(``chemsmart/cli/thermochemistry/boltzmann.py``) never passes
``filename``, so it currently raises ``AttributeError: 'NoneType'
object has no attribute 'endswith'`` on any real invocation — this
looks like a genuine bug, flagged separately. These tests construct the
job directly with an explicit ``filename`` to exercise the class's own
logic in isolation.
"""

import os
from unittest.mock import MagicMock, patch

import pytest

from chemsmart.jobs.thermochemistry.boltzmann import (
    BoltzmannAverageThermochemistryJob,
)
from chemsmart.jobs.thermochemistry.settings import (
    ThermochemistryJobSettings,
)


class TestBoltzmannJobConstruction:
    def test_label_actually_derives_from_filename_not_files(
        self, gaussian_co2_opt_outfile
    ):
        """The boltzmann-specific "common prefix of ``files``" label logic
        in ``BoltzmannAverageThermochemistryJob.__init__`` is dead code in
        practice: ``ThermochemistryJob.__init__`` (called first via
        ``super().__init__``) already derives ``self.label`` from
        ``filename`` whenever one is given, so the child's own
        ``if self.label is None`` branch never fires as long as the
        (mandatory, see docstring) ``filename`` is supplied."""
        job = BoltzmannAverageThermochemistryJob(
            filename=gaussian_co2_opt_outfile,
            files=["conformer_c1.log", "conformer_c2.log"],
            energy_type="gibbs",
        )
        assert job.label == "co2"

    def test_custom_label_preserved(self, gaussian_co2_opt_outfile):
        job = BoltzmannAverageThermochemistryJob(
            filename=gaussian_co2_opt_outfile,
            files=["a.log", "b.log"],
            label="my_boltzmann_job",
        )
        assert job.label == "my_boltzmann_job"

    def test_default_energy_type_is_gibbs(self, gaussian_co2_opt_outfile):
        job = BoltzmannAverageThermochemistryJob(
            filename=gaussian_co2_opt_outfile,
            files=["a.log", "b.log"],
            label="x",
        )
        assert job.energy_type == "gibbs"

    def test_settings_class_returns_expected_type(self):
        assert (
            BoltzmannAverageThermochemistryJob.settings_class()
            is ThermochemistryJobSettings
        )


class TestBoltzmannJobPaths:
    def test_inputfile_is_absolute_path(self, gaussian_co2_opt_outfile):
        job = BoltzmannAverageThermochemistryJob(
            filename=gaussian_co2_opt_outfile,
            files=["a.log", "b.log"],
            label="x",
        )
        assert job.inputfile == os.path.abspath(gaussian_co2_opt_outfile)

    def test_inputfile_none_when_no_filename(self):
        job = BoltzmannAverageThermochemistryJob.__new__(
            BoltzmannAverageThermochemistryJob
        )
        job.filename = None
        assert job.inputfile is None

    def test_outputfile_and_errfile_use_label(self, gaussian_co2_opt_outfile):
        job = BoltzmannAverageThermochemistryJob(
            filename=gaussian_co2_opt_outfile,
            files=["a.log", "b.log"],
            label="my_boltzmann_job",
        )
        assert job.outputfile == os.path.join(
            job.folder, "my_boltzmann_job_boltzmann.dat"
        )
        assert job.errfile == os.path.join(
            job.folder, "my_boltzmann_job_boltzmann.err"
        )


class TestBoltzmannJobFromFiles:
    def test_from_files_requires_filename_workaround(self):
        """Demonstrates the underlying bug: ``from_files`` alone (as used
        by the CLI) crashes because no ``filename`` is ever supplied."""
        with pytest.raises(AttributeError):
            BoltzmannAverageThermochemistryJob.from_files(
                files=["a.log", "b.log"]
            )

    def test_from_files_builds_job_with_explicit_filename(
        self, gaussian_co2_opt_outfile
    ):
        mock_jobrunner = MagicMock()
        job = BoltzmannAverageThermochemistryJob.from_files(
            files=["a.log", "b.log"],
            filename=gaussian_co2_opt_outfile,
            jobrunner=mock_jobrunner,
        )
        assert isinstance(job, BoltzmannAverageThermochemistryJob)
        assert job.files == ["a.log", "b.log"]
        assert job.jobrunner is mock_jobrunner

    def test_from_files_creates_jobrunner_when_not_given(
        self, gaussian_co2_opt_outfile
    ):
        mock_runner_instance = MagicMock()
        with patch(
            "chemsmart.jobs.thermochemistry.boltzmann.JobRunner.from_job",
            return_value=mock_runner_instance,
        ) as mock_from_job:
            job = BoltzmannAverageThermochemistryJob.from_files(
                files=["a.log", "b.log"],
                filename=gaussian_co2_opt_outfile,
            )

        mock_from_job.assert_called_once()
        assert job.jobrunner is mock_runner_instance


class TestBoltzmannJobCompute:
    def test_compute_requires_files(self, gaussian_co2_opt_outfile):
        job = BoltzmannAverageThermochemistryJob(
            filename=gaussian_co2_opt_outfile, files=None, label="x"
        )
        with pytest.raises(ValueError, match="No input file provided"):
            job.compute_boltzmann_averages()

    def test_compute_calls_analysis_and_logs_results(
        self, gaussian_co2_opt_outfile, tmp_path
    ):
        job = BoltzmannAverageThermochemistryJob(
            filename=gaussian_co2_opt_outfile,
            files=["a.log", "b.log"],
            label="x",
        )
        job.settings.outputfile = str(tmp_path / "out.dat")

        mock_thermo = MagicMock()
        mock_thermo.compute_boltzmann_averages.return_value = (
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
            "chemsmart.analysis.thermochemistry.BoltzmannAverageThermochemistry",
            return_value=mock_thermo,
        ) as mock_cls:
            job.compute_boltzmann_averages()

        mock_cls.assert_called_once()
        assert mock_cls.call_args.kwargs["files"] == ["a.log", "b.log"]
        mock_thermo.compute_boltzmann_averages.assert_called_once()
        mock_thermo.log_results_to_file.assert_called_once()

    def test_compute_reraises_on_failure(
        self, gaussian_co2_opt_outfile, tmp_path
    ):
        job = BoltzmannAverageThermochemistryJob(
            filename=gaussian_co2_opt_outfile,
            files=["a.log", "b.log"],
            label="x",
        )
        job.settings.outputfile = str(tmp_path / "out.dat")

        with patch(
            "chemsmart.analysis.thermochemistry.BoltzmannAverageThermochemistry",
            side_effect=RuntimeError("boom"),
        ):
            with pytest.raises(RuntimeError, match="boom"):
                job.compute_boltzmann_averages()
