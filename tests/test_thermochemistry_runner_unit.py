"""
Direct unit tests for :class:`ThermochemistryJobRunner`.

Exercises path setup (job-directory vs scratch), command/executable
stubs, process execution success/failure paths, and post-run file
handling — none of which are reached by the CLI-propagation tests in
``test_thermochemistry.py``.
"""

import os
from unittest.mock import MagicMock, patch

import pytest

from chemsmart.jobs.thermochemistry.job import ThermochemistryJob
from chemsmart.jobs.thermochemistry.runner import ThermochemistryJobRunner


@pytest.fixture()
def thermo_job(gaussian_co2_opt_outfile, tmp_path):
    job = ThermochemistryJob(
        filename=gaussian_co2_opt_outfile, label="co2_thermo"
    )
    job.folder = str(tmp_path)
    return job


class TestThermochemistryJobRunnerSetup:
    def test_defaults_scratch_to_false(self, pbs_server):
        runner = ThermochemistryJobRunner(server=pbs_server)
        assert runner.scratch is False

    def test_prerun_sets_up_job_directory_paths(self, pbs_server, thermo_job):
        runner = ThermochemistryJobRunner(server=pbs_server, scratch=False)
        runner._prerun(thermo_job)

        assert runner.running_directory == thermo_job.folder
        assert runner.job_inputfile == os.path.abspath(thermo_job.inputfile)
        assert runner.job_outputfile == os.path.abspath(thermo_job.outputfile)
        assert runner.job_errfile == os.path.abspath(thermo_job.errfile)

    def test_prerun_sets_up_scratch_directory_paths(
        self, pbs_server, thermo_job, tmp_path
    ):
        scratch_dir = tmp_path / "scratch"
        runner = ThermochemistryJobRunner(
            server=pbs_server, scratch=True, scratch_dir=str(scratch_dir)
        )
        runner._prerun(thermo_job)

        expected_dir = os.path.join(str(scratch_dir), thermo_job.label)
        assert runner.running_directory == expected_dir
        assert os.path.isdir(expected_dir)
        assert runner.job_inputfile.startswith(expected_dir)
        assert runner.job_outputfile.endswith("co2_thermo.dat")
        assert runner.job_errfile.endswith("co2_thermo.err")

    def test_get_command_and_executable_are_none(self, pbs_server, thermo_job):
        runner = ThermochemistryJobRunner(server=pbs_server)
        assert runner._get_command(thermo_job) is None
        assert runner._get_executable() is None

    def test_run_is_a_no_op(self, pbs_server):
        runner = ThermochemistryJobRunner(server=pbs_server)
        # Should not raise regardless of the "process" passed in.
        runner._run(process=None)


class TestThermochemistryJobRunnerWriteInput:
    def test_write_input_copies_when_scratch_and_paths_differ(
        self, pbs_server, thermo_job, tmp_path
    ):
        scratch_dir = tmp_path / "scratch"
        runner = ThermochemistryJobRunner(
            server=pbs_server, scratch=True, scratch_dir=str(scratch_dir)
        )
        runner._prerun(thermo_job)

        with patch("chemsmart.jobs.thermochemistry.runner.copy") as mock_copy:
            runner._write_input(thermo_job)

        mock_copy.assert_called_once_with(
            thermo_job.inputfile, runner.job_inputfile
        )

    def test_write_input_no_op_without_scratch(self, pbs_server, thermo_job):
        runner = ThermochemistryJobRunner(server=pbs_server, scratch=False)
        runner._prerun(thermo_job)

        with patch("chemsmart.jobs.thermochemistry.runner.copy") as mock_copy:
            runner._write_input(thermo_job)

        mock_copy.assert_not_called()


class TestThermochemistryJobRunnerCreateProcess:
    def test_create_process_returns_zero_on_success(
        self, pbs_server, thermo_job
    ):
        runner = ThermochemistryJobRunner(server=pbs_server)
        runner._prerun(thermo_job)
        thermo_job.compute_thermochemistry = MagicMock()

        exit_code = runner._create_process(thermo_job, None, None)

        assert exit_code == 0
        thermo_job.compute_thermochemistry.assert_called_once()

    def test_create_process_returns_one_and_writes_errfile_on_failure(
        self, pbs_server, thermo_job
    ):
        runner = ThermochemistryJobRunner(server=pbs_server)
        runner._prerun(thermo_job)
        thermo_job.compute_thermochemistry = MagicMock(
            side_effect=RuntimeError("boom")
        )

        exit_code = runner._create_process(thermo_job, None, None)

        assert exit_code == 1
        assert os.path.exists(runner.job_errfile)
        with open(runner.job_errfile) as f:
            assert "boom" in f.read()


class TestThermochemistryJobRunnerPostrun:
    def test_postrun_no_op_without_scratch(self, pbs_server, thermo_job):
        runner = ThermochemistryJobRunner(server=pbs_server, scratch=False)
        runner._prerun(thermo_job)
        # Should not raise even though output/err files don't exist.
        runner._postrun(thermo_job)

    def test_postrun_copies_files_still_missing_in_scratch(
        self, pbs_server, thermo_job, tmp_path
    ):
        """``_postrun`` copies whichever of outputfile/errfile is still
        absent at the scratch location (only the errfile here, since the
        outputfile was already produced there)."""
        scratch_dir = tmp_path / "scratch"
        runner = ThermochemistryJobRunner(
            server=pbs_server, scratch=True, scratch_dir=str(scratch_dir)
        )
        runner._prerun(thermo_job)

        # Simulate output produced in the scratch directory; errfile absent.
        os.makedirs(os.path.dirname(runner.job_outputfile), exist_ok=True)
        with open(runner.job_outputfile, "w") as f:
            f.write("result")

        with patch("chemsmart.jobs.thermochemistry.runner.copy") as mock_copy:
            runner._postrun(thermo_job)

        mock_copy.assert_called_once_with(
            runner.job_errfile, thermo_job.folder
        )
