"""Tests for JobRunner class and run_in_serial flag."""
import pytest

from chemsmart.jobs.runner import JobRunner
from chemsmart.settings.server import Server


class TestJobRunner:
    """Test JobRunner functionality."""

    def test_jobrunner_default_run_in_serial(self, pbs_server):
        """Test that run_in_serial defaults to False."""
        runner = JobRunner(server=pbs_server, fake=True)
        assert hasattr(runner, 'run_in_serial')
        assert runner.run_in_serial is False

    def test_jobrunner_run_in_serial_true(self, pbs_server):
        """Test that run_in_serial can be set to True."""
        runner = JobRunner(
            server=pbs_server, fake=True, run_in_serial=True
        )
        assert runner.run_in_serial is True

    def test_jobrunner_run_in_serial_false(self, pbs_server):
        """Test that run_in_serial can be set to False."""
        runner = JobRunner(
            server=pbs_server, fake=True, run_in_serial=False
        )
        assert runner.run_in_serial is False


