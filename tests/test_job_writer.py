"""Direct unit tests for the abstract :class:`chemsmart.jobs.writer.InputWriter`."""

from unittest.mock import MagicMock

import pytest

from chemsmart.jobs.writer import InputWriter


class TestInputWriter:
    def test_init_stores_job_settings_and_jobrunner(self):
        job = MagicMock()
        writer = InputWriter(job)

        assert writer.job is job
        assert writer.settings is job.settings
        assert writer.jobrunner is job.jobrunner

    def test_write_not_implemented(self):
        job = MagicMock()
        writer = InputWriter(job)

        with pytest.raises(NotImplementedError):
            writer.write()
