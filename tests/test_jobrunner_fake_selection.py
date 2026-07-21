from chemsmart.jobs.gaussian.runner import (
    FakeGaussianJobRunner,
    GaussianJobRunner,
)
from chemsmart.jobs.orca.runner import FakeORCAJobRunner, ORCAJobRunner
from chemsmart.jobs.runner import JobRunner


class _Job:
    def __init__(self, jobtype: str):
        self.TYPE = jobtype


def test_from_job_selects_fake_gaussian_runner_when_requested(pbs_server):
    runner = JobRunner.from_job(
        _Job("g16opt"),
        server=pbs_server,
        scratch=False,
        fake=True,
    )

    assert isinstance(runner, FakeGaussianJobRunner)
    assert runner.fake is True


def test_from_job_selects_regular_gaussian_runner_by_default(pbs_server):
    runner = JobRunner.from_job(
        _Job("g16opt"),
        server=pbs_server,
        scratch=False,
        fake=False,
    )

    assert isinstance(runner, GaussianJobRunner)
    assert not isinstance(runner, FakeGaussianJobRunner)
    assert runner.fake is False


def test_from_job_selects_fake_orca_runner_when_requested(pbs_server):
    runner = JobRunner.from_job(
        _Job("orcaopt"),
        server=pbs_server,
        scratch=False,
        fake=True,
    )

    assert isinstance(runner, FakeORCAJobRunner)
    assert runner.fake is True


def test_from_job_selects_regular_orca_runner_by_default(pbs_server):
    runner = JobRunner.from_job(
        _Job("orcaopt"),
        server=pbs_server,
        scratch=False,
        fake=False,
    )

    assert isinstance(runner, ORCAJobRunner)
    assert not isinstance(runner, FakeORCAJobRunner)
    assert runner.fake is False
