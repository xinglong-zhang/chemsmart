from types import SimpleNamespace

from chemsmart.jobs.gaussian.runner import (
    FakeGaussianJobRunner,
    GaussianJobRunner,
)
from chemsmart.jobs.iterate.runner import IterateJobRunner
from chemsmart.jobs.orca.runner import FakeORCAJobRunner, ORCAJobRunner
from chemsmart.jobs.runner import JobRunner


class TestJobRunnerSelection:
    def test_fake_gaussian_runner_selected_when_fake_enabled(self, pbs_server):
        job = SimpleNamespace(TYPE="g16opt")
        runner = JobRunner.from_job(
            job=job, server=pbs_server, scratch=False, fake=True
        )
        assert isinstance(runner, FakeGaussianJobRunner)
        assert runner.fake is True

    def test_fake_orca_runner_selected_when_fake_enabled(self, pbs_server):
        job = SimpleNamespace(TYPE="orcasp")
        runner = JobRunner.from_job(
            job=job, server=pbs_server, scratch=False, fake=True
        )
        assert isinstance(runner, FakeORCAJobRunner)
        assert runner.fake is True

    def test_real_runner_selected_when_fake_disabled(self, pbs_server):
        gaussian_job = SimpleNamespace(TYPE="g16opt")
        gaussian_runner = JobRunner.from_job(
            job=gaussian_job, server=pbs_server, scratch=False, fake=False
        )
        assert isinstance(gaussian_runner, GaussianJobRunner)

        orca_job = SimpleNamespace(TYPE="orcasp")
        orca_runner = JobRunner.from_job(
            job=orca_job, server=pbs_server, scratch=False, fake=False
        )
        assert isinstance(orca_runner, ORCAJobRunner)

    def test_fake_flag_propagates_when_no_fake_runner_exists(
        self, pbs_server
    ):
        job = SimpleNamespace(TYPE="iterate")
        runner = JobRunner.from_job(
            job=job, server=pbs_server, scratch=False, fake=True
        )
        assert isinstance(runner, IterateJobRunner)
        assert runner.fake is True
