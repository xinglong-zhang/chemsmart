from pathlib import Path
from types import SimpleNamespace

from chemsmart.jobs.gaussian.runner import (
    FakeGaussianJobRunner,
    GaussianJobRunner,
)
from chemsmart.jobs.iterate.runner import IterateJobRunner
from chemsmart.jobs.orca.runner import FakeORCAJobRunner, ORCAJobRunner
from chemsmart.jobs.runner import JobRunner


class DummyORCAJob:
    def __init__(self, folder, label):
        self.folder = str(folder)
        self.label = label

    @property
    def inputfile(self):
        return str(Path(self.folder) / f"{self.label}.inp")

    @property
    def gbwfile(self):
        return str(Path(self.folder) / f"{self.label}.gbw")

    @property
    def errfile(self):
        return str(Path(self.folder) / f"{self.label}.err")

    @property
    def outputfile(self):
        return str(Path(self.folder) / f"{self.label}.out")


class DummyGaussianJob:
    def __init__(self, folder, label):
        self.folder = str(folder)
        self.label = label

    @property
    def inputfile(self):
        return str(Path(self.folder) / f"{self.label}.com")

    @property
    def chkfile(self):
        return str(Path(self.folder) / f"{self.label}.chk")

    @property
    def errfile(self):
        return str(Path(self.folder) / f"{self.label}.err")


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

    def test_fake_flag_propagates_when_no_fake_runner_exists(self, pbs_server):
        job = SimpleNamespace(TYPE="iterate")
        runner = JobRunner.from_job(
            job=job, server=pbs_server, scratch=False, fake=True
        )
        assert isinstance(runner, IterateJobRunner)
        assert runner.fake is True

    def test_fake_orca_appends_fake_suffix_in_job_directory(
        self, pbs_server, tmp_path
    ):
        runner = FakeORCAJobRunner(server=pbs_server, scratch=False, fake=True)
        job = DummyORCAJob(folder=tmp_path, label="orca_opt")

        runner._set_up_variables_in_job_directory(job)

        assert job.label == "orca_opt_fake"
        assert Path(runner.job_inputfile).name == "orca_opt_fake.inp"
        assert Path(runner.job_gbwfile).name == "orca_opt_fake.gbw"
        assert Path(runner.job_errfile).name == "orca_opt_fake.err"
        assert Path(runner.job_outputfile).name == "orca_opt_fake.out"

    def test_fake_orca_appends_fake_suffix_in_scratch(
        self, pbs_server, tmp_path
    ):
        runner = FakeORCAJobRunner(
            server=pbs_server,
            scratch=True,
            scratch_dir=str(tmp_path),
            fake=True,
        )
        job = DummyORCAJob(folder=tmp_path, label="orca_opt")

        runner._set_up_variables_in_scratch(job)

        assert job.label == "orca_opt_fake"
        assert Path(runner.job_inputfile).name == "orca_opt_fake.inp"
        assert Path(runner.job_gbwfile).name == "orca_opt_fake.gbw"
        assert Path(runner.job_errfile).name == "orca_opt_fake.err"
        assert Path(runner.job_outputfile).name == "orca_opt_fake.out"

    def test_fake_gaussian_does_not_duplicate_fake_suffix(
        self, pbs_server, tmp_path
    ):
        runner = FakeGaussianJobRunner(server=pbs_server, scratch=False, fake=True)
        job = DummyGaussianJob(folder=tmp_path, label="gaussian_opt_fake")

        runner._set_up_variables_in_job_directory(job)

        assert job.label == "gaussian_opt_fake"
        assert Path(runner.job_inputfile).name == "gaussian_opt_fake.com"
        assert Path(runner.job_chkfile).name == "gaussian_opt_fake.chk"
        assert Path(runner.job_errfile).name == "gaussian_opt_fake.err"

    def test_fake_orca_does_not_duplicate_fake_suffix(
        self, pbs_server, tmp_path
    ):
        runner = FakeORCAJobRunner(server=pbs_server, scratch=False, fake=True)
        job = DummyORCAJob(folder=tmp_path, label="orca_opt_fake")

        runner._set_up_variables_in_job_directory(job)

        assert job.label == "orca_opt_fake"
        assert Path(runner.job_inputfile).name == "orca_opt_fake.inp"
        assert Path(runner.job_gbwfile).name == "orca_opt_fake.gbw"
        assert Path(runner.job_errfile).name == "orca_opt_fake.err"
        assert Path(runner.job_outputfile).name == "orca_opt_fake.out"
