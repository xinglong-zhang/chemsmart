from pathlib import Path
from types import SimpleNamespace

from click.testing import CliRunner

from chemsmart.cli.run import run
from chemsmart.cli.sub import sub
from chemsmart.jobs.chain import ChainJob
from chemsmart.jobs.chain_runner import ChainJobRunner, FakeChainJobRunner
from chemsmart.jobs.gaussian.runner import (
    FakeGaussianJobRunner,
    GaussianJobRunner,
)
from chemsmart.jobs.iterate.runner import IterateJobRunner
from chemsmart.jobs.job import Job
from chemsmart.jobs.orca.runner import FakeORCAJobRunner, ORCAJobRunner
from chemsmart.jobs.runner import JobRunner
from chemsmart.jobs.xtb.runner import FakeXTBJobRunner, XTBJobRunner
from chemsmart.settings.server import Server


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

    def test_fake_xtb_runner_selected_when_fake_enabled(self, pbs_server):
        job = SimpleNamespace(TYPE="xtbhess")
        runner = JobRunner.from_job(
            job=job, server=pbs_server, scratch=False, fake=True
        )
        assert isinstance(runner, FakeXTBJobRunner)
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

        xtb_job = SimpleNamespace(TYPE="xtbsp")
        xtb_runner = JobRunner.from_job(
            job=xtb_job, server=pbs_server, scratch=False, fake=False
        )
        assert isinstance(xtb_runner, XTBJobRunner)

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
        runner = FakeGaussianJobRunner(
            server=pbs_server, scratch=False, fake=True
        )
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


class TestPropagateRunnerRetype:
    def test_same_program_child_copies_runner(self, pbs_server):
        parent = GaussianJobRunner(
            server=pbs_server, scratch=False, num_cores=8
        )
        child = SimpleNamespace(TYPE="g16opt", jobrunner=None)
        result = Job._propagate_runner(parent, child, num_cores=2)

        assert isinstance(result, GaussianJobRunner)
        assert result is not parent
        assert child.jobrunner is result
        assert result.num_cores == 2
        assert parent.num_cores == 8

    def test_different_program_child_retypes_runner(self, pbs_server):
        parent = GaussianJobRunner(server=pbs_server, scratch=False)
        child = SimpleNamespace(TYPE="xtbopt", jobrunner=None)
        result = Job._propagate_runner(parent, child)

        assert isinstance(result, XTBJobRunner)
        assert not isinstance(result, FakeXTBJobRunner)
        assert child.jobrunner is result
        assert result.server is parent.server
        assert result.scratch is False
        assert result.fake is False

    def test_fake_parent_retypes_to_fake_child(self, pbs_server):
        parent = FakeGaussianJobRunner(
            server=pbs_server, scratch=False, fake=True
        )
        child = SimpleNamespace(TYPE="xtbopt", jobrunner=None)
        result = Job._propagate_runner(parent, child)

        assert isinstance(result, FakeXTBJobRunner)
        assert result.fake is True

    def test_resource_overrides_apply_after_retype(self, pbs_server):
        parent = GaussianJobRunner(
            server=pbs_server, scratch=False, num_cores=8, mem_gb=16
        )
        child = SimpleNamespace(TYPE="xtbopt", jobrunner=None)
        result = Job._propagate_runner(
            parent, child, num_cores=2, mem_gb=4, num_gpus=1
        )

        assert isinstance(result, XTBJobRunner)
        assert result.num_cores == 2
        assert result.mem_gb == 4
        assert result.num_gpus == 1

    def test_untyped_parent_runner_retypes_child(self, pbs_server):
        parent = JobRunner(server=pbs_server, scratch=False, fake=False)
        child = SimpleNamespace(TYPE="g16opt", jobrunner=None)
        result = Job._propagate_runner(parent, child)

        assert isinstance(result, GaussianJobRunner)

    def test_chain_parent_retypes_gaussian_child(self, tmp_path):
        yaml_path = _write_server_yaml(
            tmp_path / "g16_on.yaml",
            gaussian_scratch=True,
            orca_scratch=False,
        )
        server = Server.from_yaml(str(yaml_path))
        scratch_dir = tmp_path / "scratch"
        scratch_dir.mkdir()
        parent = ChainJobRunner(server=server, scratch_dir=str(scratch_dir))
        child = SimpleNamespace(TYPE="g16opt", jobrunner=None)
        result = Job._propagate_runner(parent, child)

        assert isinstance(result, GaussianJobRunner)
        assert child.jobrunner is result
        assert result.scratch is True

    def test_chain_parent_no_scratch_forces_child_off(self, pbs_server):
        parent = ChainJobRunner(server=pbs_server, scratch=False)
        child = SimpleNamespace(TYPE="g16opt", jobrunner=None)
        result = Job._propagate_runner(parent, child)

        assert isinstance(result, GaussianJobRunner)
        assert result.scratch is False

    def test_none_runner_returns_none(self):
        child = SimpleNamespace(TYPE="g16opt", jobrunner="unchanged")
        result = Job._propagate_runner(None, child)
        assert result is None
        assert child.jobrunner == "unchanged"


class TestChainJobRunnerSelection:
    def test_chain_type_selects_chain_runner(self, pbs_server):
        job = SimpleNamespace(TYPE="chain")
        runner = JobRunner.from_job(
            job=job, server=pbs_server, scratch=False, fake=False
        )
        assert isinstance(runner, ChainJobRunner)
        assert not isinstance(runner, FakeChainJobRunner)
        assert runner.fake is False
        assert runner.scratch is False

    def test_fake_chain_type_selects_fake_chain_runner(self, pbs_server):
        job = SimpleNamespace(TYPE="chain")
        runner = JobRunner.from_job(
            job=job, server=pbs_server, scratch=False, fake=True
        )
        assert isinstance(runner, FakeChainJobRunner)
        assert runner.fake is True

    def test_chain_job_instance_selects_chain_runner(self, pbs_server):
        job = ChainJob(molecule=None, label="chain", jobrunner=None)
        runner = JobRunner.from_job(
            job=job, server=pbs_server, scratch=None, fake=False
        )
        assert isinstance(runner, ChainJobRunner)
        assert runner.scratch is None


class TestScratchCLI:
    """Scratch CLI wiring: omit vs explicit flags."""

    def test_run_omitted_scratch_leaves_none_for_from_job(
        self,
        monkeypatch,
        single_molecule_xyz_file,
        gaussian_project_config_dir,
    ):
        monkeypatch.setenv(
            "CHEMSMART_CONFIG_DIR", str(gaussian_project_config_dir)
        )
        observed = {"scratch_arg": "unset"}

        def _from_job(cls, job, server, scratch=None, fake=False, **kwargs):
            observed["scratch_arg"] = scratch
            return type("R", (), {"scratch": scratch})()

        monkeypatch.setattr(
            "chemsmart.jobs.runner.JobRunner.from_job",
            classmethod(_from_job),
        )
        monkeypatch.setattr("chemsmart.jobs.job.Job.run", lambda self: None)

        result = CliRunner().invoke(
            run,
            [
                "--fake",
                "gaussian",
                "-p",
                "test",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "opt",
            ],
        )
        assert result.exit_code == 0, result.output
        assert observed["scratch_arg"] is None

    def test_run_explicit_scratch_reaches_from_job(
        self,
        monkeypatch,
        single_molecule_xyz_file,
        gaussian_project_config_dir,
    ):
        monkeypatch.setenv(
            "CHEMSMART_CONFIG_DIR", str(gaussian_project_config_dir)
        )
        observed = {"scratch_arg": "unset"}

        def _from_job(cls, job, server, scratch=None, fake=False, **kwargs):
            observed["scratch_arg"] = scratch
            return type("R", (), {"scratch": scratch})()

        monkeypatch.setattr(
            "chemsmart.jobs.runner.JobRunner.from_job",
            classmethod(_from_job),
        )
        monkeypatch.setattr("chemsmart.jobs.job.Job.run", lambda self: None)

        result = CliRunner().invoke(
            run,
            [
                "--fake",
                "--scratch",
                "gaussian",
                "-p",
                "test",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "opt",
            ],
        )
        assert result.exit_code == 0, result.output
        assert observed["scratch_arg"] is True

    def test_sub_omitted_scratch_does_not_reconstruct_no_scratch(
        self,
        monkeypatch,
        single_molecule_xyz_file,
        gaussian_project_config_dir,
        captured,
    ):
        monkeypatch.setenv(
            "CHEMSMART_CONFIG_DIR", str(gaussian_project_config_dir)
        )
        fake_server = Server(name="dummy")
        captured["cli_args"] = None
        fake_server.submit = (
            lambda job, test=False, cli_args=None, **kw: captured.update(
                cli_args=cli_args
            )
        )
        monkeypatch.setattr(
            "chemsmart.settings.server.Server.from_servername",
            lambda _name: fake_server,
        )

        result = CliRunner().invoke(
            sub,
            [
                "--test",
                "--server",
                "dummy",
                "gaussian",
                "-p",
                "test",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "opt",
            ],
        )
        assert result.exit_code == 0, result.output
        assert "--no-scratch" not in captured["cli_args"]
        assert "--scratch" not in captured["cli_args"]


def _write_server_yaml(
    path, *, gaussian_scratch, orca_scratch, xtb_scratch=None
):
    """Write a minimal server YAML with optional program SCRATCH keys."""
    gaussian_line = (
        f"    SCRATCH: {gaussian_scratch}\n"
        if gaussian_scratch is not None
        else ""
    )
    orca_line = (
        f"    SCRATCH: {orca_scratch}\n" if orca_scratch is not None else ""
    )
    xtb_line = (
        f"    SCRATCH: {xtb_scratch}\n" if xtb_scratch is not None else ""
    )
    path.write_text(
        "SERVER:\n"
        "    SCHEDULER: PBS\n"
        "    MEM_GB: 8\n"
        "    NUM_CORES: 2\n"
        "    NUM_GPUS: 0\n"
        "    NUM_THREADS: 2\n"
        "    SUBMIT_COMMAND: qsub\n"
        "    SCRATCH_DIR: null\n"
        "GAUSSIAN:\n"
        "    EXEFOLDER: ~/programs/g16\n"
        "    LOCAL_RUN: True\n"
        f"{gaussian_line}"
        "    ENVARS: |\n"
        "        export SCRATCH=~/scratch\n"
        "ORCA:\n"
        "    EXEFOLDER: ~/programs/orca\n"
        "    LOCAL_RUN: False\n"
        f"{orca_line}"
        "    ENVARS: |\n"
        "        export SCRATCH=~/scratch\n"
        "XTB:\n"
        "    EXEFOLDER: null\n"
        "    LOCAL_RUN: True\n"
        f"{xtb_line}"
        "    ENVARS: |\n"
        "        export SCRATCH=~/scratch\n"
    )
    return path


class TestScratchYamlOverride:
    """CLI → YAML program SCRATCH → runner class SCRATCH."""

    def test_program_scratch_from_servername_reads_yaml(self, tmp_path):
        from chemsmart.settings.executable import (
            GaussianExecutable,
            ORCAExecutable,
            XTBExecutable,
        )

        yaml_path = _write_server_yaml(
            tmp_path / "mixed.yaml",
            gaussian_scratch=True,
            orca_scratch=False,
            xtb_scratch=True,
        )
        assert (
            GaussianExecutable.program_scratch_from_servername(str(yaml_path))
            is True
        )
        assert (
            ORCAExecutable.program_scratch_from_servername(str(yaml_path))
            is False
        )
        assert (
            XTBExecutable.program_scratch_from_servername(str(yaml_path))
            is True
        )

    def test_program_scratch_from_servername_missing_key_is_none(
        self, tmp_path
    ):
        from chemsmart.settings.executable import ORCAExecutable

        yaml_path = _write_server_yaml(
            tmp_path / "no_orca_scratch.yaml",
            gaussian_scratch=True,
            orca_scratch=None,
        )
        assert (
            ORCAExecutable.program_scratch_from_servername(str(yaml_path))
            is None
        )

    def test_omit_uses_yaml_false_over_orca_class_true(self, tmp_path):
        yaml_path = _write_server_yaml(
            tmp_path / "orca_off.yaml",
            gaussian_scratch=True,
            orca_scratch=False,
        )
        server = Server.from_yaml(str(yaml_path))
        runner = JobRunner.from_job(
            job=SimpleNamespace(TYPE="orcasp"),
            server=server,
            scratch=None,
            fake=True,
        )
        assert isinstance(runner, FakeORCAJobRunner)
        assert ORCAJobRunner.SCRATCH is True
        assert runner.scratch is False

    def test_omit_uses_yaml_true_for_gaussian(self, tmp_path):
        yaml_path = _write_server_yaml(
            tmp_path / "g16_on.yaml",
            gaussian_scratch=True,
            orca_scratch=False,
        )
        server = Server.from_yaml(str(yaml_path))
        scratch_dir = tmp_path / "scratch"
        scratch_dir.mkdir()
        runner = JobRunner.from_job(
            job=SimpleNamespace(TYPE="g16opt"),
            server=server,
            scratch=None,
            fake=True,
            scratch_dir=str(scratch_dir),
        )
        assert isinstance(runner, FakeGaussianJobRunner)
        assert runner.scratch is True

    def test_omit_falls_back_to_class_when_yaml_key_absent(self, tmp_path):
        yaml_path = _write_server_yaml(
            tmp_path / "orca_absent.yaml",
            gaussian_scratch=True,
            orca_scratch=None,
        )
        server = Server.from_yaml(str(yaml_path))
        scratch_dir = tmp_path / "scratch"
        scratch_dir.mkdir()
        runner = JobRunner.from_job(
            job=SimpleNamespace(TYPE="orcasp"),
            server=server,
            scratch=None,
            fake=True,
            scratch_dir=str(scratch_dir),
        )
        assert runner.scratch is ORCAJobRunner.SCRATCH

    def test_cli_scratch_true_overrides_yaml_false(self, tmp_path):
        yaml_path = _write_server_yaml(
            tmp_path / "orca_off.yaml",
            gaussian_scratch=True,
            orca_scratch=False,
        )
        server = Server.from_yaml(str(yaml_path))
        scratch_dir = tmp_path / "scratch"
        scratch_dir.mkdir()
        runner = JobRunner.from_job(
            job=SimpleNamespace(TYPE="orcasp"),
            server=server,
            scratch=True,
            fake=True,
            scratch_dir=str(scratch_dir),
        )
        assert runner.scratch is True

    def test_cli_scratch_false_overrides_yaml_true(self, tmp_path):
        yaml_path = _write_server_yaml(
            tmp_path / "g16_on.yaml",
            gaussian_scratch=True,
            orca_scratch=False,
        )
        server = Server.from_yaml(str(yaml_path))
        runner = JobRunner.from_job(
            job=SimpleNamespace(TYPE="g16opt"),
            server=server,
            scratch=False,
            fake=True,
        )
        assert runner.scratch is False

    def test_missing_scratch_dir_falls_back_to_job_folder(
        self, pbs_server, tmp_path, monkeypatch
    ):
        runner = FakeGaussianJobRunner(
            server=pbs_server, scratch=False, fake=True
        )
        runner.scratch = True
        runner._scratch_dir = None
        runner._set_scratch.cache_clear()
        fake_exe = SimpleNamespace(
            scratch_dir=str(tmp_path / "does_not_exist"),
            local_run=None,
        )
        monkeypatch.setattr(
            type(runner),
            "executable",
            property(lambda self: fake_exe),
        )
        runner._set_scratch()
        job = DummyGaussianJob(folder=tmp_path, label="gaussian_opt")
        runner._assign_variables(job)

        assert runner.scratch is False
        assert runner.running_directory == job.folder


class TestSubResourceOverrides:
    def test_cli_resources_reach_submission_server(
        self,
        monkeypatch,
        single_molecule_xyz_file,
        gaussian_project_config_dir,
        captured,
    ):
        monkeypatch.setenv(
            "CHEMSMART_CONFIG_DIR", str(gaussian_project_config_dir)
        )
        fake_server = Server(
            name="dummy",
            NUM_CORES=2,
            NUM_GPUS=0,
            MEM_GB=8,
            NUM_HOURS=1,
            QUEUE_NAME="normal",
        )

        def _capture_script(self, job, cli_args, **kwargs):
            captured.update(
                server=self,
                num_cores=self.num_cores,
                num_gpus=self.num_gpus,
                mem_gb=self.mem_gb,
                num_hours=self.num_hours,
                queue_name=self.queue_name,
            )

        monkeypatch.setattr(
            Server,
            "_check_running_jobs",
            lambda self, job: None,
        )
        monkeypatch.setattr(
            Server,
            "_write_submission_script",
            _capture_script,
        )
        monkeypatch.setattr(
            Server,
            "current",
            classmethod(lambda cls: fake_server),
        )

        result = CliRunner().invoke(
            sub,
            [
                "--test",
                "--num-cores",
                "7",
                "--num-gpus",
                "2",
                "--mem-gb",
                "23",
                "--time-hours",
                "12.5",
                "--queue",
                "debug",
                "gaussian",
                "-p",
                "test",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "opt",
            ],
        )

        assert result.exit_code == 0, result.output
        assert captured["server"] is fake_server
        assert captured["num_cores"] == 7
        assert captured["num_gpus"] == 2
        assert captured["mem_gb"] == 23
        assert captured["num_hours"] == 12.5
        assert captured["queue_name"] == "debug"
