from pathlib import Path
from types import SimpleNamespace

import pytest
from click.testing import CliRunner

from chemsmart.cli.run import run
from chemsmart.cli.sub import sub
from chemsmart.jobs.gaussian.runner import (
    FakeGaussian,
    FakeGaussianJobRunner,
    GaussianJobRunner,
)
from chemsmart.jobs.iterate.runner import IterateJobRunner
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

    @property
    def outputfile(self):
        return str(Path(self.folder) / f"{self.label}.log")


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
    ):
        monkeypatch.setenv(
            "CHEMSMART_CONFIG_DIR", str(gaussian_project_config_dir)
        )
        fake_server = Server(name="dummy")
        captured = {"cli_args": None}
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

    def test_sub_time_hours_and_queue_override_server_settings(
        self,
        monkeypatch,
        single_molecule_xyz_file,
        gaussian_project_config_dir,
    ):
        monkeypatch.setenv(
            "CHEMSMART_CONFIG_DIR", str(gaussian_project_config_dir)
        )
        fake_server = Server(name="dummy")
        fake_server.submit = lambda job, test=False, cli_args=None, **kw: None
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
                "--time-hours",
                "48.0",
                "--queue",
                "gpu",
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
        assert fake_server.num_hours == 48.0
        assert fake_server.queue_name == "gpu"

    def test_sub_verbose_flag_enables_stream_and_debug_logging(
        self,
        monkeypatch,
        single_molecule_xyz_file,
        gaussian_project_config_dir,
    ):
        monkeypatch.setenv(
            "CHEMSMART_CONFIG_DIR", str(gaussian_project_config_dir)
        )
        fake_server = Server(name="dummy")
        fake_server.submit = lambda job, test=False, cli_args=None, **kw: None
        monkeypatch.setattr(
            "chemsmart.settings.server.Server.from_servername",
            lambda _name: fake_server,
        )
        captured = {}
        monkeypatch.setattr(
            "chemsmart.cli.sub.create_logger",
            lambda **kw: captured.update(kw),
        )

        result = CliRunner().invoke(
            sub,
            [
                "--test",
                "--server",
                "dummy",
                "--verbose",
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
        assert captured == {"stream": True, "debug": True}


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


class TestGaussianJobRunnerDirectHelpers:
    """Direct unit coverage for GaussianJobRunner's file-management and
    execution-plumbing helpers, which are otherwise only exercised
    indirectly (and incompletely) through full job-run integration
    tests elsewhere."""

    def test_init_defaults_scratch_to_class_attribute(self, pbs_server):
        runner = GaussianJobRunner(server=pbs_server, scratch=None, fake=True)
        assert runner.scratch is GaussianJobRunner.SCRATCH

    def test_executable_property_wraps_filenotfounderror(
        self, pbs_server, monkeypatch
    ):
        runner = GaussianJobRunner(server=pbs_server, scratch=False, fake=True)
        runner.server = SimpleNamespace(name="no-such-server")

        from chemsmart.jobs.gaussian import runner as runner_module

        def _raise(servername):
            raise FileNotFoundError(f"no config for {servername}")

        monkeypatch.setattr(
            runner_module.GaussianExecutable,
            "from_servername",
            staticmethod(_raise),
        )
        with pytest.raises(FileNotFoundError):
            runner.executable

    def test_assign_variables_dispatches_to_scratch_setup(
        self, pbs_server, tmp_path, monkeypatch
    ):
        runner = GaussianJobRunner(
            server=pbs_server,
            scratch=True,
            scratch_dir=str(tmp_path),
            fake=True,
        )
        calls = []
        monkeypatch.setattr(
            runner,
            "_set_up_variables_in_scratch",
            lambda job: calls.append("scratch"),
        )
        monkeypatch.setattr(
            runner,
            "_set_up_variables_in_job_directory",
            lambda job: calls.append("job_directory"),
        )
        monkeypatch.setattr(
            type(runner),
            "executable",
            property(lambda self: SimpleNamespace(local_run=None)),
        )
        job = DummyGaussianJob(folder=tmp_path, label="assign_test")
        runner._assign_variables(job)
        assert calls == ["scratch"]
        assert not hasattr(job, "local")

    def test_assign_variables_dispatches_to_job_directory_setup(
        self, pbs_server, tmp_path, monkeypatch
    ):
        runner = GaussianJobRunner(server=pbs_server, scratch=False, fake=True)
        calls = []
        monkeypatch.setattr(
            runner,
            "_set_up_variables_in_scratch",
            lambda job: calls.append("scratch"),
        )
        monkeypatch.setattr(
            runner,
            "_set_up_variables_in_job_directory",
            lambda job: calls.append("job_directory"),
        )
        monkeypatch.setattr(
            type(runner),
            "executable",
            property(lambda self: SimpleNamespace(local_run=True)),
        )
        job = DummyGaussianJob(folder=tmp_path, label="assign_test")
        runner._assign_variables(job)
        assert calls == ["job_directory"]
        assert job.local is True

    def test_set_up_variables_in_scratch_creates_fresh_directory(
        self, pbs_server, tmp_path
    ):
        runner = GaussianJobRunner(
            server=pbs_server,
            scratch=True,
            scratch_dir=str(tmp_path / "scratch_root"),
            fake=True,
        )
        job = DummyGaussianJob(
            folder=tmp_path / "jobfolder", label="mytest_job"
        )
        runner._set_up_variables_in_scratch(job)

        import os

        expected_dir = os.path.join(str(tmp_path / "scratch_root"), job.label)
        assert os.path.isdir(expected_dir)
        assert runner.running_directory == expected_dir
        assert runner.job_inputfile == os.path.abspath(
            os.path.join(expected_dir, f"{job.label}.com")
        )
        assert runner.job_chkfile == os.path.abspath(
            os.path.join(expected_dir, f"{job.label}.chk")
        )
        assert runner.job_errfile == os.path.abspath(
            os.path.join(expected_dir, f"{job.label}.err")
        )
        assert runner.job_outputfile == os.path.abspath(
            os.path.join(expected_dir, f"{job.label}.log")
        )

    def test_set_up_variables_in_scratch_reuses_existing_directory(
        self, pbs_server, tmp_path
    ):
        import os

        scratch_root = tmp_path / "scratch_root2"
        job = DummyGaussianJob(
            folder=tmp_path / "jobfolder2", label="existing_job"
        )
        os.makedirs(os.path.join(str(scratch_root), job.label))

        runner = GaussianJobRunner(
            server=pbs_server,
            scratch=True,
            scratch_dir=str(scratch_root),
            fake=True,
        )
        runner._set_up_variables_in_scratch(job)
        assert runner.running_directory == os.path.join(
            str(scratch_root), job.label
        )

    def test_set_up_variables_in_job_directory(self, pbs_server, tmp_path):
        import os

        runner = GaussianJobRunner(server=pbs_server, scratch=False, fake=True)
        job = DummyGaussianJob(folder=tmp_path / "jobdir", label="direct_job")
        runner._set_up_variables_in_job_directory(job)

        assert runner.running_directory == job.folder
        assert runner.job_inputfile == os.path.abspath(job.inputfile)
        assert runner.job_chkfile == os.path.abspath(job.chkfile)
        assert runner.job_errfile == os.path.abspath(job.errfile)
        assert runner.job_outputfile == os.path.abspath(job.outputfile)

    def test_prerun_delegates_to_assign_variables(
        self, pbs_server, tmp_path, monkeypatch
    ):
        runner = GaussianJobRunner(server=pbs_server, scratch=False, fake=True)
        calls = []
        monkeypatch.setattr(
            runner, "_assign_variables", lambda job: calls.append(job)
        )
        job = DummyGaussianJob(folder=tmp_path, label="prerun_job")
        runner._prerun(job)
        assert calls == [job]

    def test_write_input_uses_gaussian_input_writer(
        self, pbs_server, tmp_path, monkeypatch
    ):
        runner = GaussianJobRunner(server=pbs_server, scratch=False, fake=True)
        runner.running_directory = str(tmp_path)
        job = DummyGaussianJob(folder=tmp_path, label="write_input_job")

        captured = {}

        class FakeInputWriter:
            def __init__(self, job):
                captured["job"] = job

            def write(self, target_directory):
                captured["target_directory"] = target_directory

        monkeypatch.setattr(
            "chemsmart.jobs.gaussian.writer.GaussianInputWriter",
            FakeInputWriter,
        )
        runner._write_input(job)
        assert captured["job"] is job
        assert captured["target_directory"] == str(tmp_path)

    def test_get_command_combines_executable_and_inputfile(
        self, pbs_server, monkeypatch
    ):
        runner = GaussianJobRunner(server=pbs_server, scratch=False, fake=True)
        monkeypatch.setattr(runner, "_get_executable", lambda: "/opt/g16")
        runner.job_inputfile = "/scratch/job/mytest.com"
        assert (
            runner._get_command(job=None) == "/opt/g16 /scratch/job/mytest.com"
        )

    def test_get_executable_delegates_to_executable_object(
        self, pbs_server, monkeypatch
    ):
        runner = GaussianJobRunner(server=pbs_server, scratch=False, fake=True)
        monkeypatch.setattr(
            type(runner),
            "executable",
            property(
                lambda self: SimpleNamespace(
                    get_executable=lambda: "/usr/bin/g16"
                )
            ),
        )
        assert runner._get_executable() == "/usr/bin/g16"

    def test_create_process_opens_files_and_spawns_popen(
        self, pbs_server, tmp_path, monkeypatch
    ):
        runner = GaussianJobRunner(server=pbs_server, scratch=False, fake=True)
        runner.running_directory = str(tmp_path)
        runner.job_outputfile = str(tmp_path / "job.log")
        runner.job_errfile = str(tmp_path / "job.err")
        monkeypatch.setattr(
            type(runner),
            "executable",
            property(lambda self: SimpleNamespace(env={"PATH": "/bin"})),
        )

        captured = {}

        class FakeProcess:
            pass

        def fake_popen(args, stdout, stderr, env, cwd):
            captured["args"] = args
            captured["env"] = env
            captured["cwd"] = cwd
            captured["stdout_closed"] = stdout.closed
            captured["stderr_closed"] = stderr.closed
            return FakeProcess()

        from chemsmart.jobs.gaussian import runner as runner_module

        monkeypatch.setattr(runner_module.subprocess, "Popen", fake_popen)

        result = runner._create_process(
            job=None, command="g16 job.com", env={"PATH": "/bin"}
        )

        assert isinstance(result, FakeProcess)
        assert captured["args"] == ["g16", "job.com"]
        assert captured["cwd"] == str(tmp_path)
        assert (tmp_path / "job.log").exists()
        assert (tmp_path / "job.err").exists()

    def test_postrun_noop_when_scratch_disabled(self, pbs_server, tmp_path):
        runner = GaussianJobRunner(server=pbs_server, scratch=False, fake=True)
        job = DummyGaussianJob(folder=tmp_path, label="noop_job")
        runner._postrun(job)  # should not raise, nothing to copy

    def test_postrun_copies_matching_files_to_job_folder(
        self, pbs_server, tmp_path
    ):
        import os

        running_dir = tmp_path / "scratch_run"
        running_dir.mkdir()
        job_folder = tmp_path / "job_folder"
        job_folder.mkdir()

        (running_dir / "copy_job.log").write_text("fake log contents")

        runner = GaussianJobRunner(server=pbs_server, scratch=True, fake=True)
        runner.running_directory = str(running_dir)
        job = DummyGaussianJob(folder=job_folder, label="copy_job")

        runner._postrun(job)

        assert (job_folder / "copy_job.log").exists()
        assert (job_folder / "copy_job.log").read_text() == (
            "fake log contents"
        )
        assert os.path.exists(str(running_dir / "copy_job.log"))

    def test_postrun_skips_files_whose_full_path_starts_with_gau_prefix(
        self, pbs_server, tmp_path, monkeypatch
    ):
        """Defensive branch: the "except files starting with Gau-" check
        tests the full glob'd path (not the basename), so it only ever
        skips a file when the running_directory string itself begins
        with "Gau-" (an edge case, not Gaussian's own scratch files,
        which never match the job.label* glob pattern in the first
        place)."""
        monkeypatch.chdir(tmp_path)
        import os

        os.makedirs("Gau-scratch")
        with open(os.path.join("Gau-scratch", "prefixed_job.log"), "w") as f:
            f.write("should not be copied")
        job_folder = tmp_path / "job_folder"
        job_folder.mkdir()

        runner = GaussianJobRunner(server=pbs_server, scratch=True, fake=True)
        runner.running_directory = "Gau-scratch"
        job = DummyGaussianJob(folder=job_folder, label="prefixed_job")

        runner._postrun(job)

        assert not (job_folder / "prefixed_job.log").exists()

    def test_postrun_logs_and_continues_on_copy_failure(
        self, pbs_server, tmp_path, monkeypatch
    ):
        running_dir = tmp_path / "scratch_run2"
        running_dir.mkdir()
        job_folder = tmp_path / "job_folder2"
        job_folder.mkdir()
        (running_dir / "fail_job.log").write_text("contents")

        from chemsmart.jobs.gaussian import runner as runner_module

        def _raise_copy(src, dst):
            raise OSError("disk full")

        monkeypatch.setattr(runner_module, "copy", _raise_copy)

        runner = GaussianJobRunner(server=pbs_server, scratch=True, fake=True)
        runner.running_directory = str(running_dir)
        job = DummyGaussianJob(folder=job_folder, label="fail_job")

        runner._postrun(job)  # should not raise; failure is logged

        assert not (job_folder / "fail_job.log").exists()


class TestFakeGaussianJobRunnerDirectHelpers:
    def test_set_up_variables_in_scratch_reuses_existing_directory(
        self, pbs_server, tmp_path
    ):
        import os

        scratch_root = tmp_path / "fake_scratch_root"
        job = DummyGaussianJob(
            folder=tmp_path / "jobfolder", label="fake_existing_job"
        )
        os.makedirs(os.path.join(str(scratch_root), job.label))

        expected_dir = os.path.join(str(scratch_root), job.label)
        runner = FakeGaussianJobRunner(
            server=pbs_server,
            scratch=True,
            scratch_dir=str(scratch_root),
            fake=True,
        )
        # _set_up_variables_in_scratch appends "_fake" to job.label as a
        # side effect, so the expected directory must be captured first.
        runner._set_up_variables_in_scratch(job)
        assert runner.running_directory == expected_dir

    def test_set_up_variables_in_scratch_creates_fresh_directory(
        self, pbs_server, tmp_path
    ):
        import os

        scratch_root = tmp_path / "fake_scratch_root_fresh"
        job = DummyGaussianJob(
            folder=tmp_path / "jobfolder_fresh", label="fake_fresh_job"
        )
        expected_dir = os.path.join(str(scratch_root), job.label)

        runner = FakeGaussianJobRunner(
            server=pbs_server,
            scratch=True,
            scratch_dir=str(scratch_root),
            fake=True,
        )
        runner._set_up_variables_in_scratch(job)
        assert os.path.isdir(expected_dir)
        assert runner.running_directory == expected_dir

    def test_run_executes_full_fake_pipeline_in_order(
        self, pbs_server, tmp_path, monkeypatch
    ):
        runner = FakeGaussianJobRunner(
            server=pbs_server, scratch=False, fake=True
        )
        calls = []
        monkeypatch.setattr(
            runner, "_prerun", lambda job: calls.append("prerun")
        )
        monkeypatch.setattr(
            runner, "_write_input", lambda job: calls.append("write_input")
        )
        monkeypatch.setattr(
            runner, "_postrun", lambda job: calls.append("postrun")
        )
        monkeypatch.setattr(
            runner,
            "_postrun_cleanup",
            lambda job: calls.append("postrun_cleanup"),
        )
        runner.job_inputfile = "unused.com"

        class FakeFakeGaussian:
            def __init__(self, file_to_run):
                calls.append(f"init:{file_to_run}")

            def run(self):
                calls.append("fake_run")
                return 0

        from chemsmart.jobs.gaussian import runner as runner_module

        monkeypatch.setattr(runner_module, "FakeGaussian", FakeFakeGaussian)

        job = DummyGaussianJob(folder=tmp_path, label="pipeline_job")
        result = runner.run(job)

        assert result == 0
        assert calls == [
            "prerun",
            "write_input",
            "init:unused.com",
            "fake_run",
            "postrun",
            "postrun_cleanup",
        ]


def _write_gaussian_com(
    tmp_path, filename="test.com", route="# opt freq b3lyp/6-31g(d)", mult=1
):
    lines = [
        "%chk=test.chk",
        "%mem=4GB",
        route,
        "",
        "title",
        "",
        f"0 {mult}",
        "C 0.0 0.0 0.0",
        "H 0.0 0.0 1.0",
        "",
    ]
    path = tmp_path / filename
    path.write_text("\n".join(lines))
    return str(path)


class TestFakeGaussianDirect:
    """Direct unit coverage for the FakeGaussian simulator (used by
    FakeGaussianJobRunner.run) that isn't reached by the higher-level
    fake-job-run tests elsewhere: the missing-file error, odd-multiplicity
    spin, extra input blocks, and the semiempirical/ab-initio/unknown-method
    branches of the fake SCF-energy line."""

    def test_missing_input_file_raises_filenotfounderror(self):
        with pytest.raises(FileNotFoundError):
            FakeGaussian("/no/such/path/does_not_exist.com")

    def test_spin_is_u_for_non_singlet_multiplicity(self, tmp_path):
        path = _write_gaussian_com(tmp_path, mult=2)
        fake = FakeGaussian(path)
        assert fake.spin == "U"

    def test_run_writes_additional_input_blocks(self, tmp_path):
        path = _write_gaussian_com(tmp_path, route="# opt freq b3lyp/gen")
        # append a gen/genecp-style basis-set block (content_groups[3:]),
        # separated from the coordinate block by a blank line.
        with open(path, "a") as f:
            f.write("\nC 0\n6-31g(d)\n****\n\n")

        fake = FakeGaussian(path)
        assert len(fake.input_blocks) > 3
        fake.run()
        output = Path(fake.output_filepath).read_text()
        assert "6-31g(d)" in output

    def test_run_uses_semiempirical_energy_line(self, tmp_path):
        path = _write_gaussian_com(tmp_path, route="# opt freq pm6")
        fake = FakeGaussian(path)
        fake.run()
        output = Path(fake.output_filepath).read_text()
        assert "E(PM6)" in output

    def test_run_uses_ab_initio_energy_line(self, tmp_path):
        # A bare ab initio method with no "/basis" suffix is required so
        # that .functional stays None (any "method/basis" token is always
        # parsed as a functional first, regardless of the method).
        path = _write_gaussian_com(tmp_path, route="# opt freq hf")
        fake = FakeGaussian(path)
        assert fake.input_object.functional is None
        assert fake.input_object.ab_initio == "hf"
        fake.run()
        output = Path(fake.output_filepath).read_text()
        assert "E(hf)" in output

    def test_run_uses_unknown_method_energy_line(self, tmp_path):
        path = _write_gaussian_com(tmp_path, route="# opt freq uff")
        fake = FakeGaussian(path)
        assert fake.input_object.functional is None
        assert fake.input_object.semiempirical is None
        assert fake.input_object.ab_initio is None
        fake.run()
        output = Path(fake.output_filepath).read_text()
        assert "E(Unknow Method)" in output


class TestSubResourceOverrides:
    def test_cli_resources_reach_submission_server(
        self,
        monkeypatch,
        single_molecule_xyz_file,
        gaussian_project_config_dir,
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
        captured = {}

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
