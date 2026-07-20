from pathlib import Path
from types import SimpleNamespace

from click.testing import CliRunner

from chemsmart.cli.run import run
from chemsmart.cli.sub import sub
from chemsmart.jobs.gaussian.runner import (
    FakeGaussianJobRunner,
    GaussianJobRunner,
)
from chemsmart.jobs.iterate.runner import IterateJobRunner
from chemsmart.jobs.orca.runner import FakeORCAJobRunner, ORCAJobRunner
from chemsmart.jobs.runner import JobRunner
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
            job=job,
            server=pbs_server,
            scratch=False,
            scratch_from_cli=True,
            fake=True,
        )
        assert isinstance(runner, FakeGaussianJobRunner)
        assert runner.fake is True

    def test_fake_orca_runner_selected_when_fake_enabled(self, pbs_server):
        job = SimpleNamespace(TYPE="orcasp")
        runner = JobRunner.from_job(
            job=job,
            server=pbs_server,
            scratch=False,
            scratch_from_cli=True,
            fake=True,
        )
        assert isinstance(runner, FakeORCAJobRunner)
        assert runner.fake is True

    def test_real_runner_selected_when_fake_disabled(self, pbs_server):
        gaussian_job = SimpleNamespace(TYPE="g16opt")
        gaussian_runner = JobRunner.from_job(
            job=gaussian_job,
            server=pbs_server,
            scratch=False,
            scratch_from_cli=True,
            fake=False,
        )
        assert isinstance(gaussian_runner, GaussianJobRunner)

        orca_job = SimpleNamespace(TYPE="orcasp")
        orca_runner = JobRunner.from_job(
            job=orca_job,
            server=pbs_server,
            scratch=False,
            scratch_from_cli=True,
            fake=False,
        )
        assert isinstance(orca_runner, ORCAJobRunner)

    def test_fake_flag_propagates_when_no_fake_runner_exists(self, pbs_server):
        job = SimpleNamespace(TYPE="iterate")
        runner = JobRunner.from_job(
            job=job,
            server=pbs_server,
            scratch=False,
            scratch_from_cli=True,
            fake=True,
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


class TestScratchCLIDefaults:
    """CLI scratch flag and from_job scratch_from_cli policy."""

    @staticmethod
    def _write_gaussian_project(tmp_path):
        config_root = tmp_path / "chemsmart_cfg"
        gaussian_cfg = config_root / "gaussian"
        gaussian_cfg.mkdir(parents=True)
        (gaussian_cfg / "test.yaml").write_text(
            "gas:\n  functional: B3LYP\n  basis: def2-SVP\n"
            "solv:\n  functional: B3LYP\n  basis: def2-SVP\n"
            "  solvent_model: smd\n  solvent_id: water\n"
        )
        return config_root

    def test_run_omitted_scratch_passes_from_cli_false(
        self, tmp_path, monkeypatch, single_molecule_xyz_file
    ):
        monkeypatch.setenv(
            "CHEMSMART_CONFIG_DIR",
            str(self._write_gaussian_project(tmp_path)),
        )
        observed = {}

        def _from_job(
            cls,
            job,
            server,
            scratch=False,
            fake=False,
            scratch_from_cli=False,
            **kwargs,
        ):
            observed["scratch"] = scratch
            observed["scratch_from_cli"] = scratch_from_cli
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
        assert observed["scratch"] is False
        assert observed["scratch_from_cli"] is False

    def test_sub_omitted_scratch_does_not_force_no_scratch(
        self, tmp_path, monkeypatch, single_molecule_xyz_file
    ):
        monkeypatch.setenv(
            "CHEMSMART_CONFIG_DIR",
            str(self._write_gaussian_project(tmp_path)),
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

    def test_run_explicit_scratch_is_not_silently_disabled(
        self, tmp_path, monkeypatch, single_molecule_xyz_file
    ):
        """Placeholder must not clear --scratch before from_job."""
        monkeypatch.setenv(
            "CHEMSMART_CONFIG_DIR",
            str(self._write_gaussian_project(tmp_path)),
        )
        observed = {}

        def _from_job(
            cls,
            job,
            server,
            scratch=False,
            fake=False,
            scratch_from_cli=False,
            **kwargs,
        ):
            observed["scratch"] = scratch
            observed["scratch_from_cli"] = scratch_from_cli
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
        assert observed["scratch"] is True
        assert observed["scratch_from_cli"] is True

    def test_sub_explicit_scratch_keeps_jobrunner_scratch_true(
        self, tmp_path, monkeypatch, single_molecule_xyz_file
    ):
        monkeypatch.setenv(
            "CHEMSMART_CONFIG_DIR",
            str(self._write_gaussian_project(tmp_path)),
        )
        fake_server = Server(name="dummy")
        captured = {"job": None, "cli_args": None}

        def _submit(job, test=False, cli_args=None, **kw):
            captured["job"] = job
            captured["cli_args"] = cli_args

        fake_server.submit = _submit
        monkeypatch.setattr(
            "chemsmart.settings.server.Server.from_servername",
            lambda _name: fake_server,
        )

        result = CliRunner().invoke(
            sub,
            [
                "--test",
                "--scratch",
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
        assert "--scratch" in captured["cli_args"]
        assert captured["job"].jobrunner.scratch is True
        assert captured["job"].jobrunner.scratch_from_cli is True

    def test_from_job_omitted_uses_yaml_scratch(self, pbs_server, tmp_path):
        scratch_dir = tmp_path / "scratch"
        scratch_dir.mkdir()
        job = SimpleNamespace(TYPE="g16opt")
        runner = JobRunner.from_job(
            job=job,
            server=pbs_server,
            scratch=False,
            scratch_from_cli=False,
            fake=True,
            scratch_dir=str(scratch_dir),
        )
        assert runner.scratch_from_cli is False
        assert runner.scratch is True

    def test_from_job_omitted_orca_uses_yaml_scratch_false(self, pbs_server):
        job = SimpleNamespace(TYPE="orcaopt")
        runner = JobRunner.from_job(
            job=job,
            server=pbs_server,
            scratch=False,
            scratch_from_cli=False,
            fake=True,
        )
        assert runner.scratch_from_cli is False
        assert runner.scratch is False

    def test_from_job_explicit_no_scratch_ignores_program_default(
        self, pbs_server
    ):
        job = SimpleNamespace(TYPE="g16opt")
        runner = JobRunner.from_job(
            job=job,
            server=pbs_server,
            scratch=False,
            scratch_from_cli=True,
            fake=True,
        )
        assert runner.scratch is False
        assert runner.scratch_from_cli is True


class TestScratchPathFallback:
    """Invalid configured scratch paths fall back to job.folder."""

    def test_set_scratch_missing_path_disables_scratch(
        self, pbs_server, tmp_path, monkeypatch
    ):
        runner = FakeGaussianJobRunner(
            server=pbs_server, scratch=False, fake=True
        )
        runner.scratch = True
        runner._scratch_dir = None
        runner._set_scratch.cache_clear()
        fake_exe = SimpleNamespace(
            scratch_dir=str(tmp_path / "does_not_exist")
        )
        monkeypatch.setattr(
            type(runner),
            "executable",
            property(lambda self: fake_exe),
        )

        assert runner._set_scratch() is None
        assert runner.scratch is False

    def test_set_scratch_non_directory_disables_scratch(
        self, pbs_server, tmp_path, monkeypatch
    ):
        not_a_dir = tmp_path / "scratch.txt"
        not_a_dir.write_text("x")
        runner = FakeGaussianJobRunner(
            server=pbs_server, scratch=False, fake=True
        )
        runner.scratch = True
        runner._scratch_dir = None
        runner._set_scratch.cache_clear()
        fake_exe = SimpleNamespace(scratch_dir=str(not_a_dir))
        monkeypatch.setattr(
            type(runner),
            "executable",
            property(lambda self: fake_exe),
        )

        assert runner._set_scratch() is None
        assert runner.scratch is False

    def test_invalid_scratch_path_uses_job_folder(
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
