import os
from io import StringIO
from unittest.mock import MagicMock, patch

import pytest

from chemsmart.jobs.job import Job
from chemsmart.settings.executable import (
    CRESTExecutable,
    Executable,
    GaussianExecutable,
    NCIPLOTExecutable,
    ORCAExecutable,
    XTBExecutable,
)
from chemsmart.settings.server import (
    LSFServer,
    PBSServer,
    Server,
    ServerSettingsManager,
    SGE_Server,
    SLURMServer,
    YamlServerSettings,
)
from chemsmart.settings.submitters import PBSSubmitter, SLURMSubmitter


@pytest.fixture
def _clear_scheduler_caches():
    Server.detect_server_scheduler.cache_clear()
    Server.from_scheduler_type.cache_clear()
    yield
    Server.detect_server_scheduler.cache_clear()
    Server.from_scheduler_type.cache_clear()


class TestServer:
    @pytest.mark.parametrize(
        "program,label,queued_names,duplicate",
        [
            ("ORCA", "molecule_opt", ["orca_molecule_opt"], True),
            ("Gaussian", "molecule_opt", ["gaussian_molecule_opt"], True),
            ("ORCA", "molecule_opt", ["gaussian_molecule_opt"], False),
            ("ORCA", "molecule_opt", [], False),
            ("ORCA", None, [], False),
        ],
        ids=[
            "orca-duplicate",
            "gaussian-duplicate",
            "other-program",
            "empty-queue",
            "no-label",
        ],
    )
    def test_check_running_jobs(
        self, mocker, program, label, queued_names, duplicate
    ):
        """Reject queued names only when both program and label match.

        Mock the queue lookup so no scheduler is needed, and verify unlabeled
        jobs skip the lookup entirely.
        """
        helper_cls = mocker.patch("chemsmart.utils.cluster.ClusterHelper")
        query = helper_cls.return_value.get_running_jobs
        query.return_value = (list(range(len(queued_names))), queued_names)
        job = Job(molecule=None, label=label, jobrunner=None)
        job.PROGRAM = program

        if duplicate:
            with pytest.raises(
                SystemExit, match="Duplicate job NOT submitted"
            ):
                Server._check_running_jobs(job)
        else:
            Server._check_running_jobs(job)

        if label is None:
            helper_cls.assert_not_called()
        else:
            query.assert_called_once_with()
        assert job.label == label

    def test_server_yaml(self, server_yaml_file):
        assert os.path.exists(server_yaml_file)
        assert os.path.isfile(server_yaml_file)
        server = Server.from_yaml(name=server_yaml_file)
        assert server.scheduler.lower() == "pbs"
        assert server.queue_name == "normal"
        assert server.num_hours == 24
        assert server.mem_gb == 375
        assert server.num_cores == 64
        assert server.num_gpus == 0
        assert server.num_threads == 64
        assert server.submit_command == "qsub"
        assert server.scratch_dir is None
        assert server.use_hosts is True
        assert (
            server.extra_commands == """export PATH=$HOME/bin/chemsmart:$PATH
export PATH=$HOME/bin/chemsmart/chemsmart/cli:$PATH
export PATH=$HOME/bin/chemsmart/chemsmart/scripts:$PATH
export PYTHONPATH=$HOME/bin/chemsmart:$PYTHONPATH
"""
        )
        assert server.extra_scheduler_directives == "#PBS -m abe\n"

    def test_gaussian_executable(self, server_yaml_file):
        gaussian_executable = GaussianExecutable.from_servername(
            server_yaml_file
        )
        assert gaussian_executable.executable_folder == os.path.expanduser(
            "~/programs/g16"
        )
        assert gaussian_executable.local_run is True

        gaussian_conda_env = """source ~/anaconda3/etc/profile.d/conda.sh
conda activate ~/anaconda3/envs/chemsmart
"""
        assert gaussian_executable.conda_env == gaussian_conda_env

        gaussian_modules = """module purge
module load craype-x86-rome
module load libfabric/1.11.0.4.125
"""
        assert gaussian_executable.modules == gaussian_modules

        assert (
            gaussian_executable.scripts
            == 'tcsh -c "source ~/programs/g16/bsd/g16.login"\n'
        )

        gassian_envars = """export SCRATCH=~/scratch
export GAUSS_EXEDIR=~/programs/g16
export g16root=~/programs/g16

"""
        assert gaussian_executable.envars == gassian_envars

    def test_orca_executable(self, server_yaml_file):
        orca_executable = ORCAExecutable.from_servername(server_yaml_file)
        assert orca_executable.executable_folder == os.path.expanduser(
            "~/programs/orca_6_0_0"
        )
        assert orca_executable.local_run is False

        assert orca_executable.conda_env is None

        assert orca_executable.modules is None

        assert orca_executable.scripts is None

        orca_envars = """export PATH=~/programs/openmpi-4.1.6/build/bin:$PATH
export LD_LIBRARY_PATH=~/programs/openmpi-4.1.6/build/lib:$LD_LIBRARY_PATH
"""
        assert orca_executable.envars == orca_envars

    def test_xtb_executable(self, server_yaml_file):
        xtb_executable = XTBExecutable.from_servername(server_yaml_file)
        assert xtb_executable.executable_folder is None
        assert xtb_executable.get_executable() == "xtb"
        assert xtb_executable.local_run is True

        xtb_conda_env = """source ~/anaconda3/etc/profile.d/conda.sh
conda activate ~/anaconda3/envs/chemsmart
"""
        assert xtb_executable.conda_env == xtb_conda_env
        assert xtb_executable.modules is None
        assert xtb_executable.scripts is None
        assert xtb_executable.envars is None

    def test_xtb_executable_with_explicit_folder(self):
        xtb_executable = XTBExecutable(executable_folder="/opt/xtb/bin")
        assert xtb_executable.get_executable() == "/opt/xtb/bin/xtb"

    def test_crest_executable(self, server_yaml_file):
        crest_executable = CRESTExecutable.from_servername(server_yaml_file)
        assert crest_executable.executable_folder is None
        assert crest_executable.get_executable() == "crest"
        assert crest_executable.local_run is True

        crest_conda_env = """source ~/anaconda3/etc/profile.d/conda.sh
conda activate ~/anaconda3/envs/chemsmart
"""
        assert crest_executable.conda_env == crest_conda_env
        assert crest_executable.modules is None
        assert crest_executable.scripts is None
        assert crest_executable.envars is None

    def test_crest_executable_with_explicit_folder(self):
        crest_executable = CRESTExecutable(executable_folder="/opt/crest/bin")
        assert crest_executable.get_executable() == "/opt/crest/bin/crest"

    def test_slurm_submitter_writes_extra_scheduler_directives(self):
        server = Server(
            "custom-slurm",
            SCHEDULER="SLURM",
            NUM_CORES=8,
            MEM_GB=24,
            NUM_GPUS=0,
            EXTRA_SCHEDULER_DIRECTIVES="#SBATCH --reservation=xlzhang_1\n",
        )
        job = type(
            "DummyJob", (), {"label": "job1", "execution_name": "orca_job1"}
        )()
        submitter = SLURMSubmitter(job=job, server=server)

        buffer = StringIO()
        submitter._write_scheduler_options(buffer)
        assert "#SBATCH --reservation=xlzhang_1\n" in buffer.getvalue()
        assert "#SBATCH --job-name=orca_job1\n" in buffer.getvalue()

    def test_pbs_submitter_writes_extra_scheduler_directives(self):
        server = Server(
            "custom-pbs",
            SCHEDULER="PBS",
            NUM_CORES=8,
            MEM_GB=24,
            NUM_GPUS=0,
            EXTRA_SCHEDULER_DIRECTIVES="#PBS -m abe\n",
        )
        job = type(
            "DummyJob", (), {"label": "job1", "execution_name": "orca_job1"}
        )()
        submitter = PBSSubmitter(job=job, server=server)

        buffer = StringIO()
        submitter._write_scheduler_options(buffer)
        assert "#PBS -m abe\n" in buffer.getvalue()
        assert "#PBS -N orca_job1\n" in buffer.getvalue()


class TestMissingProgramSectionFallback:
    """Tests that Executable.from_servername raises ValueError when the
    program block (e.g. XTB, CREST) is absent from the server YAML, and
    that programs which ARE present still parse correctly.

    Covers users who installed CHEMSMART before xTB support was added."""

    def test_xtb_missing_section_raises(self, legacy_server_yaml):
        """XTBExecutable.from_servername raises ValueError when XTB block absent."""
        import pytest

        with pytest.raises(ValueError, match="XTB"):
            XTBExecutable.from_servername(legacy_server_yaml)

    def test_xtb_missing_section_error_message(self, legacy_server_yaml):
        """The ValueError message mentions updating the server YAML."""
        import pytest

        with pytest.raises(ValueError, match="server YAML"):
            XTBExecutable.from_servername(legacy_server_yaml)

    def test_crest_missing_section_raises(self, legacy_server_yaml):
        """CRESTExecutable.from_servername raises ValueError when CREST block absent."""
        import pytest

        with pytest.raises(ValueError, match="CREST"):
            CRESTExecutable.from_servername(legacy_server_yaml)

    def test_present_section_still_parsed_correctly(self, legacy_server_yaml):
        """Programs that *are* present in the legacy YAML are still parsed."""
        exe = GaussianExecutable.from_servername(legacy_server_yaml)
        assert exe.executable_folder == os.path.expanduser("~/programs/g16")
        assert exe.local_run is True


@pytest.mark.usefixtures("_clear_scheduler_caches")
class TestServerDunderMethods:
    def test_str(self):
        assert str(Server("myserver")) == "Server: myserver"

    def test_repr(self):
        assert repr(Server("myserver")) == "Server(name=myserver)"

    def test_eq_by_name(self):
        assert Server("a") == Server("a", NUM_CORES=99)
        assert Server("a") != Server("b")

    def test_hash_by_name(self):
        assert hash(Server("a")) == hash("a")

    def test_from_dict(self):
        server = Server.from_dict({"name": "x", "NUM_CORES": 4})
        assert server.name == "x"
        assert server.num_cores == 4


@pytest.mark.usefixtures("_clear_scheduler_caches")
class TestServerProperties:
    def test_scheduler_from_kwargs(self):
        assert Server("s", SCHEDULER="SLURM").scheduler == "SLURM"

    def test_scheduler_default_none(self):
        assert Server("s").scheduler is None

    def test_queue_name_get_set(self):
        server = Server("s", QUEUE_NAME="normal")
        assert server.queue_name == "normal"
        server.queue_name = "gpu"
        assert server.queue_name == "gpu"

    def test_num_hours_get_set(self):
        server = Server("s", NUM_HOURS=24)
        assert server.num_hours == 24
        server.num_hours = 48
        assert server.num_hours == 48

    def test_mem_gb_default(self):
        assert Server("s").mem_gb == 64

    def test_mem_gb_override(self):
        assert Server("s", MEM_GB=128).mem_gb == 128

    def test_num_cores_default(self):
        assert Server("s").num_cores == 16

    def test_num_gpus_default(self):
        assert Server("s").num_gpus == 0

    def test_num_threads_default(self):
        assert Server("s").num_threads == 16

    def test_scratch_dir_default_none(self):
        assert Server("s").scratch_dir is None

    def test_scratch_true_when_scratch_dir_set(self):
        assert Server("s", SCRATCH_DIR="/tmp/scratch").scratch is True

    def test_scratch_false_when_scratch_dir_unset(self):
        assert Server("s").scratch is False

    def test_use_hosts_default_none(self):
        assert Server("s").use_hosts is None

    def test_extra_commands_default_none(self):
        assert Server("s").extra_commands is None

    def test_extra_scheduler_directives_default_none(self):
        assert Server("s").extra_scheduler_directives is None


@pytest.mark.usefixtures("_clear_scheduler_caches")
class TestSubmitCommand:
    def test_submit_command_explicit_override(self):
        server = Server("s", SUBMIT_COMMAND="my-custom-submit")
        assert server.submit_command == "my-custom-submit"

    @pytest.mark.parametrize(
        "scheduler,expected",
        [
            ("SLURM", "sbatch"),
            ("PBS", "qsub"),
            ("LSF", "bsub < "),
            ("SGE", "qsub"),
            ("HTCondor", "condor_q"),
        ],
    )
    def test_submit_command_from_scheduler(self, scheduler, expected):
        server = Server("s", SCHEDULER=scheduler)
        assert server.submit_command == expected

    def test_submit_command_unknown_scheduler_is_none(self):
        server = Server("s", SCHEDULER="BOGUS")
        assert server.submit_command is None


@pytest.mark.usefixtures("_clear_scheduler_caches")
class TestServerRegisterBug:
    """Documents a real bug: ``RegistryMixin`` shares a single global
    ``_REGISTRY`` list across *every* registrable class hierarchy
    (``Server``, ``Executable``, ``Submitter``, ...), populated with
    classes (not instances) via the metaclass. ``Server.register()``
    does ``if self in Server._REGISTRY`` to check for a prior instance,
    but that membership test invokes ``Server.__eq__`` (``self.name ==
    other.name``) against every element it isn't identical to —
    including unrelated classes like ``Executable`` that have no
    ``.name`` attribute at all. So ``register()`` always crashes for a
    freshly constructed server. See BUGS_FOUND.md."""

    def test_register_crashes_on_unrelated_registry_entries(self):
        server = Server("unique-register-test-server")
        with pytest.raises(AttributeError, match="has no attribute 'name'"):
            server.register()


@pytest.mark.usefixtures("_clear_scheduler_caches")
class TestDetectServerScheduler:
    def test_detects_slurm_via_env_var(self, monkeypatch):
        monkeypatch.setenv("SLURM_JOB_ID", "12345")
        assert Server.detect_server_scheduler() == "SLURM"

    def test_detects_pbs_via_env_var(self, monkeypatch):
        monkeypatch.delenv("SLURM_JOB_ID", raising=False)
        monkeypatch.delenv("SLURM_CLUSTER_NAME", raising=False)
        monkeypatch.setenv("PBS_JOBID", "1.server")
        assert Server.detect_server_scheduler() == "PBS"

    def test_detects_via_command_available(self, monkeypatch):
        for var in [
            "SLURM_JOB_ID",
            "SLURM_CLUSTER_NAME",
            "PBS_JOBID",
            "PBS_QUEUE",
            "LSB_JOBID",
            "LSB_MCPU_HOSTS",
        ]:
            monkeypatch.delenv(var, raising=False)

        def fake_run(command, **kwargs):
            result = MagicMock()
            if command == ["squeue"]:
                result.stdout = b""
                return result
            raise FileNotFoundError()

        with patch(
            "chemsmart.settings.server.subprocess.run", side_effect=fake_run
        ):
            assert Server.detect_server_scheduler() == "SLURM"

    def test_unknown_scheduler_when_nothing_detected(self, monkeypatch):
        for var in [
            "SLURM_JOB_ID",
            "SLURM_CLUSTER_NAME",
            "PBS_JOBID",
            "PBS_QUEUE",
            "LSB_JOBID",
            "LSB_MCPU_HOSTS",
        ]:
            monkeypatch.delenv(var, raising=False)

        with patch(
            "chemsmart.settings.server.subprocess.run",
            side_effect=FileNotFoundError(),
        ):
            assert Server.detect_server_scheduler() == "Unknown Scheduler"


@pytest.mark.usefixtures("_clear_scheduler_caches")
class TestFromScheduerTypeAndCurrent:
    def test_falls_back_to_local_when_unknown(self):
        with (
            patch.object(
                Server,
                "detect_server_scheduler",
                return_value="Unknown Scheduler",
            ),
            patch.object(
                Server, "from_servername", return_value="local-server-sentinel"
            ) as mock_from_servername,
        ):
            result = Server.from_scheduler_type()

        mock_from_servername.assert_called_once_with(servername="local")
        assert result == "local-server-sentinel"

    def test_current_delegates_to_from_scheduler_type(self):
        with patch.object(
            Server, "from_scheduler_type", return_value="sentinel"
        ) as mock_fst:
            assert Server.current() == "sentinel"
        mock_fst.assert_called_once()


@pytest.mark.usefixtures("_clear_scheduler_caches")
class TestFromServername:
    def test_none_uses_current_server(self):
        with patch.object(
            Server, "current", return_value="current-sentinel"
        ) as mock_current:
            result = Server.from_servername(None)
        mock_current.assert_called_once()
        assert result == "current-sentinel"

    def test_missing_server_raises_with_helpful_message(self):
        with pytest.raises(ValueError, match="No server implemented"):
            Server.from_servername("totally_bogus_server_name_xyz")

    def test_existing_server_loads_successfully(self, server_yaml_file):
        assert server_yaml_file.endswith(".yaml")
        bare_name = server_yaml_file[: -len(".yaml")]
        server = Server.from_servername(bare_name)
        assert server.num_cores == 64

    def test_from_servername_accepts_name_without_yaml_suffix(
        self, server_yaml_file
    ):
        # server_yaml_file already ends with .yaml; strip it to exercise
        # the "append .yaml" branch of from_servername.
        assert server_yaml_file.endswith(".yaml")
        bare_name = server_yaml_file[: -len(".yaml")]
        executable = GaussianExecutable.from_servername(bare_name)
        assert executable.executable_folder == os.path.expanduser(
            "~/programs/g16"
        )

    def test_from_servername_accepts_name_with_yaml_suffix(
        self, server_yaml_file
    ):
        executable = GaussianExecutable.from_servername(server_yaml_file)
        assert executable.executable_folder == os.path.expanduser(
            "~/programs/g16"
        )

    def test_from_servername_program_block_without_envars(self, tmp_path):
        yaml_file = tmp_path / "minimal.yaml"
        yaml_file.write_text("GAUSSIAN:\n    EXEFOLDER: ~/programs/g16\n")
        executable = GaussianExecutable.from_servername(str(yaml_file))
        assert executable.envars is None


@pytest.mark.usefixtures("_clear_scheduler_caches")
class TestGetSubmitter:
    def test_get_submitter_dispatches_via_scheduler(self):
        server = Server("s", SCHEDULER="PBS")
        job = MagicMock()
        with patch(
            "chemsmart.settings.server.Submitter.from_scheduler_type"
        ) as mock_from_scheduler:
            mock_from_scheduler.return_value = "submitter-sentinel"
            result = server.get_submitter(job)

        mock_from_scheduler.assert_called_once_with(
            scheduler_type="PBS", job=job, server=server
        )
        assert result == "submitter-sentinel"


@pytest.mark.usefixtures("_clear_scheduler_caches")
class TestServerSubmit:
    def test_submit_writes_script_and_submits_by_default(self):
        server = Server("s", SCHEDULER="PBS")
        job = MagicMock()
        with (
            patch.object(server, "_check_running_jobs") as mock_check,
            patch.object(server, "_write_submission_script") as mock_write,
            patch.object(server, "_submit_job") as mock_submit,
        ):
            server.submit(job, cli_args=["gaussian", "opt"])

        mock_check.assert_called_once_with(job)
        mock_write.assert_called_once_with(
            job=job, cli_args=["gaussian", "opt"]
        )
        mock_submit.assert_called_once_with(job)

    def test_submit_test_mode_skips_actual_submission(self):
        server = Server("s", SCHEDULER="PBS")
        job = MagicMock()
        with (
            patch.object(server, "_check_running_jobs"),
            patch.object(server, "_write_submission_script"),
            patch.object(server, "_submit_job") as mock_submit,
        ):
            server.submit(job, test=True, cli_args=[])

        mock_submit.assert_not_called()

    def test_check_running_jobs_no_op_for_unlabeled_job(self):
        job = MagicMock()
        job.label = None
        with patch("chemsmart.utils.cluster.ClusterHelper") as mock_helper_cls:
            Server._check_running_jobs(job)
        mock_helper_cls.assert_not_called()

    def test_check_running_jobs_exits_on_duplicate(self):
        from chemsmart.jobs.gaussian.job import GaussianGeneralJob

        job = MagicMock(spec=GaussianGeneralJob)
        job.label = "duplicate_job"
        job.execution_name = "gaussian_duplicate_job"

        mock_cluster_helper = MagicMock()
        mock_cluster_helper.get_running_jobs.return_value = (
            ["123"],
            ["gaussian_duplicate_job"],
        )
        with patch(
            "chemsmart.utils.cluster.ClusterHelper",
            return_value=mock_cluster_helper,
        ):
            with pytest.raises(SystemExit, match="Duplicate job NOT"):
                Server._check_running_jobs(job)

    def test_submit_job_raises_without_submit_command(self):
        server = Server("s", SCHEDULER="BOGUS")
        job = MagicMock()
        job.folder = "."
        with patch.object(server, "get_submitter", return_value=MagicMock()):
            with pytest.raises(ValueError, match="no submit command"):
                server._submit_job(job)

    def test_submit_job_uses_split_command_without_shell_operators(self):
        server = Server("s", SCHEDULER="PBS")
        job = MagicMock()
        job.folder = "/some/folder"
        mock_submitter = MagicMock()
        mock_submitter.submit_script = "chemsmart_sub_x.sh"

        with (
            patch.object(server, "get_submitter", return_value=mock_submitter),
            patch("chemsmart.settings.server.subprocess.Popen") as mock_popen,
        ):
            mock_process = MagicMock()
            mock_process.wait.return_value = 0
            mock_popen.return_value = mock_process

            result = server._submit_job(job)

        mock_popen.assert_called_once_with(["qsub", "chemsmart_sub_x.sh"])
        assert result == 0

    def test_submit_job_uses_shell_true_for_shell_operators(self):
        server = Server("s", SCHEDULER="LSF")
        job = MagicMock()
        job.folder = "/some/folder"
        mock_submitter = MagicMock()
        mock_submitter.submit_script = "chemsmart_sub_x.sh"

        with (
            patch.object(server, "get_submitter", return_value=mock_submitter),
            patch("chemsmart.settings.server.subprocess.Popen") as mock_popen,
        ):
            mock_process = MagicMock()
            mock_process.wait.return_value = 0
            mock_popen.return_value = mock_process

            server._submit_job(job)

        assert mock_popen.call_args.kwargs.get("shell") is True


@pytest.mark.usefixtures("_clear_scheduler_caches")
class TestServerSubmitArrayJob:
    def test_no_jobs_warns_and_returns(self, caplog):
        server = Server("s", SCHEDULER="PBS")
        with caplog.at_level("WARNING"):
            server.submit_array_job(jobs=[])
        assert "No jobs to submit" in caplog.text

    def test_submits_array_job_by_default(self):
        server = Server("s", SCHEDULER="PBS")
        jobs = [MagicMock(), MagicMock()]
        mock_submitter = MagicMock()

        with (
            patch.object(server, "_check_running_jobs") as mock_check,
            patch.object(server, "get_submitter", return_value=mock_submitter),
            patch.object(server, "_submit_array_job") as mock_submit_array,
        ):
            server.submit_array_job(jobs=jobs, num_nodes=2, cli_args=[])

        assert mock_check.call_count == 2
        mock_submitter.write_array_job.assert_called_once_with(
            jobs=jobs, num_nodes=2, cli_args=[]
        )
        mock_submit_array.assert_called_once_with(jobs[0], mock_submitter)

    def test_test_mode_skips_actual_submission(self):
        server = Server("s", SCHEDULER="PBS")
        jobs = [MagicMock()]
        mock_submitter = MagicMock()

        with (
            patch.object(server, "_check_running_jobs"),
            patch.object(server, "get_submitter", return_value=mock_submitter),
            patch.object(server, "_submit_array_job") as mock_submit_array,
        ):
            server.submit_array_job(jobs=jobs, test=True)

        mock_submit_array.assert_not_called()

    def test_submit_array_job_raises_without_submit_command(self):
        server = Server("s", SCHEDULER="BOGUS")
        job = MagicMock()
        job.folder = "."
        submitter = MagicMock()
        with pytest.raises(ValueError, match="no submit command"):
            server._submit_array_job(job, submitter)

    def test_submit_array_job_uses_split_command_without_shell_operators(
        self,
    ):
        server = Server("s", SCHEDULER="PBS")
        job = MagicMock()
        job.folder = "/some/folder"
        submitter = MagicMock()
        submitter.array_submit_script = "chemsmart_array_x.sh"

        with patch("chemsmart.settings.server.subprocess.Popen") as mock_popen:
            mock_process = MagicMock()
            mock_process.wait.return_value = 0
            mock_popen.return_value = mock_process

            result = server._submit_array_job(job, submitter)

        mock_popen.assert_called_once_with(["qsub", "chemsmart_array_x.sh"])
        assert result == 0

    def test_submit_array_job_uses_shell_true_for_shell_operators(self):
        server = Server("s", SCHEDULER="LSF")
        job = MagicMock()
        job.folder = "/some/folder"
        submitter = MagicMock()
        submitter.array_submit_script = "chemsmart_array_x.sh"

        with patch("chemsmart.settings.server.subprocess.Popen") as mock_popen:
            mock_process = MagicMock()
            mock_process.wait.return_value = 0
            mock_popen.return_value = mock_process

            server._submit_array_job(job, submitter)

        assert mock_popen.call_args.kwargs.get("shell") is True


@pytest.mark.usefixtures("_clear_scheduler_caches")
class TestCheckRunningJobsNoDuplicate:
    def test_check_running_jobs_returns_normally_when_not_duplicate(self):
        """A job absent from the scheduler queue can be submitted."""
        from chemsmart.jobs.gaussian.job import GaussianGeneralJob

        job = MagicMock(spec=GaussianGeneralJob)
        job.label = "not_running_job"
        job.execution_name = "gaussian_not_running_job"

        mock_cluster_helper = MagicMock()
        mock_cluster_helper.get_running_jobs.return_value = (
            ["123"],
            ["some_other_job"],
        )
        with patch(
            "chemsmart.utils.cluster.ClusterHelper",
            return_value=mock_cluster_helper,
        ):
            assert Server._check_running_jobs(job) is None


@pytest.mark.usefixtures("_clear_scheduler_caches")
class TestWriteSubmissionScript:
    def test_write_submission_script_delegates_to_submitter(self):
        server = Server("s", SCHEDULER="PBS")
        job = MagicMock()
        mock_submitter = MagicMock()

        with patch.object(
            server, "get_submitter", return_value=mock_submitter
        ) as mock_get_submitter:
            server._write_submission_script(job, cli_args=["opt"], extra=1)

        mock_get_submitter.assert_called_once_with(job, extra=1)
        mock_submitter.write.assert_called_once_with(["opt"])


@pytest.mark.usefixtures("_clear_scheduler_caches")
class TestFromYaml:
    def test_empty_name_raises(self):
        with pytest.raises(ValueError, match="No yaml file provided"):
            Server.from_yaml(name=None)


@pytest.mark.usefixtures("_clear_scheduler_caches")
class TestFromServernameYamlSuffix:
    def test_name_already_ending_in_yaml_is_not_double_suffixed(
        self, server_yaml_file
    ):
        """Covers the `if server_name.endswith(".yaml"): server_name =
        server_name` no-op branch by passing a name that already ends
        with ".yaml"."""
        assert server_yaml_file.endswith(".yaml")
        server = Server.from_servername(server_yaml_file)
        assert server.num_cores == 64


@pytest.mark.usefixtures("_clear_scheduler_caches")
class TestFromSchedulerTypeSubclassMatch:
    def test_matching_scheduler_dispatches_to_subclass(self):
        from chemsmart.settings.server import SLURMServer

        with (
            patch.object(
                Server, "detect_server_scheduler", return_value="SLURM"
            ),
            patch.object(
                SLURMServer,
                "from_servername",
                return_value="slurm-sentinel",
            ) as mock_from_servername,
        ):
            result = Server.from_scheduler_type()

        mock_from_servername.assert_called_once_with("SLURM")
        assert result == "slurm-sentinel"

    def test_no_matching_subclass_raises(self):
        with patch.object(
            Server, "detect_server_scheduler", return_value="HTCondor"
        ):
            with pytest.raises(
                ValueError, match="No server class defined for scheduler"
            ):
                Server.from_scheduler_type()


@pytest.mark.usefixtures("_clear_scheduler_caches")
class TestDetectServerSchedulerCheckOutput:
    def test_sge_detected_via_check_output_match(self, monkeypatch):
        for var in [
            "SLURM_JOB_ID",
            "SLURM_CLUSTER_NAME",
            "PBS_JOBID",
            "PBS_QUEUE",
            "LSB_JOBID",
            "LSB_MCPU_HOSTS",
        ]:
            monkeypatch.delenv(var, raising=False)

        def fake_run(command, **kwargs):
            # PBS also tries a bare "qstat" (with no check_output), so
            # that exact command must fail here to let the loop reach
            # SGE's own second command option, "qstat -help".
            if command == ["qstat", "-help"]:
                result = MagicMock()
                result.stdout = b"Sun Grid Engine\n"
                return result
            raise FileNotFoundError()

        with patch(
            "chemsmart.settings.server.subprocess.run", side_effect=fake_run
        ):
            assert Server.detect_server_scheduler() == "SGE"

    def test_sge_check_output_no_match_falls_through(self, monkeypatch):
        for var in [
            "SLURM_JOB_ID",
            "SLURM_CLUSTER_NAME",
            "PBS_JOBID",
            "PBS_QUEUE",
            "LSB_JOBID",
            "LSB_MCPU_HOSTS",
        ]:
            monkeypatch.delenv(var, raising=False)

        def fake_run(command, **kwargs):
            # Only SGE's second command option ("qstat -help") should
            # actually run; every other scheduler's command (including
            # PBS's plain "qstat") must fail so the loop reaches SGE's
            # check_output without detecting an earlier scheduler.
            if command == ["qstat", "-help"]:
                result = MagicMock()
                result.stdout = b"unrelated output\n"
                return result
            raise FileNotFoundError()

        with patch(
            "chemsmart.settings.server.subprocess.run", side_effect=fake_run
        ):
            assert Server.detect_server_scheduler() == "Unknown Scheduler"


@pytest.mark.usefixtures("_clear_scheduler_caches")
class TestRegisterSuccessPathsIsolated:
    """Bug #7 (see TestServerRegisterBug above) means register() always
    crashes when the shared global _REGISTRY contains unrelated
    non-Server entries, as it does in the real app. These tests
    temporarily swap in a clean _REGISTRY containing only Server
    instances to exercise register()'s own intended logic in
    isolation, purely for coverage."""

    def test_register_appends_new_instance(self):
        original_registry = Server._REGISTRY
        Server._REGISTRY = []
        try:
            server = Server("isolated-register-test")
            result = server.register()
            assert result is server
            assert server in Server._REGISTRY
        finally:
            Server._REGISTRY = original_registry

    def test_register_is_idempotent_for_already_registered_instance(self):
        original_registry = Server._REGISTRY
        try:
            server = Server("isolated-register-test-2")
            Server._REGISTRY = [server]
            result = server.register()
            assert result is server
            assert Server._REGISTRY == [server]
        finally:
            Server._REGISTRY = original_registry


class TestYamlServerSettings:
    def test_from_yaml_loads_server_section(self, server_yaml_file):
        settings = YamlServerSettings.from_yaml(server_yaml_file)
        assert settings.name == server_yaml_file
        assert settings.num_cores == 64
        assert settings.mem_gb == 375

    def test_repr_and_str(self, server_yaml_file):
        settings = YamlServerSettings.from_yaml(server_yaml_file)
        assert repr(settings) == f"YamlServerSettings(name={server_yaml_file})"
        assert str(settings) == f"YamlServerSettings: {server_yaml_file}"

    def test_eq_and_hash_by_name(self):
        a = YamlServerSettings(name="x")
        b = YamlServerSettings(name="x", NUM_CORES=99)
        c = YamlServerSettings(name="y")
        assert a == b
        assert a != c
        assert hash(a) == hash("x")

    def test_call_returns_self(self):
        settings = YamlServerSettings(name="x")
        assert settings() is settings

    def test_register_returns_self_without_side_effects(self):
        """Unlike ``Server.register()`` (see BUGS_FOUND.md #7),
        ``YamlServerSettings.register()`` is overridden to be a
        trivial, safe no-op that just returns ``self``."""
        settings = YamlServerSettings(name="x")
        assert settings.register() is settings


class TestServerSettingsManager:
    def test_rejects_none_filename(self):
        with pytest.raises(ValueError, match="filename is not specified"):
            ServerSettingsManager(filename=None)

    def test_stores_absolute_path(self, server_yaml_file):
        manager = ServerSettingsManager(filename=server_yaml_file)
        assert manager.filename == os.path.abspath(server_yaml_file)

    def test_create_loads_yaml_server_settings(self, server_yaml_file):
        manager = ServerSettingsManager(filename=server_yaml_file)
        settings = manager.create()
        assert isinstance(settings, YamlServerSettings)
        assert settings.num_cores == 64

    def test_create_raises_for_missing_file(self, tmp_path):
        manager = ServerSettingsManager(filename=str(tmp_path / "nope.yaml"))
        with pytest.raises(FileNotFoundError):
            manager.create()


class TestSchedulerSpecificServers:
    def test_pbs_server_construction(self):
        server = PBSServer(NUM_CORES=8)
        assert server.name == "PBS"
        assert server.SCHEDULER_TYPE == "PBS"
        assert server.num_cores == 8

    def test_lsf_server_construction(self):
        server = LSFServer(NUM_CORES=4)
        assert server.name == "LSF"
        assert server.SCHEDULER_TYPE == "LSF"

    def test_sge_server_construction(self):
        server = SGE_Server(NUM_CORES=2)
        assert server.name == "SGE"
        assert server.SCHEDULER_TYPE == "SGE"


class TestSLURMServerBug:
    """Documents a real bug: ``SLURMServer.__init__`` calls
    ``super().__init__(filename=..., **kwargs)``, but
    ``YamlServerSettings.__init__`` (and ``Server.__init__``) take a
    positional/``name`` argument, not ``filename`` — unlike
    ``PBSServer``/``LSFServer``/``SGE_Server``, which correctly pass
    ``self.NAME`` positionally. Constructing ``SLURMServer()`` directly
    always raises ``TypeError``. In practice this path isn't hit by
    ``Server.from_scheduler_type()`` (which goes through
    ``from_servername`` → ``ServerSettingsManager`` instead of calling
    the class directly), but the class is still public API and broken.
    See BUGS_FOUND.md."""

    def test_construction_raises_type_error(self):
        with pytest.raises(TypeError, match="missing 1 required positional"):
            SLURMServer(NUM_CORES=8)


class TestExecutableScratchDir:
    def test_scratch_dir_none_without_envars(self):
        exe = Executable(executable_folder="/opt/prog")
        assert exe.scratch_dir is None

    def test_scratch_dir_none_when_not_present_in_envars(self):
        exe = Executable(envars="export FOO=bar\n")
        assert exe.scratch_dir is None

    def test_scratch_dir_extracted_from_envars(self):
        exe = Executable(envars="export SCRATCH=/tmp/scratch\n")
        assert exe.scratch_dir == "/tmp/scratch"

    def test_scratch_dir_strips_inline_comment(self):
        exe = Executable(
            envars="export SCRATCH=/tmp/scratch  # local scratch\n"
        )
        assert exe.scratch_dir == "/tmp/scratch"

    def test_scratch_dir_finds_line_among_others(self):
        exe = Executable(
            envars=(
                "export PATH=/opt/bin:$PATH\n"
                "export SCRATCH=/tmp/scratch\n"
                "export OTHER=1\n"
            )
        )
        assert exe.scratch_dir == "/tmp/scratch"


class TestExecutableEnv:
    def test_env_none_without_envars(self):
        exe = Executable()
        assert exe.env is None

    def test_env_parses_export_lines(self):
        exe = Executable(
            envars=(
                "export PATH=/opt/bin:$PATH\n" "export SCRATCH=/tmp/scratch\n"
            )
        )
        env = exe.env
        assert env == {
            "PATH": "/opt/bin:$PATH",
            "SCRATCH": "/tmp/scratch",
        }

    def test_env_ignores_non_export_lines(self):
        exe = Executable(envars=("# a comment line\n" "export FOO=bar\n" "\n"))
        assert exe.env == {"FOO": "bar"}

    def test_env_strips_inline_comment(self):
        exe = Executable(envars="export FOO=bar  # a comment\n")
        assert exe.env == {"FOO": "bar"}


class TestExecutableAvailableServers:
    def test_available_servers_delegates_to_user_settings(self):
        exe = Executable()
        # Just verify it returns whatever the (real) user settings expose,
        # without asserting specific contents (environment-dependent).
        assert isinstance(exe.available_servers, list)


class TestGaussianExecutable:
    def test_get_executable_returns_none_without_folder(self):
        exe = GaussianExecutable()
        assert exe.get_executable() is None

    def test_get_executable_joins_folder_and_g16(self):
        exe = GaussianExecutable(executable_folder="/opt/g16")
        assert exe.get_executable() == os.path.join("/opt/g16", "g16")

    def test_program_identifier(self):
        assert GaussianExecutable.PROGRAM == "GAUSSIAN"


class TestORCAExecutable:
    def test_get_executable_returns_none_without_folder(self):
        exe = ORCAExecutable()
        assert exe.get_executable() is None

    def test_get_executable_joins_folder_and_orca(self):
        exe = ORCAExecutable(executable_folder="/opt/orca")
        assert exe.get_executable() == os.path.join("/opt/orca", "orca")

    def test_program_identifier(self):
        assert ORCAExecutable.PROGRAM == "ORCA"


class TestNCIPLOTExecutable:
    def test_get_executable_returns_none_without_folder(self):
        exe = NCIPLOTExecutable()
        assert exe.get_executable() is None

    def test_get_executable_joins_folder_and_nciplot(self):
        exe = NCIPLOTExecutable(executable_folder="/opt/nciplot")
        assert exe.get_executable() == os.path.join("/opt/nciplot", "nciplot")

    def test_program_identifier(self):
        assert NCIPLOTExecutable.PROGRAM == "NCIPLOT"


class TestProgramScratchFromServername:
    def test_none_when_program_unset(self):
        assert Executable.program_scratch_from_servername("anything") is None

    def test_none_when_servername_falsy(self):
        assert GaussianExecutable.program_scratch_from_servername("") is None

    def test_appends_yaml_suffix_when_missing(self, server_yaml_file):
        bare_name = server_yaml_file[: -len(".yaml")]
        assert (
            GaussianExecutable.program_scratch_from_servername(bare_name)
            is True
        )

    def test_accepts_direct_yaml_file_path(self, server_yaml_file):
        assert (
            GaussianExecutable.program_scratch_from_servername(
                server_yaml_file
            )
            is True
        )

    def test_none_when_scratch_key_null(self, tmp_path):
        yaml_file = tmp_path / "server.yaml"
        yaml_file.write_text(
            "GAUSSIAN:\n    EXEFOLDER: ~/programs/g16\n    SCRATCH: null\n"
        )
        assert (
            GaussianExecutable.program_scratch_from_servername(str(yaml_file))
            is None
        )

    def test_none_when_scratch_key_missing(self, tmp_path):
        yaml_file = tmp_path / "server.yaml"
        yaml_file.write_text("GAUSSIAN:\n    EXEFOLDER: ~/programs/g16\n")
        assert (
            GaussianExecutable.program_scratch_from_servername(str(yaml_file))
            is None
        )

    def test_none_on_missing_file(self, tmp_path):
        missing = tmp_path / "does_not_exist.yaml"
        assert (
            GaussianExecutable.program_scratch_from_servername(str(missing))
            is None
        )
