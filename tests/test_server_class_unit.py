"""
Direct unit tests for :class:`chemsmart.settings.server.Server`, beyond
what ``test_server.py`` already covers (YAML loading, executable
resolution). Focuses on: property defaults/overrides, equality/hash,
the ``current``/``from_scheduler_type``/``detect_server_scheduler``
scheduler-autodetection chain, ``from_servername`` error paths, and the
``submit``/``submit_array_job`` orchestration (with subprocess and
cluster-checking collaborators mocked out).

``detect_server_scheduler`` and ``from_scheduler_type`` are
``lru_cache``-wrapped, so every test that exercises them clears the
cache first/after to avoid cross-test pollution.
"""

from unittest.mock import MagicMock, patch

import pytest

from chemsmart.settings.server import Server


@pytest.fixture(autouse=True)
def _clear_scheduler_caches():
    Server.detect_server_scheduler.cache_clear()
    Server.from_scheduler_type.cache_clear()
    yield
    Server.detect_server_scheduler.cache_clear()
    Server.from_scheduler_type.cache_clear()


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

    def test_check_running_jobs_no_op_for_non_gaussian_job(self):
        job = MagicMock()
        job.__class__.__name__ = "NotAGaussianJob"
        # A plain MagicMock is not an instance of GaussianJob, so this
        # should return without touching ClusterHelper at all.
        with patch("chemsmart.utils.cluster.ClusterHelper") as mock_helper_cls:
            Server._check_running_jobs(job)
        mock_helper_cls.assert_not_called()

    def test_check_running_jobs_exits_on_duplicate(self):
        from chemsmart.jobs.gaussian.job import GaussianGeneralJob

        job = MagicMock(spec=GaussianGeneralJob)
        job.label = "duplicate_job"

        mock_cluster_helper = MagicMock()
        mock_cluster_helper.get_gaussian_running_jobs.return_value = (
            ["123"],
            ["duplicate_job"],
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

        mock_popen.assert_called_once_with(
            ["qsub", "chemsmart_sub_x.sh"], cwd="/some/folder"
        )
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
