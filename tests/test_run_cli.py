"""
Direct unit tests for ``chemsmart.cli.run.process_pipeline``, the
result-callback that dispatches whatever a subcommand returns (a
single Job, a list of Jobs, or None) to local execution.
"""

import importlib
from unittest.mock import patch

import click
import pytest

import chemsmart.cli.run as run_module
from chemsmart.cli.run import process_pipeline, run
from chemsmart.jobs.job import Job
from chemsmart.jobs.runner import JobRunner


def _make_ctx(pbs_server):
    ctx = click.Context(run)
    ctx.ensure_object(dict)
    ctx.obj["jobrunner"] = JobRunner(server=pbs_server, fake=True)
    return ctx


class TestProcessPipeline:
    def test_no_args_returns_none(self, pbs_server):
        """Some subcommands (e.g. post-processing) return nothing at all."""
        ctx = _make_ctx(pbs_server)
        assert process_pipeline.__wrapped__(ctx) is None

    def test_none_job_returns_none(self, pbs_server):
        ctx = _make_ctx(pbs_server)
        assert process_pipeline.__wrapped__(ctx, None) is None

    def test_empty_job_list_returns_none(self, pbs_server):
        ctx = _make_ctx(pbs_server)
        assert process_pipeline.__wrapped__(ctx, []) is None

    def test_invalid_single_job_type_raises(self, pbs_server):
        ctx = _make_ctx(pbs_server)
        with pytest.raises(ValueError, match="Invalid job type"):
            process_pipeline.__wrapped__(ctx, "not-a-job")

    def test_rejects_non_job_batch_payload(self, pbs_server):
        """Scheduler-style batch payloads that are not Job lists stay blocked."""
        ctx = _make_ctx(pbs_server)
        with pytest.raises(
            ValueError, match="Batch job submission is not supported"
        ):
            process_pipeline.__wrapped__(ctx, ["not-a-job", "also-not-a-job"])

    def test_single_job_runs_locally(self, pbs_server, mocker):
        ctx = _make_ctx(pbs_server)
        mock_job = mocker.MagicMock(spec=Job)
        mocker.patch.object(
            ctx.obj["jobrunner"], "from_job", return_value=mocker.MagicMock()
        )
        process_pipeline.__wrapped__(ctx, mock_job)
        mock_job.run.assert_called_once()

    def test_job_list_runs_each_locally(self, pbs_server, mocker):
        ctx = _make_ctx(pbs_server)
        mock_jobs = [mocker.MagicMock(spec=Job) for _ in range(2)]
        mocker.patch.object(
            ctx.obj["jobrunner"], "from_job", return_value=mocker.MagicMock()
        )
        process_pipeline.__wrapped__(ctx, mock_jobs)
        for job in mock_jobs:
            job.run.assert_called_once()


class TestProcessPipelineBranches:
    def test_no_args_returned_skips_job_execution(self, run_ctx):
        assert process_pipeline.__wrapped__(run_ctx) is None

    def test_empty_job_list_skips_job_execution(self, run_ctx):
        assert process_pipeline.__wrapped__(run_ctx, []) is None

    def test_none_job_skips_job_execution(self, run_ctx):
        assert process_pipeline.__wrapped__(run_ctx, None) is None


class TestSetStartMethodModuleLevel:
    def test_start_method_setup_covers_all_os_branches(self):
        """Reloading the module with platform.system() patched
        re-executes the module-level Darwin/Windows/else branches.
        Since a start method is already set from the real import,
        each patched reload also exercises the RuntimeError handler."""
        try:
            with patch("platform.system", return_value="Windows"):
                importlib.reload(run_module)
            assert run_module.system_type == "Windows"

            with patch("platform.system", return_value="Linux"):
                importlib.reload(run_module)
            assert run_module.system_type == "Linux"

            with patch("platform.system", return_value="Darwin"):
                importlib.reload(run_module)
            assert run_module.system_type == "Darwin"
        finally:
            # restore the real platform for any subsequent tests/modules
            # that rely on chemsmart.cli.run reflecting the actual OS
            importlib.reload(run_module)
