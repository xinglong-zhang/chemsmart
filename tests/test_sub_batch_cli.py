from types import SimpleNamespace
from unittest.mock import Mock, patch

import click
from click.testing import CliRunner

from chemsmart.cli.run import run
from chemsmart.cli.sub import sub
from chemsmart.jobs.job import Job
from chemsmart.utils.cli import MyCommand


def _invoke_test_jobs(args):
    jobs = [
        SimpleNamespace(label=f"structure_{index}", jobrunner=None)
        for index in range(2)
    ]

    @click.command("test-batch", cls=MyCommand)
    def test_batch():
        return jobs

    sub.add_command(test_batch)
    server = SimpleNamespace(
        kwargs={},
        num_cores=8,
        num_gpus=0,
        mem_gb=24,
        submit=Mock(),
        submit_batch_job=Mock(),
    )
    try:
        with (
            patch("chemsmart.cli.sub.Server.current", return_value=server),
            patch(
                "chemsmart.cli.sub.JobRunner",
                return_value=SimpleNamespace(server=server),
            ),
        ):
            result = CliRunner().invoke(sub, args, catch_exceptions=False)
    finally:
        sub.commands.pop("test-batch", None)
    return result, server, jobs


def test_num_nodes_routes_job_list_to_one_batch_submission():
    result, server, jobs = _invoke_test_jobs(
        ["-N", "2", "--test", "test-batch"]
    )

    assert result.exit_code == 0, result.output
    server.submit_batch_job.assert_called_once()
    assert server.submit_batch_job.call_args.kwargs["jobs"] == jobs
    assert server.submit_batch_job.call_args.kwargs["num_nodes"] == 2
    assert all(
        "--batch-label" in args
        for args in server.submit_batch_job.call_args.kwargs["cli_args"]
    )
    server.submit.assert_not_called()


def test_omitting_num_nodes_keeps_individual_submissions():
    result, server, jobs = _invoke_test_jobs(["--test", "test-batch"])

    assert result.exit_code == 0, result.output
    assert server.submit.call_count == len(jobs)
    server.submit_batch_job.assert_not_called()


def test_batch_run_marks_incomplete_job_as_failed():
    class IncompleteJob(Job):
        def _run(self, **kwargs):
            pass

        def is_complete(self):
            return False

    job = IncompleteJob(molecule=None, label="structure_1", jobrunner=None)

    @click.command("test-batch-run", cls=MyCommand)
    def test_batch_run():
        return job

    run.add_command(test_batch_run)
    try:
        with (
            patch(
                "chemsmart.cli.run.JobRunner",
                return_value=SimpleNamespace(),
            ),
            patch("chemsmart.cli.run._run_single_job"),
        ):
            result = CliRunner().invoke(
                run,
                ["--batch-label", "structure_1", "test-batch-run"],
            )
    finally:
        run.commands.pop("test-batch-run", None)

    assert result.exit_code != 0
    assert "did not complete" in str(result.exception)
