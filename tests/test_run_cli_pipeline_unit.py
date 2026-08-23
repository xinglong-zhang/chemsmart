"""
Direct unit tests for chemsmart/cli/run.py's process_pipeline result
callback and its module-level, OS-conditional multiprocessing
start-method setup -- neither of which is exercised by the existing
CLI-level tests, which only ever invoke real job-returning
subcommands.
"""

import importlib
from unittest.mock import patch

import click
import pytest

import chemsmart.cli.run as run_module
from chemsmart.cli.run import process_pipeline, run
from chemsmart.jobs.runner import JobRunner
from chemsmart.settings.server import Server


@pytest.fixture()
def run_ctx():
    server = Server(name="dummy")
    ctx = click.Context(run)
    ctx.ensure_object(dict)
    ctx.obj["jobrunner"] = JobRunner(server=server, fake=True)
    return ctx


class TestProcessPipelineBranches:
    def test_no_args_returned_skips_job_execution(self, run_ctx):
        assert process_pipeline.__wrapped__(run_ctx) is None

    def test_empty_job_list_skips_job_execution(self, run_ctx):
        assert process_pipeline.__wrapped__(run_ctx, []) is None

    def test_none_job_skips_job_execution(self, run_ctx):
        assert process_pipeline.__wrapped__(run_ctx, None) is None

    def test_invalid_job_type_raises_value_error(self, run_ctx):
        with pytest.raises(ValueError, match="Invalid job type"):
            process_pipeline.__wrapped__(run_ctx, "not-a-job")


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
