"""
Submission of jobs to queuing system via cli.

This module provides command-line interface for submitting jobs to
various queuing systems and cluster schedulers.
"""

import logging

import click

from chemsmart.cli.jobrunner import click_jobrunner_options
from chemsmart.cli.logger import logger_options
from chemsmart.cli.subcommands import subcommands
from chemsmart.jobs.runner import JobRunner
from chemsmart.settings.server import Server
from chemsmart.utils.cli import CtxObjArguments, MyGroup
from chemsmart.utils.logger import create_logger

logger = logging.getLogger(__name__)


@click.group(name="sub", cls=MyGroup)
@click.pass_context
@click_jobrunner_options
@logger_options
@click.option("-t", "--time-hours", type=float, default=None)
@click.option("-q", "--queue", type=str, help="queue")
@click.option(
    "-v/",
    "--verbose/--no-verbose",
    default=False,
    help="Turns on logging to stream output and debug logging.",
)
@click.option(
    "--test/--no-test",
    default=False,
    help="If true, job will not be submitted; only run and submit "
    "scripts will be written.",
)
@click.option(
    "--print-command/--no-print-command",
    default=False,
    help="print the command generated",
)
def sub(
    ctx,
    server,
    num_cores,
    num_gpus,
    mem_gb,
    fake,
    scratch,
    delete_scratch,
    debug,
    stream,
    time_hours,
    queue,
    verbose,
    test,
    print_command,
    **kwargs,
):
    """
    Main command for submitting chemsmart jobs to queuing systems.

    This command prepares and submits jobs to cluster schedulers with
    specified resource requirements and queue parameters.
    """
    # Set up logging
    if verbose:
        create_logger(stream=True, debug=True)
    else:
        create_logger(debug=debug, stream=stream)
    logger.info("Entering main program")

    # Instantiate the jobrunner with CLI options
    if server is not None:
        server = Server.from_servername(server)
        if time_hours is not None:
            server.num_hours = time_hours
        if queue is not None:
            server.queue_name = queue

    jobrunner = JobRunner(
        server=server,
        scratch=scratch,
        delete_scratch=delete_scratch,
        fake=fake,
        num_cores=num_cores,
        num_gpus=num_gpus,
        mem_gb=mem_gb,
    )

    # Log the scratch value for debugging purposes
    logger.debug(f"Scratch value passed from CLI: {scratch}")

    # Store the jobrunner and other options in the context object
    ctx.ensure_object(dict)  # Ensure ctx.obj is initialized as a dict
    ctx.obj["jobrunner"] = jobrunner


@sub.result_callback(replace=True)
@click.pass_context
def process_pipeline(ctx, *args, **kwargs):  # noqa: PLR0915
    """
    Process the job for submission to queuing system.

    This callback function handles job submission by reconstructing
    command-line arguments and interfacing with the appropriate
    scheduler system.
    """

    def _clean_command(ctx):
        """
        Remove keywords used in sub.py but not in run.py.

        Specifically: Some keywords/options (like queue, verbose, etc.)
        are only relevant to sub.py and not applicable to run.py.
        """
        # Get "sub" command and assert that there is exactly one.
        command = next(
            (
                subcommand
                for subcommand in ctx.obj["subcommand"]
                if subcommand["name"] == "sub"
            ),
            None,
        )
        if not command:
            raise ValueError("No 'sub' command found in context.")

        # Find the keywords that are valid in sub.py
        # but should not be passed to run.py and remove those
        keywords_not_in_run = [
            "time_hours",
            "queue",
            "verbose",
            "test",
            "print_command",
        ]

        for keyword in keywords_not_in_run:
            # Remove keyword if it exists
            command["kwargs"].pop(keyword, None)
        return ctx

    def _reconstruct_cli_args(ctx, job):
        """
        Get cli args that reconstruct the command line.

        Rebuilds the command-line arguments from the context object
        for job submission purposes.
        """
        commands = ctx.obj["subcommand"]

        args = CtxObjArguments(commands, entry_point="sub")
        cli_args = args.reconstruct_command_line()[
            1:
        ]  # remove the first element 'sub'
        if kwargs.get("print_command"):
            print(cli_args)
        return cli_args

    def _process_single_job(job):
        if kwargs.get("test"):
            logger.warning('Not submitting as "test" flag specified.')

        cli_args = _reconstruct_cli_args(ctx, job)

        server = Server.from_servername(kwargs.get("server"))
        server.submit(job=job, test=kwargs.get("test"), cli_args=cli_args)

    ctx = _clean_command(ctx)
    jobrunner = ctx.obj["jobrunner"]
    job = args[0]
    job.jobrunner = jobrunner

    _process_single_job(job=job)


for subcommand in subcommands:
    sub.add_command(subcommand)
