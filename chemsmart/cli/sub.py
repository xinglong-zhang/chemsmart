"""Submission of jobs to queuing system via cli."""

import logging

import click

from chemsmart.cli.jobrunner import jobrunner_options
from chemsmart.cli.logger import logger_options
from chemsmart.cli.subcommands import subcommands
from chemsmart.utils.logger import create_logger
from chemsmart.utils.cli import MyGroup
from chemsmart.settings.server import Server
from chemsmart.jobs.runner import JobRunner

# from chemsmart.utils.cli import CtxObjArguments, MyGroup, determine_folders

logger = logging.getLogger(__name__)


# def prepare_job_from_source_folder(folder):
#     from pyatoms.analysis.results.results import NebOptimizationResults, Results
#
#     try:
#         results = Results.from_folder(folder)
#     except Exception as e:
#         raise ValueError(f'Could not create job from {folder}') from e
#
#     if not results.is_complete:
#         raise ValueError(f'Job from {folder} is not complete.')
#
#     atoms = results.transition_state_atoms if isinstance(results, NebOptimizationResults) else results.optimized_atoms
#
#     return atoms.set_calculated_magmoms_as_initial(error_if_no_magmoms=False)


@click.group(name="sub", cls=MyGroup)
@click.pass_context
@jobrunner_options
@logger_options
@click.option("-t", "--time-hours", type=float, default=None)
@click.option("-q", "--queue", type=str, help="queue")
@click.option(
    "-v",
    "--verbose/--no-verbose",
    help="Turns on logging to stream output and debug logging.",
)
@click.option(
    "--test/--no-test",
    default=False,
    help="If true, job will not be submitted; only run and submit scripts will be written.",
)
def sub(
    ctx,
    server,
    num_cores,
    num_gpus,
    mem_gb,
    fake,
    scratch,
    debug,
    stream,
    time_hours,
    queue,
    verbose,
    test,
    **kwargs,
):
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
    def _clean_command(ctx):
        """Remove keywords used in sub.py but not in run.py.
        Specifically: Some keywords/options (like queue, verbose, etc.)
        are only relevant to sub.py and not applicable to run.py."""
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
        ]

        for keyword in keywords_not_in_run:
            # Remove keyword if it exists
            command["kwargs"].pop(keyword, None)
        return ctx

    def _process_single_job(job, jobrunner):
        if kwargs.get("test"):
            logger.warning('Not submitting as "test" flag specified.')

        server = Server.from_servername(kwargs.get("server"))
        server.submit(job=job, test=kwargs.get("test"))

    ctx = _clean_command(ctx)
    jobrunner = ctx.obj["jobrunner"]
    job = args[0]

    _process_single_job(job=job, jobrunner=jobrunner)


for subcommand in subcommands:
    sub.add_command(subcommand)
