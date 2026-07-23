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
from chemsmart.jobs.batch import (
    BatchJob,
    get_nestable_array_children,
    prepare_nestable_batch_jobs,
    warn_legacy_job_list,
)
from chemsmart.jobs.gaussian.batch import GaussianBatchJob
from chemsmart.jobs.job import Job
from chemsmart.jobs.orca.batch import ORCABatchJob
from chemsmart.jobs.runner import JobRunner
from chemsmart.settings.server import SchedulerArrayPolicy, Server
from chemsmart.utils.cli import CtxObjArguments, MyGroup
from chemsmart.utils.logger import create_logger

logger = logging.getLogger(__name__)


@click.group(name="sub", cls=MyGroup)
@click.pass_context
@click_jobrunner_options(entry_point="sub")
@logger_options
@click.option(
    "-t",
    "--time-hours",
    type=float,
    default=None,
    help="Time limit in hours for the job (e.g., 48.0).",
)
@click.option("-q", "--queue", type=str, help="Queue name for job submission.")
@click.option(
    "-v/",
    "--verbose/--no-verbose",
    default=False,
    help="Turn on logging to stream output and debug logging.",
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
    help="Print the generated command.",
)
def sub(
    ctx,
    server,
    num_cores,
    num_gpus,
    mem_gb,
    array_concurrency,
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
        no_run_in_parallel=not kwargs.get("run_in_parallel", False),
        num_cores=num_cores,
        num_gpus=num_gpus,
        mem_gb=mem_gb,
        array_concurrency=array_concurrency,
    )

    # Log the scratch value for debugging purposes
    logger.debug(f"Scratch value passed from CLI: {scratch}")

    # Store the jobrunner and other options in the context object
    ctx.ensure_object(dict)  # Ensure ctx.obj is initialized as a dict
    ctx.obj["jobrunner"] = jobrunner


@sub.result_callback(replace=True)
@click.pass_context
def process_pipeline(ctx, *args, **kwargs):
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

    def _reconstruct_cli_args(ctx):
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

        cli_args = _reconstruct_cli_args(ctx)

        server = Server.from_servername(kwargs.get("server"))
        server.submit(job=job, test=kwargs.get("test"), cli_args=cli_args)

    def _process_nestable_array_job(parent_job):
        """Expand nestable children and submit as a scheduler array.

        Each array task re-runs the parent CLI; the parent selects one child
        via ``--child-index``, falling back to ``SLURM_ARRAY_TASK_ID`` /
        ``PBS_ARRAYID`` / ``LSB_JOBINDEX``.
        """
        if kwargs.get("test"):
            logger.warning('Not submitting as "test" flag specified.')

        children = get_nestable_array_children(parent_job)
        if not children:
            raise ValueError(
                f"Nestable job {parent_job} has no children to array-submit."
            )

        program = (parent_job.PROGRAM or "").lower()
        batch_cls = ORCABatchJob if program == "orca" else GaussianBatchJob

        rewrite_cli = prepare_nestable_batch_jobs(children)
        batch_job = batch_cls(
            jobs=children,
            label=f"{parent_job.label}_array",
            jobrunner=jobrunner,
            rewrite_cli=rewrite_cli,
        )
        shared_cli_args = _reconstruct_cli_args(ctx)
        logger.info(
            "Expanding nestable job %r into array of %s child task(s) "
            "(--run-in-parallel)",
            parent_job.label,
            len(children),
        )
        server = Server.from_servername(kwargs.get("server"))
        server.submit_batch(
            batch_job,
            policy=SchedulerArrayPolicy.from_jobrunner(jobrunner),
            test=kwargs.get("test"),
            cli_args=shared_cli_args,
        )

    def _process_batch_job(batch_job):
        """Submit a top-level BatchJob via ``Server.submit_batch``."""
        if kwargs.get("test"):
            logger.warning('Not submitting as "test" flag specified.')

        shared_cli_args = _reconstruct_cli_args(ctx)
        server = Server.from_servername(kwargs.get("server"))
        server.submit_batch(
            batch_job,
            policy=SchedulerArrayPolicy.from_jobrunner(jobrunner),
            test=kwargs.get("test"),
            cli_args=shared_cli_args,
        )

    ctx = _clean_command(ctx)
    jobrunner = ctx.obj["jobrunner"]
    job = args[0]

    # Handle list of jobs (legacy path; prefer BatchJob from CLI fan-out)
    if isinstance(job, list):
        if not job:
            logger.debug("Empty job list. Skipping job submission.")
            return
        if not all(isinstance(single_job, Job) for single_job in job):
            raise ValueError("Expected a list of Job instances.")
        warn_legacy_job_list(stacklevel=2)
        logger.info(
            "Processing %s jobs individually (legacy list path; prefer BatchJob)",
            len(job),
        )
        for single_job in job:
            single_job.jobrunner = jobrunner
            _process_single_job(job=single_job)
    elif isinstance(job, BatchJob):
        if not job.jobs:
            raise ValueError(f"BatchJob {job} has no child jobs to submit.")
        job.jobrunner = jobrunner
        _process_batch_job(job)
    else:
        job.jobrunner = jobrunner
        nestable_children = get_nestable_array_children(job)
        # Default / --no-run-in-parallel: one parent job with nested serial
        # children. --run-in-parallel: expand nestable parents to an array.
        if nestable_children is not None and not jobrunner.no_run_in_parallel:
            _process_nestable_array_job(job)
        else:
            _process_single_job(job=job)


for subcommand in subcommands:
    sub.add_command(subcommand)
