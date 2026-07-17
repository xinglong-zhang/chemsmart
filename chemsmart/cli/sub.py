"""
Submission of jobs to queuing system via cli.

This module provides command-line interface for submitting jobs to
various queuing systems and cluster schedulers.
"""

import logging
from pathlib import Path

import click

from chemsmart.cli.jobrunner import click_jobrunner_options
from chemsmart.cli.logger import logger_options
from chemsmart.cli.pka import rewrite_pka_batch_cli_args
from chemsmart.cli.subcommands import subcommands
from chemsmart.jobs.batch import BatchJob
from chemsmart.jobs.batch_manifest import (
    build_manifest_children,
    get_job_batch_entry,
    resolve_array_cli_args,
    write_batch_manifest,
)
from chemsmart.jobs.runner import JobRunner, get_configured_max_submitters
from chemsmart.settings.server import Server
from chemsmart.utils.cli import CtxObjArguments, MyGroup
from chemsmart.utils.logger import create_logger

logger = logging.getLogger(__name__)


@click.group(name="sub", cls=MyGroup)
@click.pass_context
@click_jobrunner_options
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
    num_nodes,
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
        if num_nodes is not None:
            server.num_nodes = num_nodes
        if queue is not None:
            server.queue_name = queue

    jobrunner = JobRunner(
        server=server,
        scratch=scratch,
        delete_scratch=delete_scratch,
        fake=fake,
        no_run_in_parallel=not kwargs.get("run_in_parallel", True),
        num_cores=num_cores,
        num_gpus=num_gpus,
        mem_gb=mem_gb,
        num_nodes=num_nodes,
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

    def _array_throttle(jobrunner, num_jobs):
        """Return SLURM ``--array=1-N%M`` concurrency throttle ``M``.

        ``--no-run-in-parallel`` forces ``M=1``. Otherwise prefer ``-N`` /
        ``num_nodes``, then ``CHEMSMART_MAX_SUBMITTERS`` / cores policy,
        then ``num_jobs``.
        """
        if jobrunner.no_run_in_parallel:
            return 1
        num_nodes = jobrunner.num_nodes
        if num_nodes is not None and num_nodes > 0:
            return int(num_nodes)
        return min(num_jobs, get_configured_max_submitters(jobrunner))

    def _process_single_job(job):
        if kwargs.get("test"):
            logger.warning('Not submitting as "test" flag specified.')

        cli_args = _reconstruct_cli_args(ctx)

        server = Server.from_servername(kwargs.get("server"))
        server.submit(job=job, test=kwargs.get("test"), cli_args=cli_args)

    def _process_batch_job(batch_job):
        """Submit a top-level BatchJob as a scheduler array (one task/child)."""
        if kwargs.get("test"):
            logger.warning('Not submitting as "test" flag specified.')

        shared_cli_args = _reconstruct_cli_args(ctx)
        # Homogeneous batches keep a shared CLI; BatchJob.run() selects the
        # child via array env. Heterogeneous pKa rows use per-task CLI lists
        # (one submit command per child) from job.batch_entry.
        has_batch_entries = any(
            get_job_batch_entry(job) is not None for job in batch_job.jobs
        )
        # Currently only pKa attaches batch_entry; use its CLI rewriter.
        rewrite_cli = rewrite_pka_batch_cli_args if has_batch_entries else None
        array_cli_args = resolve_array_cli_args(
            batch_job.jobs, shared_cli_args, rewrite_cli=rewrite_cli
        )
        if has_batch_entries:
            program = batch_job.PROGRAM or (
                batch_job.jobs[0].PROGRAM if batch_job.jobs else None
            )
            first_folder = getattr(batch_job.jobs[0], "folder", None)
            try:
                manifest_dir = (
                    Path(first_folder) if first_folder else Path(".")
                )
            except TypeError:
                manifest_dir = Path(".")
            write_batch_manifest(
                batch_label=batch_job.label,
                program=str(program).lower() if program else "unknown",
                children=build_manifest_children(
                    batch_job.jobs,
                    shared_cli_args,
                    rewrite_cli=rewrite_cli,
                ),
                directory=manifest_dir,
            )
        throttle = _array_throttle(jobrunner, len(batch_job.jobs))
        logger.info(
            "Submitting BatchJob %r as array with %s task(s), "
            "concurrency throttle %%s=%s",
            batch_job.label,
            len(batch_job.jobs),
            throttle,
        )
        server = Server.from_servername(kwargs.get("server"))
        server.submit_array_job(
            jobs=batch_job.jobs,
            num_nodes=throttle,
            test=kwargs.get("test"),
            cli_args=array_cli_args,
            batch_label=batch_job.label,
        )

    ctx = _clean_command(ctx)
    jobrunner = ctx.obj["jobrunner"]
    job = args[0]

    # Handle list of jobs (when multiple molecules are specified with --index)
    if isinstance(job, list):
        logger.info(f"Processing {len(job)} jobs")
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
        _process_single_job(job=job)


for subcommand in subcommands:
    sub.add_command(subcommand)
