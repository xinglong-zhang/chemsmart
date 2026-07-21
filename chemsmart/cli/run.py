import logging
import platform
from multiprocessing import set_start_method

import click

from chemsmart.cli.jobrunner import (
    click_jobrunner_options,
    scratch_for_jobrunner,
)
from chemsmart.cli.logger import logger_options
from chemsmart.cli.subcommands import subcommands
from chemsmart.jobs.job import Job
from chemsmart.jobs.runner import JobRunner
from chemsmart.utils.logger import create_logger

logger = logging.getLogger(__name__)


system_type = platform.system()

if system_type == "Darwin":
    try:
        set_start_method("fork")
    except RuntimeError as e:
        logger.error(f"Failed to set start method to 'fork' on Darwin: {e}")
elif system_type == "Windows":
    try:
        set_start_method("spawn")
    except RuntimeError as e:
        logger.error(f"Failed to set start method to 'spawn' on Windows: {e}")
else:
    pass


@click.group(name="run")
@click.pass_context
@click_jobrunner_options
@logger_options
def run(
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
):
    """
    Main command for running chemsmart jobs.

    This command sets up the job runner with specified computational
    resources and executes jobs directly in the current environment.
    """
    # Set up logging
    create_logger(debug=debug, stream=stream)
    logger.info("Entering main program")

    resolved_scratch = scratch_for_jobrunner(ctx, scratch)

    # Instantiate the jobrunner with CLI options
    jobrunner = JobRunner(
        server=server,
        scratch=resolved_scratch,
        delete_scratch=delete_scratch,
        fake=fake,
        num_cores=num_cores,
        num_gpus=num_gpus,
        mem_gb=mem_gb,
    )

    logger.debug(f"Scratch value passed from CLI: {resolved_scratch}")

    # Store the jobrunner and other options in the context object
    ctx.ensure_object(dict)  # Ensure ctx.obj is initialized as a dict
    ctx.obj["jobrunner"] = jobrunner


def _run_single_job(job, jobrunner):
    """Attach a typed jobrunner and execute one job locally."""
    typed_runner = jobrunner.from_job(
        job=job,
        server=jobrunner.server,
        scratch=jobrunner.scratch,
        fake=jobrunner.fake,
        delete_scratch=jobrunner.delete_scratch,
        num_cores=jobrunner.num_cores,
        num_gpus=jobrunner.num_gpus,
        mem_gb=jobrunner.mem_gb,
    )
    job.jobrunner = typed_runner
    job.run()


@run.result_callback()
@click.pass_context
def process_pipeline(ctx, *args, **kwargs):
    """
    Process the job returned by subcommands.
    """
    logger.debug(f"Processing pipeline with args: {args}, kwargs: {kwargs}")
    # will give the following error if without **kwargs:
    # TypeError: process_pipeline() got an unexpected keyword argument
    # 'stream'

    # Retrieve the jobrunner from context
    # jobrunner at this stage is an instance of JobRunner class
    jobrunner = ctx.obj["jobrunner"]

    # Get the job (some command paths intentionally return no job object)
    if not args:
        logger.debug("No pipeline output returned. Skipping job execution.")
        return None
    job = args[0]
    logger.debug(f"Job to be run: {job}")

    # Handle None return (e.g., from post-processing subcommands like
    # boltzmann)
    if job is None:
        logger.debug(
            "No job to process (None returned). Skipping job execution."
        )
        return None

    if isinstance(job, list):
        if not job:
            logger.debug("Empty job list. Skipping job execution.")
            return None
        if not all(isinstance(single_job, Job) for single_job in job):
            raise ValueError(
                "Batch job submission is not supported in this branch."
            )
        logger.info(f"Running {len(job)} jobs locally")
        for single_job in job:
            _run_single_job(single_job, jobrunner)
        return None

    if isinstance(job, Job):
        _run_single_job(job, jobrunner)
    else:
        raise ValueError(f"Invalid job type: {type(job)}.")


for subcommand in subcommands:
    run.add_command(subcommand)


if __name__ == "__main__":
    obj: dict[str, str] = {}
    try:
        run(obj=obj)
    except KeyboardInterrupt as e:
        raise e
