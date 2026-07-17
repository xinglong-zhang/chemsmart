import logging
import platform
from multiprocessing import set_start_method

import click

from chemsmart.cli.jobrunner import click_jobrunner_options
from chemsmart.cli.logger import logger_options
from chemsmart.cli.subcommands import subcommands
from chemsmart.jobs.batch import (
    BatchExecutionError,
    BatchJob,
    warn_legacy_job_list,
)
from chemsmart.jobs.job import Job
from chemsmart.jobs.runner import JobRunner
from chemsmart.settings.server import Server
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
@click_jobrunner_options(entry_point="run")
@logger_options
def run(
    ctx,
    server,
    num_cores,
    num_gpus,
    mem_gb,
    num_nodes,
    fake,
    scratch,
    delete_scratch,
    run_in_parallel,
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

    # Local execution does not use scratch unless explicitly requested.
    if scratch is None:
        scratch = False

    # Instantiate the jobrunner with CLI options
    if server is not None:
        server = Server.from_servername(server)
    jobrunner = JobRunner(
        server=server,
        scratch=scratch,
        delete_scratch=delete_scratch,
        fake=fake,
        no_run_in_parallel=not run_in_parallel,
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

    def _prepare_runner(single_job):
        return jobrunner.from_job(
            job=single_job,
            server=jobrunner.server,
            scratch=jobrunner.scratch,
            fake=jobrunner.fake,
            delete_scratch=jobrunner.delete_scratch,
            no_run_in_parallel=jobrunner.no_run_in_parallel,
            num_cores=jobrunner.num_cores,
            num_gpus=jobrunner.num_gpus,
            mem_gb=jobrunner.mem_gb,
            num_nodes=jobrunner.num_nodes,
        )

    # Handle list of jobs
    if isinstance(job, list):
        if not job:
            logger.debug("Empty job list. Skipping job execution.")
            return None
        if not all(isinstance(single_job, Job) for single_job in job):
            raise ValueError("Expected a list of Job instances.")
        warn_legacy_job_list(stacklevel=2)
        logger.info(
            "Running %s jobs serially (legacy list path; prefer BatchJob)",
            len(job),
        )

        failures = []
        for single_job in job:
            logger.info(f"Running job: {single_job.label}")
            single_job.jobrunner = _prepare_runner(single_job)
            try:
                single_job.run()
            except Exception as exc:
                logger.error(
                    f"Job {single_job.label} failed in list execution: {exc}",
                    exc_info=True,
                )
                failures.append({"label": single_job.label, "error": str(exc)})
        if failures:
            lines = [
                f"- {item['label']}: {item['error']}" for item in failures
            ]
            raise BatchExecutionError(
                f"{len(failures)} of {len(job)} list job(s) failed:\n"
                + "\n".join(lines)
            )
        return None

    # Bind an engine runner from the first child so BatchJob can propagate
    # copies onto each child job.
    if isinstance(job, BatchJob):
        if not job.jobs:
            raise ValueError(f"BatchJob {job} has no child jobs to run.")
        job.jobrunner = _prepare_runner(job.jobs[0])
        job.enable_serial_local_execution()
        job.run()
        return None

    # Instantiate a specific jobrunner based on job type
    # jobrunner at this stage is an instance of specific JobRunner subclass
    # to run the job
    if isinstance(job, Job):
        # Attach jobrunner to job and run the job with the jobrunner
        job.jobrunner = _prepare_runner(job)
        job.run()
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
