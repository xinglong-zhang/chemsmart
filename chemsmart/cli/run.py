import logging
import platform
from multiprocessing import set_start_method

import click

from chemsmart.cli.jobrunner import click_jobrunner_options
from chemsmart.cli.logger import logger_options
from chemsmart.cli.subcommands import subcommands
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
    run_in_serial,
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

    # Instantiate the jobrunner with CLI options
    if server is not None:
        server = Server.from_servername(server)
    jobrunner = JobRunner(
        server=server,
        scratch=scratch,
        delete_scratch=delete_scratch,
        fake=fake,
        run_in_serial=run_in_serial,
        num_cores=num_cores,
        num_gpus=num_gpus,
        mem_gb=mem_gb,
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

    # Get the job
    job = args[0]
    logger.debug(f"Job to be run: {job}")

    # Handle None return (e.g., from post-processing subcommands like
    # boltzmann)
    if job is None:
        logger.debug(
            "No job to process (None returned). Skipping job execution."
        )
        return None

    # Handle list of jobs (when multiple molecules are specified with --index)
    if isinstance(job, list):
        logger.info(f"Running {len(job)} jobs")
        # Check if jobs should be run in serial
        if jobrunner.run_in_serial:
            logger.info("Running jobs in serial mode (one after another)")
            for single_job in job:
                logger.info(f"Running job: {single_job.label}")
                # Instantiate a specific jobrunner based on job type
                job_specific_runner = jobrunner.from_job(
                    job=single_job,
                    server=jobrunner.server,
                    scratch=jobrunner.scratch,
                    fake=jobrunner.fake,
                    delete_scratch=jobrunner.delete_scratch,
                    run_in_serial=jobrunner.run_in_serial,
                    num_cores=jobrunner.num_cores,
                    num_gpus=jobrunner.num_gpus,
                    mem_gb=jobrunner.mem_gb,
                )
                # Attach jobrunner to job and run the job with the jobrunner
                single_job.jobrunner = job_specific_runner
                single_job.run()
        else:
            logger.info("Running jobs using default behavior")
            for single_job in job:
                logger.info(f"Running job: {single_job.label}")
                # Instantiate a specific jobrunner based on job type
                job_specific_runner = jobrunner.from_job(
                    job=single_job,
                    server=jobrunner.server,
                    scratch=jobrunner.scratch,
                    fake=jobrunner.fake,
                    delete_scratch=jobrunner.delete_scratch,
                    run_in_serial=jobrunner.run_in_serial,
                    num_cores=jobrunner.num_cores,
                    num_gpus=jobrunner.num_gpus,
                    mem_gb=jobrunner.mem_gb,
                )
                # Attach jobrunner to job and run the job with the jobrunner
                single_job.jobrunner = job_specific_runner
                single_job.run()
        return None

    # Instantiate a specific jobrunner based on job type
    # jobrunner at this stage is an instance of specific JobRunner subclass
    # to run the job
    if isinstance(job, Job):
        jobrunner = jobrunner.from_job(
            job=job,
            server=jobrunner.server,
            scratch=jobrunner.scratch,
            fake=jobrunner.fake,
            delete_scratch=jobrunner.delete_scratch,
            run_in_serial=jobrunner.run_in_serial,
            num_cores=jobrunner.num_cores,
            num_gpus=jobrunner.num_gpus,
            mem_gb=jobrunner.mem_gb,
        )

        # Attach jobrunner to job and run the job with the jobrunner
        job.jobrunner = jobrunner
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
