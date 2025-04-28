import contextlib
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

system_type = platform.system()

if system_type in ("Darwin", "Windows"):
    with contextlib.suppress(RuntimeError):
        set_start_method("fork")

logger = logging.getLogger(__name__)


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
    debug,
    stream,
):
    # Set up logging
    create_logger(debug=debug, stream=stream)
    logger.info("Entering main program")

    # Instantiate the jobrunner with CLI options
    if server is not None:
        server = Server.from_servername(server)
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


@run.result_callback()
@click.pass_context
def process_pipeline(ctx, *args, **kwargs):
    """Process the job returned by subcommands."""
    # will give the following error if without **kwargs:
    # TypeError: process_pipeline() got an unexpected keyword argument 'stream'

    # Retrieve the jobrunner from context
    # jobrunner at this stage is an instance of JobRunner class
    jobrunner = ctx.obj["jobrunner"]

    # Get the job
    job = args[0]
    logger.debug(f"Job to be run: {job}")
    if isinstance(job, list) and len(job) == 1:
        job = job[0]

    # Instantiate a specific jobrunner based on job type
    # jobrunner at this stage is an instance of specific JobRunner subclass to run the job
    if isinstance(job, Job):
        jobrunner = jobrunner.from_job(
            job=job,
            server=jobrunner.server,
            scratch=jobrunner.scratch,
            fake=jobrunner.fake,
            num_cores=jobrunner.num_cores,
            num_gpus=jobrunner.num_gpus,
            mem_gb=jobrunner.mem_gb,
        )
        job.jobrunner = jobrunner
        # Run the job with the jobrunner
        print(f"Running job: {job}")
        print(f"Jobrunner: {jobrunner}")
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
