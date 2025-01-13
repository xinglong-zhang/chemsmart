import contextlib
import logging
import platform
from multiprocessing import set_start_method
import click
from chemsmart.utils.logger import create_logger
from chemsmart.cli.subcommands import subcommands
from chemsmart.jobs.runner import JobRunner
from chemsmart.cli.jobrunner import jobrunner_options
from chemsmart.cli.logger import logger_options

system_type = platform.system()

if system_type in ("Darwin", "Windows"):
    with contextlib.suppress(RuntimeError):
        set_start_method("fork")

logger = logging.getLogger(__name__)


@click.group(name="run")
@click.pass_context
@jobrunner_options
@logger_options
def run(
    ctx,
    server,
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
    jobrunner = JobRunner(server=server, scratch=scratch, fake=fake)

    # Log the scratch value for debugging purposes
    logger.debug(f"Scratch value passed from CLI: {scratch}")

    # Store the jobrunner and other options in the context object
    ctx.ensure_object(dict)  # Ensure ctx.obj is initialized as a dict
    ctx.obj["mem_gb"] = mem_gb
    ctx.obj["jobrunner"] = jobrunner


@run.result_callback()
@click.pass_context
def process_pipeline(ctx, *args, **kwargs):
    # will give the following error if without **kwargs:
    # TypeError: process_pipeline() got an unexpected keyword argument 'stream'

    # Retrieve the jobrunner from context
    # jobrunner at this stage is an instance of JobRunner class
    jobrunner = ctx.obj["jobrunner"]
    mem_gb = ctx.obj["mem_gb"]

    # Get the job
    job = args[0]
    if isinstance(job, list) and len(job) == 1:
        job = job[0]

    # Instantiate a specific jobrunner based on job type
    # jobrunner at this stage is an instance of specific JobRunner subclass to run the job
    jobrunner = jobrunner.from_jobtype(
        job=job,
        server=jobrunner.server,
        scratch=jobrunner.scratch,  # Propagate scratch
        fake=jobrunner.fake,
        mem_gb=mem_gb,
        **kwargs
    )

    # Run the job with the jobrunner
    job.run(jobrunner=jobrunner)


for subcommand in subcommands:
    run.add_command(subcommand)

if __name__ == "__main__":
    obj = {}
    try:
        run(obj=obj)
    except KeyboardInterrupt as e:
        raise e
