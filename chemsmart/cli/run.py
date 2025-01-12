import contextlib
import logging
import platform
from multiprocessing import set_start_method
import click
from chemsmart.utils.utils import create_logger
from chemsmart.cli.subcommands import subcommands
from chemsmart.jobs.runner import JobRunner

system_type = platform.system()

if system_type in ("Darwin", "Windows"):
    with contextlib.suppress(RuntimeError):
        set_start_method("fork")

logger = logging.getLogger(__name__)


@click.group(name="run")
@click.pass_context
@click.option(
    "-s",
    "--server",
    type=str,
    default=None,
    help="Server. If not specified, will try to automatically determine and use the current server.",
)
@click.option(
    "-n",
    "--num-processes",
    type=int,
    default=1,
    help="Runs jobs in parallel with specified number of processes.",
)
@click.option(
    "-d", "--debug/--no-debug", default=False, help="Turns on debug logging."
)
@click.option(
    "--fake/--no-fake",
    default=False,
    help="If true, fake jobrunners will be used.",
)
@click.option(
    "--scratch/--no-scratch",
    default=True,  # Default behavior is to use scratch
    help="Run in scratch mode or without scratch folder.",
)
@click.option(
    "--stream/--no-stream",
    default=None,
    help="Turns on logging to stdout.",
)
def run(
    ctx,
    server,
    num_processes,
    debug,
    fake,
    scratch,
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
    ctx.obj["num_processes"] = num_processes
    ctx.obj["jobrunner"] = jobrunner


@run.result_callback()
@click.pass_context
def process_pipeline(ctx, *args, **kwargs):
    # Retrieve the jobrunner from context
    jobrunner = ctx.obj["jobrunner"]

    # Log for debugging
    logger.debug(f"JobRunner before from_jobtype: {jobrunner}")

    # Get the job
    job = args[0]
    if isinstance(job, list) and len(job) == 1:
        job = job[0]

    # Instantiate a specific jobrunner based on job type
    jobrunner = jobrunner.from_jobtype(
        job=job,
        server=jobrunner.server,
        scratch=jobrunner.scratch,  # Propagate scratch
        fake=jobrunner.fake
    )

    # Log for debugging
    logger.debug(f"JobRunner after from_jobtype: {jobrunner}")

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
