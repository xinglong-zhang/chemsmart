import contextlib
import logging
import platform
from multiprocessing import set_start_method
import click
from chemsmart.cli.subcommands import subcommands

system_type = platform.system()

if system_type in ("Darwin", "Windows"):
    with contextlib.suppress(RuntimeError):
        set_start_method("fork")

logger = logging.getLogger(__name__)


@click.group(name="run")
@click.pass_context
@click.option(
    "-n",
    "--num-processes",
    type=int,
    default=1,
    help="runs jobs in parallel with specified number of processes",
)
@click.option(
    "-d", "--debug/--no-debug", default=False, help="turns on debug logging"
)
@click.option(
    "--test/--no-test",
    default=False,
    help="If true, fake jobrunners will be used",
)
# @click.option(
#     "--stream/--no-stream",
#     default=None,
#     help="Turns on logging to stdout",
# )
def run(
    ctx,
    num_processes,
    debug,
    test,
):

    logger.info("Entering main program")

    jobrunner = create_jobrunner(
        servername=servername,
        num_nodes=num_nodes,
        num_cpus=num_cpus,
        num_gpus=num_gpus,
        num_gpus_per_node_per_job=num_gpus_per_node_per_job,
        exclude_localhost=exclude_localhost,
        fake=test,
        use_host_queues=use_host_queues,
    )

    ctx.obj["num_processes"] = num_processes
    ctx.obj['jobrunner'] = jobrunner


@run.result_callback()
@click.pass_context
def process_pipeline(ctx, *args, **kwargs):
    invoked_subcommand = ctx.invoked_subcommand
    jobrunner = ctx.obj["jobrunner"]
    queue_manager = ctx.obj["queue_manager"]


for subcommand in subcommands:
    run.add_command(subcommand)

if __name__ == "__main__":
    obj = {}
    try:
        run(obj=obj)
    except KeyboardInterrupt as e:
        raise e
