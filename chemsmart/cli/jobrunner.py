import functools
import logging

import click

logger = logging.getLogger(__name__)

PARALLEL_HELP_SUB = (
    "Max concurrent SLURM array tasks when submitting a batch, and expand "
    "nestable jobs (crest/QRC/dias/traj) into one array task per child. "
    "Default is serial: top-level batches use one array task at a time (%1), "
    "and nestable jobs submit as a single parent with nested serial children. "
    "Pass --run-in-parallel to enable concurrent array tasks / nestable "
    "expansion; cap concurrency with -N (see -N help)."
)

NUM_NODES_HELP_SUB = (
    "With chemsmart sub --run-in-parallel: max concurrent SLURM array tasks "
    "(%M in --array=1-N%M). Ignored under the serial default. Each array "
    "task still uses one node (--nodes=1); -n sets cores per task."
)

NUM_NODES_HELP_RUN = (
    "Reserved for multi-node batch execution; not used for SLURM array "
    "throttling on chemsmart run."
)

PARALLEL_HELP_RUN = (
    "Batch children always run serially with full resources; "
    "use chemsmart sub for cluster concurrency. "
    "--no-run-in-parallel is the default behavior."
)


def click_jobrunner_options(f=None, *, entry_point=None):
    """Common job runner configuration options.

    Pass ``entry_point="run"`` or ``entry_point="sub"`` so
    ``--run-in-parallel/--no-run-in-parallel`` help matches that command.
    """
    if entry_point == "sub":
        parallel_help = PARALLEL_HELP_SUB
        num_nodes_help = NUM_NODES_HELP_SUB
    elif entry_point == "run":
        parallel_help = PARALLEL_HELP_RUN
        num_nodes_help = NUM_NODES_HELP_RUN
    else:
        parallel_help = f"{PARALLEL_HELP_SUB} {PARALLEL_HELP_RUN}"
        num_nodes_help = f"{NUM_NODES_HELP_SUB} {NUM_NODES_HELP_RUN}"

    def decorator(func):
        @click.option(
            "-s",
            "--server",
            type=str,
            default=None,
            help="Server name. If not specified, will try to automatically "
            "determine and use the current server.",
        )
        @click.option(
            "-n",
            "--num-cores",
            type=int,
            help="Number of cores for each job.",
        )
        @click.option(
            "-g",
            "--num-gpus",
            type=int,
            default=None,
            help="Number of GPUs per node. Defaults to number of GPUs on "
            "the specified server if None.",
        )
        @click.option(
            "-m", "--mem-gb", type=int, default=None, help="Memory in GB."
        )
        @click.option(
            "-N",
            "--num-nodes",
            "--number-of-nodes",
            type=int,
            default=None,
            help=num_nodes_help,
        )
        @click.option(
            "--fake/--no-fake",
            default=False,
            type=bool,
            help="If true, fake job runners will be used.",
        )
        @click.option(
            "--scratch/--no-scratch",
            default=None,
            type=bool,
            help="Run in scratch mode or without a scratch folder.",
        )
        @click.option(
            "--delete-scratch/--no-delete-scratch",
            default=False,
            type=bool,
            help="If job was run in scratch, delete the scratch folder after the job "
            "is completed successfully.",
        )
        @click.option(
            "--run-in-parallel/--no-run-in-parallel",
            default=False,
            type=bool,
            help=parallel_help,
        )
        @functools.wraps(func)
        def wrapper_common_options(*args, **kwargs):
            return func(*args, **kwargs)

        return wrapper_common_options

    if f is not None:
        return decorator(f)
    return decorator
