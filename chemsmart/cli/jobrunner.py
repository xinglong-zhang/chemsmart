import functools
import logging

import click

logger = logging.getLogger(__name__)

PARALLEL_HELP_SUB = (
    "Max concurrent SLURM array tasks when submitting a batch. "
    "Use --no-run-in-parallel to run one array task at a time."
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
    elif entry_point == "run":
        parallel_help = PARALLEL_HELP_RUN
    else:
        parallel_help = f"{PARALLEL_HELP_SUB} {PARALLEL_HELP_RUN}"

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
            help="Number of nodes to request for each job.",
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
            default=True,
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
