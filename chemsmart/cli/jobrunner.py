import functools
import logging

import click

logger = logging.getLogger(__name__)


def click_jobrunner_options(f):
    """Common job runner configuration options."""

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
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options
