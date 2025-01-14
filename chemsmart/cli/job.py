"""CLI options for all jobs that can be run in this package."""

import functools

import click


def click_job_options(f):
    @click.option(
        "-S/-R",
        "--skip-completed/--no-skip-completed",
        is_flag=True,
        default=True,
        type=bool,
        help="To run completed job again. Use -R to rerun completed job.",
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options


def click_gaussian_link_jobs_options(f):
    """Common click options for Gaussian link/crest jobs."""

    @click.option(
        "-j",
        "--jobtype",
        type=str,
        default=None,
        help='Gaussian job type. Options: ["opt", "ts", "modred", "scan", "sp"]',
    )
    @click.option(
        "-c",
        "--coordinates",
        default=None,
        help="list of coordinates to be fixed for modred or scan job",
    )
    @click.option(
        "-s",
        "--step-size",
        default=None,
        help="step size of coordinates to scan",
    )
    @click.option(
        "-n",
        "--num-steps",
        default=None,
        help="step size of coordinates to scan",
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options
