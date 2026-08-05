"""Reusable Click options for Iterate input commands."""

import functools

import click


def click_iterate_combination_options(f):
    """Common combination-generation options for Iterate inputs."""

    @click.option(
        "-ms",
        "--max-substituted-sites",
        default=None,
        type=click.IntRange(min=0),
        help="Maximum number of substituted skeleton attachment sites per "
        "combination. Unset or 0 means unlimited.",
    )
    @functools.wraps(f)
    def wrapper_iterate_combination_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_iterate_combination_options
