"""CLI options for ALL jobs that can be run in this package."""

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


def click_pubchem_options(f):

    @click.option(
        "-P",
        "--pubchem",
        type=str,
        default=None,
        help="Queries structure from PubChem using name, smiles, cid and conformer information.",
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options
