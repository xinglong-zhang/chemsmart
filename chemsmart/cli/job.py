"""
CLI options for ALL jobs that can be run in this package.

This module provides common command-line interface decorators and options
that are shared across different job types in the chemsmart package.
"""

import functools

import click


def click_job_options(f):
    """
    Common job control options for all job types.
    
    Provides standard options for job execution control, including
    the ability to skip completed jobs or force re-execution.
    """
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
    """
    PubChem integration options for molecular structure retrieval.
    
    Provides command-line options for querying molecular structures
    from the PubChem database using various identifiers.
    """

    @click.option(
        "-P",
        "--pubchem",
        type=str,
        default=None,
        help="Queries structure from PubChem using name, smiles, cid and "
        "conformer information.",
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options
