"""
CLI options for ALL jobs that can be run in this package.

This module provides common command-line interface decorators and options
that are shared across different job types in the chemsmart package.
"""

import functools

import click
import numpy as np


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


def click_folder_options(f):
    """
    Common click options for Thermochemistry.
    """

    @click.option(
        "-d",
        "--directory",
        default=None,
        help="Directory in which to run specific jobs for all files.",
    )
    @click.option(
        "-t",
        "--filetype",
        default=None,
        help="Type of file to run specific jobs for, if directory "
        "is specified.",
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options
      
    
def click_molecule_vibrational_displacement_options(f):
    """
    CLI options for vibrationally_displaced() method of Molecule object.
    """

    @click.option(
        "-m",
        "--mode-idx",
        type=int,
        default=1,
        help="Mode number of vibrational mode. 1-indexed.\n"
        "Defaults to 1, first vibrational mode.",
    )
    @click.option(
        "-a",
        "--amp",
        type=float,
        default=0.5,
        help="Amplitude for displacement for chosen vibrational mode.\n"
        "Defaults to 0.5 Å.",
    )
    @click.option(
        "-N",
        "--nframes",
        type=int,
        default=None,
        help="Number of frames to create molecules along the vibrational mode.\n"
        "Defaults to None, for which only one molecule will be created.",
    )
    @click.option(
        "-p",
        "--phase",
        type=float,
        default=np.pi / 2,
        help="Phase angle (radians) for displacement.\n"
        "Defaults to pi/2 so that sin(phase) = 1.",
    )
    @click.option(
        "--normalize/--no-normalize",
        type=bool,
        default=False,
        help="If True, scale the (un-weighted) mode so its largest "
        "per-atom displacement is 1.0, making `amp` the max "
        "displacement. \nDefaults to False.",
    )
    @click.option(
        "--return-xyz/--no-return-xyz",
        type=bool,
        default=False,
        help="If True and `nframes` is set, return a multi-frame XYZ string.",
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options
