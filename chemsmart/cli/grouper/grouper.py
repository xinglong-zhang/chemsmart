"""
CLI command group for molecular structure grouping/clustering.

This module provides the main `grouper` command that serves as a parent
for various grouping strategy subcommands (irmsd, tfd, tanimoto, etc.).
"""

import functools
import logging
import os

import click

from chemsmart.cli.job import (
    click_file_label_and_index_options,
    click_filenames_options,
)
from chemsmart.io.molecules.structure import Molecule
from chemsmart.utils.cli import MyGroup
from chemsmart.utils.io import clean_label
from chemsmart.utils.utils import get_list_from_string_range

logger = logging.getLogger(__name__)


def click_grouper_common_options(f):
    """Common click options for all grouper subcommands."""

    @click.option(
        "-ih",
        "--ignore-hydrogens",
        is_flag=True,
        default=False,
        help="Whether to ignore hydrogens in grouping.",
    )
    @click.option(
        "-np",
        "--num-procs",
        type=int,
        default=1,
        help="Number of processors to use for the grouper.",
    )
    @click.option(
        "-T",
        "--threshold",
        type=float,
        default=None,
        help="Threshold for grouping. If not specified, uses strategy-specific "
        "defaults: RMSD=0.5, Tanimoto=0.9, TFD=0.1, Connectivity=0.0.",
    )
    @click.option(
        "-N",
        "--num-groups",
        type=int,
        default=None,
        help="Target number of groups to return. This uses adaptive threshold "
        "finding to return approximately this many unique structures.",
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options


@click.group(name="grouper", cls=MyGroup)
@click_filenames_options
@click_file_label_and_index_options
@click_grouper_common_options
@click.pass_context
def grouper(
    ctx,
    filenames,
    label,
    append_label,
    index,
    ignore_hydrogens,
    num_procs,
    threshold,
    num_groups,
):
    """
    Group/cluster molecular structures based on various similarity metrics.

    This command provides multiple grouping strategies as subcommands:
    - irmsd: Invariant RMSD (considers molecular symmetry)
    - tfd: Torsion Fingerprint Deviation
    - tanimoto: Tanimoto similarity with various fingerprints
    - rmsd: Simple Kabsch RMSD
    - hrmsd: Hungarian RMSD
    - spyrmsd: spyrmsd-based RMSD
    - pymolrmsd: PyMOL-based RMSD alignment

    Common options (-f, -T, -N, -np, -ih) go BEFORE the subcommand.
    Strategy-specific options go AFTER the subcommand.

    Examples:
        chemsmart run grouper -f conformers.xyz irmsd
        chemsmart run grouper -f conformers.xyz -T 0.5 irmsd --check-stereo on
        chemsmart run grouper -f conformers.xyz -T 0.1 tfd --use-weights
        chemsmart run grouper -f conformers.xyz -N 10 tanimoto --fingerprint-type morgan
    """
    ctx.ensure_object(dict)

    # Store common options in context for subcommands
    ctx.obj["ignore_hydrogens"] = ignore_hydrogens
    ctx.obj["num_procs"] = num_procs
    ctx.obj["threshold"] = threshold
    ctx.obj["num_groups"] = num_groups

    # Initialize molecules
    molecules = None

    if filenames is None:
        logger.warning("[filename] has not been specified!")
        ctx.obj["molecules"] = None
        ctx.obj["grouper_label"] = None
        return

    # Handle single file input (like mol.py)
    if len(filenames) == 1:
        filename = filenames[0]
        # Load all structures from file
        molecules = Molecule.from_filepath(
            filepath=filename, index=":", return_list=True
        )
        assert (
            molecules is not None
        ), f"Could not obtain molecules from {filename}!"
        logger.debug(f"Obtained {len(molecules)} molecules from {filename}")
    else:
        # For now, only support single file
        raise click.BadParameter(
            "Currently only single file input is supported for grouper. "
            "Please provide one file with multiple structures."
        )

    # Update labels
    if label is not None and append_label is not None:
        raise ValueError("Only give label or append_label, but not both!")

    filename = filenames[0]
    if append_label is not None:
        grouper_label = os.path.splitext(os.path.basename(filename))[0]
        grouper_label = f"{grouper_label}_{append_label}"
    elif label is not None:
        grouper_label = label
    else:
        grouper_label = os.path.splitext(os.path.basename(filename))[0]

    grouper_label = clean_label(grouper_label)

    logger.info(
        f"Obtained {len(molecules)} molecules before applying indices, "
        f"with label: {grouper_label}"
    )

    # If user specified index, filter structures
    if index is not None:
        logger.debug(f"Using molecules with index: {index}")
        try:
            from chemsmart.utils.utils import string2index_1based

            idx = string2index_1based(index)
            molecules = molecules[idx]
            if not isinstance(molecules, list):
                molecules = [molecules]
        except ValueError:
            # Handle user defined indices like '[1-3,28-31,34-41]' or '1-3,28-31,34-41'
            idx_list = get_list_from_string_range(index)
            molecules = [molecules[i - 1] for i in idx_list]
    # For grouper, we need all structures by default (no else clause to take last one)

    logger.debug(f"Final molecules count: {len(molecules)}")

    if len(molecules) < 2:
        # Don't raise error here - let subcommand handle it
        # This allows --help to work without files
        logger.warning("Less than 2 molecules loaded, grouping may not work.")

    ctx.obj["molecules"] = molecules
    ctx.obj["grouper_label"] = grouper_label

    logger.info(f"Loaded {len(molecules)} molecules for grouping")
