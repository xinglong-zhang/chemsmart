"""
CLI subcommand for PyMOL RMSD-based molecular grouping.
"""

import logging

import click

from chemsmart.cli.grouper.grouper import grouper
from chemsmart.utils.cli import MyCommand

logger = logging.getLogger(__name__)


@grouper.command("pymolrmsd", cls=MyCommand)
@click.pass_context
def pymolrmsd(ctx):
    """
    Group structures using PyMOL-based RMSD alignment.

    Uses PyMOL's align command for structure alignment and RMSD calculation.
    Requires PyMOL to be installed.

    Default threshold: 0.5 Å

    Examples:
        chemsmart run grouper -f conformers.xyz pymolrmsd
        chemsmart run grouper -f conformers.xyz -T 0.3 pymolrmsd
        chemsmart run grouper -f conformers.xyz -N 10 pymolrmsd
    """
    molecules = ctx.obj["molecules"]
    ignore_hydrogens = ctx.obj["ignore_hydrogens"]
    num_procs = ctx.obj["num_procs"]
    label = ctx.obj["grouper_label"]
    num_groups = ctx.obj["num_groups"]

    # Use threshold from parent command, with strategy-specific default
    threshold = ctx.obj["threshold"]
    if threshold is None and num_groups is None:
        threshold = 0.5  # Strategy-specific default

    logger.info(
        f"Running PyMOL RMSD grouping with threshold={threshold}, "
        f"ignore_hydrogens={ignore_hydrogens}"
    )

    from chemsmart.jobs.grouper import GrouperJob

    return GrouperJob(
        molecules=molecules,
        grouping_strategy="pymolrmsd",
        threshold=threshold,
        num_groups=num_groups,
        ignore_hydrogens=ignore_hydrogens,
        num_procs=num_procs,
        label=f"{label}_pymolrmsd",
    )
