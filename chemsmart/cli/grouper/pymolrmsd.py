"""
CLI subcommand for PyMOL RMSD-based molecular grouping.
"""

import logging

import click

from chemsmart.cli.grouper.grouper import (
    create_grouper_job_from_context,
    grouper,
)
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
    logger.info("Running PyMOL RMSD grouping")
    return create_grouper_job_from_context(
        ctx, strategy="pymolrmsd", default_threshold=0.5
    )
