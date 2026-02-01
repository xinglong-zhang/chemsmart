"""
CLI subcommand for spyrmsd-based molecular grouping.
"""

import logging

import click

from chemsmart.cli.grouper.grouper import (
    create_grouper_job_from_context,
    grouper,
)
from chemsmart.utils.cli import MyCommand

logger = logging.getLogger(__name__)


@grouper.command("spyrmsd", cls=MyCommand)
@click.pass_context
def spyrmsd(ctx):
    """
    Group structures using spyrmsd library.

    Uses the spyrmsd library for RMSD calculation with support for
    molecular symmetry.

    Default threshold: 0.5 Å

    Examples:
        chemsmart run grouper -f conformers.xyz spyrmsd
        chemsmart run grouper -f conformers.xyz -T 0.3 spyrmsd
        chemsmart run grouper -f conformers.xyz -N 10 spyrmsd
    """
    logger.info("Running spyrmsd grouping")
    return create_grouper_job_from_context(
        ctx, strategy="spyrmsd", default_threshold=0.5
    )
