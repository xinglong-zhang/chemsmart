"""
CLI subcommand for basic Kabsch RMSD-based molecular grouping.
"""

import logging

import click

from chemsmart.cli.grouper.grouper import (
    create_grouper_job_from_context,
    grouper,
)
from chemsmart.utils.cli import MyCommand

logger = logging.getLogger(__name__)


@grouper.command("rmsd", cls=MyCommand)
@click.pass_context
def rmsd(ctx):
    """
    Group structures using simple Kabsch RMSD.

    This method uses the Kabsch algorithm to optimally align structures
    and compute RMSD. Assumes atoms are already in corresponding order.

    Default threshold: 0.5 Å

    Examples:
        chemsmart run grouper -f conformers.xyz rmsd
        chemsmart run grouper -f conformers.xyz -T 0.3 rmsd
        chemsmart run grouper -f conformers.xyz -N 10 rmsd
    """
    logger.info("Running Kabsch RMSD grouping")
    return create_grouper_job_from_context(ctx, strategy="rmsd")
