"""
CLI subcommand for Hungarian RMSD-based molecular grouping.
"""

import logging

import click

from chemsmart.cli.grouper.grouper import (
    create_grouper_job_from_context,
    grouper,
)
from chemsmart.utils.cli import MyCommand

logger = logging.getLogger(__name__)


@grouper.command("hrmsd", cls=MyCommand)
@click.pass_context
def hrmsd(ctx):
    """
    Group structures using Hungarian RMSD.

    Hungarian RMSD uses the Hungarian algorithm to find optimal atom-to-atom
    correspondence before computing RMSD. Useful when atom ordering differs
    between structures.

    Default threshold: 0.5 Å

    Examples:
        chemsmart run grouper -f conformers.xyz hrmsd
        chemsmart run grouper -f conformers.xyz -T 0.3 hrmsd
        chemsmart run grouper -f conformers.xyz -N 10 hrmsd
    """
    logger.info("Running Hungarian RMSD grouping")
    return create_grouper_job_from_context(ctx, strategy="hrmsd")
