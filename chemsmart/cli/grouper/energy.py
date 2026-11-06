"""
CLI subcommand for energy-based molecular grouping.
"""

import logging

import click

from chemsmart.cli.grouper.grouper import (
    create_grouper_job_from_context,
    grouper,
)
from chemsmart.utils.cli import MyCommand

logger = logging.getLogger(__name__)


@grouper.command("energy", cls=MyCommand)
@click.pass_context
def energy(ctx):
    """
    Group structures using energy differences.

    This method groups molecules based on their energy differences.
    Molecules with energy differences below the threshold are grouped
    together using complete linkage clustering.

    Default threshold: 1.0 kcal/mol

    Note: All molecules must have energy information for this method.

    Examples:
        chemsmart run grouper -f conformers.xyz energy
        chemsmart run grouper -f conformers.xyz -T 0.5 energy
        chemsmart run grouper -f conformers.xyz -N 5 energy
    """
    logger.info("Running Energy-based grouping")
    return create_grouper_job_from_context(ctx, strategy="energy")
