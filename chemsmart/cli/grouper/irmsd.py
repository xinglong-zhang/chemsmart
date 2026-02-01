"""
CLI subcommand for iRMSD-based molecular grouping.
"""

import logging

import click

from chemsmart.cli.grouper.grouper import (
    create_grouper_job_from_context,
    grouper,
)
from chemsmart.utils.cli import MyCommand

logger = logging.getLogger(__name__)


@grouper.command("irmsd", cls=MyCommand)
@click.option(
    "--check-stereo",
    type=click.Choice(["auto", "on", "off"], case_sensitive=False),
    default="auto",
    help="Control stereochemistry/inversion checking: "
    "'auto' (default): automatically detect, 'on': force check, 'off': disable.",
)
@click.pass_context
def irmsd(ctx, check_stereo):
    """
    Group structures using invariant RMSD (iRMSD).

    iRMSD considers molecular symmetry and equivalent atom permutations
    when calculating RMSD values. This is particularly useful for molecules
    with symmetric groups (e.g., methyl groups, phenyl rings).

    Default threshold: 0.5 Å

    Examples:
        chemsmart run grouper -f conformers.xyz irmsd
        chemsmart run grouper -f conformers.xyz -T 0.3 irmsd --check-stereo on
        chemsmart run grouper -f conformers.xyz -N 10 irmsd
    """
    logger.info(f"Running iRMSD grouping with check_stereo={check_stereo}")
    return create_grouper_job_from_context(
        ctx, strategy="irmsd", default_threshold=0.5, check_stereo=check_stereo
    )
