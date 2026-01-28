"""
CLI subcommand for iRMSD-based molecular grouping.
"""

import logging

import click

from chemsmart.cli.grouper.grouper import grouper
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
        f"Running iRMSD grouping with threshold={threshold}, "
        f"check_stereo={check_stereo}, ignore_hydrogens={ignore_hydrogens}"
    )

    from chemsmart.jobs.grouper import GrouperJob

    return GrouperJob(
        molecules=molecules,
        grouping_strategy="irmsd",
        threshold=threshold,
        num_groups=num_groups,
        ignore_hydrogens=ignore_hydrogens,
        num_procs=num_procs,
        check_stereo=check_stereo,
        label=f"{label}_irmsd",
    )
