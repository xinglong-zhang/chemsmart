"""
CLI subcommand for TFD (Torsion Fingerprint Deviation) based molecular grouping.
"""

import logging

import click

from chemsmart.cli.grouper.grouper import grouper
from chemsmart.utils.cli import MyCommand

logger = logging.getLogger(__name__)


@grouper.command("tfd", cls=MyCommand)
@click.option(
    "--use-weights/--no-use-weights",
    type=bool,
    default=True,
    help="Whether to use torsion weights in TFD calculation. Default=True.",
)
@click.pass_context
def tfd(ctx, use_weights):
    """
    Group structures using Torsion Fingerprint Deviation (TFD).

    TFD compares the torsion angles of rotatable bonds between conformers.
    This method is particularly effective for flexible molecules where
    conformational differences are primarily due to bond rotations.

    Default threshold: 0.1

    Examples:
        chemsmart run grouper -f conformers.xyz tfd
        chemsmart run grouper -f conformers.xyz -T 0.2 tfd --no-use-weights
        chemsmart run grouper -f conformers.xyz -N 10 tfd
    """
    molecules = ctx.obj["molecules"]
    num_procs = ctx.obj["num_procs"]
    label = ctx.obj["grouper_label"]
    num_groups = ctx.obj["num_groups"]

    # Use threshold from parent command, with strategy-specific default
    threshold = ctx.obj["threshold"]
    if threshold is None and num_groups is None:
        threshold = 0.1  # Strategy-specific default

    logger.info(
        f"Running TFD grouping with threshold={threshold}, "
        f"use_weights={use_weights}"
    )

    from chemsmart.jobs.grouper import GrouperJob

    return GrouperJob(
        molecules=molecules,
        grouping_strategy="torsion",
        threshold=threshold,
        num_groups=num_groups,
        ignore_hydrogens=False,
        num_procs=num_procs,
        use_weights=use_weights,
        label=f"{label}_tfd",
    )
