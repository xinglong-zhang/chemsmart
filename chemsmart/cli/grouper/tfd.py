"""
CLI subcommand for TFD (Torsion Fingerprint Deviation) based molecular grouping.
"""

import logging

import click

from chemsmart.cli.grouper.grouper import (
    create_grouper_job_from_context,
    grouper,
)
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
    logger.info(f"Running TFD grouping with use_weights={use_weights}")
    return create_grouper_job_from_context(
        ctx, strategy="torsion", default_threshold=0.1, use_weights=use_weights
    )
