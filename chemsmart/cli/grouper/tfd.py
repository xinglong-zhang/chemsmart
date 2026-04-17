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
@click.option(
    "--max-dev",
    type=click.Choice(["equal", "spec"], case_sensitive=False),
    default="equal",
    help="Normalization method for torsion deviations. "
    "'equal': all torsions normalized using 180.0 degrees (default). "
    "'spec': each torsion normalized using its specific maximal deviation.",
)
@click.pass_context
def tfd(ctx, use_weights, max_dev):
    """
    Group structures using Torsion Fingerprint Deviation (TFD).

    TFD compares the torsion angles of rotatable bonds between conformers.
    This method is particularly effective for flexible molecules where
    conformational differences are primarily due to bond rotations.

    Default threshold: 0.1

    Examples:
        chemsmart run grouper -f conformers.xyz tfd
        chemsmart run grouper -f conformers.xyz -T 0.2 tfd --no-use-weights
        chemsmart run grouper -f conformers.xyz -N 10 tfd --max-dev spec
    """
    logger.info(
        f"Running TFD grouping with use_weights={use_weights}, max_dev={max_dev}"
    )
    return create_grouper_job_from_context(
        ctx, strategy="torsion", use_weights=use_weights, max_dev=max_dev
    )
