"""
CLI subcommand for Tanimoto similarity based molecular grouping.
"""

import logging

import click

from chemsmart.cli.grouper.grouper import (
    create_grouper_job_from_context,
    grouper,
)
from chemsmart.utils.cli import MyCommand

logger = logging.getLogger(__name__)


@grouper.command("tanimoto", cls=MyCommand)
@click.option(
    "--fingerprint-type",
    "-ft",
    type=click.Choice(
        [
            "rdkit",
            "rdk",
            "morgan",
            "maccs",
            "atompair",
            "torsion",
            "usr",
            "usrcat",
        ],
        case_sensitive=False,
    ),
    default="rdkit",
    help="Fingerprint type to use. Options: rdkit (default), rdk, morgan, maccs, atompair, torsion, usr, usrcat.",
)
@click.pass_context
def tanimoto(ctx, fingerprint_type):
    """
    Group structures using Tanimoto similarity with molecular fingerprints.

    This method computes molecular fingerprints and uses Tanimoto coefficient
    to measure similarity. Different fingerprint types capture different
    aspects of molecular structure.

    Default threshold: 0.9

    Examples:
        chemsmart run grouper -f conformers.xyz tanimoto
        chemsmart run grouper -f conformers.xyz -T 0.85 tanimoto --fingerprint-type morgan
        chemsmart run grouper -f conformers.xyz -N 10 tanimoto
    """
    logger.info(
        f"Running Tanimoto grouping with fingerprint_type={fingerprint_type}"
    )
    return create_grouper_job_from_context(
        ctx,
        strategy="tanimoto",
        default_threshold=0.9,
        fingerprint_type=fingerprint_type,
    )
