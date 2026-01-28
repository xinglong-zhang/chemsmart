"""
CLI subcommand for Tanimoto similarity based molecular grouping.
"""

import logging

import click

from chemsmart.cli.grouper.grouper import grouper
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
    molecules = ctx.obj["molecules"]
    num_procs = ctx.obj["num_procs"]
    label = ctx.obj["grouper_label"]
    num_groups = ctx.obj["num_groups"]

    # Use threshold from parent command, with strategy-specific default
    threshold = ctx.obj["threshold"]
    if threshold is None and num_groups is None:
        threshold = 0.9  # Strategy-specific default

    logger.info(
        f"Running Tanimoto grouping with threshold={threshold}, "
        f"fingerprint_type={fingerprint_type}"
    )

    from chemsmart.jobs.grouper import GrouperJob

    return GrouperJob(
        molecules=molecules,
        grouping_strategy="tanimoto",
        threshold=threshold,
        num_groups=num_groups,
        ignore_hydrogens=False,
        num_procs=num_procs,
        fingerprint_type=fingerprint_type,
        label=f"{label}_tanimoto",
    )
