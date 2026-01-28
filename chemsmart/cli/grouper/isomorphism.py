"""
CLI subcommand for RDKit isomorphism-based molecular grouping.
"""

import logging

import click

from chemsmart.cli.grouper.grouper import grouper
from chemsmart.utils.cli import MyCommand

logger = logging.getLogger(__name__)


@grouper.command("isomorphism", cls=MyCommand)
@click.pass_context
def isomorphism(ctx):
    """
    Group structures using RDKit graph isomorphism.

    Groups molecules based on their molecular graph structure using
    RDKit's isomorphism detection. This grouper does not use a threshold.

    Examples:
        chemsmart run grouper -f conformers.xyz isomorphism
    """
    molecules = ctx.obj["molecules"]
    num_procs = ctx.obj["num_procs"]
    label = ctx.obj["grouper_label"]
    num_groups = ctx.obj["num_groups"]
    threshold = ctx.obj["threshold"]

    # Isomorphism grouper does not support threshold
    if threshold is not None:
        raise click.UsageError(
            "Isomorphism grouper does not support threshold (-T). "
            "Molecules are grouped by RDKit molecular hash and isomorphism checks."
        )

    # Isomorphism grouper does not support num_groups
    if num_groups is not None:
        raise click.UsageError(
            "Isomorphism grouper does not support num_groups (-N). "
            "Molecules are grouped by RDKit molecular hash and isomorphism checks."
        )

    logger.info("Running isomorphism grouping")

    from chemsmart.jobs.grouper import GrouperJob

    return GrouperJob(
        molecules=molecules,
        grouping_strategy="isomorphism",
        threshold=None,
        num_groups=None,
        ignore_hydrogens=False,
        num_procs=num_procs,
        label=f"{label}_isomorphism",
    )
