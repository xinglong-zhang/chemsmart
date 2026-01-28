"""
CLI subcommand for molecular connectivity-based grouping.
"""

import logging

import click

from chemsmart.cli.grouper.grouper import grouper
from chemsmart.utils.cli import MyCommand

logger = logging.getLogger(__name__)


@grouper.command("connectivity", cls=MyCommand)
@click.pass_context
def connectivity(ctx):
    """
    Group structures by molecular connectivity.

    Groups molecules based on their bond connectivity patterns (graph isomorphism).
    This grouper does not use a threshold - molecules are either isomorphic or not.

    Examples:
        chemsmart run grouper -f molecules.xyz connectivity
    """
    molecules = ctx.obj["molecules"]
    num_procs = ctx.obj["num_procs"]
    label = ctx.obj["grouper_label"]
    num_groups = ctx.obj["num_groups"]
    threshold = ctx.obj["threshold"]

    # Connectivity grouper does not support threshold
    if threshold is not None:
        raise click.UsageError(
            "Connectivity grouper does not support threshold (-T). "
            "Molecules are grouped by graph isomorphism (same/different connectivity)."
        )

    # Connectivity grouper does not support num_groups
    if num_groups is not None:
        raise click.UsageError(
            "Connectivity grouper does not support num_groups (-N). "
            "Molecules are grouped by graph isomorphism (same/different connectivity)."
        )

    logger.info("Running connectivity grouping (graph isomorphism)")

    from chemsmart.jobs.grouper import GrouperJob

    return GrouperJob(
        molecules=molecules,
        grouping_strategy="connectivity",
        threshold=None,
        num_groups=None,
        ignore_hydrogens=False,
        num_procs=num_procs,
        label=f"{label}_connectivity",
    )
