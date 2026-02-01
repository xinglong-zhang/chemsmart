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
    This grouper does not use threshold or num_groups - molecules are either
    isomorphic (same connectivity) or not.

    Examples:
        chemsmart run grouper -f molecules.xyz connectivity
    """
    # Validate: connectivity grouper does not support threshold or num_groups
    if ctx.obj["threshold"] is not None:
        raise click.UsageError(
            "Connectivity grouper does not support threshold (-T). "
            "Molecules are grouped by graph isomorphism (same/different connectivity)."
        )
    if ctx.obj["num_groups"] is not None:
        raise click.UsageError(
            "Connectivity grouper does not support num_groups (-N). "
            "Molecules are grouped by graph isomorphism (same/different connectivity)."
        )

    logger.info("Running connectivity grouping (graph isomorphism)")

    from chemsmart.jobs.grouper import GrouperJob

    return GrouperJob(
        molecules=ctx.obj["molecules"],
        grouping_strategy="connectivity",
        threshold=None,
        num_groups=None,
        ignore_hydrogens=ctx.obj["ignore_hydrogens"],
        num_procs=ctx.obj["num_procs"],
        label=f"{ctx.obj['grouper_label']}_connectivity",
        conformer_ids=ctx.obj.get("conformer_ids"),
    )
