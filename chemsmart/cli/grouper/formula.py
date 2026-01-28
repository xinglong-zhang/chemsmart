"""
CLI subcommand for molecular formula-based grouping.
"""

import logging

import click

from chemsmart.cli.grouper.grouper import grouper
from chemsmart.utils.cli import MyCommand

logger = logging.getLogger(__name__)


@grouper.command("formula", cls=MyCommand)
@click.pass_context
def formula(ctx):
    """
    Group structures by molecular formula.

    Groups molecules that have the same molecular formula together.
    This grouper does not use a threshold - molecules are grouped by identical formula.

    Examples:
        chemsmart run grouper -f molecules.xyz formula
    """
    molecules = ctx.obj["molecules"]
    num_procs = ctx.obj["num_procs"]
    label = ctx.obj["grouper_label"]
    num_groups = ctx.obj["num_groups"]
    threshold = ctx.obj["threshold"]

    # Formula grouper does not support threshold
    if threshold is not None:
        raise click.UsageError(
            "Formula grouper does not support threshold (-T). "
            "Molecules are grouped by identical chemical formula."
        )

    # Formula grouper does not support num_groups
    if num_groups is not None:
        raise click.UsageError(
            "Formula grouper does not support num_groups (-N). "
            "Molecules are grouped by identical chemical formula."
        )

    logger.info("Running formula grouping")

    from chemsmart.jobs.grouper import GrouperJob

    return GrouperJob(
        molecules=molecules,
        grouping_strategy="formula",
        threshold=None,
        num_groups=None,
        ignore_hydrogens=False,
        num_procs=num_procs,
        label=f"{label}_formula",
    )
