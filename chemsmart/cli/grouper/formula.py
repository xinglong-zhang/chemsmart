"""
CLI subcommand for molecular formula-based grouping.
"""

import logging

import click

from chemsmart.cli.grouper.grouper import (
    create_grouper_job_from_context,
    grouper,
)
from chemsmart.utils.cli import MyCommand

logger = logging.getLogger(__name__)


@grouper.command("formula", cls=MyCommand)
@click.pass_context
def formula(ctx):
    """
    Group structures by molecular formula.

    Groups molecules that have the same molecular formula together.
    This grouper does not use threshold or num_groups - molecules are grouped
    by identical chemical formula.

    Examples:
        chemsmart run grouper -f molecules.xyz formula
    """
    # Validate: formula grouper does not support threshold or num_groups
    if ctx.obj["threshold"] is not None:
        raise click.UsageError(
            "Formula grouper does not support threshold (-T). "
            "Molecules are grouped by identical chemical formula."
        )
    if ctx.obj["num_groups"] is not None:
        raise click.UsageError(
            "Formula grouper does not support num_groups (-N). "
            "Molecules are grouped by identical chemical formula."
        )

    logger.info("Running formula grouping")
    return create_grouper_job_from_context(ctx, strategy="formula")
