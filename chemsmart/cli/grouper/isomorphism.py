"""
CLI subcommand for RDKit isomorphism-based molecular grouping.
"""

import logging

import click

from chemsmart.cli.grouper.grouper import (
    create_grouper_job_from_context,
    grouper,
)
from chemsmart.utils.cli import MyCommand

logger = logging.getLogger(__name__)


@grouper.command("isomorphism", cls=MyCommand)
@click.pass_context
def isomorphism(ctx):
    """
    Group structures using RDKit graph isomorphism.

    Groups molecules based on their molecular graph structure using
    RDKit's isomorphism detection. This grouper does not use threshold
    or num_groups - molecules are grouped by RDKit molecular hash.

    Examples:
        chemsmart run grouper -f conformers.xyz isomorphism
    """
    # Validate: isomorphism grouper does not support threshold or num_groups
    if ctx.obj["threshold"] is not None:
        raise click.UsageError(
            "Isomorphism grouper does not support threshold (-T). "
            "Molecules are grouped by RDKit molecular hash and isomorphism checks."
        )
    if ctx.obj["num_groups"] is not None:
        raise click.UsageError(
            "Isomorphism grouper does not support num_groups (-N). "
            "Molecules are grouped by RDKit molecular hash and isomorphism checks."
        )

    logger.info("Running isomorphism grouping")
    return create_grouper_job_from_context(ctx, strategy="isomorphism")
