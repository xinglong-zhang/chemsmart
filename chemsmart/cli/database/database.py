import logging

import click

from chemsmart.utils.cli import MyGroup

logger = logging.getLogger(__name__)


@click.group(name="database", cls=MyGroup)
@click.pass_context
def database(ctx):
    """CLI for database operations.

    This command group provides subcommands for managing the chemsmart
    database, including assembling calculation data from output files.
    """
    pass
