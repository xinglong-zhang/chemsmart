import logging

import click

from chemsmart.utils.cli import MyGroup

logger = logging.getLogger(__name__)


@click.group(name="database", cls=MyGroup)
@click.pass_context
def database(ctx):
    """CLI for chemsmart database operations."""
    pass
