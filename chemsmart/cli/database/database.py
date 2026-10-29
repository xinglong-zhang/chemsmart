import functools
import logging

import click

from chemsmart.utils.cli import MyGroup

logger = logging.getLogger(__name__)


@click.group(name="database", cls=MyGroup)
@click.pass_context
def database(ctx):
    """CLI for chemsmart database operations."""
    pass


def click_database_id_options(f):
    """CLI options for selecting records/molecules/structures from a chemsmart database
    by ID or index. Only relevant when the input file (-f) is a chemsmart database (.db).
    """

    @click.option(
        "--ri",
        "--record-index",
        "record_index",
        type=int,
        default=None,
        help="Record index (1-based) within a chemsmart database.",
    )
    @click.option(
        "--rid",
        "--record-id",
        "record_id",
        type=str,
        default=None,
        help="Record ID (or unique prefix) within a chemsmart database.",
    )
    @click.option(
        "--si",
        "--structure-index",
        "structure_index",
        type=str,
        default=None,
        help="Structure index (1-based) within the selected record.",
    )
    @click.option(
        "--sid",
        "--structure-id",
        "structure_id",
        type=str,
        default=None,
        help="Global structure ID (or unique prefix) within a chemsmart database.",
    )
    @click.option(
        "--mid",
        "--molecule-id",
        "molecule_id",
        type=str,
        default=None,
        help="Molecule ID (or unique prefix); selects every conformer of that molecule.",
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options
