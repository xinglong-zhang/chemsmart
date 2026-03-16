import functools
import logging
import os

import click

from chemsmart.assembler.inspect import DatabaseInspector
from chemsmart.utils.cli import MyCommand

from .database import database

logger = logging.getLogger(__name__)


def click_inspect_options(f):
    """Common click options for database inspect."""

    @click.option(
        "-f",
        "--file",
        type=str,
        required=True,
        help="Path to the input database file (.db).",
    )
    @click.option(
        "-i",
        "--index",
        type=int,
        default=None,
        help="Record index (1-based) to inspect.",
    )
    @click.option(
        "--id",
        "record_id",
        type=str,
        default=None,
        help="Record ID (or prefix, at least 12 chars) to inspect.",
    )
    @click.option(
        "-m",
        "--molecule",
        type=int,
        default=None,
        help="Molecule index (1-based) within the record.",
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options


@database.command(cls=MyCommand)
@click_inspect_options
@click.pass_context
def inspect(ctx, file, index, record_id, molecule):
    """Inspect a chemsmart database, record, or molecule.

    Without -i/--id, show a database overview (metadata and statistics).
    With -i or --id, show detailed information for one record.
    With -i/--id and -m, show detailed information for one molecule.

    \b
    Examples:
        chemsmart run database inspect -f my.db
        chemsmart run database inspect -f my.db -i 3
        chemsmart run database inspect -f my.db --id a1b2c3d4e5f6
        chemsmart run database inspect -f my.db -i 3 -m 1
    """
    # Validate input database
    file = os.path.abspath(file)
    if not os.path.isfile(file):
        raise click.UsageError(f"Database file not found: {file}")

    # Mutual exclusivity: -i and --id
    if index is not None and record_id is not None:
        raise click.UsageError(
            "Options -i/--index and --id are mutually exclusive."
        )

    # -m requires -i or --id
    if molecule is not None and index is None and record_id is None:
        raise click.UsageError(
            "Option -m/--molecule requires -i/--index or --id."
        )

    inspector = DatabaseInspector(
        file, index=index, record_id=record_id, molecule=molecule
    )

    if index is None and record_id is None:
        # Database overview
        logger.info(
            f"Displaying database overview for {os.path.basename(file)}."
        )
        print(inspector.format_overview())
    elif molecule is None:
        # Record detail
        if index is not None:
            logger.info(
                f"Displaying record at index {index} from {os.path.basename(file)}."
            )
        else:
            logger.info(
                f"Displaying record with ID {record_id} from {os.path.basename(file)}."
            )
        print(inspector.format_record_detail())
    else:
        # Molecule detail
        if index is not None:
            logger.info(
                f"Displaying molecule {molecule} from record at index {index} in {os.path.basename(file)}."
            )
        else:
            logger.info(
                f"Displaying molecule {molecule} from record with ID {record_id} in {os.path.basename(file)}."
            )
        print(inspector.format_molecule_detail())

    return None
