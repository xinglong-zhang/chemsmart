import functools
import logging
import os

import click

from chemsmart.cli.job import click_database_id_options
from chemsmart.database.inspect import DatabaseInspector
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
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options


@database.command(cls=MyCommand)
@click_inspect_options
@click_database_id_options
@click.pass_context
def inspect(
    ctx,
    file,
    record_index,
    record_id,
    structure_index,
    molecule_id,
    structure_id,
):
    """Inspect a chemsmart database, record, molecule, or structure.

    Without entity options, show a database overview (metadata and statistics).
    With --ri or --rid, show detailed information for one record.
    With --ri/--rid and --si, show detailed information for one structure within a record.
    With --mid, show detailed information for a molecule (chemical species).
    With --sid, show detailed information for a structure (conformer).

    Options --ri, --rid, --mid, and --sid are mutually exclusive.
    Option --si requires --ri or --rid.

    \b
    Examples:
        chemsmart run database inspect -f my.db
        chemsmart run database inspect -f my.db --ri 3
        chemsmart run database inspect -f my.db --rid a1b2c3d4e5f6
        chemsmart run database inspect -f my.db --ri 3 --si 1
        chemsmart run database inspect -f my.db --mid CURLTUGMZLYLDI-U
        chemsmart run database inspect -f my.db --sid c4d5e6f78a9b
    """
    logger.info(f"Validating database: {file}")
    file = os.path.abspath(file)
    if not os.path.isfile(file):
        raise click.UsageError(f"Database file not found: {file}")

    from chemsmart.database.utils import (
        check_schema_version,
        is_chemsmart_database,
    )

    if not is_chemsmart_database(file):
        raise click.UsageError(
            f"File {file} is not a valid chemsmart database file."
        )
    try:
        check_schema_version(file)
    except RuntimeError as e:
        raise click.UsageError(str(e))

    # Mutual exclusivity: --ri / --rid / --mid / --sid
    entity_options = [
        ("--ri/--record-index", record_index),
        ("--rid/--record-id", record_id),
        ("--mid/--molecule-id", molecule_id),
        ("--sid/--structure-id", structure_id),
    ]
    specified = [
        (name, val) for name, val in entity_options if val is not None
    ]
    if len(specified) > 1:
        names = " and ".join(name for name, _ in specified)
        raise click.UsageError(f"Options {names} are mutually exclusive.")

    # --si requires --ri or --rid
    if (
        structure_index is not None
        and record_index is None
        and record_id is None
    ):
        raise click.UsageError(
            "Option --si/--structure-index requires --ri/--record-index or --rid/--record-id."
        )

    # Convert structure_index to int for DatabaseInspector (shared decorator
    # uses str type to support slice syntax; inspect accepts only single index)
    if structure_index is not None:
        try:
            structure_index = int(structure_index)
        except ValueError:
            raise click.UsageError(
                f"--si/--structure-index must be an integer index (1-based) for inspect, "
                f"but got '{structure_index}'."
            )

    inspector = DatabaseInspector(
        file,
        index=record_index,
        record_id=record_id,
        structure_index=structure_index,
        molecule_id=molecule_id,
        structure_id=structure_id,
    )

    if molecule_id is not None:
        logger.info(
            f"Displaying molecule with ID {molecule_id} from {os.path.basename(file)}."
        )
        print(inspector.format_molecule_detail())
    elif structure_id is not None:
        logger.info(
            f"Displaying structure with ID {structure_id} from {os.path.basename(file)}."
        )
        print(inspector.format_standalone_structure_detail())
    elif record_index is None and record_id is None:
        # Database overview
        logger.info(
            f"Displaying database overview for {os.path.basename(file)}."
        )
        print(inspector.format_overview())
    elif structure_index is None:
        # Record detail
        if record_index is not None:
            logger.info(
                f"Displaying record at index {record_index} from {os.path.basename(file)}."
            )
        else:
            logger.info(
                f"Displaying record with ID {record_id} from {os.path.basename(file)}."
            )
        print(inspector.format_record_detail())
    else:
        # Structure detail
        if record_index is not None:
            logger.info(
                f"Displaying structure {structure_index} from record at index {record_index} in {os.path.basename(file)}."
            )
        else:
            logger.info(
                f"Displaying structure {structure_index} from record with ID {record_id} in {os.path.basename(file)}."
            )
        print(inspector.format_structure_detail())

    return None
