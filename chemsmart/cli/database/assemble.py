import functools
import logging
import os

import click

from chemsmart.assembler.assemble import SingleFileAssembler
from chemsmart.assembler.database import Database
from chemsmart.io.folder import BaseFolder
from chemsmart.utils.cli import MyCommand

from .database import database

logger = logging.getLogger(__name__)


def click_assemble_options(f):
    """
    Common click options for database assemble.
    """

    @click.option(
        "-d",
        "--directory",
        default=".",
        show_default=True,
        help="Directory containing calculation output files.",
    )
    @click.option(
        "-t",
        "--filetype",
        type=click.Choice(["gaussian", "orca"], case_sensitive=False),
        required=True,
        help="Type of calculation output files to assemble.",
    )
    @click.option(
        "-i",
        "--index",
        default=":",
        show_default=True,
        help="Index (1-based) of the molecule to extract from multi-molecule files.",
    )
    @click.option(
        "-o",
        "--output",
        type=str,
        default="database.db",
        show_default=True,
        help="Output database file (.db extension). ",
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options


@database.command(cls=MyCommand)
@click_assemble_options
@click.pass_context
def assemble(
    ctx,
    directory,
    filetype,
    index,
    output,
):
    """Assemble calculation output files into a SQLite database.

    This command collects calculation data from output files in the specified
    directory and assembles them into a unified SQLite database (.db).

    Example usage:
    chemsmart run database assemble -d results/ -t gaussian -o database.db
    """

    # Ensure the output filename ends with .db
    if not output.endswith(".db"):
        output = output + ".db"

    directory = os.path.abspath(directory)
    if not os.path.isdir(directory):
        raise FileNotFoundError(f"Directory does not exist: {directory}")

    # Collect available output files
    if filetype.lower() in {"gaussian", "orca"}:
        folder = BaseFolder(folder=directory)
        files = folder.get_all_output_files_in_current_folder_and_subfolders_by_program(
            program=filetype.lower()
        )
    else:
        raise ValueError(
            f"Unsupported filetype '{filetype}'. Use 'gaussian' or 'orca'."
        )

    if not files:
        logger.error(
            f"No {filetype} output files found in directory: {directory}"
        )
        return None

    logger.info(f"Found {len(files)} {filetype} files, assembling...")

    # Parse all collected files
    rows = []
    for file in files:
        try:
            assembler = SingleFileAssembler(filename=file, index=index)
            data = assembler.assemble_data
            if data:
                rows.append(data)
        except Exception as e:
            logger.error(f"Failed to parse {file}: {e}")

    if not rows:
        logger.error("No valid data parsed. Aborting.")
        return None

    # Write to SQLite database
    db = Database(db_file=output)
    db.create()
    count = db.insert_records(rows, program=filetype.lower())
    logger.info(f"Assembled {count} record(s) into database: {output}")
    return None
