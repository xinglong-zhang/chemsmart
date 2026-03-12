import logging
import os

import click

from chemsmart.assembler.assemble import SingleFileAssembler
from chemsmart.assembler.database import Database
from chemsmart.cli.job import click_output_folder_options
from chemsmart.io.folder import BaseFolder
from chemsmart.utils.cli import MyCommand

from .database import database

logger = logging.getLogger(__name__)


@database.command(cls=MyCommand)
@click_output_folder_options
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
@click.pass_context
def assemble(
    ctx,
    directory,
    program,
    index,
    output,
):
    """Assemble calculation output files into a SQLite database.

    This command collects calculation data from output files in the specified
    directory and assembles them into a unified SQLite database (.db).

    Example usage:
    chemsmart run database assemble -d results/ -p gaussian -o database.db
    """

    # Ensure the output filename ends with .db
    if not output.lower().endswith(".db"):
        output = output + ".db"

    if directory is None:
        directory = "./"
    directory = os.path.abspath(directory)
    if not os.path.isdir(directory):
        raise FileNotFoundError(f"Directory does not exist: {directory}")

    if program is None:
        raise ValueError(
            "Program must be specified with -p or --program option"
        )

    # Collect available output files
    if program.lower() in {"gaussian", "orca"}:
        folder = BaseFolder(folder=directory)
        files = folder.get_all_output_files_in_current_folder_and_subfolders_by_program(
            program=program.lower()
        )
    else:
        raise ValueError(
            f"Unsupported program '{program}'. Use 'gaussian' or 'orca'."
        )

    if not files:
        logger.error(
            f"No {program} output files found in directory: {directory}"
        )
        return None

    logger.info(f"Found {len(files)} {program} files, assembling...")

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
    count = db.insert_records(rows, program=program.lower())
    logger.info(f"Assembled {count} record(s) into database: {output}")
    return None
