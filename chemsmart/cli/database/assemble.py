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
    """Assemble calculation output files into a chemsmart database.

    This command collects calculation data from output files in the specified
    directory and assembles them into a unified chemsmart database (.db).

    \b
    Examples:
        chemsmart run database assemble -d results/ -o my.db
        chemsmart run database assemble -d ./ -p gaussian -o database.db
    """

    # Ensure the output filename ends with .db
    if not output.lower().endswith(".db"):
        output = output + ".db"

    if directory is None:
        directory = "./"
    directory = os.path.abspath(directory)
    if not os.path.isdir(directory):
        raise FileNotFoundError(f"Directory does not exist: {directory}")

    # Collect available output files
    supported_programs = {"gaussian", "orca"}
    if program is None:
        programs = supported_programs
    elif program.lower() in supported_programs:
        programs = {program.lower()}
    else:
        raise ValueError(
            f"Unsupported program '{program}'. Use 'gaussian' or 'orca'."
        )
    folder = BaseFolder(folder=directory)
    files = []
    for prog in programs:
        files.extend(
            folder.get_all_output_files_in_current_folder_and_subfolders_by_program(
                program=prog
            )
        )

    if not files:
        logger.error(
            f"No {', '.join(programs)} output files found in directory: {directory}"
        )
        return None

    if program is None:
        logger.info(f"Found {len(files)} output files, assembling...")
    else:
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

    # Write to chemsmart database
    db = Database(db_file=output)
    db.create()
    attempted = db.insert_records(rows)
    actual = db.count_records()
    if attempted != actual:
        logger.warning(
            f"Processed {attempted} file(s), but only {actual} unique record(s) "
            f"were stored ({attempted - actual} duplicates were replaced)."
        )
    logger.info(f"Assembled {actual} record(s) into database: {output}")
    return None
