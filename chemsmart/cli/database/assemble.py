import logging
import os

import click

from chemsmart.cli.job import click_folder_options
from chemsmart.database.assemble import (
    SingleFileAssembler,
    SingleFolderAssembler,
)
from chemsmart.database.database import Database
from chemsmart.io.folder import BaseFolder
from chemsmart.io.xtb.folder import XTBFolder
from chemsmart.utils.cli import MyCommand

from .database import database

logger = logging.getLogger(__name__)


@database.command(cls=MyCommand)
@click_folder_options
@click.option(
    "-i",
    "--index",
    default=":",
    show_default=True,
    help="Index (1-based) of the structures to extract from multi-structure files.",
)
@click.option(
    "-o",
    "--output",
    type=str,
    default="database.db",
    show_default=True,
    help="Output database file (.db extension). ",
)
@click.option(
    "--include-failed",
    is_flag=True,
    default=False,
    help="Include records from calculations that did not terminate normally. "
    "Partial data (e.g. intermediate geometries) will be assembled when available.",
)
@click.pass_context
def assemble(
    ctx,
    directory,
    filetype,
    program,
    index,
    output,
    include_failed,
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

    supported_programs = {"gaussian", "orca", "xtb"}
    if program is None:
        programs = supported_programs
    elif program.lower() in supported_programs:
        programs = {program.lower()}
    else:
        raise ValueError(
            f"Unsupported program '{program}'. "
            "Use 'gaussian', 'orca', or 'xtb'."
        )

    files = []
    file_programs = programs & {"gaussian", "orca"}
    if file_programs:
        base_folder = BaseFolder(folder=directory)
        for prog in sorted(file_programs):
            files.extend(
                base_folder.get_all_output_files_in_current_folder_and_subfolders_by_program(
                    program=prog
                )
            )
        files = sorted(set(files))

    folders = []
    folder_programs = programs & {"xtb"}
    if folder_programs:
        base_folder = BaseFolder(folder=directory)
        for prog in sorted(folder_programs):
            output_files = base_folder.get_all_output_files_in_current_folder_and_subfolders_by_program(
                program=prog
            )
            candidate_folders = {
                os.path.dirname(os.path.abspath(output_file))
                for output_file in output_files
            }
            if prog == "xtb":
                folders.extend(
                    folder
                    for folder in candidate_folders
                    if XTBFolder(folder=folder).is_xtb_calculation_directory
                )
        folders = sorted(set(folders))

    if not files and not folders:
        program_label = program or "gaussian, orca, or xtb"
        logger.error(
            f"No {program_label} output found in directory: {directory}"
        )
        return None

    if program is None:
        logger.info(
            f"Found {len(files)} output file(s) and "
            f"{len(folders)} xTB folder(s), assembling..."
        )
    elif program == "xtb":
        logger.info(f"Found {len(folders)} xTB folder(s), assembling...")
    else:
        logger.info(
            f"Found {len(files)} {program} output files, assembling..."
        )

    # Parse all collected files or folders
    rows = []
    for file in files:
        try:
            assembler = SingleFileAssembler(
                filename=file, index=index, include_failed=include_failed
            )
            data = assembler.assemble_data
            if data:
                rows.append(data)
        except Exception as e:
            logger.error(f"Failed to parse {file}: {e}")

    for folder in folders:
        try:
            assembler = SingleFolderAssembler(
                folder=folder, index=index, include_failed=include_failed
            )
            data = assembler.assemble_data
            if data:
                rows.append(data)
        except Exception as e:
            logger.error(f"Failed to parse {folder}: {e}")

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
            f"Processed {attempted} source(s), but only {actual} unique "
            f"record(s) were stored ({attempted - actual} duplicates were "
            f"replaced)."
        )
    logger.info(f"Assembled {actual} record(s) into database: {output}")
    return None
