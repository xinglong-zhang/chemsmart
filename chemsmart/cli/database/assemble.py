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
from chemsmart.jobs.job import Job
from chemsmart.utils.cli import MyCommand

from .database import database

logger = logging.getLogger(__name__)


class DatabaseAssemblyJob(Job):
    """Submission metadata for a database assembly command."""

    PROGRAM = "database"

    def __init__(self, output, jobrunner):
        super().__init__(
            molecule=None,
            label=os.path.splitext(os.path.basename(output))[0],
            jobrunner=jobrunner,
            skip_completed=False,
        )

    def is_complete(self):
        return False


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

    # A queued invocation must not discover or parse files on the submit host.
    if ctx.find_root().invoked_subcommand == "sub":
        output = os.path.abspath(output)
        for command in ctx.obj["subcommand"]:
            if command["name"] == "assemble":
                command["kwargs"]["directory"]["value"] = directory
                command["kwargs"]["output"]["value"] = output
        return DatabaseAssemblyJob(output, ctx.obj["jobrunner"])

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
                for folder in sorted(candidate_folders):
                    try:
                        is_calculation_directory = XTBFolder(
                            folder=folder
                        ).is_xtb_calculation_directory
                    except ValueError as e:
                        logger.error(
                            f"Skipping invalid xTB calculation directory "
                            f"'{folder}': {e}"
                        )
                        continue
                    if is_calculation_directory:
                        folders.append(folder)
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

    # Commit bounded batches so long jobs make progress and retain completed work.
    db = Database(db_file=output)
    db.create()
    batch = []
    inserted = 0
    processed = 0
    total = len(files) + len(folders)
    for target, assembler_class, keyword in [
        (file, SingleFileAssembler, "filename") for file in files
    ] + [(folder, SingleFolderAssembler, "folder") for folder in folders]:
        try:
            assembler = assembler_class(
                **{keyword: target}, index=index, include_failed=include_failed
            )
            data = assembler.assemble_data
            if data:
                batch.append(data)
        except Exception:
            logger.exception("Failed to parse %s", target)
        processed += 1
        if len(batch) >= 100:
            inserted += db.insert_records(batch)
            batch.clear()
        if processed % 100 == 0 or processed == total:
            logger.info(
                "Processed %d/%d sources; wrote %d records",
                processed,
                total,
                inserted,
            )

    if batch:
        inserted += db.insert_records(batch)
    if not inserted:
        logger.error("No valid data inserted.")
        return None
    logger.info("Inserted %d records into database: %s", inserted, output)
    return None
