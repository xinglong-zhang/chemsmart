import functools
import logging
import os
from pathlib import Path

import click

from chemsmart.assembler.export import DataExporter
from chemsmart.assembler.single import SingleFileAssembler
from chemsmart.utils.cli import MyGroup

logger = logging.getLogger(__name__)


def click_assemble_options(f):
    """
    Common click options for DataAssembler.
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
        "-o",
        "--output",
        type=str,
        required=True,
        help="Output filename; extension determines format (.json / .csv / .sqlite). "
        "If no extension is provided, defaults to JSON format.",
    )
    @click.option(
        "-i",
        "--index",
        default="-1",
        show_default=True,
        help="Index (1-based) of the molecule to extract from multi-molecule files.",
    )
    @click.option(
        "-k",
        "--keys",
        type=str,
        show_default=True,
        help="Comma-separated list of export fields.",
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options


# use MyGroup to allow potential subcommands in the future
@click.group(cls=MyGroup, invoke_without_command=True)
@click_assemble_options
@click.pass_context
def assemble(
    ctx,
    directory,
    filetype,
    output,
    index,
    keys,
    **kwargs,
):
    """CLI for running assemble jobs using the chemsmart framework.

    This command collects calculation data from output files in the specified directory and assembles them into a
    unified dataset, which can then be exported in the desired format.

    Example usage:
    chemsmart run assemble -d results/ -t gaussian -o output.json
    """

    directory = os.path.abspath(directory)
    if not os.path.isdir(directory):
        raise FileNotFoundError(f"Directory does not exist: {directory}")

    # Collect available output files
    if filetype.lower() == "gaussian":
        from chemsmart.io.gaussian.folder import GaussianLogFolder

        folder = GaussianLogFolder(directory)
        files = folder.all_logfiles
    elif filetype.lower() == "orca":
        from chemsmart.io.orca.folder import ORCAOutFolder

        folder = ORCAOutFolder(directory)
        files = folder.all_outfiles
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
        logger.error("No valid data parsed. Aborting export.")
        return None

    # Prepare output base name (DataExporter will append extension)
    fmt = _deduce_format(output)
    base, ext = os.path.splitext(output)
    outputfile = base if ext else output

    if keys:
        keys_list = [k.strip() for k in keys.split(",") if k.strip()]
    else:
        keys_list = None

    exporter = DataExporter(rows, outputfile=outputfile, keys=keys_list)
    if fmt == "json":
        exporter.to_json()
    elif fmt == "csv":
        exporter.to_csv()
    elif fmt == "sqlite":
        exporter.to_sqlite()
    logger.info(f"Export completed: {outputfile}.{fmt}")
    return None


def _deduce_format(output):
    ext = Path(output).suffix.lower()
    if ext in {".json", ".csv", ".sqlite"}:
        return ext[1:]
    elif ext == "":
        logger.warning("No extension specified. Defaulting to JSON format.")
        return "json"
    raise ValueError(
        f"Unsupported output format: {ext}. Only .json .csv .sqlite are supported."
    )
