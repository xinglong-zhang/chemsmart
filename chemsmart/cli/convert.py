"""
File format conversion CLI.

Registered as ``chemsmart run convert``. Converts molecular structure
files between formats using ``Molecule`` as the intermediate
representation. Supports both single-file and batch directory modes.
"""

import logging

import click

from chemsmart.cli.logger import logger_options
from chemsmart.io.converter import FileConverter
from chemsmart.utils.cli import MyCommand
from chemsmart.utils.logger import create_logger

logger = logging.getLogger(__name__)


@click.command(name="convert", cls=MyCommand)
@click.option(
    "-i",
    "--input",
    "input_file",
    type=click.Path(exists=True),
    default=None,
    help="Input molecular structure file for single-file conversion.",
)
@click.option(
    "-o",
    "--output",
    "output_file",
    type=click.Path(),
    default=None,
    help="Output file path (format is inferred from the extension). "
    "May be an extension only (e.g. .xyz or xyz) to keep the input "
    "basename. When omitted with --input, writes {input_stem}.xyz.",
)
@click.option(
    "-d",
    "--directory",
    type=click.Path(exists=True, file_okay=False),
    default=None,
    help="Directory containing files to convert in batch mode.",
)
@click.option(
    "-t",
    "--filetype",
    "filetype",
    type=str,
    default=None,
    help="Input file type for batch directory conversion "
    "(e.g. log, com, gjf, out, inp, xyz, sdf, pdb, cdxml, cdx).",
)
@click.option(
    "-p",
    "--program",
    type=str,
    default=None,
    help="Computational program (gaussian or orca). Required when "
    "--filetype is 'out', because both Gaussian and ORCA use this extension.",
)
@click.option(
    "--output-filetype",
    "output_filetype",
    type=str,
    default="xyz",
    show_default=True,
    help="Output file type when --output is not specified "
    "(e.g. xyz, com, pdb). Used for single-file and batch conversion.",
)
@click.option(
    "-z/--no-z",
    "--include-intermediate-structures/"
    "--no-include-intermediate-structures",
    "include_intermediate_structures",
    default=False,
    show_default=True,
    help="Include all intermediate structures from multi-structure files.",
)
@logger_options
@click.pass_context
def convert(
    ctx,
    input_file,
    output_file,
    directory,
    filetype,
    program,
    output_filetype,
    include_intermediate_structures,
    debug,
    stream,
):
    """
    Convert molecular structure files between formats.

    Single-file conversion (defaults to {input_stem}.xyz):

    \b
        chemsmart run convert -i molecule.pdb
        chemsmart run convert -i molecule.pdb -o .mol2
        chemsmart run convert --input a.pdb --output a.xyz

    Batch directory conversion:

    \b
        chemsmart run convert --directory /path/to/dir --filetype log \\
            --output-filetype xyz
    """
    create_logger(debug=debug, stream=stream)

    if input_file is not None and directory is not None:
        raise click.UsageError(
            "Provide either --input (single-file) or "
            "--directory/--filetype (batch), not both."
        )

    if input_file is not None:
        output_path = FileConverter.resolve_output_path(
            input_file,
            output=output_file,
            output_filetype=output_filetype,
        )
        FileConverter.convert_file(
            input_file,
            output_path,
            include_intermediate_structures=include_intermediate_structures,
        )
        return None

    if directory is not None:
        if filetype is None:
            raise click.UsageError(
                "--filetype is required when --directory is specified."
            )
        file_converter = FileConverter(
            directory=directory,
            type=filetype,
            program=program,
            output_filetype=output_filetype,
            include_intermediate_structures=include_intermediate_structures,
        )
        file_converter.convert_files()
        return None

    raise click.UsageError(
        "Provide either --input (single-file) or "
        "--directory/--filetype (batch)."
    )
