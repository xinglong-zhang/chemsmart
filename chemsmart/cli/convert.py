"""
File format conversion CLI.

Registered as ``chemsmart run convert``. Converts a single molecular
structure file to another format using ``Molecule`` as the intermediate
representation.
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
    required=True,
    type=click.Path(exists=True),
    help="Input molecular structure file.",
)
@click.option(
    "-o",
    "--output",
    "output_file",
    required=True,
    type=click.Path(),
    help="Output file path (format is inferred from the extension).",
)
@logger_options
@click.pass_context
def convert(ctx, input_file, output_file, debug, stream):
    """
    Convert a molecular structure file to another format.

    The input format is detected from the input file extension, and the
    output format is detected from the output file extension. Conversion is
    performed via the ``Molecule`` class, so any format pair that can be
    read and written by ``Molecule`` is supported.

    \b
    Example:
        chemsmart run convert --input a.pdb --output a.xyz
    """
    create_logger(debug=debug, stream=stream)
    FileConverter.convert_file(input_file, output_file)
    return None
