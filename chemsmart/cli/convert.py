import logging

import click

from chemsmart.io.converter import FileConverter

logger = logging.getLogger(__name__)


@click.command(name="convert")
@click.option(
    "-i",
    "--input",
    "input_filepath",
    type=str,
    default=None,
    help="Input file path to convert (e.g. abc.pdb).",
)
@click.option(
    "-o",
    "--output",
    "output_filepath",
    type=str,
    default=None,
    help="Output file path. The format is inferred from the file extension "
    "(e.g. xyz.xyz). When omitted, the output filename is derived from the "
    "input filename using --output-filetype.",
)
@click.option(
    "-d",
    "--directory",
    type=str,
    default=None,
    help="Directory containing files to convert in batch mode.",
)
@click.option(
    "-t",
    "--filetype",
    "type",
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
    help="Output file type used when --output is not specified "
    "(e.g. xyz, com, pdb).",
)
@click.option(
    "-z/",
    "--include-intermediate-structures/--no-include-intermediate-structures",
    default=False,
    show_default=True,
    help="Include all intermediate structures from multi-structure files.",
)
def convert(
    input_filepath,
    output_filepath,
    directory,
    type,
    program,
    output_filetype,
    include_intermediate_structures,
):
    """Convert molecular files between supported formats via Molecule object.

    Single-file conversion with explicit output path:

    \b
        chemsmart run convert --input abc.pdb --output xyz.xyz
        chemsmart run convert -i molecule.log -o molecule.xyz

    Batch directory conversion:

    \b
        chemsmart run convert --directory /path/to/dir --filetype log --output-filetype xyz
    """
    file_converter = FileConverter(
        directory=directory,
        type=type,
        program=program,
        filename=input_filepath,
        output_filetype=output_filetype,
        output_filepath=output_filepath,
        include_intermediate_structures=include_intermediate_structures,
    )
    file_converter.convert_files()
    return None
