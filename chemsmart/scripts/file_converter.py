#!/usr/bin/env python
import logging
import os

import click

from chemsmart.io.converter import FileConverter
from chemsmart.utils.logger import create_logger

logger = logging.getLogger(__name__)
os.environ["OMP_NUM_THREADS"] = "1"


@click.command()
@click.option(
    "-d",
    "--directory",
    default=None,
    help="Directory in which to convert files.",
)
@click.option(
    "-t",
    "--type",
    default=None,
    help="Type of file to be converted, if directory is specified.",
)
@click.option(
    "-f", "--filename", default=None, help="Input filename to be converted."
)
@click.option(
    "-o",
    "--output-filetype",
    default="xyz",
    help="Type of files to convert to, defaults to .xyz",
)
@click.option(
    "-i/",
    "--include-intermediate-structures/--no-include-intermediate-structures",
    is_flag=True,
    type=bool,
    default=False,
    help="Include intermediate structures in the conversion.",
)
def entry_point(
    directory, type, filename, output_filetype, include_intermediate_structures
):
    """
    Script for converting structures in different formats.
    """
    create_logger()
    if directory is not None:
        logger.info(f"Converting files in directory: {directory}")
        assert (
            type is not None
        ), "Type of file to be converted must be specified."
        assert filename is None, "Filename cannot be specified with directory."
        converter = FileConverter(
            directory=directory,
            type=type,
            output_filetype=output_filetype,
            include_intermediate_structures=include_intermediate_structures,
        )
        converter.convert_files()
    else:
        assert (
            filename is not None
        ), "Filename must be specified, since directory is not specified."
        logger.info(f"Converting file: {filename}")
        converter = FileConverter(
            filename=filename,
            output_filetype=output_filetype,
            include_intermediate_structures=include_intermediate_structures,
        )
        converter.convert_files()


if __name__ == "__main__":
    entry_point()
