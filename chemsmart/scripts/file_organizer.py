#!/usr/bin/env python
"""
File organization script.

This script provides command-line functionality to organize
computational chemistry files into structured directories.
"""

import logging
import os

import click

from chemsmart.io.organizer import FileOrganizer
from chemsmart.utils.logger import create_logger

logger = logging.getLogger(__name__)
os.environ["OMP_NUM_THREADS"] = "1"


@click.command()
@click.option(
    "-d",
    "--directory",
    default=".",
    help="Directory in which to organize the files.",
)
@click.option(
    "-f",
    "--filename",
    default=None,
    required=True,
    help="Filename of Excel file to use for organizing.",
)
@click.option(
    "-n",
    "--name",
    default=None,
    required=True,
    help="Sheet name of Excel file to use for organizing.",
)
@click.option(
    "-t",
    "--type",
    default="log",
    help="Type of file to be organized.",
)
@click.option(
    "-c",
    "--cols",
    default=None,
    help="Columns to be used.",
)
@click.option(
    "-s",
    "--skip",
    type=int,
    default=2,
    help="Number of rows to skip.",
)
@click.option(
    "-r",
    "--row",
    type=int,
    default=100,
    help="Number of rows to organize.",
)
@click.option(
    "--keep-default-na/--no-keep-default-na",
    type=bool,
    default=False,
    help="Keep default NA values in Excel file.",
)
def entry_point(
    directory, filename, name, type, cols, skip, row, keep_default_na
):
    """
    Script for organizing files for supporting information.

    Example usage:
        file_organizer.py -f jq.xlsx -n co2 -c B:D -s 2 -r 100
    """

    create_logger()
    logger.info(f"Organizing files in directory: {directory}")

    organizer = FileOrganizer(
        directory=directory,
        filename=filename,
        sheetname=name,
        type=type,
        cols=cols,
        skip=skip,
        row=row,
        keep_default_na=keep_default_na,
    )
    organizer.organize_files()


if __name__ == "__main__":
    entry_point()
