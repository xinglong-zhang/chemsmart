#!/usr/bin/env python
"""
Cube file operations script.

This script provides command-line functionality for performing mathematical
operations on Gaussian cube files, such as addition, subtraction, and
other mathematical transformations.
"""

import logging
import os

import click

from chemsmart.io.gaussian.cube import CubeFileOperator

logger = logging.getLogger(__name__)
os.environ["OMP_NUM_THREADS"] = "1"


@click.command()
@click.option(
    "-c1",
    "--cube1",
    required=True,
    default=None,
    type=str,
    help="Cube file 1.",
)
@click.option(
    "-c2",
    "--cube2",
    required=True,
    default=None,
    type=str,
    help="Cube file 2.",
)
@click.option(
    "-x",
    "--operation",
    default="subtract",
    type=click.Choice(["subtract", "add"]),
    help="Operation type. Defaults to subtract.",
)
@click.option(
    "-o",
    "--outputname",
    default=None,
    type=str,
    help="outputname of the operated cube file.",
)
def entry_point(cube1, cube2, operation, outputname):
    """
    Perform mathematical operations on Gaussian cube files.
    """
    cube_operator = CubeFileOperator(
        cubefile1=cube1,
        cubefile2=cube2,
        operation=operation,
        output_cubefile=outputname,
    )
    logger.info(f"Operation done: {operation}")
    logger.info(f"Writing results to {outputname}")
    cube_operator.write_results()


if __name__ == "__main__":
    entry_point()
