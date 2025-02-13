#!/usr/bin/env python
import logging
import os
import shlex
import subprocess

import click

from chemsmart.utils.logger import create_logger

logger = logging.getLogger(__name__)
os.environ["OMP_NUM_THREADS"] = "1"


@click.command()
@click.option(
    "-f",
    "--filename",
    required=True,
    type=str,
    help=".txt file containing the names of files to be submitted.",
)
@click.option(
    "-c",
    "--command",
    type=str,
    required=True,
    default=None,
    help="submission command for all jobs.",
)
def entry_point(filename, command):
    create_logger()
    logger.info(f"Reading filenames from {filename}")
    with open(filename, "r") as f:
        filenames = f.readlines()
    filenames = [filename.strip() for filename in filenames]

    for filename in filenames:
        logger.info(f"Submitting job for {filename}")
        cmd = command.replace("file", filename)
        logger.info(f"Command executed: {cmd}")
        # use subprocess to run the command
        p = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE)
        p.communicate()
        logger.info(f"Job submitted for {filename}")


if __name__ == "__main__":
    entry_point()
