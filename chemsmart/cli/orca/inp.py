"""
ORCA Input File Job CLI Module

This module provides the command-line interface for running ORCA jobs
directly from existing input files. It allows users to execute ORCA
calculations using pre-prepared .inp files without modification.

The module defines:
    - inp: Command for executing ORCA input files as-is

This is useful when you have already prepared ORCA input files and want
to run them through the chemsmart job management system without any
additional processing or modification.
"""

import logging

import click

from chemsmart.cli.job import click_job_options
from chemsmart.cli.orca.orca import orca
from chemsmart.utils.cli import MyCommand

logger = logging.getLogger(__name__)


@orca.command("inp", cls=MyCommand)
@click_job_options
@click.pass_context
def inp(ctx, skip_completed, **kwargs):
    """
    Run an ORCA input job as it is.
    Only requires the file that is to be run.
    """
    filename = ctx.obj["filename"]
    logger.info(f"Creating ORCA input job from file: {filename}")

    from chemsmart.jobs.orca.job import ORCAInpJob

    job = ORCAInpJob.from_filename(
        filename=filename, skip_completed=skip_completed, **kwargs
    )
    logger.debug(f"Created ORCA input job: {job}")
    return job
