"""
xTB Single Point Calculation CLI Module

This module provides the command-line interface for xTB single point
energy calculations.
"""

import logging

import click

from chemsmart.cli.job import click_job_options
from chemsmart.cli.xtb.xtb import xtb
from chemsmart.utils.cli import MyCommand
from chemsmart.utils.utils import check_charge_and_multiplicity

logger = logging.getLogger(__name__)


@xtb.command("sp", cls=MyCommand)
@click_job_options
@click.pass_context
def singlepoint(ctx, skip_completed, **kwargs):
    """
    Run xTB single point energy calculations.

    This command performs single point energy calculations using xTB.

    The calculation uses settings from the project configuration merged
    with any command-line overrides. Charge and multiplicity are validated
    before running the calculation.
    """
    # get settings from project
    project_settings = ctx.obj["project_settings"]
    sp_settings = project_settings.sp_settings()

    # job setting from filename or default, with updates from user in cli specified in keywords
    # e.g., `sub.py xtb -c <user_charge> -m <user_multiplicity>`
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]

    # merge project sp settings with job settings from cli keywords from cli.xtb.py subcommands
    sp_settings = sp_settings.merge(job_settings, keywords=keywords)

    check_charge_and_multiplicity(sp_settings)

    # get molecule
    molecules = ctx.obj["molecules"]
    molecule = molecules[-1]
    logger.info(f"Running single point calculation for molecule: {molecule}.")

    # get label for the job
    label = ctx.obj["label"]

    logger.info(
        f"Single point job settings from project: {sp_settings.__dict__}"
    )

    from chemsmart.jobs.xtb.singlepoint import XTBSinglePointJob

    return XTBSinglePointJob(
        molecule=molecule,
        settings=sp_settings,
        label=label,
        skip_completed=skip_completed,
        **kwargs,
    )
