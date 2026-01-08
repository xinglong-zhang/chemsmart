"""
xTB Molecular Dynamics CLI Module

This module provides the command-line interface for xTB molecular
dynamics simulations.
"""

import logging

import click

from chemsmart.cli.job import click_job_options
from chemsmart.cli.xtb.xtb import xtb
from chemsmart.utils.cli import MyCommand
from chemsmart.utils.utils import check_charge_and_multiplicity

logger = logging.getLogger(__name__)


@xtb.command("md", cls=MyCommand)
@click_job_options
@click.pass_context
def md(ctx, skip_completed, **kwargs):
    """
    Run xTB molecular dynamics simulations.

    This command performs molecular dynamics simulations using xTB.

    The simulation uses settings from the project configuration merged
    with any command-line overrides. Charge and multiplicity are validated
    before running the calculation.
    """
    # get settings from project
    project_settings = ctx.obj["project_settings"]
    md_settings = project_settings.md_settings()

    # job setting from filename or default, with updates from user in cli specified in keywords
    # e.g., `sub.py xtb -c <user_charge> -m <user_multiplicity>`
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]

    # merge project md settings with job settings from cli keywords from cli.xtb.py subcommands
    md_settings = md_settings.merge(job_settings, keywords=keywords)

    check_charge_and_multiplicity(md_settings)

    # get molecule
    molecules = ctx.obj["molecules"]
    molecule = molecules[-1]
    logger.info(f"Running molecular dynamics for molecule: {molecule}.")

    # get label for the job
    label = ctx.obj["label"]

    logger.info(f"MD job settings from project: {md_settings.__dict__}")

    from chemsmart.jobs.xtb.md import XTBMDJob

    return XTBMDJob(
        molecule=molecule,
        settings=md_settings,
        label=label,
        skip_completed=skip_completed,
        **kwargs,
    )
