"""
ORCA Modified Redundant Coordinate (Modred) CLI Module

This module provides the command-line interface for ORCA modified redundant
coordinate calculations. Modred calculations allow users to constrain or
modify specific internal coordinates during geometry optimizations, enabling
controlled structural changes and conformational searches.
"""

import logging

import click

from chemsmart.cli.job import click_job_options
from chemsmart.cli.orca.orca import click_orca_jobtype_options, orca
from chemsmart.utils.cli import MyCommand, get_setting_from_jobtype_for_orca
from chemsmart.utils.utils import check_charge_and_multiplicity

logger = logging.getLogger(__name__)


@orca.command("modred", cls=MyCommand)
@click_job_options
@click_orca_jobtype_options
@click.pass_context
def modred(
    ctx, jobtype, coordinates, dist_start, dist_end, num_steps, **kwargs
):
    """
    Run ORCA modified redundant coordinate (modred) calculations.

    This command performs geometry optimizations with modified redundant
    coordinates, allowing users to constrain or fix specific internal
    coordinates during the optimization process. This is useful for
    conformational searches, reaction coordinate analysis, and studying
    specific structural changes.
    """
    # set default job type if not specified
    if jobtype is None:
        jobtype = "modred"
        logger.debug("Using default jobtype: modred")

    # get modred settings from project configuration based on job type
    project_settings = ctx.obj["project_settings"]
    modred_settings = get_setting_from_jobtype_for_orca(
        project_settings, jobtype, coordinates, dist_start, dist_end, num_steps
    )
    logger.debug(
        f"Loaded modred settings for jobtype {jobtype}: {modred_settings}"
    )

    # job setting from filename or default, with updates from user in cli
    # specified in keywords
    # e.g., `chemsmart sub orca -c <user_charge> -m <user_multiplicity> modred`
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]

    # merge project modred settings with job settings from cli keywords
    modred_settings = modred_settings.merge(job_settings, keywords=keywords)
    logger.info(f"Final modred settings: {modred_settings.__dict__}")

    # validate charge and multiplicity consistency
    check_charge_and_multiplicity(modred_settings)

    # get molecule from context (use the last molecule if multiple)
    molecules = ctx.obj["molecules"]
    molecule = molecules[-1]
    logger.info(f"Running modred calculation on molecule: {molecule}")

    # get label for the job output files
    label = ctx.obj["label"]
    logger.debug(f"Job label: {label}")

    from chemsmart.jobs.orca.modred import ORCAModredJob

    job = ORCAModredJob(
        molecule=molecule, settings=modred_settings, label=label, **kwargs
    )
    logger.debug(f"Created ORCA modred job: {job}")
    return job
