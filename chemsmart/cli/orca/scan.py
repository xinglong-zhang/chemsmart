"""
ORCA Coordinate Scan CLI Module

This module provides the command-line interface for ORCA coordinate scan
calculations. Scan calculations systematically vary specific internal
coordinates (bonds, angles, dihedrals) while optimizing the rest of the
structure, creating potential energy surfaces along reaction coordinates.
"""

import logging

import click

from chemsmart.cli.job import click_job_options
from chemsmart.cli.orca.orca import click_orca_jobtype_options, orca
from chemsmart.utils.cli import MyCommand, get_setting_from_jobtype_for_orca
from chemsmart.utils.utils import check_charge_and_multiplicity

logger = logging.getLogger(__name__)


@orca.command("scan", cls=MyCommand)
@click_job_options
@click_orca_jobtype_options
@click.pass_context
def scan(ctx, jobtype, coordinates, dist_start, dist_end, num_steps, **kwargs):
    """
    Run ORCA coordinate scan calculations.

    This command performs systematic coordinate scans by varying specific
    internal coordinates (bonds, angles, dihedrals) while optimizing the
    rest of the molecular structure. Scan calculations are essential for
    generating potential energy surfaces, studying conformational changes,
    and analyzing reaction pathways.

    The calculation uses settings from the project configuration merged
    with command-line overrides. Coordinate specifications and scan
    parameters are inherited from the jobtype options decorator.
    """
    # set default job type if not specified
    if jobtype is None:
        jobtype = "scan"
        logger.debug("Using default jobtype: scan")

    # get scan settings from project configuration based on job type
    project_settings = ctx.obj["project_settings"]
    scan_settings = get_setting_from_jobtype_for_orca(
        project_settings, jobtype, coordinates, dist_start, dist_end, num_steps
    )
    logger.debug(
        f"Loaded scan settings for jobtype {jobtype}: {scan_settings}"
    )

    # job setting from filename or default, with updates from user in cli
    # specified in keywords
    # e.g., `chemsmart sub orca -c <user_charge> -m <user_multiplicity> scan`
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]

    # merge project scan settings with job settings from cli keywords
    scan_settings = scan_settings.merge(job_settings, keywords=keywords)
    logger.info(f"Final scan settings: {scan_settings.__dict__}")

    # validate charge and multiplicity consistency
    check_charge_and_multiplicity(scan_settings)

    # get molecule from context (use the last molecule if multiple)
    molecules = ctx.obj["molecules"]
    molecule = molecules[-1]
    logger.info(f"Running coordinate scan on molecule: {molecule}")

    # get label for the job output files
    label = ctx.obj["label"]
    logger.debug(f"Job label: {label}")

    from chemsmart.jobs.orca.scan import ORCAScanJob

    job = ORCAScanJob(
        molecule=molecule, settings=scan_settings, label=label, **kwargs
    )
    logger.debug(f"Created ORCA scan job: {job}")
    return job
