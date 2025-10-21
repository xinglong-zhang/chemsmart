"""
ORCA Modified Redundant Coordinate (Modred) CLI Module

This module provides the command-line interface for ORCA modified redundant
coordinate calculations. Modred calculations allow users to constrain or
modify specific internal coordinates during geometry optimizations, enabling
controlled structural changes and conformational searches.
"""

import logging

import click

from chemsmart.cli.job import (
    click_job_options,
    click_molecule_vibrational_displacement_options,
)
from chemsmart.cli.orca.orca import click_orca_jobtype_options, orca
from chemsmart.utils.cli import MyCommand, get_setting_from_jobtype_for_orca
from chemsmart.utils.utils import check_charge_and_multiplicity

logger = logging.getLogger(__name__)


@orca.command("qrc", cls=MyCommand)
@click_job_options
@click_orca_jobtype_options
@click_molecule_vibrational_displacement_options
@click.pass_context
def qrc(
    ctx,
    jobtype,
    coordinates,
    dist_start,
    dist_end,
    num_steps,
    mode_idx,
    amp,
    nframes,
    phase,
    normalize,
    return_xyz,
    skip_completed,
    **kwargs,
):
    """
     Run ORCA quick reaction coordinate (qrc) calculations.

    `chemsmart sub gaussian -p proj -f ts.log qrc` will run QRC optimization jobs
     after displacing along vibrational mode 1 (default) by ± 0.5 Å (default).
     `chemsmart sub gaussian -p proj -f ts.log qrc -m 2 -a 1.2` will run QRC opt jobs
     after displacing along vibrational mode 2 (-m 2) by ± 1.2 Å (-a 1.2).
     `chemsmart sub gaussian -p proj -f ts.log qrc -j ts -m 2 -a 1.5` will run QRC TS jobs
     after displacing along vibrational mode 2 (-m 2) by ± 1.5 Å (-a 1.5).
     `chemsmart sub gaussian -p proj -f ts.log qrc -j modred -c [1,2] -m 2 -a 1.5`
     will run QRC modred jobs after displacing along vibrational mode 2 by ± 1.5 Å.
    """

    # get jobrunner for running Gaussian crest jobs
    jobrunner = ctx.obj["jobrunner"]

    # get settings from project
    project_settings = ctx.obj["project_settings"]

    # get job settings from job type
    if jobtype is None:
        jobtype = "opt"  # default using opt job for QRC, if jobtype not given

    qrc_settings = get_setting_from_jobtype_for_orca(
        project_settings, jobtype, coordinates, dist_start, dist_end, num_steps
    )
    logger.debug(f"Loaded QRC settings for jobtype {jobtype}: {qrc_settings}")

    # job setting from filename or default, with updates from user in cli
    # specified in keywords
    # e.g., `chemsmart sub orca -c <user_charge> -m <user_multiplicity> modred`
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]

    # merge project modred settings with job settings from cli keywords
    qrc_settings = qrc_settings.merge(job_settings, keywords=keywords)
    logger.info(f"Final QRC settings: {qrc_settings.__dict__}")

    # validate charge and multiplicity consistency
    check_charge_and_multiplicity(qrc_settings)

    # get molecule from context (use the last molecule if multiple)
    molecules = ctx.obj["molecules"]
    molecule = molecules[-1]
    logger.info(f"Running QRC calculation on molecule: {molecule}")

    # get label for the job output files
    label = ctx.obj["label"]
    logger.debug(f"Job label: {label}")

    from chemsmart.jobs.orca.qrc import ORCAQRCJob

    return ORCAQRCJob(
        molecule=molecule,
        settings=qrc_settings,
        label=label,
        jobrunner=jobrunner,
        mode_idx=mode_idx,
        amp=amp,
        nframes=nframes,
        phase=phase,
        normalize=normalize,
        return_xyz=return_xyz,
        skip_completed=skip_completed,
        **kwargs,
    )
