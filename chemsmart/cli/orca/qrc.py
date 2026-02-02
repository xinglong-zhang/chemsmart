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
from chemsmart.cli.orca.qmmm import create_orca_qmmm_subcommand
from chemsmart.utils.cli import MyGroup, get_setting_from_jobtype_for_orca
from chemsmart.utils.utils import check_charge_and_multiplicity

logger = logging.getLogger(__name__)


@orca.group("qrc", cls=MyGroup, invoke_without_command=True)
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
    CLI subcommand for running ORCA quick reaction coordinate (QRC) calculations.

    Examples:
        `chemsmart sub orca -p proj -f ts.log qrc` runs QRC optimization jobs
        after displacing along vibrational mode 1 (default) by +/- 0.5 Angstroms (default).

        `chemsmart sub orca -p proj -f ts.log qrc -m 2 -a 1.2` runs QRC optimization jobs
        after displacing along vibrational mode 2 (-m 2) by +/- 1.2 Angstroms (-a 1.2).

        `chemsmart sub orca -p proj -f ts.log qrc -j ts -m 2 -a 1.5` runs QRC TS jobs
        after displacing along vibrational mode 2 (-m 2) by +/- 1.5 Angstroms (-a 1.5).

        `chemsmart sub orca -p proj -f ts.log qrc -j modred -c [1,2] -m 2 -a 1.5`
        runs QRC modred jobs after displacing along vibrational mode 2 by +/- 1.5 Angstroms.
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

    ctx.obj["parent_skip_completed"] = skip_completed
    ctx.obj["parent_kwargs"] = kwargs
    ctx.obj["parent_settings"] = qrc_settings
    ctx.obj["parent_jobtype"] = jobtype

    if ctx.invoked_subcommand is not None:
        return

    # validate charge and multiplicity consistency only for direct qrc jobs
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


create_orca_qmmm_subcommand(qrc)
