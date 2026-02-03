"""
ORCA Single Point Calculation CLI Module

This module provides the command-line interface for ORCA single point
energy calculations. Single point calculations compute the energy and
properties of a molecule at a fixed geometry without optimization.
"""

import logging

import click

from chemsmart.cli.job import click_job_options
from chemsmart.cli.orca.orca import orca
from chemsmart.cli.orca.qmmm import create_orca_qmmm_subcommand
from chemsmart.utils.cli import MyGroup
from chemsmart.utils.utils import check_charge_and_multiplicity

logger = logging.getLogger(__name__)


@orca.group("sp", cls=MyGroup, invoke_without_command=True)
@click_job_options
@click.pass_context
def sp(ctx, skip_completed, **kwargs):
    """
    Run ORCA single point energy calculations.

    This command performs single point energy calculations using ORCA
    at the current molecular geometry without any optimization. It's
    useful for computing energies, electronic properties, and other
    molecular characteristics at a fixed structure.

    The calculation uses settings from the project configuration merged
    with any command-line overrides. Charge and multiplicity are validated
    before running the calculation.
    """
    # get single point settings from project configuration
    project_settings = ctx.obj["project_settings"]
    sp_settings = project_settings.sp_settings()
    logger.debug(f"Loaded single point settings from project: {sp_settings}")

    # job setting from filename or default, with updates from user in cli
    # specified in keywords
    # e.g., `chemsmart sub orca -c <user_charge> -m <user_multiplicity> sp`
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]

    # merge project sp settings with job settings from cli keywords
    sp_settings = sp_settings.merge(job_settings, keywords=keywords)
    logger.info(f"Final single point settings: {sp_settings.__dict__}")

    ctx.obj["parent_skip_completed"] = skip_completed
    ctx.obj["parent_kwargs"] = kwargs
    ctx.obj["parent_settings"] = sp_settings
    ctx.obj["parent_jobtype"] = "sp"

    if ctx.invoked_subcommand is not None:
        return

    # validate charge and multiplicity consistency only for direct sp jobs
    check_charge_and_multiplicity(sp_settings)

    # get molecules from context
    molecules = ctx.obj["molecules"]

    # get label for the job output files
    label = ctx.obj["label"]

    from chemsmart.jobs.orca.singlepoint import ORCASinglePointJob
    from chemsmart.utils.cli import create_sp_label

    # Get the original molecule indices from context
    molecule_indices = ctx.obj["molecule_indices"]

    # Handle multiple molecules: create one job per molecule
    if len(molecules) > 1 and molecule_indices is not None:
        logger.info(f"Creating {len(molecules)} ORCA single point jobs")
        jobs = []
        for molecule, idx in zip(molecules, molecule_indices):
            molecule_label = f"{label}_idx{idx}"
            final_label = create_sp_label(molecule_label, sp_settings)
            logger.info(
                f"Running single point for molecule {idx}: {molecule} with label {final_label}"
            )

            job = ORCASinglePointJob(
                molecule=molecule,
                settings=sp_settings,
                label=final_label,
                skip_completed=skip_completed,
                **kwargs,
            )
            jobs.append(job)
        logger.debug(f"Created {len(jobs)} ORCA single point jobs")
        return jobs
    else:
        # Single molecule case
        molecule = molecules[-1]
        label = create_sp_label(label, sp_settings)
        logger.info(
            f"Running single point calculation on molecule: {molecule} with label: {label}"
        )

        job = ORCASinglePointJob(
            molecule=molecule,
            settings=sp_settings,
            label=label,
            skip_completed=skip_completed,
            **kwargs,
        )
        logger.debug(f"Created ORCA single point job: {job}")
        return job


create_orca_qmmm_subcommand(sp)
