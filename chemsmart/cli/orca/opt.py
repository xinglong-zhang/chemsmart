"""
ORCA Geometry Optimization CLI Module

This module provides the command-line interface for ORCA geometry
optimization calculations. It supports both unconstrained and constrained
optimizations with options for freezing specific atoms during the
optimization process.
"""

import logging

import click

from chemsmart.cli.job import click_job_options
from chemsmart.cli.orca.orca import orca
from chemsmart.utils.cli import MyCommand
from chemsmart.utils.utils import check_charge_and_multiplicity

logger = logging.getLogger(__name__)


@orca.command("opt", cls=MyCommand)
@click_job_options
@click.option(
    "-f",
    "--freeze-atoms",
    type=str,
    help="Indices of atoms to freeze for constrained optimization. "
    "1-indexed.",
)
@click.option(
    "-i",
    "--invert-constraints/--no-invert-constraints",
    default=False,
    type=bool,
    help="Invert the constraints for frozen atoms in optimization.",
)
@click.pass_context
def opt(ctx, freeze_atoms, invert_constraints, skip_completed, **kwargs):
    """
    Run ORCA geometry optimization calculations.

    This command performs geometry optimization using ORCA with support
    for both unconstrained and constrained optimizations. Users can
    specify atoms to freeze during optimization and configure various
    optimization parameters.

    The optimization uses settings from the project configuration merged
    with any command-line overrides. Charge and multiplicity are validated
    before running the calculation.
    """
    # get optimization settings from project configuration
    project_settings = ctx.obj["project_settings"]
    opt_settings = project_settings.opt_settings()
    logger.debug(f"Loaded optimization settings from project: {opt_settings}")

    # job setting from filename or default, with updates from user in cli
    # specified in keywords
    # e.g., `chemsmart sub orca -c <user_charge> -m <user_multiplicity> opt`
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]

    # merge project opt settings with job settings from cli keywords
    opt_settings = opt_settings.merge(job_settings, keywords=keywords)
    opt_settings.invert_constraints = invert_constraints
    logger.info(f"Final optimization settings: {opt_settings.__dict__}")

    # validate charge and multiplicity consistency
    check_charge_and_multiplicity(opt_settings)

    # get molecules from context
    molecules = ctx.obj["molecules"]

    # get label for the job output files
    label = ctx.obj["label"]

    # Set atoms to freeze for constrained optimization
    from chemsmart.jobs.orca.opt import ORCAOptJob
    from chemsmart.utils.utils import (
        convert_list_to_gaussian_frozen_list,
        get_list_from_string_range,
    )

    # Handle multiple molecules: create one job per molecule
    if len(molecules) > 1:
        logger.info(f"Creating {len(molecules)} ORCA optimization jobs")
        jobs = []
        for idx, molecule in enumerate(molecules, start=1):
            # Create a copy to avoid side effects from mutation
            molecule = molecule.copy()
            molecule_label = f"{label}_idx{idx}"
            logger.info(
                f"Optimizing molecule {idx}: {molecule} with label {molecule_label}"
            )

            # Apply frozen atoms if specified
            if freeze_atoms is not None:
                frozen_atoms_list = get_list_from_string_range(freeze_atoms)
                logger.debug(f"Freezing atoms: {frozen_atoms_list}")
                molecule.frozen_atoms = convert_list_to_gaussian_frozen_list(
                    frozen_atoms_list, molecule
                )
            else:
                logger.debug("No atoms will be frozen during optimization")

            job = ORCAOptJob(
                molecule=molecule,
                settings=opt_settings,
                label=molecule_label,
                skip_completed=skip_completed,
                **kwargs,
            )
            jobs.append(job)
        logger.debug(f"Created {len(jobs)} ORCA optimization jobs")
        return jobs
    else:
        # Single molecule case
        molecule = molecules[-1]
        molecule = molecule.copy()
        logger.info(f"Optimizing molecule: {molecule}")

        if freeze_atoms is not None:
            frozen_atoms_list = get_list_from_string_range(freeze_atoms)
            logger.debug(f"Freezing atoms: {frozen_atoms_list}")
            molecule.frozen_atoms = convert_list_to_gaussian_frozen_list(
                frozen_atoms_list, molecule
            )
        else:
            logger.debug("No atoms will be frozen during optimization")

        job = ORCAOptJob(
            molecule=molecule,
            settings=opt_settings,
            label=label,
            skip_completed=skip_completed,
            **kwargs,
        )
        logger.debug(f"Created ORCA optimization job: {job}")
        return job
