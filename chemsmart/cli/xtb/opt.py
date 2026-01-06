"""
xTB Geometry Optimization CLI Module

This module provides the command-line interface for xTB geometry
optimization calculations. It supports both unconstrained and constrained
optimizations with options for freezing specific atoms during the
optimization process.
"""

import logging

import click

from chemsmart.cli.job import click_job_options
from chemsmart.cli.xtb.xtb import xtb
from chemsmart.utils.cli import MyCommand
from chemsmart.utils.utils import check_charge_and_multiplicity

logger = logging.getLogger(__name__)


@xtb.command("opt", cls=MyCommand)
@click_job_options
@click.option(
    "-f",
    "--freeze-atoms",
    type=str,
    help="Indices of atoms to freeze for constrained optimization.",
)
@click.pass_context
def opt(ctx, freeze_atoms, skip_completed, **kwargs):
    """
    Run xTB geometry optimization calculations.

    This command performs geometry optimization using xTB with support
    for both unconstrained and constrained optimizations. Users can
    specify atoms to freeze during optimization and configure various
    optimization parameters.

    The optimization uses settings from the project configuration merged
    with any command-line overrides. Charge and multiplicity are validated
    before running the calculation.
    """
    # get settings from project
    project_settings = ctx.obj["project_settings"]
    opt_settings = project_settings.opt_settings()

    # job setting from filename or default, with updates from user in cli specified in keywords
    # e.g., `sub.py xtb -c <user_charge> -m <user_multiplicity>`
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]

    # merge project opt settings with job settings from cli keywords from cli.xtb.py subcommands
    opt_settings = opt_settings.merge(job_settings, keywords=keywords)

    check_charge_and_multiplicity(opt_settings)

    # get molecule
    molecules = ctx.obj["molecules"]
    molecule = molecules[-1]
    logger.info(f"Optimizing molecule: {molecule}.")

    # get label for the job
    label = ctx.obj["label"]

    # Set atoms to freeze

    from chemsmart.utils.utils import (
        convert_list_to_gaussian_frozen_list,
        get_list_from_string_range,
    )

    if freeze_atoms is not None:
        frozen_atoms_list = get_list_from_string_range(freeze_atoms)
        logger.debug(f"Freezing atoms: {frozen_atoms_list}")
        molecule.frozen_atoms = convert_list_to_gaussian_frozen_list(
            frozen_atoms_list, molecule
        )

    logger.info(f"Opt job settings from project: {opt_settings.__dict__}")

    from chemsmart.jobs.xtb.opt import XTBOptJob

    return XTBOptJob(
        molecule=molecule,
        settings=opt_settings,
        label=label,
        skip_completed=skip_completed,
        **kwargs,
    )
