import logging

import click

from chemsmart.cli.gaussian.gaussian import gaussian
from chemsmart.cli.gaussian.qmmm import create_qmmm_subcommand
from chemsmart.cli.job import click_job_options
from chemsmart.utils.cli import MyGroup
from chemsmart.utils.utils import check_charge_and_multiplicity

logger = logging.getLogger(__name__)


@gaussian.group("ts", cls=MyGroup, invoke_without_command=True)
@click_job_options
@click.option(
    "-f",
    "--freeze-atoms",
    type=str,
    help="Indices of atoms to freeze for constrained optimization. 1-indexed.",
)
@click.pass_context
def ts(ctx, freeze_atoms, skip_completed, **kwargs):
    """CLI subcommand for running Gaussian transition state calculation.

    Can be used standalone for regular TS search or with the 'qmmm'
    subcommand for QM/MM TS calculations.

    Examples:
        chemsmart sub gaussian ts              # Regular TS search
        chemsmart sub gaussian ts qmmm         # QM/MM TS search
    """

    # get jobrunner for transition state calculation
    jobrunner = ctx.obj["jobrunner"]
    # get settings from project
    project_settings = ctx.obj["project_settings"]
    ts_settings = project_settings.ts_settings()

    # job setting from filename or default, with updates from user in cli
    # specified in keywords
    # e.g., `sub.py gaussian -c <user_charge> -m <user_multiplicity>`
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]

    # merge project opt settings with job settings from cli keywords from
    # cli.gaussian.py subcommands
    ts_settings = ts_settings.merge(job_settings, keywords=keywords)

    if ctx.invoked_subcommand is not None:
        return

    check_charge_and_multiplicity(ts_settings)

    # get molecule
    molecules = ctx.obj["molecules"]
    molecule = molecules[-1]

    # get label for the job
    label = ctx.obj["label"]
    logger.debug(f"Label for job: {label}")

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

    logger.info(f"TS job settings from project: {ts_settings.__dict__}")

    # Store parent context for potential qmmm subcommand
    ctx.obj["parent_skip_completed"] = skip_completed
    ctx.obj["parent_freeze_atoms"] = freeze_atoms
    ctx.obj["parent_kwargs"] = kwargs
    ctx.obj["parent_settings"] = ts_settings

    # If no subcommand invoked, run regular TS search
    if ctx.invoked_subcommand is None:
        from chemsmart.jobs.gaussian.ts import GaussianTSJob

        return GaussianTSJob(
            molecule=molecule,
            settings=ts_settings,
            label=label,
            jobrunner=jobrunner,
            skip_completed=skip_completed,
            **kwargs,
        )


create_qmmm_subcommand(ts)
