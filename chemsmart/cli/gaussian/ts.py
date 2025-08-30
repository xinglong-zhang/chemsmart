import logging

import click

from chemsmart.cli.gaussian.gaussian import gaussian
from chemsmart.cli.job import click_job_options
from chemsmart.utils.cli import MyCommand
from chemsmart.utils.utils import check_charge_and_multiplicity

logger = logging.getLogger(__name__)


@gaussian.command("ts", cls=MyCommand)
@click_job_options
@click.option(
    "-f",
    "--freeze-atoms",
    type=str,
    help="Indices of atoms to freeze for constrained optimization.",
)
@click.pass_context
def ts(ctx, freeze_atoms, skip_completed, **kwargs):
    """CLI for transition state calculation for Gaussian."""

    # get jobrunner for transition state calculation
    jobrunner = ctx.obj["jobrunner"]
    # get settings from project
    project_settings = ctx.obj["project_settings"]
    ts_settings = project_settings.ts_settings()

    # job setting from filename or default, with updates from user in cli specified in keywords
    # e.g., `sub.py gaussian -c <user_charge> -m <user_multiplicity>`
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]

    # merge project opt settings with job settings from cli keywords from cli.gaussian.py subcommands
    ts_settings = ts_settings.merge(job_settings, keywords=keywords)

    check_charge_and_multiplicity(ts_settings)

    # get molecule
    molecules = ctx.obj["molecules"]
    molecule = molecules[
        -1
    ]  # get last molecule from list of molecules from cli.gaussian.py subcommands
    # index = '-1' would access the right structure from the list of molecule returned from cli.gaussian.py subcommands
    # user specified index was used there to return the right molecule and store it as a list of single element/itself

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

    from chemsmart.jobs.gaussian.ts import GaussianTSJob

    return GaussianTSJob(
        molecule=molecule,
        settings=ts_settings,
        label=label,
        jobrunner=jobrunner,
        skip_completed=skip_completed,
        **kwargs,
    )
