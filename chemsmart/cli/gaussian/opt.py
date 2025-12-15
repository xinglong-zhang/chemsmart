import logging

import click

from chemsmart.cli.gaussian.gaussian import gaussian
from chemsmart.cli.job import click_job_options
from chemsmart.utils.cli import MyCommand
from chemsmart.utils.utils import check_charge_and_multiplicity

logger = logging.getLogger(__name__)


@gaussian.command("opt", cls=MyCommand)
@click_job_options
@click.option(
    "-f",
    "--freeze-atoms",
    type=str,
    help="Indices of atoms to freeze for constrained optimization. 1-indexed.",
)
@click.pass_context
def opt(ctx, freeze_atoms, skip_completed, **kwargs):
    """CLI subcommand for running Gaussian optimization calculation."""

    # get jobrunner for optimization
    jobrunner = ctx.obj["jobrunner"]

    # get settings from project
    project_settings = ctx.obj["project_settings"]
    opt_settings = project_settings.opt_settings()

    # job setting from filename or default, with updates from user in cli
    # specified in keywords
    # e.g., `sub.py gaussian -c <user_charge> -m <user_multiplicity>`
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]

    # merge project opt settings with job settings from cli keywords from
    # cli.gaussian.py subcommands
    opt_settings = opt_settings.merge(job_settings, keywords=keywords)

    check_charge_and_multiplicity(opt_settings)

    # get molecules
    molecules = ctx.obj["molecules"]

    # get label for the job
    label = ctx.obj["label"]

    # Set atoms to freeze

    from chemsmart.utils.utils import (
        convert_list_to_gaussian_frozen_list,
        get_list_from_string_range,
    )

    logger.info(f"Opt job settings from project: {opt_settings.__dict__}")

    from chemsmart.jobs.gaussian.opt import GaussianOptJob

    # Get the original molecule indices from context
    molecule_indices = ctx.obj.get("molecule_indices", list(range(1, len(molecules) + 1)))

    # Handle multiple molecules: create one job per molecule
    if len(molecules) > 1:
        logger.info(f"Creating {len(molecules)} optimization jobs")
        jobs = []
        for molecule, idx in zip(molecules, molecule_indices):
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

            job = GaussianOptJob(
                molecule=molecule,
                settings=opt_settings,
                label=molecule_label,
                jobrunner=jobrunner,
                skip_completed=skip_completed,
                **kwargs,
            )
            jobs.append(job)
        return jobs
    else:
        # Single molecule case
        molecule = molecules[-1].copy()
        logger.info(f"Optimizing molecule: {molecule}.")

        if freeze_atoms is not None:
            frozen_atoms_list = get_list_from_string_range(freeze_atoms)
            logger.debug(f"Freezing atoms: {frozen_atoms_list}")
            molecule.frozen_atoms = convert_list_to_gaussian_frozen_list(
                frozen_atoms_list, molecule
            )
        else:
            logger.debug("No atoms will be frozen during optimization")

        return GaussianOptJob(
            molecule=molecule,
            settings=opt_settings,
            label=label,
            jobrunner=jobrunner,
            skip_completed=skip_completed,
            **kwargs,
        )
