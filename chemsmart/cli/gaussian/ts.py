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
    """CLI subcommand for running Gaussian transition state calculation."""

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

    # get molecules
    molecules = ctx.obj["molecules"]

    # get label for the job
    label = ctx.obj["label"]
    logger.debug(f"Label for job: {label}")

    # Set atoms to freeze

    from chemsmart.utils.utils import (
        convert_list_to_gaussian_frozen_list,
        get_list_from_string_range,
    )

    logger.info(f"TS job settings from project: {ts_settings.__dict__}")

    from chemsmart.jobs.gaussian.ts import GaussianTSJob

    # Get the original molecule indices from context
    molecule_indices = ctx.obj["molecule_indices"]

    # Store parent context for potential qmmm subcommand
    ctx.obj["parent_skip_completed"] = skip_completed
    ctx.obj["parent_freeze_atoms"] = freeze_atoms
    ctx.obj["parent_kwargs"] = kwargs
    ctx.obj["parent_settings"] = ts_settings
    ctx.obj["parent_jobtype"] = "ts"

    if ctx.invoked_subcommand is None:
        check_charge_and_multiplicity(ts_settings)

        # Handle multiple molecules: create one job per molecule
        if len(molecules) > 1 and molecule_indices is not None:
            logger.info(f"Creating {len(molecules)} TS jobs")
            jobs = []
            for molecule, idx in zip(molecules, molecule_indices):
                # Create a copy to avoid side effects from mutation
                molecule = molecule.copy()
                molecule_label = f"{label}_idx{idx}"
                logger.info(
                    f"Running TS search for molecule {idx}: {molecule} with label {molecule_label}"
                )

                # Apply frozen atoms if specified
                if freeze_atoms is not None:
                    frozen_atoms_list = get_list_from_string_range(
                        freeze_atoms
                    )
                    logger.debug(f"Freezing atoms: {frozen_atoms_list}")
                    molecule.frozen_atoms = (
                        convert_list_to_gaussian_frozen_list(
                            frozen_atoms_list, molecule
                        )
                    )
                else:
                    logger.debug("No atoms will be frozen during TS search")

                job = GaussianTSJob(
                    molecule=molecule,
                    settings=ts_settings,
                    label=molecule_label,
                    jobrunner=jobrunner,
                    skip_completed=skip_completed,
                    **kwargs,
                )
                jobs.append(job)
            return jobs
        else:
            # Single molecule case
            molecule = molecules[-1]
            molecule = molecule.copy()

            if freeze_atoms is not None:
                frozen_atoms_list = get_list_from_string_range(freeze_atoms)
                logger.debug(f"Freezing atoms: {frozen_atoms_list}")
                molecule.frozen_atoms = convert_list_to_gaussian_frozen_list(
                    frozen_atoms_list, molecule
                )
            else:
                logger.debug("No atoms will be frozen during TS search")

            return GaussianTSJob(
                molecule=molecule,
                settings=ts_settings,
                label=label,
                jobrunner=jobrunner,
                skip_completed=skip_completed,
                **kwargs,
            )


create_qmmm_subcommand(ts)
