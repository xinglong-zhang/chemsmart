import logging

import click

from chemsmart.cli.gaussian.gaussian import (
    click_gaussian_solvent_options,
    gaussian,
)
from chemsmart.cli.gaussian.qmmm import create_qmmm_subcommand
from chemsmart.cli.job import click_job_options
from chemsmart.utils.cli import MyGroup
from chemsmart.utils.utils import check_charge_and_multiplicity

logger = logging.getLogger(__name__)


@gaussian.group("sp", cls=MyGroup, invoke_without_command=True)
@click_job_options
@click_gaussian_solvent_options
@click.pass_context
def sp(
    ctx,
    remove_solvent,
    solvent_model,
    solvent_id,
    solvent_options,
    skip_completed,
    **kwargs,
):
    """CLI subcommand for running Gaussian single point calculation."""

    # get jobrunner for single point
    jobrunner = ctx.obj["jobrunner"]

    # get settings from project
    project_settings = ctx.obj["project_settings"]
    sp_settings = project_settings.sp_settings()

    # job setting from filename or default, with updates from user in cli
    # specified in keywords
    # e.g., `sub.py gaussian -c <user_charge> -m <user_multiplicity>`
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]

    # merge project settings with job settings from cli keywords from
    # cli.gaussian.py subcommands
    sp_settings = sp_settings.merge(job_settings, keywords=keywords)

    # get molecules
    molecules = ctx.obj["molecules"]

    # get label for the job
    label = ctx.obj["label"]

    # cli-supplied solvent model and solvent id
    sp_settings.modify_solvent(
        remove_solvent=remove_solvent,
        solvent_model=solvent_model,
        solvent_id=solvent_id,
    )

    if solvent_options is not None:
        sp_settings.additional_solvent_options = solvent_options

    logger.info(
        f"Single point job settings from project: {sp_settings.__dict__}"
    )

    from chemsmart.jobs.gaussian.singlepoint import GaussianSinglePointJob
    from chemsmart.utils.cli import create_sp_label

    # Get the original molecule indices from context
    molecule_indices = ctx.obj["molecule_indices"]

    # Handle multiple molecules: create one job per molecule
    # Store parent context for potential qmmm subcommand
    ctx.obj["parent_skip_completed"] = skip_completed
    ctx.obj["parent_freeze_atoms"] = None  # sp doesn't have freeze_atoms
    ctx.obj["parent_kwargs"] = kwargs
    ctx.obj["parent_settings"] = sp_settings
    ctx.obj["parent_jobtype"] = "sp"

    # If no subcommand invoked, run regular single point
    if ctx.invoked_subcommand is None:
        check_charge_and_multiplicity(sp_settings)

        if len(molecules) > 1 and molecule_indices is not None:
            logger.info(f"Creating {len(molecules)} single point jobs")
            jobs = []
            for molecule, idx in zip(molecules, molecule_indices):
                molecule_label = f"{label}_idx{idx}"
                final_label = create_sp_label(molecule_label, sp_settings)
                logger.info(
                    f"Running single point for molecule {idx}: {molecule} with label {final_label}"
                )

                job = GaussianSinglePointJob(
                    molecule=molecule,
                    settings=sp_settings,
                    label=final_label,
                    jobrunner=jobrunner,
                    **kwargs,
                )
                jobs.append(job)
            return jobs
        else:
            # Single molecule case
            molecule = molecules[-1]
            label = create_sp_label(label, sp_settings)
            logger.info(
                f"Running single point calculation on molecule: {molecule} with label: {label}"
            )

            return GaussianSinglePointJob(
                molecule=molecule,
                settings=sp_settings,
                label=label,
                jobrunner=jobrunner,
                **kwargs,
            )


# Register qmmm subcommand
create_qmmm_subcommand(sp)
