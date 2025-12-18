import logging

import click

from chemsmart.cli.gaussian.gaussian import (
    click_gaussian_jobtype_options,
    gaussian,
)
from chemsmart.cli.job import click_job_options
from chemsmart.utils.cli import (
    MyCommand,
    get_setting_from_jobtype_for_gaussian,
)
from chemsmart.utils.utils import check_charge_and_multiplicity

logger = logging.getLogger(__name__)


@gaussian.command("modred", cls=MyCommand)
@click_job_options
@click_gaussian_jobtype_options
@click.pass_context
def modred(ctx, jobtype, coordinates, step_size, num_steps, **kwargs):
    """CLI subcommand for running Gaussian modred jobs."""

    # get jobrunner for running Gaussian modred jobs
    jobrunner = ctx.obj["jobrunner"]

    if jobtype is None:
        jobtype = "modred"

    # get settings from project
    project_settings = ctx.obj["project_settings"]
    modred_settings = get_setting_from_jobtype_for_gaussian(
        project_settings, jobtype, coordinates, step_size, num_steps
    )

    # job setting from filename or default, with updates from user in cli
    # specified in keywords
    # e.g., `sub.py gaussian -c <user_charge> -m <user_multiplicity>`
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]

    # merge project settings with job settings from cli keywords from
    # cli.gaussian.py subcommands
    modred_settings = modred_settings.merge(job_settings, keywords=keywords)
    check_charge_and_multiplicity(modred_settings)

    # get molecules
    molecules = ctx.obj["molecules"]

    # get label for the job
    label = ctx.obj["label"]
    logger.debug(f"Label for job: {label}")

    logger.info(f"Modred settings from project: {modred_settings.__dict__}")

    from chemsmart.jobs.gaussian.modred import GaussianModredJob

    # Get the original molecule indices from context
    molecule_indices = ctx.obj.get(
        "molecule_indices", list(range(1, len(molecules) + 1))
    )

    # Handle multiple molecules: create one job per molecule
    if len(molecules) > 1 and molecule_indices is not None:
        logger.info(f"Creating {len(molecules)} modred jobs")
        jobs = []
        for molecule, idx in zip(molecules, molecule_indices):
            molecule_label = f"{label}_idx{idx}"
            logger.info(
                f"Running modred for molecule {idx}: {molecule} with label {molecule_label}"
            )

            job = GaussianModredJob(
                molecule=molecule,
                settings=modred_settings,
                label=molecule_label,
                jobrunner=jobrunner,
                **kwargs,
            )
            jobs.append(job)
        return jobs
    else:
        # Single molecule case
        molecule = molecules[-1]
        return GaussianModredJob(
            molecule=molecule,
            settings=modred_settings,
            label=label,
            jobrunner=jobrunner,
            **kwargs,
        )
