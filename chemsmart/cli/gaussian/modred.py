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
    if jobtype is None:
        jobtype = "modred"

    # get settings from project
    project_settings = ctx.obj["project_settings"]
    modred_settings = get_setting_from_jobtype_for_gaussian(
        project_settings, jobtype, coordinates, step_size, num_steps
    )

    # job setting from filename or default, with updates from user in cli specified in keywords
    # e.g., `sub.py gaussian -c <user_charge> -m <user_multiplicity>`
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]

    # merge project settings with job settings from cli keywords from cli.gaussian.py subcommands
    modred_settings = modred_settings.merge(job_settings, keywords=keywords)
    check_charge_and_multiplicity(modred_settings)

    # get molecule
    molecules = ctx.obj["molecules"]
    molecule = molecules[-1]

    # get label for the job
    label = ctx.obj["label"]
    logger.debug(f"Label for job: {label}")

    logger.info(f"Modred settings from project: {modred_settings.__dict__}")

    from chemsmart.jobs.gaussian.modred import GaussianModredJob

    return GaussianModredJob(
        molecule=molecule, settings=modred_settings, label=label, **kwargs
    )
