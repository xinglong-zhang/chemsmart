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


@gaussian.command(cls=MyCommand)
@click_job_options
@click_gaussian_jobtype_options
@click.option(
    "-N",
    "--num-confs-to-run",
    type=int,
    default=None,
    help="Number of conformers to optimize.",
)
@click.pass_context
def crest(
    ctx,
    jobtype,
    coordinates,
    step_size,
    num_steps,
    num_confs_to_run,
    skip_completed,
    **kwargs,
):
    """CLI for running Gaussian CREST jobs."""

    # get jobrunner for running Gaussian crest jobs
    jobrunner = ctx.obj["jobrunner"]

    # get settings from project
    project_settings = ctx.obj["project_settings"]
    crest_settings = get_setting_from_jobtype_for_gaussian(
        project_settings, jobtype, coordinates, step_size, num_steps
    )

    # job setting from filename or default, with updates from user in cli specified in keywords
    # e.g., `sub.py gaussian -c <user_charge> -m <user_multiplicity>`
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]

    crest_settings = crest_settings.merge(job_settings, keywords=keywords)

    check_charge_and_multiplicity(crest_settings)

    # get molecule
    molecules = ctx.obj[
        "molecules"
    ]  # use all molecules as a list for crest jobs

    # get label for the job
    label = ctx.obj["label"]
    label = f"{label}_{jobtype}"
    logger.debug(f"Label for job: {label}")

    logger.info(
        f"Crest {type} settings from project: {crest_settings.__dict__}"
    )

    from chemsmart.jobs.gaussian.crest import GaussianCrestJob

    return GaussianCrestJob(
        molecules=molecules,
        settings=crest_settings,
        label=label,
        jobrunner=jobrunner,
        num_confs_to_run=num_confs_to_run,
        skip_completed=skip_completed,
        **kwargs,
    )
