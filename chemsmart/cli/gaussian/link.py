import logging

import click

from chemsmart.cli.gaussian.gaussian import (
    click_gaussian_jobtype_options,
    click_gaussian_solvent_options,
    gaussian,
)
from chemsmart.cli.job import click_job_options
from chemsmart.utils.cli import (
    MyCommand,
    get_setting_from_jobtype_for_gaussian,
)
from chemsmart.utils.utils import check_charge_and_multiplicity

logger = logging.getLogger(__name__)


@gaussian.command("link", cls=MyCommand)
@click_job_options
@click_gaussian_jobtype_options
@click_gaussian_solvent_options
@click.option(
    "-st",
    "--stable",
    type=str,
    default="opt",
    help='Gaussian stability test. See https://gaussian.com/stable/ for options. Defaults to "stable=opt".',
)
@click.option(
    "-g",
    "--guess",
    type=str,
    default="mix",
    help='Gaussian guess options. See https://gaussian.com/guess/ for options. Defaults to "guess=mix".',
)
@click.option(
    "-r", "--route", type=str, default=None, help="Route for link section."
)
@click.pass_context
def link(
    ctx,
    stable,
    guess,
    jobtype,
    coordinates,
    step_size,
    num_steps,
    remove_solvent,
    solvent_model,
    solvent_id,
    solvent_options,
    route,
    **kwargs,
):
    """CLI for running Gaussian link jobs."""

    # get jobrunner for running Gaussian link jobs
    jobrunner = ctx.obj["jobrunner"]

    # get settings from project
    from chemsmart.jobs.gaussian.settings import GaussianLinkJobSettings

    project_settings = ctx.obj["project_settings"]
    link_settings = get_setting_from_jobtype_for_gaussian(
        project_settings, jobtype, coordinates, step_size, num_steps
    )

    # job setting from filename or default, with updates from user in cli specified in keywords
    # e.g., `sub.py gaussian -c <user_charge> -m <user_multiplicity>`
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]

    # merge project settings with job settings from cli keywords from cli.gaussian.py subcommands
    link_settings = link_settings.merge(job_settings, keywords=keywords)
    check_charge_and_multiplicity(link_settings)

    # convert from GaussianJobSettings instance to GaussianLinkJobSettings instance
    link_settings = GaussianLinkJobSettings(**link_settings.__dict__)

    # populate GaussianLinkJobSettings
    link_settings.stable = stable
    link_settings.guess = guess
    link_settings.remove_solvent = remove_solvent
    if solvent_model is not None:
        link_settings.solvent_model = solvent_model
    if solvent_id is not None:
        link_settings.solvent_id = solvent_id
    if solvent_options is not None:
        link_settings.additional_solvent_options = solvent_options

    if route is not None:
        link_settings.link_route = route

    # get molecule
    molecules = ctx.obj["molecules"]
    molecule = molecules[-1]

    # get label for the job
    label = ctx.obj["label"]

    if jobtype is None:
        label = label
    else:
        label = f"{label}_{jobtype}_link"

    logger.debug(f"Label for job: {label}")

    # automatically use unrestricted dft if link job
    if not link_settings.functional.lower().startswith("u"):
        link_settings.functional = "u" + link_settings.functional

    logger.info(
        f"Link job {jobtype} settings from project: {link_settings.__dict__}"
    )

    from chemsmart.jobs.gaussian.link import GaussianLinkJob

    return GaussianLinkJob(
        molecule=molecule,
        settings=link_settings,
        label=label,
        jobrunner=jobrunner,
        **kwargs,
    )
