import logging

import click

from chemsmart.cli.gaussian.gaussian import click_gaussian_td_options, gaussian
from chemsmart.cli.job import click_job_options
from chemsmart.utils.cli import MyCommand
from chemsmart.utils.utils import check_charge_and_multiplicity

logger = logging.getLogger(__name__)


@gaussian.command("td", cls=MyCommand)
@click_job_options
@click_gaussian_td_options
@click.pass_context
def td(ctx, states, root, nstates, eqsolv, **kwargs):
    """Run Gaussian TDDFT jobs."""

    # get jobrunner for running Gaussian TDDFT jobs
    jobrunner = ctx.obj["jobrunner"]

    from chemsmart.jobs.gaussian.settings import GaussianTDDFTJobSettings

    # get settings from project
    project_settings = ctx.obj["project_settings"]
    td_settings = project_settings.td_settings()

    # job setting from filename or default, with updates from user in cli
    # specified in keywords
    # e.g., `sub.py gaussian -c <user_charge> -m <user_multiplicity>`
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]

    # merge project settings with job settings from cli keywords from
    # cli.gaussian.py subcommands
    td_settings = td_settings.merge(job_settings, keywords=keywords)
    check_charge_and_multiplicity(td_settings)

    # get molecule
    molecules = ctx.obj["molecules"]
    molecule = molecules[-1]

    # get label for the job
    label = ctx.obj["label"]
    logger.debug(f"Label for job: {label}")

    # convert from GaussianJobSettings instance to GaussianTDJobSettings
    # instance
    td_settings = GaussianTDDFTJobSettings(**td_settings.__dict__)

    # populate cli options
    td_settings.states = states
    td_settings.root = root
    td_settings.nstates = nstates
    td_settings.eqsolv = eqsolv

    logger.info(f"TDDFT job settings from project: {td_settings.__dict__}")

    from chemsmart.jobs.gaussian.tddft import GaussianTDDFTJob

    return GaussianTDDFTJob(
        molecule=molecule,
        settings=td_settings,
        label=label,
        jobrunner=jobrunner,
        **kwargs,
    )
