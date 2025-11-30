import logging

import click

from chemsmart.cli.gaussian.gaussian import gaussian
from chemsmart.cli.job import click_job_options
from chemsmart.utils.cli import MyCommand
from chemsmart.utils.utils import check_charge_and_multiplicity

logger = logging.getLogger(__name__)


@gaussian.command("nci", cls=MyCommand)
@click_job_options
@click.pass_context
def nci(ctx, **kwargs):
    """CLI subcommand for Gaussian NCI jobs."""

    # get jobrunner for running Gaussian NCI jobs
    jobrunner = ctx.obj["jobrunner"]

    # get settings from project
    project_settings = ctx.obj["project_settings"]
    nci_settings = project_settings.nci_settings()

    # job setting from filename or default, with updates from user in cli
    # specified in keywords
    # e.g., `sub.py gaussian -c <user_charge> -m <user_multiplicity>`
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]

    # merge project settings with job settings from cli keywords from
    # cli.gaussian.py subcommands
    nci_settings = nci_settings.merge(job_settings, keywords=keywords)
    check_charge_and_multiplicity(nci_settings)

    # get molecule
    molecules = ctx.obj["molecules"]
    molecule = molecules[-1]

    # get label for the job
    label = ctx.obj["label"]
    logger.debug(f"Label for job: {label}")

    logger.info(f"IRC job settings from project: {nci_settings.__dict__}")

    from chemsmart.jobs.gaussian.nci import GaussianNCIJob

    return GaussianNCIJob(
        molecule=molecule,
        settings=nci_settings,
        label=label,
        jobrunner=jobrunner,
        **kwargs,
    )
