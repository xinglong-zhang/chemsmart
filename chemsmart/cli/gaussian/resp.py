import logging

import click

from chemsmart.cli.gaussian.gaussian import gaussian
from chemsmart.cli.job import click_job_options
from chemsmart.utils.cli import MyCommand
from chemsmart.utils.utils import check_charge_and_multiplicity

logger = logging.getLogger(__name__)


@gaussian.command("resp", cls=MyCommand)
@click_job_options
@click.pass_context
def resp(ctx, **kwargs):
    """CLI for running Gaussian RESP jobs."""

    # get jobrunner for running Gaussian RESP jobs
    jobrunner = ctx.obj["jobrunner"]

    # get settings from project
    project_settings = ctx.obj["project_settings"]
    resp_settings = project_settings.sp_settings()

    # job setting from filename or default, with updates from user in cli
    # specified in keywords
    # e.g., `sub.py gaussian -c <user_charge> -m <user_multiplicity>`
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]

    # merge project settings with job settings from cli keywords from
    # cli.gaussian.py subcommands
    resp_settings = resp_settings.merge(job_settings, keywords=keywords)
    check_charge_and_multiplicity(resp_settings)

    # get molecule
    molecules = ctx.obj["molecules"]
    molecule = molecules[-1]

    # get label for the job
    label = ctx.obj["label"]
    logger.debug(f"Label for job: {label}")

    # fixed route for resp job
    resp_settings.route_to_be_written = (
        "HF/6-31+G(d) SCF=Tight Pop=MK IOp(6/33=2,6/41=10,6/42=17,6/50=1)"
    )

    logger.info(f"RESP settings from project: {resp_settings.__dict__}")

    from chemsmart.jobs.gaussian.resp import GaussianRESPJob

    return GaussianRESPJob(
        molecule=molecule,
        settings=resp_settings,
        label=label,
        jobrunner=jobrunner,
        **kwargs,
    )
