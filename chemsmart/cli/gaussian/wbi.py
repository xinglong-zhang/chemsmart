import logging
import click

from chemsmart.cli.job import click_job_options
from chemsmart.cli.gaussian import gaussian
from chemsmart.utils.cli import MyCommand
from chemsmart.utils.utils import check_charge_and_multiplicity

logger = logging.getLogger(__name__)


@gaussian.command("wbi", cls=MyCommand)
@click_job_options
@click.pass_context
def wbi(ctx, **kwargs):

    # get settings from project
    project_settings = ctx.obj["project_settings"]
    wbi_settings = project_settings.wbi_settings()

    # job setting from filename or default, with updates from user in cli specified in keywords
    # e.g., `sub.py gaussian -c <user_charge> -m <user_multiplicity>`
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]

    # merge project settings with job settings from cli keywords from cli.gaussian.py subcommands
    wbi_settings = wbi_settings.merge(job_settings, keywords=keywords)
    check_charge_and_multiplicity(wbi_settings)

    # get molecule
    molecules = ctx.obj["molecules"]
    molecule = molecules[-1]

    # get label for the job
    label = ctx.obj["label"]
    logger.debug(f"Label for job: {label}")

    logger.info(f"Weiberg bond index job settings from project: {wbi_settings.__dict__}")

    from chemsmart.jobs.gaussian.wbi import GaussianWBIJob

    return GaussianWBIJob(
        molecule=molecule, settings=wbi_settings, label=label, **kwargs
    )
