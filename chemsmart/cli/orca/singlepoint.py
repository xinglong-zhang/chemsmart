import logging

import click

from chemsmart.cli.job import click_job_options
from chemsmart.cli.orca.orca import orca
from chemsmart.utils.cli import MyCommand
from chemsmart.utils.utils import check_charge_and_multiplicity

logger = logging.getLogger(__name__)


@orca.command("sp", cls=MyCommand)
@click_job_options
@click.pass_context
def sp(ctx, **kwargs):
    # get settings from project
    project_settings = ctx.obj["project_settings"]
    sp_settings = project_settings.sp_settings()

    # job setting from filename or default, with updates from user in cli specified in keywords
    # e.g., `sub.py gaussian -c <user_charge> -m <user_multiplicity>`
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]

    # merge project sp settings with job settings from cli keywords from cli.orca.py subcommands
    sp_settings = sp_settings.merge(job_settings, keywords=keywords)

    check_charge_and_multiplicity(sp_settings)

    # get molecule
    molecules = ctx.obj["molecules"]
    molecule = molecules[-1]
    logger.info(f"Single point calculation on molecule: {molecule}.")

    # get label for the job
    label = ctx.obj["label"]

    logger.info(f"Opt job settings from project: {sp_settings.__dict__}")

    from chemsmart.jobs.orca.singlepoint import ORCASinglePointJob

    return ORCASinglePointJob(
        molecule=molecule,
        settings=sp_settings,
        label=label,
        **kwargs,
    )
