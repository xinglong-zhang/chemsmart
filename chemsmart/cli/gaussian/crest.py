import click
import logging

from chemsmart.utils.utils import check_charge_and_multiplicity
from chemsmart.cli.job import click_job_options
from chemsmart.utils.cli import MyCommand

logger = logging.getLogger(__name__)

@click.command(cls=MyCommand)
@click_job_options
@click.option(
    "-t",
    "--type",
    type=str,
    required=True,
    help="Type of job to run for crest.",
)
@click.pass_context
def crest(ctx, type, skip_completed, **kwargs):
    # get settings from project
    project_settings = ctx.obj["project_settings"]
    if type.lower() == "opt":
        settings = project_settings.opt_settings()
    elif type.lower() == "ts":
        settings = project_settings.ts_settings()
    elif type.lower() == "modred":
        settings = project_settings.modred_settings()
    elif type.lower() == "irc":
        settings = project_settings.irc_settings()
    elif type.lower() == "scan":
        settings = project_settings.scan_settings()
    elif type.lower() == "sp":
        settings = project_settings.sp_settings()
    elif type.lower() == "td":
        settings = project_settings.td_settings()
    elif type.lower() == "wbi":
        settings = project_settings.wbi_settings()
    elif type.lower() == "nci":
        settings = project_settings.nci_settings()
    else:
        raise ValueError(f"Invalid crest type: {type}")

    # job setting from filename or default, with updates from user in cli specified in keywords
    # e.g., `sub.py gaussian -c <user_charge> -m <user_multiplicity>`
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]

    # merge project opt settings with job settings from cli keywords from cli.gaussian.py subcommands
    crest_settings = settings.merge(job_settings, keywords=keywords)

    check_charge_and_multiplicity(crest_settings)

    # get molecule
    molecules = ctx.obj["molecules"]  # use all molecules as a list for crest jobs

    # get label for the job
    label = ctx.obj["label"]

    logger.info(f"Crest {type} settings from project: {crest_settings.__dict__}")

    from chemsmart.jobs.gaussian.job import GaussianJob

    return GaussianJob.from_jobtype(
        jobtype=type,
        molecules=molecules,
        settings=crest_settings,
        label=label,
        skip_completed=skip_completed,
        **kwargs,
    )
