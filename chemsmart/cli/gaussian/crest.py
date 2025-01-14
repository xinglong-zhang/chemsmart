import click
import logging
from chemsmart.utils.utils import check_charge_and_multiplicity
from chemsmart.cli.job import click_job_options
from chemsmart.utils.cli import MyCommand
from chemsmart.cli.gaussian.gaussian import gaussian

logger = logging.getLogger(__name__)


@gaussian.command(cls=MyCommand)
@click_job_options
@click.option(
    "-j",
    "--jobtype",
    type=str,
    required=True,
    help="Type of job to run for crest.",
)
@click.option(
    "-n",
    "--num-confs-to-run",
    type=int,
    default=None,
    help="Number of conformers to optimize.",
)
@click.pass_context
def crest(ctx, jobtype, num_confs_to_run, skip_completed, **kwargs):
    # get settings from project
    project_settings = ctx.obj["project_settings"]

    if jobtype.lower() == "opt":
        settings = project_settings.opt_settings()
    elif jobtype.lower() == "ts":
        settings = project_settings.ts_settings()
    elif jobtype.lower() == "modred":
        settings = project_settings.modred_settings()
    elif jobtype.lower() == "irc":
        settings = project_settings.irc_settings()
    elif jobtype.lower() == "scan":
        settings = project_settings.scan_settings()
    elif jobtype.lower() == "sp":
        settings = project_settings.sp_settings()
    elif jobtype.lower() == "td":
        settings = project_settings.td_settings()
    elif jobtype.lower() == "wbi":
        settings = project_settings.wbi_settings()
    elif jobtype.lower() == "nci":
        settings = project_settings.nci_settings()

    # job setting from filename or default, with updates from user in cli specified in keywords
    # e.g., `sub.py gaussian -c <user_charge> -m <user_multiplicity>`
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]

    settings = settings.merge(job_settings, keywords=keywords)

    check_charge_and_multiplicity(settings)

    # get molecule
    molecules = ctx.obj[
        "molecules"
    ]  # use all molecules as a list for crest jobs

    # get label for the job
    label = ctx.obj["label"]

    logger.info(f"Crest {type} settings from project: {settings.__dict__}")

    from chemsmart.jobs.gaussian import GaussianCrestJob

    return GaussianCrestJob(
        molecules=molecules,
        settings=settings,
        label=label,
        num_confs_to_run=num_confs_to_run,
        skip_completed=skip_completed,
        **kwargs,
    )
