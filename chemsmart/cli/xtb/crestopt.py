import logging

import click

from chemsmart.cli.xtb.xtb import xtb
from chemsmart.cli.job import click_job_options
from chemsmart.utils.cli import MyCommand
from chemsmart.utils.utils import check_charge_and_multiplicity

logger = logging.getLogger(__name__)


@xtb.command("crestopt", cls=MyCommand)
@click_job_options
@click.option(
    "-n",
    "--num-confs-to-opt",
    type=int,
    default=None,
    help="Number of conformers to optimize.",
)
@click.pass_context
def crestopt(ctx, num_confs_to_opt, skip_completed, **kwargs):
    """CLI for running xTB CREST optimization jobs."""

    # get jobrunner for running xTB crest optimization jobs
    jobrunner = ctx.obj["jobrunner"]

    # get settings from project
    project_settings = ctx.obj["project_settings"]
    opt_settings = project_settings.opt_settings()

    # job setting from filename or default, with updates from user in cli specified in keywords
    # e.g., `sub.py xtb -c <user_charge> -m <user_multiplicity>`
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]

    # merge project opt settings with job settings from cli keywords from cli.xtb.py subcommands
    opt_settings = opt_settings.merge(job_settings, keywords=keywords)

    check_charge_and_multiplicity(opt_settings)

    # get molecule
    molecules = ctx.obj[
        "molecules"
    ]  # use all molecules as a list for crestopt

    # get label for the job
    label = ctx.obj["label"]

    logger.info(f"Crest opt settings from project: {opt_settings.__dict__}")

    from chemsmart.jobs.xtb.crestopt import XTBCrestOptJob

    return XTBCrestOptJob(
        molecules=molecules,
        settings=opt_settings,
        label=label,
        jobrunner=jobrunner,
        num_confs_to_opt=num_confs_to_opt,
        skip_completed=skip_completed,
        **kwargs,
    )