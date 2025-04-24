import logging

import click

from chemsmart.cli.gaussian.gaussian import gaussian
from chemsmart.cli.job import click_job_options
from chemsmart.utils.cli import MyCommand
from chemsmart.utils.utils import check_charge_and_multiplicity

logger = logging.getLogger(__name__)


@gaussian.command(cls=MyCommand)
@click_job_options
@click.option(
    "-i",
    "--fragment-indices",
    required=True,
    type=str,
    help="Indices of one fragment for DI-AS analysis.",
)
@click.option(
    "-n",
    "--every-n-points",
    type=int,
    default=3,
    help="Every nth points along the IRC file to prepare for DI-AS analysis.",
)
@click.option(
    "-s",
    "--solv/--no-solv",
    type=bool,
    default=False,
    help="Turn on/off solvent for DI-AS job calculations.",
)
@click.option(
    "-m",
    "--mode",
    type=click.Choice(["irc", "ts"], case_sensitive=False),
    default="irc",
    help="Mode of DI-AS analysis. Default is IRC.",
)
@click.option(
    "-c1",
    "--charge-of-fragment1",
    type=int,
    default=None,
    help="Charge of fragment 1.",
)
@click.option(
    "-m1",
    "--multiplicity-of-fragment1",
    type=int,
    default=None,
    help="Multiplicity of fragment 1.",
)
@click.option(
    "-c2",
    "--charge-of-fragment2",
    type=int,
    default=None,
    help="Charge of fragment 2.",
)
@click.option(
    "-m2",
    "--multiplicity-of-fragment2",
    type=int,
    default=None,
    help="Multiplicity of fragment 2.",
)
@click.pass_context
def dias(
    ctx,
    fragment_indices,
    every_n_points,
    skip_completed,
    solv=False,
    mode="irc",
    charge_of_fragment1=None,
    multiplicity_of_fragment1=None,
    charge_of_fragment2=None,
    multiplicity_of_fragment2=None,
    **kwargs,
):
    # get settings from project
    project_settings = ctx.obj["project_settings"]
    sp_settings = project_settings.sp_settings()

    # job setting from filename or default, with updates from user in cli specified in keywords
    # e.g., `sub.py gaussian -c <user_charge> -m <user_multiplicity>`
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]

    sp_settings = sp_settings.merge(job_settings, keywords=keywords)

    check_charge_and_multiplicity(sp_settings)

    # get molecule
    molecules = ctx.obj[
        "molecules"
    ]  # use all molecules as a list for crest jobs

    # get label for the job
    label = ctx.obj["label"]

    if not solv:
        sp_settings.solvent_model = None
        sp_settings.solvent_id = None

    logger.info(f"DI-AS settings from project: {sp_settings.__dict__}")

    from chemsmart.jobs.gaussian.dias import GaussianDIASJob

    return GaussianDIASJob(
        molecules=molecules,
        settings=sp_settings,
        label=label,
        fragment_indices=fragment_indices,
        every_n_points=every_n_points,
        skip_completed=skip_completed,
        **kwargs,
    )
