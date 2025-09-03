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
    "-fl/",
    "--flat-irc/--no-flat-irc",
    type=bool,
    default=False,
    help="whether to run flat irc or not",
)
@click.option(
    "-pt",
    "--predictor",
    type=click.Choice(
        ["LQA", "HPC", "EulerPC", "DVV", "Euler"], case_sensitive=False
    ),
    default=None,
    help="Type of predictors used for IRC. Examples include[HPC, EulerPC, "
    "LQA, DVV, Euler].",
)
@click.option(
    "-rc",
    "--recorrect",
    type=click.Choice(["Never", "Always", "Test"], case_sensitive=False),
    default=None,
    help='Recorrection step of HPC and EulerPC IRCs. options are: '
    '["Never", "Always", "Test"].',
)
@click.option(
    "-rs",
    "--recalc-step",
    type=int,
    default=6,
    help="Compute the Hessian analytically every N predictor steps or every "
    "|N| corrector steps if N<0. ",
)
@click.option(
    "-p",
    "--maxpoints",
    type=int,
    default=512,
    help="Number of points along reaction path to examine.",
)
@click.option(
    "-c",
    "--maxcycles",
    type=int,
    default=128,
    help="Maximum number of steps along IRC to run.",
)
@click.option(
    "-s",
    "--stepsize",
    type=int,
    default=20,
    help="Step size along reaction path, in units of 0.01 Bohr.",
)
@click.pass_context
def irc(
    ctx,
    flat_irc,
    predictor,
    recorrect,
    recalc_step,
    maxpoints,
    maxcycles,
    stepsize,
    skip_completed,
    **kwargs,
):
    """CLI for running Gaussian IRC jobs."""

    # get jobrunner for running Gaussian IRC jobs
    jobrunner = ctx.obj["jobrunner"]

    # get settings from project
    project_settings = ctx.obj["project_settings"]
    irc_project_settings = project_settings.irc_settings()

    # job setting from filename or default, with updates from user in cli
    # specified in keywords
    # e.g., `sub.py gaussian -c <user_charge> -m <user_multiplicity>`
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]

    # merge project irc settings with job settings from cli keywords from
    # cli.orca.py subcommands
    irc_settings = irc_project_settings.merge(job_settings, keywords=keywords)

    # update irc_settings if any attribute is specified in cli options
    # suppose project has a non None value, and user does not specify a
    # value (None),
    # then the project value should be used and unmodified, ie, should not be
    # merged.
    # update value only if user specifies a value for the attribute:
    if predictor is not None:
        irc_settings.predictor = predictor
    if recorrect is not None:
        irc_settings.recorrect = recorrect
    if recalc_step is not None:
        irc_settings.recalc_step = recalc_step
    if maxpoints is not None:
        irc_settings.maxpoints = maxpoints
    if maxcycles is not None:
        irc_settings.maxcycles = maxcycles
    if stepsize is not None:
        irc_settings.stepsize = stepsize
    if flat_irc is not None:
        irc_settings.flat_irc = flat_irc

    check_charge_and_multiplicity(irc_settings)

    # get molecule
    molecules = ctx.obj[
        "molecules"
    ]  # use all molecules as a list for crest jobs
    molecule = molecules[-1]  # get last molecule from list of molecules

    # get label for the job
    label = ctx.obj["label"]
    logger.debug(f"Label for job: {label}")

    logger.info(f"IRC job settings from project: {irc_settings.__dict__}")

    from chemsmart.jobs.gaussian.irc import GaussianIRCJob

    return GaussianIRCJob(
        molecule=molecule,
        settings=irc_settings,
        label=label,
        jobrunner=jobrunner,
        skip_completed=skip_completed,
        **kwargs,
    )
