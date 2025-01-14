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
    help="Type of predictors used for IRC. Examples include[HPC, EulerPC, LQA, DVV, Euler].",
)
@click.option(
    "-rc",
    "--recorrect",
    type=click.Choice(["Never", "Always", "Test"], case_sensitive=False),
    default=None,
    help='Recorrection step of HPC and EulerPC IRCs. options are: ["Never", "Always", "Test"].',
)
@click.option(
    "-rs",
    "--recalc-step",
    type=int,
    default=6,
    help="Compute the Hessian analytically every N predictor steps or every |N| corrector steps if N<0. ",
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
    # get settings from project
    project_settings = ctx.obj["project_settings"]
    irc_project_settings = project_settings.irc_settings()

    # get settings from cli for irc
    from chemsmart.jobs.gaussian.settings import GaussianIRCJobSettings

    irc_settings = GaussianIRCJobSettings(
        predictor=predictor,
        recorrect=recorrect,
        recalc_step=recalc_step,
        maxpoints=maxpoints,
        maxcycles=maxcycles,
        stepsize=stepsize,
        flat_irc=flat_irc,
        **kwargs,
    )
    irc_keywords = [
        "predictor",
        "recorrect",
        "recalc_step",
        "direction",
        "maxpoints",
        "maxcycles",
        "stepsize",
        "flat_irc",
    ]

    # merge project settings with user settings from cli
    irc_settings = irc_project_settings.merge(
        irc_settings, keywords=irc_keywords
    )

    # job setting from filename or default, with updates from user in cli specified in keywords
    # e.g., `sub.py gaussian -c <user_charge> -m <user_multiplicity>`
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]

    # merge project irc settings with job settings from cli keywords from cli.gaussian.py subcommands
    irc_settings = irc_settings.merge(job_settings, keywords=keywords)
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

    from chemsmart.jobs.gaussian import GaussianIRCJob

    return GaussianIRCJob(
        molecule=molecule,
        settings=irc_settings,
        label=label,
        skip_completed=skip_completed,
        **kwargs,
    )
