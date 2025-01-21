import logging

import click

from chemsmart.cli.orca.orca import orca
from chemsmart.cli.job import click_job_options
from chemsmart.utils.cli import MyCommand
from chemsmart.utils.utils import check_charge_and_multiplicity

logger = logging.getLogger(__name__)


@orca.command("irc",cls=MyCommand)
@click_job_options
@click.option(
    "--maxiter",
    type=int,
    default=None,
    help="Maximum number of iterations.",
)
@click.option(
    "-p", "--printlevel", type=int, default=None, help="Print level."
)
@click.option(
    "-d",
    "--direction",
    type=click.Choice(
        ["both", "forward", "backward", "down"], case_sensitive=False
    ),
    default=None,
    help="IRC drirection. Available options: both, forward, backward, down.",
)
@click.option(
    "-i",
    "--inithess",
    type=click.Choice(
        ["read", "calc_anfreq", "calc_numfreq"], case_sensitive=False
    ),
    default=None,
    help="Initial Hessian. Available options: read, calc_anfreq, calc_numfreq.",
)
@click.option(
    "-f",
    "--hess-filename",
    type=str,
    default=None,
    help="Filename of initial Hessian.",
)
@click.option(
    "-m",
    "--hessmode",
    type=int,
    default=None,
    help="Hessian mode used for the initial displacement.Default 0.",
)
@click.option(
    "-M",
    "--monitor-internals/--no-monitor-internals",
    type=bool,
    default=False,
    help="Monitor internals to print out up to three internal coordinates",
)
@click.option(
    "--init-displ",
    type=click.Choice(["DE", "length"], case_sensitive=False),
    default=None,
    help="Initial displacement. Available options: DE, length. DE for energy difference, length for step size.",
)
@click.option(
    "--scale-init-displ",
    type=float,
    default=None,
    help="Step size for initial displacement from TS. Default 0.1 a.u.",
)
@click.option(
    "--de-init-displ",
    type=float,
    default=None,
    help="Energy difference for initial displacement based on provided Hessian (Default: 2 mEh)",
)
@click.option(
    "--follow-coordtype",
    type=str,
    default=None,
    help="Follow coordinate type. Default cartesian. The only option.",
)
@click.option(
    "--scale-displ-sd",
    type=float,
    default=None,
    help="Scaling factor for scaling the 1st SD step.Default to 0.15.",
)
@click.option(
    "--adapt-scale-displ/--no-adapt-scale-displ",
    type=bool,
    default=False,
    help="Modify Scale_Displ_SD when the step size becomes smaller or larger.",
)
@click.option(
    "--sd-parabolicfit/--no-sd-parabolicfit",
    type=bool,
    default=False,
    help="Do a parabolic fit for finding an optimal SD step length.",
)
@click.option(
    "--interpolate-only/--no-interpolate-only",
    type=bool,
    default=False,
    help="Only allow interpolation for parabolic fit, not extrapolation.",
)
@click.option(
    "--do-sd-corr/--no-do-sd-corr",
    type=bool,
    default=False,
    help="Do SD correction to 1st step.",
)
@click.option(
    "--scale-displ-sd-corr",
    type=float,
    default=None,
    help="Scaling factor for scaling the correction step to the SD step. "
    "It is multiplied by the length of the final 1st SD step.",
)
@click.option(
    "--sd-corr-parabolicfit/--no-sd-corr-parabolicfit",
    type=bool,
    default=False,
    help="Do a parabolic fit for finding an optimal correction step length.",
)
@click.option(
    "--tolrmsg",
    type=float,
    default=None,
    help="Tolerance for RMS gradient (a.u.). Default 5.e-4",
)
@click.option(
    "--tolmaxg",
    type=float,
    default=None,
    help="Tolerance for maximum gradient (a.u.). Default 2.e-3",
)
@click.option(
    "-I",
    "--internal-modred",
    default=None,
    help="Internal modred. Up to three internal coordinates can be defined and values printed.",
)
@click.pass_context
def irc(
    ctx,
    maxiter,
    printlevel,
    direction,
    inithess,
    hess_filename,
    hessmode,
    init_displ,
    scale_init_displ,
    de_init_displ,
    follow_coordtype,
    scale_displ_sd,
    adapt_scale_displ,
    sd_parabolicfit,
    interpolate_only,
    do_sd_corr,
    scale_displ_sd_corr,
    sd_corr_parabolicfit,
    tolrmsg,
    tolmaxg,
    internal_modred,
    skip_completed,
    **kwargs,
):
    # get settings from project
    project_settings = ctx.obj["project_settings"]
    irc_project_settings = project_settings.irc_settings()

    # job setting from filename or default, with updates from user in cli specified in keywords
    # e.g., `sub.py orca -c <user_charge> -m <user_multiplicity>`
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]

    # merge project irc settings with job settings from cli keywords from cli.orca.py subcommands
    irc_settings = irc_project_settings.merge(job_settings, keywords=keywords)

    # update irc_settings if any attribute is specified in cli options
    # suppose project has a non None value, and user does not specify a value (None),
    # then the project value should be used and unmodified, ie, should not be merged.
    # update value only if user specifies a value for the attribute:
    if maxiter is not None:
        irc_settings.maxiter = maxiter
    if printlevel is not None:
        irc_settings.printlevel = printlevel
    if direction is not None:
        irc_settings.direction = direction
    if inithess is not None:
        irc_settings.inithess = inithess
    if hess_filename is not None:
        irc_settings.hess_filename = hess_filename
    if hessmode is not None:
        irc_settings.hessmode = hessmode
    if init_displ is not None:
        irc_settings.init_displ = init_displ
    if scale_init_displ is not None:
        irc_settings.scale_init_displ = scale_init_displ
    if de_init_displ is not None:
        irc_settings.de_init_displ = de_init_displ
    if follow_coordtype is not None:
        irc_settings.follow_coordtype = follow_coordtype
    if scale_displ_sd is not None:
        irc_settings.scale_displ_sd = scale_displ_sd
    if adapt_scale_displ is not None:
        irc_settings.adapt_scale_displ = adapt_scale_displ
    if sd_parabolicfit is not None:
        irc_settings.sd_parabolicfit = sd_parabolicfit
    if interpolate_only is not None:
        irc_settings.interpolate_only = interpolate_only
    if do_sd_corr is not None:
        irc_settings.do_sd_corr = do_sd_corr
    if scale_displ_sd_corr is not None:
        irc_settings.scale_displ_sd_corr = scale_displ_sd_corr
    if sd_corr_parabolicfit is not None:
        irc_settings.sd_corr_parabolicfit = sd_corr_parabolicfit
    if tolrmsg is not None:
        irc_settings.tolrmsg = tolrmsg
    if tolmaxg is not None:
        irc_settings.tolmaxg = tolmaxg
    if internal_modred is not None:
        modred_info = eval(internal_modred)
        irc_settings.internal_modred = modred_info

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

    from chemsmart.jobs.orca.irc import ORCAIRCJob

    return ORCAIRCJob(
        molecule=molecule,
        settings=irc_settings,
        label=label,
        skip_completed=skip_completed,
        **kwargs,
    )
