import logging

import click

from pyatoms.cli.job import click_job_options
from pyatoms.jobs.orca.job import ORCAJob
from pyatoms.utils.cli import MyCommand

logger = logging.getLogger(__name__)


class ORCAIRCJob(ORCAJob):
    TYPE = "orcairc"

    def __init__(self, folder, atoms, settings, **kwargs):
        super().__init__(
            folder=folder, atoms=atoms, settings=settings, **kwargs
        )

    @staticmethod
    @click.command("irc", cls=MyCommand)
    @click.option(
        "-n",
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
        "-mi",
        "--monitor-internals",
        type=click.Choice(["true", "false"], case_sensitive=False),
        default=None,
        help="Monitor internals to print out up to three internal coordinates",
    )
    @click.option(
        "-c",
        "--coordinates",
        default=None,
        help="Internal modred. Up to three internal coordinates can be defined and values printed.",
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
        "--adapt-scale-displ",
        type=click.Choice(["true", "false"], case_sensitive=False),
        default=None,
        help="Modify Scale_Displ_SD when the step size becomes smaller or larger.",
    )
    @click.option(
        "--sd-parabolicfit",
        type=click.Choice(["true", "false"], case_sensitive=False),
        default=None,
        help="Do a parabolic fit for finding an optimal SD step length.",
    )
    @click.option(
        "--interpolate-only",
        type=click.Choice(["true", "false"], case_sensitive=False),
        default=None,
        help="Only allow interpolation for parabolic fit, not extrapolation.",
    )
    @click.option(
        "--do-sd-corr",
        type=click.Choice(["true", "false"], case_sensitive=False),
        default=None,
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
        "--sd-corr-parabolicfit",
        type=click.Choice(["true", "false"], case_sensitive=False),
        default=None,
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
    @click_job_options
    @click.pass_context
    def cli_creation(
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
        monitor_internals,
        coordinates,
        **kwargs,
    ):
        folder = ctx.obj["folder"]

        # get settings from project
        project_settings = ctx.obj["project_settings"]
        irc_settings = project_settings.irc_settings()

        # job setting from filename or default, with updates from user in cli specified in keywords
        # e.g., `sub.py orca -c <user_charge> -m <user_multiplicity>`
        job_settings = ctx.obj["job_settings"]
        keywords = ctx.obj["keywords"]

        # merge project opt settings with job settings from cli keywords from cli.orca.py subcommands
        irc_settings = irc_settings.merge(job_settings, keywords=keywords)

        # update settings from cli
        irc_settings.maxiter = maxiter
        irc_settings.printlevel = printlevel
        irc_settings.direction = direction
        irc_settings.inithess = inithess
        irc_settings.hess_filename = hess_filename
        irc_settings.hessmode = hessmode
        irc_settings.init_displ = init_displ
        irc_settings.scale_init_displ = scale_init_displ
        irc_settings.de_init_displ = de_init_displ
        irc_settings.follow_coordtype = follow_coordtype
        irc_settings.scale_displ_sd = scale_displ_sd
        irc_settings.adapt_scale_displ = adapt_scale_displ
        irc_settings.sd_parabolicfit = sd_parabolicfit
        irc_settings.interpolate_only = interpolate_only
        irc_settings.do_sd_corr = do_sd_corr
        irc_settings.scale_displ_sd_corr = scale_displ_sd_corr
        irc_settings.sd_corr_parabolicfit = sd_corr_parabolicfit
        irc_settings.tolrmsg = tolrmsg
        irc_settings.tolmaxg = tolmaxg
        irc_settings.monitor_internals = monitor_internals
        if coordinates is not None:
            modred_info = eval(coordinates)
            irc_settings.internal_modred = modred_info

        # get molecule
        atoms = ctx.obj["molecule"]

        # get label for the job
        label = ctx.obj["label"]

        # update irc settings
        # irc_settings.constrained_atoms = frozen_atoms
        # irc_settings.invert_constraints = invert_atoms

        logger.info(f"IRC settings from project: {irc_settings.__dict__}")

        if len(atoms) < 1:
            raise ValueError(
                "Atoms object must contain at least one image for optimization. "
            )

        if len(atoms) == 1:
            atoms = atoms[-1]
            return ORCAIRCJob(
                folder=folder,
                atoms=atoms,
                settings=irc_settings,
                label=label,
                **kwargs,
            )

        return [
            ORCAIRCJob(
                folder=folder,
                atoms=image,
                settings=irc_settings,
                label=f"{label}_opt_c{idx + 1}",
                **kwargs,
            )
            for idx, image in enumerate(atoms)
        ]
