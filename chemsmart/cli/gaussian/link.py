import logging

import click

from chemsmart.cli.gaussian.gaussian import (
    click_gaussian_irc_options,
    click_gaussian_jobtype_options,
    click_gaussian_solvent_options,
    gaussian,
)
from chemsmart.cli.job import click_job_options
from chemsmart.utils.cli import (
    MyCommand,
    get_setting_from_jobtype_for_gaussian,
    update_irc_label,
)
from chemsmart.utils.utils import check_charge_and_multiplicity

logger = logging.getLogger(__name__)


@gaussian.command("link", cls=MyCommand)
@click_job_options
@click_gaussian_jobtype_options
@click_gaussian_solvent_options
@click_gaussian_irc_options
@click.option(
    "-st",
    "--stable",
    type=str,
    default="opt",
    help="Gaussian stability test. See https://gaussian.com/stable/ for "
    'options. Defaults to "stable=opt".',
)
@click.option(
    "-g",
    "--guess",
    type=str,
    default="mix",
    help="Gaussian guess options. See https://gaussian.com/guess/ for "
    'options. Defaults to "guess=mix".',
)
@click.option(
    "--route", type=str, default=None, help="Route for the link section."
)
# MECP-specific options (used when --jobtype mecp)
@click.option(
    "--multiplicity-a",
    type=int,
    default=None,
    help="[MECP] Spin multiplicity for state A.",
)
@click.option(
    "--multiplicity-b",
    type=int,
    default=None,
    help="[MECP] Spin multiplicity for state B. Defaults to multiplicity-a + 2.",
)
@click.option(
    "--charge-a",
    type=int,
    default=None,
    help="[MECP] Charge for state A.",
)
@click.option(
    "--charge-b",
    type=int,
    default=None,
    help="[MECP] Charge for state B. Defaults to charge-a.",
)
@click.option(
    "--max-steps",
    type=int,
    default=None,
    help="[MECP] Maximum number of MECP optimization steps.",
)
@click.option(
    "--step-size",
    type=float,
    default=None,
    help="[MECP] Step size (Bohr^2/Hartree) for MECP optimization.",
)
@click.option(
    "--trust-radius",
    type=float,
    default=None,
    help="[MECP] Maximum Cartesian displacement per atom (Bohr) per step.",
)
@click.option(
    "--energy-diff-tol",
    type=float,
    default=None,
    help="[MECP] Convergence threshold for |E(A)-E(B)| in Hartree.",
)
@click.option(
    "--force-max-tol",
    type=float,
    default=None,
    help="[MECP] Convergence threshold for max effective gradient (Hartree/Bohr).",
)
@click.option(
    "--force-rms-tol",
    type=float,
    default=None,
    help="[MECP] Convergence threshold for RMS effective gradient (Hartree/Bohr).",
)
@click.option(
    "--disp-max-tol",
    type=float,
    default=None,
    help="[MECP] Convergence threshold for max displacement (Bohr).",
)
@click.option(
    "--disp-rms-tol",
    type=float,
    default=None,
    help="[MECP] Convergence threshold for RMS displacement (Bohr).",
)
@click.option(
    "--adaptive-step-size/--no-adaptive-step-size",
    default=None,
    help="[MECP] Enable adaptive step size scaling (default: enabled).",
)
@click.option(
    "--step-size-method",
    type=click.Choice(["bb", "grow_shrink"], case_sensitive=False),
    default=None,
    help=(
        "[MECP] Step size adaptation algorithm when adaptive is enabled. "
        "'bb' uses the Barzilai-Borwein secant rule (default). "
        "'grow_shrink' uses a merit-based grow/shrink rule."
    ),
)
@click.option(
    "--step-size-grow",
    type=float,
    default=None,
    help="[MECP] Factor to grow step size when making progress (default: 1.2).",
)
@click.option(
    "--step-size-shrink",
    type=float,
    default=None,
    help="[MECP] Factor to shrink step size on overshoot/oscillation (default: 0.7).",
)
@click.option(
    "--step-size-min",
    type=float,
    default=None,
    help="[MECP] Minimum allowed adaptive step size in Bohr^2/Hartree (default: 1e-4).",
)
@click.option(
    "--step-size-max",
    type=float,
    default=None,
    help="[MECP] Maximum allowed adaptive step size in Bohr^2/Hartree (default: 1.0).",
)
@click.pass_context
def link(
    ctx,
    stable,
    guess,
    jobtype,
    coordinates,
    step_size,
    num_steps,
    remove_solvent,
    solvent_model,
    solvent_id,
    solvent_options,
    route,
    flat_irc,
    predictor,
    recorrect,
    recalc_step,
    maxpoints,
    maxcycles,
    stepsize,
    direction,
    # MECP options
    multiplicity_a,
    multiplicity_b,
    charge_a,
    charge_b,
    max_steps,
    energy_diff_tol,
    force_max_tol,
    force_rms_tol,
    disp_max_tol,
    disp_rms_tol,
    adaptive_step_size,
    step_size_method,
    step_size_grow,
    step_size_shrink,
    step_size_min,
    step_size_max,
    **kwargs,
):
    """CLI subcommand for running Gaussian link jobs."""

    # Dispatch to MECP handler when jobtype is 'mecp'
    if jobtype == "mecp":
        return _link_mecp(
            ctx=ctx,
            stable=stable,
            guess=guess,
            remove_solvent=remove_solvent,
            solvent_model=solvent_model,
            solvent_id=solvent_id,
            solvent_options=solvent_options,
            multiplicity_a=multiplicity_a,
            multiplicity_b=multiplicity_b,
            charge_a=charge_a,
            charge_b=charge_b,
            max_steps=max_steps,
            step_size=step_size,
            energy_diff_tol=energy_diff_tol,
            force_max_tol=force_max_tol,
            force_rms_tol=force_rms_tol,
            disp_max_tol=disp_max_tol,
            disp_rms_tol=disp_rms_tol,
            adaptive_step_size=adaptive_step_size,
            step_size_method=step_size_method,
            step_size_grow=step_size_grow,
            step_size_shrink=step_size_shrink,
            step_size_min=step_size_min,
            step_size_max=step_size_max,
            **kwargs,
        )

    # get jobrunner for running Gaussian link jobs
    jobrunner = ctx.obj["jobrunner"]

    # get settings from project
    from chemsmart.jobs.gaussian.settings import GaussianLinkJobSettings

    project_settings = ctx.obj["project_settings"]
    link_settings = get_setting_from_jobtype_for_gaussian(
        project_settings, jobtype, coordinates, step_size, num_steps
    )

    # job setting from filename or default, with updates from user in cli
    # specified in keywords
    # e.g., `sub.py gaussian -c <user_charge> -m <user_multiplicity>`
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]

    # merge project settings with job settings from cli keywords from
    # cli.gaussian.py subcommands
    link_settings = link_settings.merge(job_settings, keywords=keywords)
    check_charge_and_multiplicity(link_settings)

    # convert from GaussianJobSettings instance to GaussianLinkJobSettings
    # instance with IRC parameters
    link_kwargs = link_settings.__dict__.copy()

    # Add IRC-specific parameters with defaults if this is an IRC job
    if jobtype in ["irc", "ircf", "ircr"]:
        irc_params = {
            "predictor": predictor,
            "recorrect": recorrect,
            "recalc_step": recalc_step if recalc_step is not None else 6,
            "direction": direction,  # Will be set based on jobtype
            "maxpoints": maxpoints if maxpoints is not None else 512,
            "maxcycles": maxcycles if maxcycles is not None else 128,
            "stepsize": stepsize if stepsize is not None else 20,
            "flat_irc": flat_irc if flat_irc is not None else False,
        }
        link_kwargs.update(irc_params)
        logger.info(f"Adding IRC parameters to link job: {irc_params}")

    link_settings = GaussianLinkJobSettings(**link_kwargs)

    # populate GaussianLinkJobSettings
    link_settings.stable = stable
    link_settings.guess = guess
    link_settings.remove_solvent = remove_solvent
    if solvent_model is not None:
        link_settings.solvent_model = solvent_model
    if solvent_id is not None:
        link_settings.solvent_id = solvent_id
    if solvent_options is not None:
        link_settings.additional_solvent_options = solvent_options

    if route is not None:
        link_settings.link_route = route

    # get molecule
    molecules = ctx.obj["molecules"]
    molecule = molecules[-1]

    # get label for the job
    label = ctx.obj["label"]

    if jobtype is None:
        label = label
    else:
        label += f"_{jobtype}"
        if jobtype.lower() == "irc":
            label = update_irc_label(
                label=label,
                direction=link_settings.direction,
                flat_irc=link_settings.flat_irc,
            )
        else:
            label += "_link"

    logger.debug(f"Label for job: {label}")

    # automatically use unrestricted dft if link job
    if not link_settings.functional.lower().startswith("u"):
        link_settings.functional = "u" + link_settings.functional

    logger.info(
        f"Link job {jobtype} settings from project: {link_settings.__dict__}"
    )

    from chemsmart.jobs.gaussian.link import GaussianLinkJob

    return GaussianLinkJob(
        molecule=molecule,
        settings=link_settings,
        label=label,
        jobrunner=jobrunner,
        **kwargs,
    )


def _link_mecp(
    ctx,
    stable,
    guess,
    remove_solvent,
    solvent_model,
    solvent_id,
    solvent_options,
    multiplicity_a,
    multiplicity_b,
    charge_a,
    charge_b,
    max_steps,
    step_size,
    energy_diff_tol,
    force_max_tol,
    force_rms_tol,
    disp_max_tol,
    disp_rms_tol,
    adaptive_step_size,
    step_size_method,
    step_size_grow,
    step_size_shrink,
    step_size_min,
    step_size_max,
    **kwargs,
):
    """
    Handle ``link -j mecp``: broken-symmetry MECP via link sub-jobs.

    Each MECP iteration step runs two ``GaussianLinkJob`` sub-jobs (one per
    spin state) that use ``stable=opt`` to converge the broken-symmetry
    wavefunction and then compute forces on the stable solution.

    The number of α and β electrons is determined by the charge/multiplicity
    line, not by Guess options.  Use ``guess=mix`` (the default) to break
    α/β spatial symmetry and ``stable=opt`` to verify/optimise the
    wavefunction stability.
    """
    from chemsmart.jobs.gaussian.mecp import GaussianMECPJob
    from chemsmart.jobs.gaussian.settings import GaussianMECPJobSettings

    jobrunner = ctx.obj["jobrunner"]
    project_settings = ctx.obj["project_settings"]
    mecp_project_settings = project_settings.opt_settings()

    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]

    mecp_project_settings = mecp_project_settings.merge(
        job_settings, keywords=keywords
    )
    check_charge_and_multiplicity(mecp_project_settings)

    mecp_settings = GaussianMECPJobSettings.from_settings(mecp_project_settings)

    # --- state A charge / multiplicity ---
    if multiplicity_a is None:
        mecp_settings.multiplicity_a = mecp_project_settings.multiplicity
    else:
        mecp_settings.multiplicity_a = multiplicity_a
    if mecp_settings.multiplicity_a is None:
        raise ValueError(
            "State A multiplicity is not set. "
            "Use gaussian -m/--multiplicity or link -j mecp --multiplicity-a."
        )

    if multiplicity_b is None:
        mecp_settings.multiplicity_b = mecp_settings.multiplicity_a + 2
    else:
        mecp_settings.multiplicity_b = multiplicity_b

    if charge_a is None:
        mecp_settings.charge_a = mecp_settings.charge
    else:
        mecp_settings.charge_a = charge_a
    if mecp_settings.charge_a is None:
        raise ValueError(
            "State A charge is not set. "
            "Use gaussian -c/--charge or link -j mecp --charge-a."
        )

    if charge_b is None:
        mecp_settings.charge_b = mecp_settings.charge_a
    else:
        mecp_settings.charge_b = charge_b

    # --- broken-symmetry (link) settings ---
    mecp_settings.use_link = True
    mecp_settings.stable = stable
    mecp_settings.guess = guess

    # --- solvent ---
    if remove_solvent:
        mecp_settings.remove_solvent()
    if solvent_model is not None:
        mecp_settings.solvent_model = solvent_model
    if solvent_id is not None:
        mecp_settings.solvent_id = solvent_id
    if solvent_options is not None:
        mecp_settings.additional_solvent_options = solvent_options

    # --- MECP optimizer parameters ---
    if max_steps is not None:
        mecp_settings.max_steps = max_steps
    if step_size is not None:
        mecp_settings.step_size = step_size
    if energy_diff_tol is not None:
        mecp_settings.energy_diff_tol = energy_diff_tol
    if force_max_tol is not None:
        mecp_settings.force_max_tol = force_max_tol
    if force_rms_tol is not None:
        mecp_settings.force_rms_tol = force_rms_tol
    if disp_max_tol is not None:
        mecp_settings.disp_max_tol = disp_max_tol
    if disp_rms_tol is not None:
        mecp_settings.disp_rms_tol = disp_rms_tol
    if adaptive_step_size is not None:
        mecp_settings.adaptive_step_size = adaptive_step_size
    if step_size_method is not None:
        mecp_settings.step_size_method = step_size_method
    if step_size_grow is not None:
        mecp_settings.step_size_grow = step_size_grow
    if step_size_shrink is not None:
        mecp_settings.step_size_shrink = step_size_shrink
    if step_size_min is not None:
        mecp_settings.step_size_min = step_size_min
    if step_size_max is not None:
        mecp_settings.step_size_max = step_size_max

    # automatically use unrestricted DFT for broken-symmetry calculations
    if not mecp_settings.functional.lower().startswith("u"):
        mecp_settings.functional = "u" + mecp_settings.functional

    molecule = ctx.obj["molecules"][-1]
    label = ctx.obj["label"] + "_mecp_link"

    logger.info(
        f"Link MECP job settings from project: {mecp_settings.__dict__}"
    )

    return GaussianMECPJob(
        molecule=molecule,
        settings=mecp_settings,
        label=label,
        jobrunner=jobrunner,
        **kwargs,
    )
