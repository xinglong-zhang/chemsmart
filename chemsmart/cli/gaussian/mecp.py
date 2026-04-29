import logging

import click

from chemsmart.cli.gaussian.gaussian import gaussian
from chemsmart.cli.job import click_job_options
from chemsmart.jobs.gaussian.settings import GaussianMECPJobSettings
from chemsmart.utils.cli import MyCommand
from chemsmart.utils.utils import check_charge_and_multiplicity

logger = logging.getLogger(__name__)


@gaussian.command("mecp", cls=MyCommand)
@click_job_options
@click.option(
    "--multiplicity-a",
    type=int,
    default=None,
    help="Spin multiplicity for state A.",
)
@click.option(
    "--multiplicity-b",
    type=int,
    default=None,
    help="Spin multiplicity for state B. Defaults to multiplicity-a + 2.",
)
@click.option(
    "--charge-a",
    type=int,
    default=None,
    help="Charge for state A.",
)
@click.option(
    "--charge-b",
    type=int,
    default=None,
    help="Charge for state B. Defaults to charge-a.",
)
@click.option(
    "--title-a",
    type=str,
    default="First",
    help="Title for state A.",
)
@click.option(
    "--title-b",
    type=str,
    default="Second",
    help="Title for state B.",
)
@click.option(
    "--max-steps",
    type=int,
    default=None,
    help="Maximum number of MECP optimization steps.",
)
@click.option(
    "--step-size",
    type=float,
    default=None,
    help="Step size (Bohr^2/Hartree) for internal MECP optimization.",
)
@click.option(
    "--trust-radius",
    type=float,
    default=None,
    help="Maximum Cartesian displacement per atom (Bohr) per step.",
)
@click.option(
    "--energy-diff-tol",
    type=float,
    default=None,
    help="Convergence threshold for |E(A)-E(B)| in Hartree.",
)
@click.option(
    "--force-max-tol",
    type=float,
    default=None,
    help="Convergence threshold for max effective gradient (Hartree/Bohr).",
)
@click.option(
    "--force-rms-tol",
    type=float,
    default=None,
    help="Convergence threshold for RMS effective gradient (Hartree/Bohr).",
)
@click.option(
    "--disp-max-tol",
    type=float,
    default=None,
    help="Convergence threshold for max displacement (Bohr).",
)
@click.option(
    "--disp-rms-tol",
    type=float,
    default=None,
    help="Convergence threshold for RMS displacement (Bohr).",
)
@click.pass_context
def mecp(
    ctx,
    multiplicity_a,
    multiplicity_b,
    charge_a,
    charge_b,
    title_a,
    title_b,
    max_steps,
    step_size,
    trust_radius,
    energy_diff_tol,
    force_max_tol,
    force_rms_tol,
    disp_max_tol,
    disp_rms_tol,
    skip_completed,
    **kwargs,
):
    """
    CLI subcommand for running Gaussian Minimum Energy Cross Point (MECP) jobs.
    """

    # get jobrunner for running Gaussian MECP jobs
    jobrunner = ctx.obj["jobrunner"]

    # get settings from project
    project_settings = ctx.obj["project_settings"]
    mecp_project_settings = project_settings.opt_settings()

    # job setting from filename or default, with updates from user in cli
    # specified in keywords
    # e.g., `sub.py gaussian -c <user_charge> -m <user_multiplicity>`
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]

    # merge project mecp settings with job settings from cli keywords from
    # cli.orca.py subcommands
    mecp_project_settings = mecp_project_settings.merge(
        job_settings, keywords=keywords
    )
    check_charge_and_multiplicity(mecp_project_settings)

    # convert project mecp settings to GaussianMECPJobSettings
    mecp_settings = GaussianMECPJobSettings.from_settings(
        mecp_project_settings
    )

    # update mecp_settings if any attribute is specified in cli options
    if multiplicity_a is None:
        # if not given, then takes value from gaussian project settings or input file
        mecp_settings.multiplicity_a = mecp_project_settings.multiplicity
    else:
        mecp_settings.multiplicity_a = multiplicity_a
    if mecp_settings.multiplicity_a is None:
        raise ValueError(
            "State A multiplicity is not set. "
            "Use gaussian -m/--multiplicity or mecp --multiplicity-a."
        )

    if multiplicity_b is None:
        # if not given, then defaults to multiplicity_a + 2
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
            "Use gaussian -c/--charge or mecp --charge-a."
        )
    if charge_b is None:
        mecp_settings.charge_b = mecp_settings.charge_a
    else:
        mecp_settings.charge_b = charge_b

    if title_a is not None:
        mecp_settings.title_a = title_a
    if title_b is not None:
        mecp_settings.title_b = title_b

    if max_steps is not None:
        mecp_settings.max_steps = max_steps
    if step_size is not None:
        mecp_settings.step_size = step_size
    if trust_radius is not None:
        mecp_settings.trust_radius = trust_radius
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

    # get molecule
    molecules = ctx.obj[
        "molecules"
    ]  # use all molecules as a list for crest jobs
    molecule = molecules[-1]  # get last molecule from list of molecules

    # get label for the job
    label = ctx.obj["label"]

    logger.debug(f"Label for job: {label}")

    logger.info(f"MECP job settings from project: {mecp_settings.__dict__}")

    from chemsmart.jobs.gaussian.mecp import GaussianMECPJob

    return GaussianMECPJob(
        molecule=molecule,
        settings=mecp_settings,
        label=label,
        jobrunner=jobrunner,
        skip_completed=skip_completed,
        **kwargs,
    )
