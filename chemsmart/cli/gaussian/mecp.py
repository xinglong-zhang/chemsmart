import logging

import click

from chemsmart.cli.gaussian.gaussian import gaussian
from chemsmart.cli.job import click_job_options
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
    "--state-a-title",
    type=str,
    default="First",
    help="Title for state A.",
)
@click.option(
    "--state-b-title",
    type=str,
    default="Second",
    help="Title for state B.",
)
@click.option(
    "--max-steps",
    type=int,
    default=50,
    show_default=True,
    help="Maximum number of MECP optimization steps.",
)
@click.option(
    "--step-size",
    type=float,
    default=0.05,
    show_default=True,
    help="Step size (Bohr/Hartree) for internal MECP optimization.",
)
@click.option(
    "--trust-radius",
    type=float,
    default=0.1,
    show_default=True,
    help="Maximum Cartesian displacement per atom (Bohr) per step.",
)
@click.option(
    "--energy-diff-tol",
    type=float,
    default=1.0e-4,
    show_default=True,
    help="Convergence threshold for |E(A)-E(B)| in Hartree.",
)
@click.option(
    "--force-max-tol",
    type=float,
    default=7.0e-4,
    show_default=True,
    help="Convergence threshold for max effective gradient (Hartree/Bohr).",
)
@click.option(
    "--force-rms-tol",
    type=float,
    default=5.0e-4,
    show_default=True,
    help="Convergence threshold for RMS effective gradient (Hartree/Bohr).",
)
@click.option(
    "--disp-max-tol",
    type=float,
    default=1.8e-3,
    show_default=True,
    help="Convergence threshold for max displacement (Bohr).",
)
@click.option(
    "--disp-rms-tol",
    type=float,
    default=1.2e-3,
    show_default=True,
    help="Convergence threshold for RMS displacement (Bohr).",
)
@click.pass_context
def mecp(
    ctx,
    multiplicity_a,
    multiplicity_b,
    charge_a,
    charge_b,
    state_a_title,
    state_b_title,
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

    # merge project irc settings with job settings from cli keywords from
    # cli.orca.py subcommands
    mecp_settings = mecp_project_settings.merge(
        job_settings, keywords=keywords
    )

    # update mecp_settings if any attribute is specified in cli options
    if multiplicity_a is None:
        multiplicity_a = mecp_settings.multiplicity
    if multiplicity_a is None:
        raise ValueError(
            "State A multiplicity is not set. "
            "Use gaussian -m/--multiplicity or mecp --multiplicity-a."
        )
    if multiplicity_b is None:
        multiplicity_b = multiplicity_a + 2

    if charge_a is None:
        charge_a = mecp_settings.charge
    if charge_a is None:
        raise ValueError(
            "State A charge is not set. "
            "Use gaussian -c/--charge or mecp --charge-a."
        )
    if charge_b is None:
        charge_b = charge_a

    check_charge_and_multiplicity(mecp_settings)

    # get molecule
    molecules = ctx.obj[
        "molecules"
    ]  # use all molecules as a list for crest jobs
    molecule = molecules[-1]  # get last molecule from list of molecules

    # get label for the job
    label = ctx.obj["label"]

    logger.debug(f"Label for job: {label}")

    logger.info(f"IRC job settings from project: {mecp_settings.__dict__}")

    from chemsmart.jobs.gaussian.mecp import GaussianMECPJob

    return GaussianMECPJob(
        molecule=molecule,
        settings=mecp_settings,
        label=label,
        jobrunner=jobrunner,
        state_a_multiplicity=multiplicity_a,
        state_b_multiplicity=multiplicity_b,
        state_a_charge=charge_a,
        state_b_charge=charge_b,
        state_a_title=state_a_title,
        state_b_title=state_b_title,
        max_steps=max_steps,
        step_size=step_size,
        trust_radius=trust_radius,
        energy_diff_tol=energy_diff_tol,
        force_max_tol=force_max_tol,
        force_rms_tol=force_rms_tol,
        disp_max_tol=disp_max_tol,
        disp_rms_tol=disp_rms_tol,
        skip_completed=skip_completed,
        **kwargs,
    )
