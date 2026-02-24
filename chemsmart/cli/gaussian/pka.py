"""
CLI for Gaussian pKa calculations.

This module provides the CLI interface for running pKa calculations
using Gaussian. It creates two sequential optimization jobs: one for
the protonated form (HA) and one for the conjugate base (A-).
"""

import functools
import logging

import click

from chemsmart.cli.gaussian.gaussian import (
    click_gaussian_solvent_options,
    gaussian,
)
from chemsmart.cli.job import click_job_options
from chemsmart.utils.cli import MyCommand

logger = logging.getLogger(__name__)


def click_pka_options(f):
    """Click options specific to pKa calculations."""

    @click.option(
        "-pi",
        "--proton-index",
        type=int,
        required=True,
        help="1-based index of the proton to remove for deprotonation.",
    )
    @click.option(
        "-r",
        "--reference",
        type=str,
        default="water",
        help="Reference acid/base for pKa calculation (default: water).",
    )
    @click.option(
        "-t",
        "--thermodynamic-cycle",
        type=click.Choice(["direct", "isodesmic"]),
        default="direct",
        help="Thermodynamic cycle type (default: direct).",
    )
    @click.option(
        "-cc",
        "--conjugate-base-charge",
        type=int,
        default=None,
        help="Charge of the conjugate base (A-). Defaults to (charge - 1).",
    )
    @click.option(
        "-cm",
        "--conjugate-base-multiplicity",
        type=int,
        default=None,
        help="Multiplicity of the conjugate base (A-). Defaults to multiplicity.",
    )
    @functools.wraps(f)
    def wrapper(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper


@gaussian.command("pka", cls=MyCommand)
@click_job_options
@click_gaussian_solvent_options
@click_pka_options
@click.pass_context
def pka(
    ctx,
    proton_index,
    reference,
    thermodynamic_cycle,
    conjugate_base_charge,
    conjugate_base_multiplicity,
    remove_solvent,
    solvent_model,
    solvent_id,
    solvent_options,
    skip_completed,
    **kwargs,
):
    """
    CLI subcommand for running Gaussian pKa calculations.

    Creates two sequential optimization jobs with frequency calculations:
    one for the protonated form (HA) and one for the conjugate base (A-).

    The proton to be removed is specified by --proton-index (1-based).
    The protonated form uses charge and multiplicity from -c and -m options.
    The conjugate base charge defaults to (charge - 1).

    Example:
        chemsmart run gaussian -f acetic_acid.xyz -c 0 -m 1 pka -pi 10 -sm SMD -si water
    """
    from chemsmart.jobs.gaussian.pka import GaussianpKaJob
    from chemsmart.jobs.gaussian.settings import GaussianpKaJobSettings

    # Get jobrunner
    jobrunner = ctx.obj["jobrunner"]

    # Get settings from project
    project_settings = ctx.obj["project_settings"]

    # Get base opt settings from project (pKa uses opt+freq)
    opt_settings = project_settings.opt_settings()

    # Job settings from CLI keywords
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]

    # Merge project opt settings with CLI job settings
    opt_settings = opt_settings.merge(job_settings, keywords=keywords)

    # Get molecules
    molecules = ctx.obj["molecules"]
    molecule = molecules[-1]

    # Get label for the job
    label = ctx.obj["label"]

    # Create pKa settings - charge and multiplicity are inherited from opt_settings
    pka_settings = GaussianpKaJobSettings(
        proton_index=proton_index,
        reference=reference,
        thermodynamic_cycle=thermodynamic_cycle,
        conjugate_base_charge=conjugate_base_charge,
        conjugate_base_multiplicity=conjugate_base_multiplicity,
        # Inherit charge and multiplicity from opt_settings (protonated form)
        charge=opt_settings.charge,
        multiplicity=opt_settings.multiplicity,
        # Inherit other settings from opt_settings
        functional=opt_settings.functional,
        basis=opt_settings.basis,
        ab_initio=opt_settings.ab_initio,
        semiempirical=opt_settings.semiempirical,
        solvent_id=opt_settings.solvent_id,
        additional_route_parameters=opt_settings.additional_route_parameters,
        gen_genecp_file=opt_settings.gen_genecp_file,
        heavy_elements=opt_settings.heavy_elements,
        heavy_elements_basis=opt_settings.heavy_elements_basis,
        light_elements_basis=opt_settings.light_elements_basis,
    )

    # Apply CLI-supplied solvent options
    if remove_solvent:
        pka_settings.solvation_model = None
        pka_settings.solvent_id = None
    else:
        if solvent_model is not None:
            pka_settings.solvation_model = solvent_model
        if solvent_id is not None:
            pka_settings.solvent_id = solvent_id

    if solvent_options is not None:
        pka_settings.additional_solvent_options = solvent_options

    # Validate charge and multiplicity for protonated form
    if pka_settings.charge is None:
        raise click.UsageError(
            "Charge must be specified via -c/--charge option"
        )
    if pka_settings.multiplicity is None:
        raise click.UsageError(
            "Multiplicity must be specified via -m/--multiplicity option"
        )

    logger.info(f"pKa job settings: {pka_settings.__dict__}")
    logger.info(f"Proton index to remove: {proton_index}")
    logger.info(
        f"Protonated form: charge={pka_settings.charge}, "
        f"mult={pka_settings.multiplicity}"
    )

    # Calculate conjugate base charge/mult for logging
    cb_charge = (
        conjugate_base_charge
        if conjugate_base_charge is not None
        else pka_settings.charge - 1
    )
    cb_mult = (
        conjugate_base_multiplicity
        if conjugate_base_multiplicity is not None
        else pka_settings.multiplicity
    )
    logger.info(f"Conjugate base: charge={cb_charge}, mult={cb_mult}")

    # Create and return pKa job
    return GaussianpKaJob(
        molecule=molecule,
        settings=pka_settings,
        label=label,
        jobrunner=jobrunner,
        skip_completed=skip_completed,
        **kwargs,
    )
