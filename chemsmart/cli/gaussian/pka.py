"""
CLI for Gaussian pKa calculations.

This module provides the CLI interface for running pKa calculations
using Gaussian with a proper thermodynamic cycle:
1. Gas phase optimization + frequency for both HA and A-
2. Solution phase single point for both HA and A- at the same level of theory

Two thermodynamic cycles are supported:
- **proton exchange**: Uses a reference acid to cancel systematic errors (default)
- **direct**: Uses the absolute free energy of a proton in water

Using the same level of theory ensures proper error cancellation for
solvation free energy calculations.
"""

import functools
import logging

import click

from chemsmart.cli.gaussian.gaussian import (
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
        "-t",
        "--thermodynamic-cycle",
        type=click.Choice(["direct", "proton exchange"]),
        default="proton exchange",
        help="Thermodynamic cycle type. 'proton exchange' uses a reference acid "
        "(default). 'direct' uses absolute free energy of H+ in water.",
    )
    @click.option(
        "-r",
        "--reference",
        type=str,
        default="water",
        help="Reference acid for proton exchange cycle (default: water). "
        "Ignored when using direct cycle.",
    )
    @click.option(
        "-dG",
        "--delta-g-proton",
        type=float,
        default=-265.9,
        help="Absolute free energy of H+ in water (kcal/mol) for direct cycle. "
        "Default: -265.9 kcal/mol (Tissandier et al., 1998).",
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
    @click.option(
        "-sm",
        "--solvent-model",
        type=str,
        default="SMD",
        help="Solvation model for solution phase SP (default: SMD).",
    )
    @click.option(
        "-si",
        "--solvent-id",
        type=str,
        default="water",
        help="Solvent ID for solution phase SP (default: water).",
    )
    @functools.wraps(f)
    def wrapper(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper


@gaussian.command("pka", cls=MyCommand)
@click_job_options
@click_pka_options
@click.pass_context
def pka(
    ctx,
    proton_index,
    thermodynamic_cycle,
    reference,
    delta_g_proton,
    conjugate_base_charge,
    conjugate_base_multiplicity,
    solvent_model,
    solvent_id,
    skip_completed,
    **kwargs,
):
    """
    CLI subcommand for running Gaussian pKa calculations.

    Performs pKa calculations using a proper thermodynamic cycle:
    1. Gas phase optimization + frequency for HA and A-
    2. Solution phase single point for HA and A- at the SAME level of theory

    Two thermodynamic cycles are available:

    \b
    - **proton exchange** (default): Uses a reference acid to cancel errors.
      HA + Ref- → A- + HRef
      pKa(HA) = pKa(HRef) + ΔG_exchange / (2.303 * R * T)

    \b
    - **direct**: Uses absolute free energy of H+ in water.
      pKa = [G(A-)_aq - G(HA)_aq + ΔG°(H+)_aq] / (2.303 * R * T)
      Default ΔG°(H+)_aq = -265.9 kcal/mol (Tissandier et al., 1998)

    The proton to be removed is specified by --proton-index (1-based).
    The protonated form uses charge and multiplicity from -c and -m options.
    The conjugate base charge defaults to (charge - 1).

    \b
    Examples:
        # Proton exchange cycle with water as reference
        chemsmart run gaussian -f acetic_acid.xyz -c 0 -m 1 pka -pi 10 -tc "proton exchange" -ref water

        # Direct cycle (no reference needed)
        chemsmart run gaussian -f acetic_acid.xyz -c 0 -m 1 pka -pi 10 -tc direct
    """
    from chemsmart.jobs.gaussian.pka import GaussianpKaJob
    from chemsmart.jobs.gaussian.settings import GaussianpKaJobSettings

    # Get jobrunner
    jobrunner = ctx.obj["jobrunner"]

    # Get settings from project
    project_settings = ctx.obj["project_settings"]

    # Get base opt settings from project (pKa uses opt+freq in gas phase)
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

    # Create pKa settings
    # - Gas phase opt uses functional/basis from opt_settings
    # - Solution phase SP uses SAME functional/basis for error cancellation
    pka_settings = GaussianpKaJobSettings(
        proton_index=proton_index,
        thermodynamic_cycle=thermodynamic_cycle,
        reference=reference,
        delta_G_proton=delta_g_proton,
        conjugate_base_charge=conjugate_base_charge,
        conjugate_base_multiplicity=conjugate_base_multiplicity,
        solvent_model=solvent_model,
        solvent_id=solvent_id,
        # Inherit charge and multiplicity from opt_settings (protonated form)
        charge=opt_settings.charge,
        multiplicity=opt_settings.multiplicity,
        # Inherit level of theory from opt_settings (used for BOTH gas and solution)
        functional=opt_settings.functional,
        basis=opt_settings.basis,
        ab_initio=opt_settings.ab_initio,
        semiempirical=opt_settings.semiempirical,
        additional_route_parameters=opt_settings.additional_route_parameters,
        gen_genecp_file=opt_settings.gen_genecp_file,
        heavy_elements=opt_settings.heavy_elements,
        heavy_elements_basis=opt_settings.heavy_elements_basis,
        light_elements_basis=opt_settings.light_elements_basis,
    )

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
    logger.info(f"Thermodynamic cycle: {thermodynamic_cycle}")
    logger.info(
        f"Protonated form (HA): charge={pka_settings.charge}, "
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
    logger.info(f"Conjugate base (A-): charge={cb_charge}, mult={cb_mult}")

    if thermodynamic_cycle == "proton exchange":
        logger.info(f"Reference acid: {reference}")
    else:
        logger.info(f"ΔG°(H+)_aq: {delta_g_proton} kcal/mol")

    logger.info(
        f"Gas phase optimization: {pka_settings.functional}/{pka_settings.basis}"
    )
    logger.info(
        f"Solution phase SP: {pka_settings.functional}/{pka_settings.basis} "
        f"with {solvent_model}({solvent_id})"
    )

    # Create and return pKa job
    return GaussianpKaJob(
        molecule=molecule,
        settings=pka_settings,
        label=label,
        jobrunner=jobrunner,
        skip_completed=skip_completed,
        **kwargs,
    )
