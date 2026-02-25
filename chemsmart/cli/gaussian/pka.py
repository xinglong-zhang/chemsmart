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
from chemsmart.utils.cli import MyCommand, MyGroup

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
        type=click.Path(exists=True),
        default=None,
        help="Path to geometry file for reference acid (HB) for proton exchange cycle. "
        "When provided, optimization and SP calculations will also be run for HB and B-.",
    )
    @click.option(
        "-rpi",
        "--reference-proton-index",
        type=int,
        default=None,
        help="1-based index of the proton to remove from reference acid (HB). "
        "Required when --reference is provided.",
    )
    @click.option(
        "-rc",
        "--reference-charge",
        type=int,
        default=None,
        help="Charge of the reference acid (HB). Required when --reference is provided.",
    )
    @click.option(
        "-rm",
        "--reference-multiplicity",
        type=int,
        default=None,
        help="Multiplicity of the reference acid (HB). Required when --reference is provided.",
    )
    @click.option(
        "-rcc",
        "--reference-conjugate-base-charge",
        type=int,
        default=None,
        help="Charge of the reference conjugate base (B-). Defaults to (reference_charge - 1).",
    )
    @click.option(
        "-rcm",
        "--reference-conjugate-base-multiplicity",
        type=int,
        default=None,
        help="Multiplicity of the reference conjugate base (B-). Defaults to reference_multiplicity.",
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
    @click.option(
        "-T",
        "--temperature",
        type=float,
        default=298.15,
        help="Temperature in Kelvin for thermochemistry calculation (default: 298.15 K).",
    )
    @click.option(
        "-conc",
        "--concentration",
        type=float,
        default=1.0,
        help="Concentration in mol/L for thermochemistry (default: 1.0 mol/L).",
    )
    @click.option(
        "-csg",
        "--cutoff-entropy-grimme",
        type=float,
        default=100.0,
        help="Cutoff frequency for entropy (cm^-1) using Grimme's quasi-RRHO method (default: 100).",
    )
    @click.option(
        "-ch",
        "--cutoff-enthalpy",
        type=float,
        default=100.0,
        help="Cutoff frequency for enthalpy (cm^-1) using Head-Gordon's method (default: 100).",
    )
    @functools.wraps(f)
    def wrapper(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper


@gaussian.group("pka", cls=MyGroup, invoke_without_command=True)
@click_job_options
@click_pka_options
@click.pass_context
def pka(
    ctx,
    proton_index,
    thermodynamic_cycle,
    reference,
    reference_proton_index,
    reference_charge,
    reference_multiplicity,
    reference_conjugate_base_charge,
    reference_conjugate_base_multiplicity,
    delta_g_proton,
    conjugate_base_charge,
    conjugate_base_multiplicity,
    solvent_model,
    solvent_id,
    temperature,
    concentration,
    cutoff_entropy_grimme,
    cutoff_enthalpy,
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
      HA + B- → A- + HB
      pKa(HA) = pKa(HB) + ΔG_exchange / (2.303 * R * T)

      When using proton exchange, provide a reference acid geometry file
      with --reference (-r) option. This will run optimization and SP
      calculations for both HB and B- alongside HA and A-.

    \b
    - **direct**: Uses absolute free energy of H+ in water.
      pKa = [G(A-)_aq - G(HA)_aq + ΔG°(H+)_aq] / (2.303 * R * T)
      Default ΔG°(H+)_aq = -265.9 kcal/mol (Tissandier et al., 1998)

    The proton to be removed is specified by --proton-index (1-based).
    The protonated form uses charge and multiplicity from -c and -m options.
    The conjugate base charge defaults to (charge - 1).

    \b
    Subcommands:
        thermo    Extract thermochemistry (E, qh-G) from completed output files

    \b
    Examples:
        # Run pKa job with proton exchange cycle
        chemsmart run gaussian -f acetic_acid.xyz -c 0 -m 1 pka -pi 10 \\
            -t "proton exchange" -r water.xyz -rpi 1 -rc 0 -rm 1

        # Run pKa job with direct cycle
        chemsmart run gaussian -f acetic_acid.xyz -c 0 -m 1 pka -pi 10 -t direct

        # Extract thermochemistry from completed jobs
        chemsmart run gaussian pka thermo -ha acid_opt.log -a base_opt.log -T 298.15
    """
    # If a subcommand is invoked, don't run the pKa job
    if ctx.invoked_subcommand is not None:
        # Store settings in context for subcommands
        ctx.obj["pka_settings"] = {
            "temperature": temperature,
            "concentration": concentration,
            "cutoff_entropy_grimme": cutoff_entropy_grimme,
            "cutoff_enthalpy": cutoff_enthalpy,
        }
        return

    from chemsmart.jobs.gaussian.pka import GaussianpKaJob
    from chemsmart.jobs.gaussian.settings import GaussianpKaJobSettings

    # Validate reference acid settings if provided
    if reference is not None:
        if thermodynamic_cycle != "proton exchange":
            raise click.UsageError(
                "Reference acid file can only be used with 'proton exchange' cycle. "
                "Use -t 'proton exchange' or remove the -r option."
            )
        missing = []
        if reference_proton_index is None:
            missing.append("-rpi/--reference-proton-index")
        if reference_charge is None:
            missing.append("-rc/--reference-charge")
        if reference_multiplicity is None:
            missing.append("-rm/--reference-multiplicity")
        if missing:
            raise click.UsageError(
                f"When --reference is provided, the following options are required: "
                f"{', '.join(missing)}"
            )

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
        reference_file=reference,
        reference_proton_index=reference_proton_index,
        reference_charge=reference_charge,
        reference_multiplicity=reference_multiplicity,
        reference_conjugate_base_charge=reference_conjugate_base_charge,
        reference_conjugate_base_multiplicity=reference_conjugate_base_multiplicity,
        delta_G_proton=delta_g_proton,
        conjugate_base_charge=conjugate_base_charge,
        conjugate_base_multiplicity=conjugate_base_multiplicity,
        solvent_model=solvent_model,
        solvent_id=solvent_id,
        # Thermochemistry settings
        temperature=temperature,
        concentration=concentration,
        cutoff_entropy_grimme=cutoff_entropy_grimme,
        cutoff_enthalpy=cutoff_enthalpy,
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
        if reference is not None:
            logger.info(f"Reference acid file: {reference}")
            logger.info(f"Reference proton index: {reference_proton_index}")
            logger.info(
                f"Reference acid (HB): charge={reference_charge}, "
                f"mult={reference_multiplicity}"
            )
            ref_cb_charge = (
                reference_conjugate_base_charge
                if reference_conjugate_base_charge is not None
                else reference_charge - 1
            )
            ref_cb_mult = (
                reference_conjugate_base_multiplicity
                if reference_conjugate_base_multiplicity is not None
                else reference_multiplicity
            )
            logger.info(
                f"Reference conjugate base (B-): charge={ref_cb_charge}, "
                f"mult={ref_cb_mult}"
            )
        else:
            logger.info("No reference acid file provided")
    else:
        logger.info(f"ΔG°(H+)_aq: {delta_g_proton} kcal/mol")

    logger.info(
        f"Gas phase optimization: {pka_settings.functional}/{pka_settings.basis}"
    )
    logger.info(
        f"Solution phase SP: {pka_settings.functional}/{pka_settings.basis} "
        f"with {solvent_model}({solvent_id})"
    )
    logger.info(
        f"Thermochemistry: T={temperature}K, c={concentration}mol/L, "
        f"csg={cutoff_entropy_grimme}cm^-1, ch={cutoff_enthalpy}cm^-1"
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


def click_pka_thermo_options(f):
    """Click options for pKa thermochemistry extraction."""

    @click.option(
        "-ha",
        "--ha-file",
        type=click.Path(exists=True),
        default=None,
        help="Path to HA (protonated acid) optimization output file.",
    )
    @click.option(
        "-a",
        "--a-file",
        type=click.Path(exists=True),
        default=None,
        help="Path to A- (conjugate base) optimization output file.",
    )
    @click.option(
        "-hb",
        "--hb-file",
        type=click.Path(exists=True),
        default=None,
        help="Path to HB (reference acid) optimization output file.",
    )
    @click.option(
        "-b",
        "--b-file",
        type=click.Path(exists=True),
        default=None,
        help="Path to B- (reference conjugate base) optimization output file.",
    )
    @click.option(
        "-T",
        "--temperature",
        type=float,
        default=298.15,
        help="Temperature in Kelvin for thermochemistry calculation (default: 298.15 K).",
    )
    @click.option(
        "-conc",
        "--concentration",
        type=float,
        default=1.0,
        help="Concentration in mol/L for thermochemistry (default: 1.0 mol/L).",
    )
    @click.option(
        "-csg",
        "--cutoff-entropy-grimme",
        type=float,
        default=100.0,
        help="Cutoff frequency for entropy (cm^-1) using Grimme's quasi-RRHO method (default: 100).",
    )
    @click.option(
        "-ch",
        "--cutoff-enthalpy",
        type=float,
        default=100.0,
        help="Cutoff frequency for enthalpy (cm^-1) using Head-Gordon's method (default: 100).",
    )
    @click.option(
        "-u",
        "--energy-units",
        type=click.Choice(
            ["hartree", "eV", "kcal/mol", "kJ/mol"], case_sensitive=False
        ),
        default="hartree",
        help="Energy units for output (default: hartree).",
    )
    @click.option(
        "-o",
        "--output",
        type=click.Path(),
        default=None,
        help="Output file to save results (default: print to stdout).",
    )
    @functools.wraps(f)
    def wrapper(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper


@pka.command("thermo", cls=MyCommand)
@click_pka_thermo_options
@click.pass_context
def thermo(
    ctx,
    ha_file,
    a_file,
    hb_file,
    b_file,
    temperature,
    concentration,
    cutoff_entropy_grimme,
    cutoff_enthalpy,
    energy_units,
    output,
    **kwargs,
):
    """
    Extract thermochemistry (E and qh-G(T)) from pKa optimization output files.

    This command extracts electronic energies (E) and quasi-harmonic Gibbs free
    energies (qh-G(T)) from Gaussian optimization output files for pKa calculations.

    The quasi-harmonic corrections use Grimme's quasi-RRHO method for entropy
    and Head-Gordon's quasi-RRHO method for enthalpy, equivalent to running:
        chemsmart run thermochemistry -f <file> -T <temp> -c <conc> -csg <cutoff> -ch <cutoff>

    Species files:
        -ha: HA (protonated acid) output file
        -a:  A- (conjugate base) output file
        -hb: HB (reference acid) output file (for proton exchange cycle)
        -b:  B- (reference conjugate base) output file (for proton exchange cycle)

    \b
    Examples:
        # Direct cycle (HA and A- only)
        chemsmart run gaussian pka thermo -ha acid_opt.log -a base_opt.log -T 298.15

        # Proton exchange cycle (all four species)
        chemsmart run gaussian pka thermo \\
            -ha acetic_acid_opt.log -a acetate_opt.log \\
            -hb water_opt.log -b hydroxide_opt.log \\
            -T 298.15 -csg 100 -ch 100

        # With custom temperature and units
        chemsmart run gaussian pka thermo \\
            -ha acid_opt.log -a base_opt.log \\
            -T 333.15 -u kcal/mol -o results.txt
    """
    from chemsmart.io.gaussian.output import Gaussian16pKaOutput

    # Validate that at least one file is provided
    if all(f is None for f in [ha_file, a_file, hb_file, b_file]):
        raise click.UsageError(
            "At least one output file must be provided. Use -ha, -a, -hb, or -b options."
        )

    # Compute and display results
    if output is not None:
        # Save to file
        import sys
        from io import StringIO

        # Capture stdout
        old_stdout = sys.stdout
        sys.stdout = StringIO()
        Gaussian16pKaOutput.print_pka_summary(
            ha_file=ha_file,
            a_file=a_file,
            hb_file=hb_file,
            b_file=b_file,
            temperature=temperature,
            concentration=concentration,
            cutoff_entropy_grimme=cutoff_entropy_grimme,
            cutoff_enthalpy=cutoff_enthalpy,
            energy_units=energy_units,
        )
        result = sys.stdout.getvalue()
        sys.stdout = old_stdout

        with open(output, "w") as f:
            f.write(result)
        logger.info(f"Results saved to {output}")
        click.echo(f"Results saved to {output}")
    else:
        # Print to stdout
        Gaussian16pKaOutput.print_pka_summary(
            ha_file=ha_file,
            a_file=a_file,
            hb_file=hb_file,
            b_file=b_file,
            temperature=temperature,
            concentration=concentration,
            cutoff_entropy_grimme=cutoff_entropy_grimme,
            cutoff_enthalpy=cutoff_enthalpy,
            energy_units=energy_units,
        )

    # Return the results dictionary for programmatic access
    return Gaussian16pKaOutput.compute_pka_thermochemistry(
        ha_file=ha_file,
        a_file=a_file,
        hb_file=hb_file,
        b_file=b_file,
        temperature=temperature,
        concentration=concentration,
        cutoff_entropy_grimme=cutoff_entropy_grimme,
        cutoff_enthalpy=cutoff_enthalpy,
        energy_units=energy_units,
    )
