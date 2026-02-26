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


def click_pka_output_options(f):
    """Click options for parsing pKa output files.

    These options allow users to provide completed Gaussian output files
    for extracting thermochemistry and computing pKa values using the
    Dual-level Proton Exchange scheme.

    All energies are in Hartree (au) except ΔG_soln which is converted
    to kcal/mol for the pKa formula.
    """
    # Apply options in reverse order (last applied = first in help)
    # Reference pKa value
    f = click.option(
        "-rp",
        "--reference-pka",
        type=float,
        default=None,
        help="Experimental pKa of reference acid HB for proton exchange cycle. "
        "Required when using output file parsing mode with reference acid.",
    )(f)
    # Solvent single-point output files
    f = click.option(
        "-bs",
        "--b-solv-output",
        type=click.Path(exists=True),
        default=None,
        help="Path to B⁻ solvent single-point output file.",
    )(f)
    f = click.option(
        "-hbs",
        "--hb-solv-output",
        type=click.Path(exists=True),
        default=None,
        help="Path to HB solvent single-point output file.",
    )(f)
    f = click.option(
        "-as",
        "--a-solv-output",
        type=click.Path(exists=True),
        default=None,
        help="Path to A⁻ solvent single-point output file.",
    )(f)
    f = click.option(
        "-has",
        "--ha-solv-output",
        type=click.Path(exists=True),
        default=None,
        help="Path to HA solvent single-point output file.",
    )(f)
    # Gas-phase optimization+frequency output files
    f = click.option(
        "-b",
        "--b-output",
        type=click.Path(exists=True),
        default=None,
        help="Path to B⁻ (reference conjugate base) gas-phase opt+freq output file.",
    )(f)
    f = click.option(
        "-hb",
        "--hb-output",
        type=click.Path(exists=True),
        default=None,
        help="Path to HB (reference acid) gas-phase opt+freq output file.",
    )(f)
    f = click.option(
        "-a",
        "--a-output",
        type=click.Path(exists=True),
        default=None,
        help="Path to A⁻ (conjugate base) gas-phase opt+freq output file.",
    )(f)
    f = click.option(
        "-ha",
        "--ha-output",
        type=click.Path(exists=True),
        default=None,
        help="Path to HA (protonated acid) gas-phase opt+freq output file.",
    )(f)

    return f


@gaussian.group("pka", cls=MyGroup, invoke_without_command=True)
@click_job_options
@click_pka_options
@click_pka_output_options
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
    # Output file parsing options
    ha_output,
    a_output,
    hb_output,
    b_output,
    ha_solv_output,
    a_solv_output,
    hb_solv_output,
    b_solv_output,
    reference_pka,
    **kwargs,
):
    """
    CLI subcommand for running Gaussian pKa calculations.

    Performs pKa calculations using a proper thermodynamic cycle:
    1. Gas phase optimization + frequency for HA and A-
    2. Solution phase single point for HA and A- at the SAME level of theory

    Two modes are available:

    \b
    **Mode 1: Run new calculations**
    Provide input geometry files and run optimization + SP calculations.

    \b
    **Mode 2: Parse existing output files**
    Provide completed Gaussian output files to compute pKa directly using
    the Dual-level Proton Exchange scheme. All energies in Hartree (au),
    except ΔG_soln which is converted to kcal/mol for the pKa formula.

    Two thermodynamic cycles are available:

    \b
    - **proton exchange** (default): Uses a reference acid to cancel errors.
      HA + B⁻ → A⁻ + HB
      pKa(HA) = pKa(HB) + ΔG_soln / (2.303 × R × T)

    \b
    - **direct**: Uses absolute free energy of H⁺ in water.
      pKa = [G(A⁻)_soln - G(HA)_soln + ΔG°(H⁺)_aq] / (2.303 × R × T)
      Default ΔG°(H⁺)_aq = -265.9 kcal/mol (Tissandier et al., 1998)

    \b
    Output file options (for parsing mode):
        -ha, --ha-output        HA gas-phase opt+freq output
        -a,  --a-output         A⁻ gas-phase opt+freq output
        -hb, --hb-output        HB gas-phase opt+freq output
        -b,  --b-output         B⁻ gas-phase opt+freq output
        -has, --ha-solv-output  HA solvent SP output
        -as,  --a-solv-output   A⁻ solvent SP output
        -hbs, --hb-solv-output  HB solvent SP output
        -bs,  --b-solv-output   B⁻ solvent SP output
        -rpka, --reference-pka  Experimental pKa of reference acid HB

    \b
    Examples:
        # Run pKa job with proton exchange cycle
        chemsmart run gaussian -f acetic_acid.xyz -c 0 -m 1 pka -pi 10 \\
            -t "proton exchange" -r water.xyz -rpi 1 -rc 0 -rm 1

        # Run pKa job with direct cycle
        chemsmart run gaussian -f acetic_acid.xyz -c 0 -m 1 pka -pi 10 -t direct

        # Compute pKa from existing output files (Dual-level Proton Exchange)
        chemsmart run gaussian pka -pi 1 \\
            -ha 5PQ_Me_ts1_no_pd_opt.log \\
            -a 5PQ_Me_ts1_b_no_pd_opt.log \\
            -hb collidine-H_opt.log \\
            -b collidine_opt.log \\
            -has 5PQ_Me_ts1_no_pd_opt_sp_smd.log \\
            -as 5PQ_Me_ts1_b_no_pd_opt_sp_smd.log \\
            -hbs collidine-H_opt_sp_smd.log \\
            -bs collidine_opt_sp_smd.log \\
            -rpka 6.75 -T 298.15
    """
    # Check if we're in output file parsing mode
    output_files_provided = any(
        [
            ha_output,
            a_output,
            hb_output,
            b_output,
            ha_solv_output,
            a_solv_output,
            hb_solv_output,
            b_solv_output,
        ]
    )

    if output_files_provided:
        # Output file parsing mode - compute pKa from existing files
        from chemsmart.io.gaussian.output import Gaussian16pKaOutput

        # Validate required files for proton exchange cycle
        required_gas = [ha_output, a_output, hb_output, b_output]
        required_solv = [
            ha_solv_output,
            a_solv_output,
            hb_solv_output,
            b_solv_output,
        ]

        missing_gas = []
        missing_solv = []
        file_names = ["ha", "a", "hb", "b"]

        for i, (gas, solv, name) in enumerate(
            zip(required_gas, required_solv, file_names)
        ):
            if gas is None:
                missing_gas.append(f"-{name}/--{name}-output")
            if solv is None:
                missing_solv.append(f"-{name}s/--{name}-solv-output")

        if missing_gas or missing_solv:
            raise click.UsageError(
                f"For pKa calculation from output files, all 8 files are required.\n"
                f"Missing gas-phase files: {', '.join(missing_gas) if missing_gas else 'none'}\n"
                f"Missing solvent SP files: {', '.join(missing_solv) if missing_solv else 'none'}"
            )

        if reference_pka is None:
            raise click.UsageError(
                "When using output file parsing mode, -rpka/--reference-pka is required "
                "to specify the experimental pKa of the reference acid (HB)."
            )

        # Compute pKa using Dual-level Proton Exchange scheme
        logger.info(
            "Computing pKa from output files using Dual-level Proton Exchange scheme..."
        )
        logger.info(f"  Temperature: {temperature} K")
        logger.info(f"  Reference pKa (HB): {reference_pka}")

        # Print summary
        Gaussian16pKaOutput.print_pka_summary(
            ha_gas_file=ha_output,
            a_gas_file=a_output,
            hb_gas_file=hb_output,
            b_gas_file=b_output,
            ha_solv_file=ha_solv_output,
            a_solv_file=a_solv_output,
            hb_solv_file=hb_solv_output,
            b_solv_file=b_solv_output,
            pka_reference=reference_pka,
            temperature=temperature,
            concentration=concentration,
            cutoff_entropy_grimme=cutoff_entropy_grimme,
            cutoff_enthalpy=cutoff_enthalpy,
        )

        return

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

    # Compute thermochemistry results
    results = Gaussian16pKaOutput.compute_pka_thermochemistry(
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

    # Format output
    def format_results(results, energy_units):
        """Format thermochemistry results for display."""
        lines = []
        lines.append("=" * 78)
        lines.append("pKa Thermochemistry Extraction")
        lines.append("=" * 78)
        lines.append(f"Temperature: {results['settings']['temperature']} K")
        lines.append(
            f"Concentration: {results['settings']['concentration']} mol/L"
        )
        lines.append(
            f"Entropy cutoff (Grimme): {results['settings']['cutoff_entropy_grimme']} cm⁻¹"
        )
        lines.append(
            f"Enthalpy cutoff (Head-Gordon): {results['settings']['cutoff_enthalpy']} cm⁻¹"
        )
        lines.append(f"Energy units: {energy_units}")
        lines.append("-" * 78)
        lines.append("")
        lines.append(
            f"{'Species':<10} {'E':<20} {'qh-G(T)':<20} {'G_corr':<20}"
        )
        lines.append("-" * 78)

        for species_key in ["HA", "A", "HB", "B"]:
            if species_key in results:
                species = results[species_key]
                E = species["E"]
                qh_G = species["qh_G"]
                G_corr = qh_G - E
                lines.append(
                    f"{species['name']:<10} {E:<20.10f} {qh_G:<20.10f} {G_corr:<20.10f}"
                )

        lines.append("=" * 78)
        return "\n".join(lines)

    output_text = format_results(results, energy_units)

    # Output results
    if output is not None:
        with open(output, "w") as f:
            f.write(output_text)
        logger.info(f"Results saved to {output}")
        click.echo(f"Results saved to {output}")
    else:
        click.echo(output_text)

    return results
