"""
CLI for ORCA pKa calculations.

This module provides the CLI interface for running pKa calculations
using ORCA with a proper thermodynamic cycle:
1. Gas phase optimization + frequency for both HA and A-
2. Solution phase single point for both HA and A- at the same level of theory

Two thermodynamic cycles are supported:
- **proton exchange**: Uses a reference acid to cancel systematic errors (default)
- **direct**: Uses the absolute free energy of a proton in water

Mirrors the Gaussian pKa CLI interface for a unified user experience.
"""

import functools
import logging

import click

from chemsmart.cli.job import click_job_options
from chemsmart.cli.orca.orca import orca
from chemsmart.utils.cli import MyCommand, MyGroup

logger = logging.getLogger(__name__)


def click_orca_pka_options(f):
    """Click options specific to ORCA pKa calculations."""

    @click.option(
        "-pi",
        "--proton-index",
        type=int,
        required=False,
        help="1-based index of the proton to remove for deprotonation. Required when launching new calculations.",
    )
    @click.option(
        "-cl",
        "--color-code",
        type=int,
        default=None,
        help="CDXML colour-table index identifying the proton to remove. "
        "When a .cdx/.cdxml file is used, the proton can be identified by "
        "its colour instead of --proton-index. If omitted, the uniquely "
        "coloured hydrogen is auto-detected.",
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
        help="Path to geometry file for reference acid (HB) for proton exchange cycle.",
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
        default="CPCM",
        help="Solvation model for solution phase SP (default: CPCM).",
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
    @click.option(
        "--parallel/--no-parallel",
        default=False,
        help="Run per-species opt->SP pipelines in parallel (default: sequential).",
    )
    @functools.wraps(f)
    def wrapper(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper


def click_orca_pka_output_options(f):
    """Click options for parsing ORCA pKa output files."""
    f = click.option(
        "-rp",
        "--reference-pka",
        type=float,
        default=None,
        help="Experimental pKa of reference acid HB for proton exchange cycle.",
    )(f)
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


@orca.group("pka", cls=MyGroup, invoke_without_command=True)
@click_job_options
@click_orca_pka_options
@click_orca_pka_output_options
@click.pass_context
def pka(
    ctx,
    proton_index,
    color_code,
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
    Run ORCA pKa calculations using the dual-level proton exchange scheme.

    \b
    **Mode 1: Run new calculations**
    Provide input geometry and run gas-phase opt+freq and solution-phase SP.

    \b
    **Mode 2: Parse existing output files**
    Provide completed ORCA output files to compute pKa directly.

    \b
    Thermodynamic cycles:
    - proton exchange (default): HA + B⁻ → A⁻ + HB
    - direct: Uses absolute free energy of H⁺ in water

    \b
    Examples:
        # Run pKa with proton exchange cycle
        chemsmart run orca -f acid.xyz -c 0 -m 1 pka -pi 10 \\
            -r ref_acid.xyz -rpi 1 -rc 0 -rm 1

        # Compute pKa from existing ORCA output files
        chemsmart run orca pka -pi 1 \\
            -ha acid_opt.out -a base_opt.out \\
            -hb ref_opt.out -b ref_base_opt.out \\
            -has acid_sp.out -as base_sp.out \\
            -hbs ref_sp.out -bs ref_base_sp.out \\
            -rp 6.75 -T 298.15
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
        from chemsmart.io.orca.output import ORCApKaOutput

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

        for gas, solv, name in zip(required_gas, required_solv, file_names):
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
                "When using output file parsing mode, -rp/--reference-pka is required."
            )

        logger.info("Computing pKa from ORCA output files...")
        ORCApKaOutput.print_pka_summary(
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

    if ctx.invoked_subcommand is not None:
        ctx.obj["pka_settings"] = {
            "temperature": temperature,
            "concentration": concentration,
            "cutoff_entropy_grimme": cutoff_entropy_grimme,
            "cutoff_enthalpy": cutoff_enthalpy,
        }
        return

    if proton_index is None:
        # Try to auto-detect proton index from CDXML colour
        filename = ctx.obj.get("filename")
        if filename and filename.endswith((".cdx", ".cdxml")):
            from chemsmart.io.file import CDXFile

            cdx_file = CDXFile(filename=filename)
            try:
                proton_index = cdx_file.get_colored_proton_index(
                    color_code=color_code
                )
                logger.info(
                    f"Detected proton index {proton_index} from CDXML "
                    f"colour in {filename}."
                )
            except ValueError as exc:
                raise click.UsageError(
                    f"Could not auto-detect proton from CDXML colour: {exc}\n"
                    "Use -pi/--proton-index to specify the proton explicitly."
                )
        elif color_code is not None:
            raise click.UsageError(
                "-cl/--color-code can only be used with .cdx/.cdxml files."
            )
        else:
            raise click.UsageError(
                "-pi/--proton-index is required when launching new ORCA pKa "
                "calculations (or use a .cdxml file with a coloured proton)."
            )

    from chemsmart.jobs.orca.pka import ORCApKaJob
    from chemsmart.jobs.orca.settings import ORCApKaJobSettings

    # Validate reference acid settings
    if reference is not None:
        if thermodynamic_cycle != "proton exchange":
            raise click.UsageError(
                "Reference acid file can only be used with 'proton exchange' cycle."
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
                f"When --reference is provided, these options are required: "
                f"{', '.join(missing)}"
            )

    # Get project settings
    project_settings = ctx.obj["project_settings"]
    opt_settings = project_settings.opt_settings()

    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]
    opt_settings = opt_settings.merge(job_settings, keywords=keywords)

    molecules = ctx.obj["molecules"]
    molecule = molecules[-1]
    label = ctx.obj["label"]

    # Create ORCA pKa settings
    pka_settings = ORCApKaJobSettings(
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
        temperature=temperature,
        concentration=concentration,
        cutoff_entropy_grimme=cutoff_entropy_grimme,
        cutoff_enthalpy=cutoff_enthalpy,
        # Inherit level of theory from project opt settings
        charge=opt_settings.charge,
        multiplicity=opt_settings.multiplicity,
        functional=opt_settings.functional,
        basis=opt_settings.basis,
        ab_initio=opt_settings.ab_initio,
        dispersion=opt_settings.dispersion,
        aux_basis=getattr(opt_settings, "aux_basis", None),
        defgrid=opt_settings.defgrid,
        semiempirical=opt_settings.semiempirical,
        additional_route_parameters=opt_settings.additional_route_parameters,
        gen_genecp_file=opt_settings.gen_genecp_file,
        heavy_elements=opt_settings.heavy_elements,
        heavy_elements_basis=opt_settings.heavy_elements_basis,
        light_elements_basis=opt_settings.light_elements_basis,
    )

    if pka_settings.charge is None:
        raise click.UsageError(
            "Charge must be specified via -c/--charge option"
        )
    if pka_settings.multiplicity is None:
        raise click.UsageError(
            "Multiplicity must be specified via -m/--multiplicity option"
        )

    logger.info(f"ORCA pKa job settings: {pka_settings.__dict__}")
    logger.info(f"Proton index to remove: {proton_index}")
    logger.info(f"Thermodynamic cycle: {thermodynamic_cycle}")
    logger.info(
        f"Protonated form (HA): charge={pka_settings.charge}, "
        f"mult={pka_settings.multiplicity}"
    )

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

    if thermodynamic_cycle == "proton exchange" and reference is not None:
        logger.info(f"Reference acid file: {reference}")
        logger.info(f"Reference proton index: {reference_proton_index}")
        logger.info(
            f"Reference acid (HB): charge={reference_charge}, "
            f"mult={reference_multiplicity}"
        )

    logger.info(
        f"Gas phase optimization: {pka_settings.functional}/{pka_settings.basis}"
    )
    logger.info(
        f"Solution phase SP: {pka_settings.functional}/{pka_settings.basis} "
        f"with {solvent_model}({solvent_id})"
    )

    # Avoid duplicating the _pka suffix if the base label already has it
    job_label = label if label.endswith("_pka") else f"{label}_pka"

    return ORCApKaJob(
        molecule=molecule,
        settings=pka_settings,
        label=job_label,
        skip_completed=skip_completed,
        parallel=kwargs.get("parallel", False),
        **kwargs,
    )


def click_orca_pka_thermo_options(f):
    """Click options for ORCA pKa thermochemistry extraction."""

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
        help="Temperature in Kelvin (default: 298.15 K).",
    )
    @click.option(
        "-conc",
        "--concentration",
        type=float,
        default=1.0,
        help="Concentration in mol/L (default: 1.0 mol/L).",
    )
    @click.option(
        "-csg",
        "--cutoff-entropy-grimme",
        type=float,
        default=100.0,
        help="Cutoff for entropy (cm^-1) using Grimme's method (default: 100).",
    )
    @click.option(
        "-ch",
        "--cutoff-enthalpy",
        type=float,
        default=100.0,
        help="Cutoff for enthalpy (cm^-1) using Head-Gordon's method (default: 100).",
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
@click_orca_pka_thermo_options
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
    Extract thermochemistry (E and qh-G(T)) from ORCA optimization output files.

    \b
    Examples:
        chemsmart run orca pka thermo -ha acid_opt.out -a base_opt.out -T 298.15

        chemsmart run orca pka thermo \\
            -ha acid.out -a base.out -hb ref.out -b ref_base.out \\
            -T 298.15 -csg 100 -ch 100
    """
    from chemsmart.io.orca.output import ORCApKaOutput

    if all(f is None for f in [ha_file, a_file, hb_file, b_file]):
        raise click.UsageError(
            "At least one output file must be provided. Use -ha, -a, -hb, or -b."
        )

    results = ORCApKaOutput.compute_pka_thermochemistry(
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

    def format_results(results, energy_units):
        lines = []
        lines.append("=" * 78)
        lines.append("pKa Thermochemistry Extraction (ORCA)")
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

    if output is not None:
        with open(output, "w") as f:
            f.write(output_text)
        logger.info(f"Results saved to {output}")
        click.echo(f"Results saved to {output}")
    else:
        click.echo(output_text)

    return results
