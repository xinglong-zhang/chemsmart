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
        "-i",
        "--input-table",
        is_flag=True,
        default=False,
        required=False,
        help="Use table-driven batch mode. In this mode, the file provided by "
        "the parent Gaussian -f/--filename option is interpreted as the input "
        "table (.txt/.csv) with columns: filepath, proton_index, charge, multiplicity.",
    )
    @click.option(
        "-O",
        "--output-table",
        type=click.Path(exists=True),
        default=None,
        required=False,
        help="Compute pKa from a table of precomputed output files. "
        "The table must contain columns: basename, ha_gas, a_gas, hb_gas, "
        "b_gas, ha_sp, a_sp, hb_sp, b_sp, pka_ref. "
        "Blank reference-acid cells are filled from the previous row.",
    )
    @click.option(
        "--output-results",
        type=click.Path(),
        default=None,
        help="Path to write computed pKa results table (.csv or .txt). "
        "Used with -O/--output-table. If omitted, results are printed to stdout.",
    )
    @click.option(
        "-pi",
        "--proton-index",
        type=int,
        required=False,
        help="1-based index of the proton to remove for deprotonation. Required for new calculations only.",
    )
    @click.option(
        "-cl",
        "--color-code",
        type=int,
        default=None,
        help="CDXML colour-table index identifying the proton to remove. "
        "When a .cdxml file is used, the proton can be identified by its "
        "colour instead of --proton-index. If omitted, the uniquely "
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
        help="Path to geometry file for reference acid (HB) for proton exchange cycle. "
        "When provided, optimization and SP calculations will also be run for HB and B-.",
    )
    @click.option(
        "-rpi",
        "--reference-proton-index",
        type=int,
        default=None,
        help="1-based index of the proton to remove from reference acid (HB). "
        "Required when --reference is provided (unless the reference file "
        "is a .cdxml with a coloured proton).",
    )
    @click.option(
        "-rcl",
        "--reference-color-code",
        type=int,
        default=None,
        help="CDXML colour-table index identifying the proton in the "
        "reference acid file. Only used when --reference is a .cdx/.cdxml.",
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
    @click.option(
        "--parallel/--no-parallel",
        default=False,
        help="Run per-species opt->SP pipelines in parallel (default: sequential).",
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
    input_table,
    output_table,
    output_results,
    proton_index,
    color_code,
    thermodynamic_cycle,
    reference,
    reference_proton_index,
    reference_color_code,
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
    parallel,
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

    Four execution modes are available:

    \b
    **Mode 1: Single molecule calculation**
    Provide input geometry file with -f and run optimization + SP calculations.

    \b
    **Mode 2: Table-driven batch calculation (job submission)**
    Enable with -o/--input-table. The table file is taken from the parent
    Gaussian -f/--filename option.
    Table format (4 columns, whitespace or comma-delimited):
        filepath    proton_index    charge    multiplicity
    In this mode, molecule geometry options are read from table rows. For
    proton exchange cycle, reference acid options (-r, -rpi, -rc, -rm)
    are still mandatory.

    \b
    **Mode 3: Parse existing output files (single system)**
    Provide completed Gaussian output files to compute pKa directly using
    the Dual-level Proton Exchange scheme. All energies in Hartree (au),
    except ΔG_soln which is converted to kcal/mol for the pKa formula.

    \b
    **Mode 4: Batch pKa from output table (post-processing)**
    Provide -O/--output-table with a table of precomputed output file paths.
    Table columns: basename, ha_gas, a_gas, hb_gas, b_gas, ha_sp, a_sp,
    hb_sp, b_sp, pka_ref.
    Blank reference-acid cells are carried forward from the previous row.
    Use --output-results to write a results table with appended pka column.

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
        # Run pKa job with proton exchange cycle (single molecule)
        chemsmart run gaussian -f acetic_acid.xyz -c 0 -m 1 pka -pi 10 \\
            -t "proton exchange" -r water.xyz -rpi 1 -rc 0 -rm 1

        # Run pKa job with direct cycle
        chemsmart run gaussian -f acetic_acid.xyz -c 0 -m 1 pka -pi 10 -t direct

        # Table-driven batch mode (proton exchange)
        chemsmart run gaussian -p myproject -f molecules.txt pka -i \
            -t "proton exchange" -r ref_acid.xyz -rpi 5 -rc 0 -rm 1

        # Table-driven batch mode (direct cycle)
        chemsmart run gaussian -p myproject -f molecules.csv pka -i -t direct

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

        # Batch pKa from output table (post-processing)
        chemsmart run gaussian pka -O outputs.csv \\
            --output-results results.csv -T 298.15
    """
    # =========================================================================
    # Output-table post-processing mode (Mode 4)
    # =========================================================================
    if output_table is not None:
        return _run_pka_from_output_table(
            output_table=output_table,
            output_results=output_results,
            temperature=temperature,
            concentration=concentration,
            cutoff_entropy_grimme=cutoff_entropy_grimme,
            cutoff_enthalpy=cutoff_enthalpy,
            program="gaussian",
        )

    # =========================================================================
    # Table-driven execution mode
    # =========================================================================
    if input_table:
        input_table_path = ctx.obj.get("filename")
        if not input_table_path:
            raise click.UsageError(
                "Table mode (-i/--input-table) requires parent Gaussian -f/--filename "
                "to specify the table file path."
            )
        return _run_pka_from_table(
            ctx=ctx,
            input_table=input_table_path,
            thermodynamic_cycle=thermodynamic_cycle,
            reference=reference,
            reference_proton_index=reference_proton_index,
            reference_color_code=reference_color_code,
            reference_charge=reference_charge,
            reference_multiplicity=reference_multiplicity,
            reference_conjugate_base_charge=reference_conjugate_base_charge,
            reference_conjugate_base_multiplicity=reference_conjugate_base_multiplicity,
            delta_g_proton=delta_g_proton,
            solvent_model=solvent_model,
            solvent_id=solvent_id,
            temperature=temperature,
            concentration=concentration,
            cutoff_entropy_grimme=cutoff_entropy_grimme,
            cutoff_enthalpy=cutoff_enthalpy,
            skip_completed=skip_completed,
            parallel=parallel,
            **kwargs,
        )

    # =========================================================================
    # Output file parsing mode
    # =========================================================================
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

    if proton_index is None:
        # Try to auto-detect proton index from CDXML colour
        filename = ctx.obj.get("filename")
        if filename and filename.endswith((".cdx", ".cdxml")):
            from chemsmart.io.file import PKaCDXFile

            cdx_file = PKaCDXFile(filename=filename)
            try:
                pka_mols = cdx_file.get_pka_molecules(color_code=color_code)
            except ValueError as exc:
                raise click.UsageError(
                    f"Could not auto-detect proton from CDXML colour: {exc}\n"
                    "Use -pi/--proton-index to specify the proton explicitly."
                )

            if len(pka_mols) > 1:
                # Multi-molecule ChemDraw file → one pKa job per molecule
                logger.info(
                    f"Detected {len(pka_mols)} molecules with per-fragment "
                    f"proton auto-detection in {filename}."
                )
                return _run_pka_from_cdxml_molecules(
                    ctx=ctx,
                    pka_molecules=pka_mols,
                    thermodynamic_cycle=thermodynamic_cycle,
                    reference=reference,
                    reference_proton_index=reference_proton_index,
                    reference_color_code=reference_color_code,
                    reference_charge=reference_charge,
                    reference_multiplicity=reference_multiplicity,
                    reference_conjugate_base_charge=reference_conjugate_base_charge,
                    reference_conjugate_base_multiplicity=reference_conjugate_base_multiplicity,
                    delta_g_proton=delta_g_proton,
                    conjugate_base_charge=conjugate_base_charge,
                    conjugate_base_multiplicity=conjugate_base_multiplicity,
                    solvent_model=solvent_model,
                    solvent_id=solvent_id,
                    temperature=temperature,
                    concentration=concentration,
                    cutoff_entropy_grimme=cutoff_entropy_grimme,
                    cutoff_enthalpy=cutoff_enthalpy,
                    skip_completed=skip_completed,
                    parallel=parallel,
                    **kwargs,
                )

            # Single molecule – use its proton_index and fall through
            proton_index = pka_mols[0].proton_index
            logger.info(
                f"Detected proton index {proton_index} from CDXML "
                f"colour in {filename}."
            )
        elif color_code is not None:
            raise click.UsageError(
                "-cl/--color-code can only be used with .cdx/.cdxml files."
            )
        else:
            raise click.UsageError(
                "-pi/--proton-index is required when launching new pKa "
                "calculations (or use a .cdxml file with a coloured proton)."
            )

    from chemsmart.jobs.gaussian.pka import GaussianpKaJob
    from chemsmart.jobs.gaussian.settings import GaussianpKaJobSettings

    # Validate reference acid settings if provided
    if reference is not None:
        if thermodynamic_cycle != "proton exchange":
            raise click.UsageError(
                "Reference acid file can only be used with 'proton exchange' cycle. "
                "Use -t 'proton exchange' or remove the -r option."
            )

        # Auto-detect reference proton index from CDXML colour
        if reference_proton_index is None:
            if reference.endswith((".cdx", ".cdxml")):
                from chemsmart.io.file import PKaCDXFile

                ref_cdx = PKaCDXFile(filename=reference)
                try:
                    ref_pka_mol = ref_cdx.get_pka_molecule(
                        color_code=reference_color_code
                    )
                    reference_proton_index = ref_pka_mol.proton_index
                    logger.info(
                        f"Detected reference proton index "
                        f"{reference_proton_index} from CDXML colour "
                        f"in {reference}."
                    )
                except ValueError as exc:
                    raise click.UsageError(
                        f"Could not auto-detect reference proton from "
                        f"CDXML colour: {exc}\n"
                        "Use -rpi/--reference-proton-index to specify "
                        "the proton explicitly."
                    )
            elif reference_color_code is not None:
                raise click.UsageError(
                    "-rcl/--reference-color-code can only be used when "
                    "--reference is a .cdx/.cdxml file."
                )

        missing = []
        if reference_proton_index is None:
            missing.append(
                "-rpi/--reference-proton-index (or use a .cdxml reference "
                "with a coloured proton)"
            )
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
    molecule_indices = ctx.obj.get("molecule_indices")

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

    # Create and return pKa job(s)
    if len(molecules) > 1 and molecule_indices:
        logger.info(f"Creating {len(molecules)} pKa jobs")
        jobs = []
        for molecule, idx in zip(molecules, molecule_indices):
            molecule_label = f"{label}_idx{idx}"
            jobs.append(
                GaussianpKaJob(
                    molecule=molecule,
                    settings=pka_settings,
                    label=molecule_label,
                    jobrunner=jobrunner,
                    skip_completed=skip_completed,
                    parallel=parallel,
                    **kwargs,
                )
            )
        return jobs

    molecule = molecules[-1]
    return GaussianpKaJob(
        molecule=molecule,
        settings=pka_settings,
        label=label,
        jobrunner=jobrunner,
        skip_completed=skip_completed,
        parallel=parallel,
        **kwargs,
    )


def _run_pka_from_cdxml_molecules(
    ctx,
    pka_molecules,
    thermodynamic_cycle,
    reference,
    reference_proton_index,
    reference_color_code,
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
    parallel,
    **kwargs,
):
    """Create one GaussianpKaJob per PKaMolecule from a multi-fragment CDXML.

    Each :class:`PKaMolecule` carries its own ``proton_index`` resolved
    during per-fragment colour detection.  This helper builds the
    corresponding :class:`GaussianpKaJobSettings` for each molecule and
    returns the list of jobs.

    Args:
        ctx: Click context (must contain ``jobrunner``, ``project_settings``,
            ``job_settings``, ``keywords``, ``label``).
        pka_molecules: List of :class:`PKaMolecule` instances.
        (remaining args): Same as the ``pka()`` CLI handler.

    Returns:
        list[GaussianpKaJob]: One job per molecule.
    """
    import os

    from chemsmart.jobs.gaussian.pka import GaussianpKaJob
    from chemsmart.jobs.gaussian.settings import GaussianpKaJobSettings

    # ---- resolve reference acid (shared across all molecules) ----
    if thermodynamic_cycle == "proton exchange" and reference is not None:
        if reference_proton_index is None:
            if reference.endswith((".cdx", ".cdxml")):
                from chemsmart.io.file import PKaCDXFile

                try:
                    ref_cdx = PKaCDXFile(filename=reference)
                    ref_pka_mol = ref_cdx.get_pka_molecule(
                        color_code=reference_color_code
                    )
                    reference_proton_index = ref_pka_mol.proton_index
                    logger.info(
                        f"Detected reference proton index "
                        f"{reference_proton_index} from CDXML colour."
                    )
                except ValueError as exc:
                    raise click.UsageError(
                        f"Could not auto-detect reference proton: {exc}\n"
                        "Use -rpi/--reference-proton-index explicitly."
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
                "When --reference is provided, the following options "
                "are required: " + ", ".join(missing)
            )

    # ---- project / job settings ----
    jobrunner = ctx.obj["jobrunner"]
    project_settings = ctx.obj["project_settings"]
    opt_settings = project_settings.opt_settings()
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]
    opt_settings = opt_settings.merge(job_settings, keywords=keywords)

    filename = ctx.obj.get("filename", "")
    base_name = os.path.splitext(os.path.basename(filename))[0]

    # ---- one job per molecule ----
    jobs = []
    for idx, pka_mol in enumerate(pka_molecules, start=1):
        mol_label = f"{base_name}_frag{idx}_pka"

        pka_settings = GaussianpKaJobSettings(
            proton_index=pka_mol.proton_index,
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
            charge=opt_settings.charge,
            multiplicity=opt_settings.multiplicity,
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

        logger.info(
            f"Creating pKa job for fragment {idx}: "
            f"proton_index={pka_mol.proton_index}, "
            f"label={mol_label}"
        )

        job = GaussianpKaJob(
            molecule=pka_mol,
            settings=pka_settings,
            label=mol_label,
            jobrunner=jobrunner,
            skip_completed=skip_completed,
            parallel=parallel,
            **kwargs,
        )
        jobs.append(job)

    logger.info(f"Created {len(jobs)} pKa jobs from multi-fragment CDXML file")
    return jobs


def _run_pka_from_table(
    ctx,
    input_table,
    thermodynamic_cycle,
    reference,
    reference_proton_index,
    reference_color_code,
    reference_charge,
    reference_multiplicity,
    reference_conjugate_base_charge,
    reference_conjugate_base_multiplicity,
    delta_g_proton,
    solvent_model,
    solvent_id,
    temperature,
    concentration,
    cutoff_entropy_grimme,
    cutoff_enthalpy,
    skip_completed,
    parallel,
    **kwargs,
):
    """Run pKa jobs from a table file.

    Table format (4 columns):
        filepath    proton_index    charge    multiplicity

    For proton exchange cycle, reference acid options are still required.
    """
    from pathlib import Path

    from chemsmart.io.molecules.structure import Molecule
    from chemsmart.jobs.gaussian.pka import GaussianpKaJob
    from chemsmart.jobs.gaussian.settings import GaussianpKaJobSettings
    from chemsmart.utils.utils import (
        parse_pka_table,
        validate_pka_table_entries,
    )

    # Parse and validate table
    logger.info(f"Reading pKa jobs from table: {input_table}")
    try:
        entries = parse_pka_table(input_table)
        validate_pka_table_entries(entries, check_file_exists=True)
    except (FileNotFoundError, ValueError) as e:
        raise click.UsageError(str(e))

    logger.info(f"Found {len(entries)} entries in table")

    # Validate reference acid options for proton exchange cycle
    if thermodynamic_cycle == "proton exchange":
        missing = []
        if reference is None:
            missing.append("-r/--reference")
        if reference_proton_index is None and reference is not None:
            # Try auto-detect from CDXML
            if reference.endswith((".cdx", ".cdxml")):
                from chemsmart.io.file import PKaCDXFile

                try:
                    ref_cdx = PKaCDXFile(filename=reference)
                    ref_pka_mol = ref_cdx.get_pka_molecule(
                        color_code=reference_color_code
                    )
                    reference_proton_index = ref_pka_mol.proton_index
                    logger.info(
                        f"Detected reference proton index "
                        f"{reference_proton_index} from CDXML colour."
                    )
                except ValueError as exc:
                    missing.append(
                        f"-rpi/--reference-proton-index (auto-detect failed: {exc})"
                    )
            else:
                missing.append("-rpi/--reference-proton-index")
        if reference_charge is None:
            missing.append("-rc/--reference-charge")
        if reference_multiplicity is None:
            missing.append("-rm/--reference-multiplicity")

        if missing:
            raise click.UsageError(
                "For proton exchange cycle with table input, the following "
                "reference acid options are required:\n  "
                + "\n  ".join(missing)
            )

    # Get jobrunner and project settings
    jobrunner = ctx.obj["jobrunner"]
    project_settings = ctx.obj["project_settings"]
    opt_settings = project_settings.opt_settings()

    # Merge with CLI keywords
    job_settings = ctx.obj.get("job_settings")
    keywords = ctx.obj.get("keywords", {})
    if job_settings:
        opt_settings = opt_settings.merge(job_settings, keywords=keywords)

    # Create jobs for each table entry
    jobs = []
    for entry in entries:
        row = entry.to_kwargs(drop_none=False)

        # Load molecule from file (supports aliasing via entry class)
        filepath = entry.get("filepath") or entry.get("path") or entry.filepath
        molecule = Molecule.from_filepath(filepath)

        # Derive label from filename
        label = Path(filepath).stem

        # Create pKa settings for this entry
        # Charge/multiplicity from table, level of theory from project
        pka_settings = GaussianpKaJobSettings(
            proton_index=int(entry.proton_index),
            thermodynamic_cycle=thermodynamic_cycle,
            reference_file=reference,
            reference_proton_index=reference_proton_index,
            reference_charge=reference_charge,
            reference_multiplicity=reference_multiplicity,
            reference_conjugate_base_charge=reference_conjugate_base_charge,
            reference_conjugate_base_multiplicity=reference_conjugate_base_multiplicity,
            delta_G_proton=delta_g_proton,
            conjugate_base_charge=None,
            conjugate_base_multiplicity=None,
            solvent_model=solvent_model,
            solvent_id=solvent_id,
            temperature=temperature,
            concentration=concentration,
            cutoff_entropy_grimme=cutoff_entropy_grimme,
            cutoff_enthalpy=cutoff_enthalpy,
            charge=int(entry.charge),
            multiplicity=int(entry.multiplicity),
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

        logger.info(
            f"Creating pKa job for {filepath}: "
            f"proton_index={entry.proton_index}, "
            f"charge={entry.charge}, mult={entry.multiplicity}, "
            f"row_fields={sorted(row.keys())}"
        )

        job = GaussianpKaJob(
            molecule=molecule,
            settings=pka_settings,
            label=label,
            jobrunner=jobrunner,
            skip_completed=skip_completed,
            parallel=parallel,
            **kwargs,
        )
        jobs.append(job)

    logger.info(f"Created {len(jobs)} pKa jobs from table")
    return jobs


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


def _run_pka_from_output_table(
    output_table,
    output_results,
    temperature,
    concentration,
    cutoff_entropy_grimme,
    cutoff_enthalpy,
    program="gaussian",
):
    """Compute pKa values from a table of precomputed output files.

    This shared helper implements Mode 4 (batch pKa post-processing).
    It is called from both Gaussian and ORCA ``pka`` CLI handlers.

    Args:
        output_table: Path to the output-table file.
        output_results: Optional path to write results table. If *None*,
            results are printed to stdout.
        temperature: Temperature in Kelvin.
        concentration: Concentration in mol/L.
        cutoff_entropy_grimme: Grimme quasi-RRHO entropy cutoff (cm⁻¹).
        cutoff_enthalpy: Head-Gordon enthalpy cutoff (cm⁻¹).
        program: ``"gaussian"`` or ``"orca"`` – selects the output class.
    """
    from chemsmart.utils.utils import (
        compute_pka_from_output_table,
        export_pka_results_table,
        parse_pka_output_table,
        resolve_pka_output_references,
    )

    # 1. Parse
    logger.info(f"Reading pKa output table: {output_table}")
    try:
        entries = parse_pka_output_table(output_table)
    except (FileNotFoundError, ValueError) as e:
        raise click.UsageError(str(e))

    logger.info(f"Found {len(entries)} entries in output table")

    # 2. Resolve blank reference-acid cells
    try:
        resolve_pka_output_references(entries)
    except ValueError as e:
        raise click.UsageError(str(e))

    # 3. Validate all file paths exist
    all_errors = []
    for entry in entries:
        try:
            entry.validate(check_file_exists=True)
        except ValueError as e:
            all_errors.append(str(e))
    if all_errors:
        raise click.UsageError(
            "Output table validation failed:\n" + "\n".join(all_errors)
        )

    # 4. Select output class
    if program == "gaussian":
        from chemsmart.io.gaussian.output import (
            Gaussian16pKaOutput as OutputCls,
        )
    elif program == "orca":
        from chemsmart.io.orca.output import ORCApKaOutput as OutputCls
    else:
        raise ValueError(f"Unknown program: {program}")

    # 5. Compute
    logger.info(
        f"Computing pKa for {len(entries)} systems "
        f"(T={temperature}K, program={program})"
    )
    results = compute_pka_from_output_table(
        entries=entries,
        output_cls=OutputCls,
        temperature=temperature,
        concentration=concentration,
        cutoff_entropy_grimme=cutoff_entropy_grimme,
        cutoff_enthalpy=cutoff_enthalpy,
    )

    # 6. Export / display
    if output_results is not None:
        export_pka_results_table(entries, results, output_results)
        click.echo(f"pKa results written to {output_results}")
    else:
        # Print summary to stdout
        click.echo("=" * 78)
        click.echo("Batch pKa Results (Dual-level Proton Exchange)")
        click.echo("=" * 78)
        click.echo(f"Temperature: {temperature} K")
        click.echo(f"{'basename':<30} {'pKa':>10} {'ΔG_soln (kcal/mol)':>20}")
        click.echo("-" * 78)
        for entry, result in zip(entries, results):
            click.echo(
                f"{entry['basename']:<30} "
                f"{result['pKa']:>10.2f} "
                f"{result['delta_G_soln_kcal_mol']:>20.4f}"
            )
        click.echo("=" * 78)

    return results
