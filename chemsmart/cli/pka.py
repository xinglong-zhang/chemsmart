"""
Backend-independent pKa output analysis.

Registered as ``chemsmart run pka`` — strictly post-processing,
never invokes Gaussian, ORCA, or any other QC backend.

Subcommands
-----------
analyze        Compute pKa from 8 existing output files.
batch-analyze  Batch pKa from a table of output file paths.
"""

import functools
import logging
import os

import click

from chemsmart.utils.cli import MyCommand, MyGroup
from chemsmart.utils.io import get_program_type_from_file

logger = logging.getLogger(__name__)

_EXT_MAP = {
    "gaussian": [".log", ".out"],
    "orca": [".out", ".log"],
}


# ── helpers moved from cli/pka_helpers.py ──────────────────────────────


def click_pka_shared_options(f):
    @click.option(
        "-t",
        "--thermodynamic-cycle",
        type=click.Choice(["direct", "proton exchange"]),
        default="proton exchange",
        help=(
            "Thermodynamic cycle type. 'proton exchange' uses a reference acid (default). 'direct' uses absolute free energy of H+ in water."
        ),
    )
    @click.option(
        "-r",
        "--reference",
        type=click.Path(exists=True),
        default=None,
        help=(
            "Path to geometry file for reference acid (HRef) for proton exchange cycle."
        ),
    )
    @click.option(
        "-rpi",
        "--reference-proton-index",
        type=int,
        default=None,
        help=(
            "1-based index of the proton to remove from reference acid (HRef). Required when --reference is provided."
        ),
    )
    @click.option(
        "-rcc",
        "--reference-color-code",
        type=int,
        default=None,
        help=(
            "CDXML colour-table index identifying the proton in the reference acid file."
        ),
    )
    @click.option(
        "-rc",
        "--reference-charge",
        type=int,
        default=None,
        help="Charge of the reference acid (HRef).",
    )
    @click.option(
        "-rm",
        "--reference-multiplicity",
        type=int,
        default=None,
        help="Multiplicity of the reference acid (HRef).",
    )
    @click.option(
        "--reference-conjugate-base-charge",
        type=int,
        default=None,
        help=(
            "Charge of the reference conjugate base (Ref-). Defaults to (reference_charge - 1)."
        ),
    )
    @click.option(
        "--reference-conjugate-base-multiplicity",
        type=int,
        default=None,
        help=(
            "Multiplicity of the reference conjugate base (Ref-). Defaults to reference_multiplicity."
        ),
    )
    @click.option(
        "-dG",
        "--delta-g-proton",
        type=float,
        default=-265.9,
        help=(
            "Absolute free energy of H+ in water (kcal/mol) for direct cycle. Default: -265.9 kcal/mol (Tissandier et al., 1998)."
        ),
    )
    @click.option(
        "--conjugate-base-charge",
        type=int,
        default=None,
        help="Charge of the conjugate base (A-). Defaults to (charge - 1).",
    )
    @click.option(
        "--conjugate-base-multiplicity",
        type=int,
        default=None,
        help=(
            "Multiplicity of the conjugate base (A-). Defaults to multiplicity."
        ),
    )
    @click.option(
        "-sm",
        "--solvent-model",
        type=str,
        default=None,
        help="Solvation model for solution phase SP.",
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
        help=(
            "Cutoff frequency for entropy (cm^-1) using Grimme's quasi-RRHO method (default: 100)."
        ),
    )
    @click.option(
        "-ch",
        "--cutoff-enthalpy",
        type=float,
        default=100.0,
        help=(
            "Cutoff frequency for enthalpy (cm^-1) using Head-Gordon's method (default: 100)."
        ),
    )
    @functools.wraps(f)
    def wrapper(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper


def click_pka_submit_options(f):
    @click.option(
        "-pi",
        "--proton-index",
        type=int,
        required=False,
        help=(
            "1-based index of the proton to remove for deprotonation. Required unless a .cdxml file with a coloured proton is used."
        ),
    )
    @click.option(
        "-cc",
        "--color-code",
        type=int,
        default=None,
        help=(
            "CDXML colour-table index identifying the proton to remove. If omitted, the uniquely coloured hydrogen is auto-detected."
        ),
    )
    @click.option(
        "--parallel/--no-parallel",
        default=False,
        help="Run per-species opt→SP pipelines in parallel.",
    )
    @functools.wraps(f)
    def wrapper(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper


def click_pka_analyze_options(f):
    f = click.option(
        "-rp",
        "--reference-pka",
        type=float,
        default=None,
        help=(
            "Experimental pKa of reference acid HRef. Required for proton exchange analysis."
        ),
    )(f)
    f = click.option(
        "-rs",
        "--ref-solv",
        "ref_solv",
        type=click.Path(exists=True),
        default=None,
        help="Ref- solvent SP output file.",
    )(f)
    f = click.option(
        "-hrs",
        "--href-solv",
        "href_solv",
        type=click.Path(exists=True),
        default=None,
        help="HRef solvent SP output file.",
    )(f)
    f = click.option(
        "-as",
        "--a-solv",
        "a_solv",
        type=click.Path(exists=True),
        default=None,
        help="A- solvent SP output file.",
    )(f)
    f = click.option(
        "-has",
        "--ha-solv",
        "ha_solv",
        type=click.Path(exists=True),
        default=None,
        help="HA solvent SP output file.",
    )(f)
    f = click.option(
        "-r",
        "--ref",
        "ref",
        type=click.Path(exists=True),
        default=None,
        help="Ref- gas-phase opt+freq output file.",
    )(f)
    f = click.option(
        "-hr",
        "--href",
        "href",
        type=click.Path(exists=True),
        default=None,
        help="HRef gas-phase opt+freq output file.",
    )(f)
    f = click.option(
        "-a",
        "--a",
        "a",
        type=click.Path(exists=True),
        default=None,
        help="A- gas-phase opt+freq output file.",
    )(f)
    f = click.option(
        "-ha",
        "--ha",
        "ha",
        type=click.Path(exists=True),
        default=None,
        help="HA gas-phase opt+freq output file.",
    )(f)
    return f


def resolve_proton_index(filename, proton_index, color_code):
    if proton_index is not None:
        return proton_index, None

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
            logger.info(
                f"Detected {len(pka_mols)} molecules with per-fragment proton auto-detection in {filename}."
            )
            return None, pka_mols

        proton_index = pka_mols[0].proton_index
        logger.info(
            f"Detected proton index {proton_index} from CDXML colour in {filename}."
        )
        return proton_index, None

    if color_code is not None:
        raise click.UsageError(
            "-cc/--color-code can only be used with .cdx/.cdxml files."
        )

    raise click.UsageError(
        "-pi/--proton-index is required when launching new pKa "
        "calculations (or use a .cdxml file with a coloured proton)."
    )


def resolve_reference_proton(
    reference, reference_proton_index, reference_color_code
):
    if reference_proton_index is not None:
        return reference_proton_index
    if reference is None:
        return None

    if reference.endswith((".cdx", ".cdxml")):
        from chemsmart.io.file import PKaCDXFile

        ref_cdx = PKaCDXFile(filename=reference)
        try:
            ref_pka_mol = ref_cdx.get_pka_molecule(
                color_code=reference_color_code
            )
            reference_proton_index = ref_pka_mol.proton_index
            logger.info(
                f"Detected reference proton index {reference_proton_index} from CDXML colour in {reference}."
            )
            return reference_proton_index
        except ValueError as exc:
            raise click.UsageError(
                f"Could not auto-detect reference proton from CDXML colour: {exc}\n"
                "Use -rpi/--reference-proton-index to specify the proton explicitly."
            )
    elif reference_color_code is not None:
        raise click.UsageError(
            "-rcc/--reference-color-code can only be used when --reference is a .cdx/.cdxml file."
        )

    return None


def validate_reference_options(shared):
    reference = shared["reference"]
    if reference is None:
        return

    if shared["thermodynamic_cycle"] != "proton exchange":
        raise click.UsageError(
            "Reference acid file can only be used with 'proton exchange' "
            "cycle. Use -t 'proton exchange' or remove the -r option."
        )

    shared["reference_proton_index"] = resolve_reference_proton(
        reference,
        shared["reference_proton_index"],
        shared["reference_color_code"],
    )

    missing = []
    if shared["reference_proton_index"] is None:
        missing.append(
            "-rpi/--reference-proton-index (or use a .cdxml reference with a coloured proton)"
        )
    if shared["reference_charge"] is None:
        missing.append("-rc/--reference-charge")
    if shared["reference_multiplicity"] is None:
        missing.append("-rm/--reference-multiplicity")
    if missing:
        raise click.UsageError(
            "When --reference is provided, the following options are "
            "required: " + ", ".join(missing)
        )


def run_pka_from_output_table(
    output_table,
    output_results,
    temperature,
    concentration,
    cutoff_entropy_grimme,
    cutoff_enthalpy,
    program="gaussian",
):
    from chemsmart.utils.utils import (
        compute_pka_from_output_table,
        export_pka_results_table,
        parse_pka_output_table,
        resolve_pka_output_references,
    )

    logger.info(f"Reading pKa output table: {output_table}")
    try:
        entries = parse_pka_output_table(output_table)
    except (FileNotFoundError, ValueError) as e:
        raise click.UsageError(str(e))

    logger.info(f"Found {len(entries)} entries in output table")

    try:
        resolve_pka_output_references(entries)
    except ValueError as e:
        raise click.UsageError(str(e))

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

    if program == "gaussian":
        from chemsmart.io.gaussian.output import (
            Gaussian16pKaOutput as OutputCls,
        )
    elif program == "orca":
        from chemsmart.io.orca.output import ORCApKaOutput as OutputCls
    else:
        raise ValueError(f"Unknown program: {program}")

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

    if output_results is not None:
        export_pka_results_table(entries, results, output_results)
        click.echo(f"pKa results written to {output_results}")
    else:
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


def validate_analyze_files(
    ha, a, href, ref, ha_solv, a_solv, href_solv, ref_solv, reference_pka
):
    required_gas = [ha, a, href, ref]
    required_solv = [ha_solv, a_solv, href_solv, ref_solv]
    file_names = ["ha", "a", "href", "ref"]

    missing_gas = []
    missing_solv = []
    for gas, solv, name in zip(required_gas, required_solv, file_names):
        if gas is None:
            missing_gas.append(f"--{name}")
        if solv is None:
            missing_solv.append(f"--{name}-solv")

    if missing_gas or missing_solv:
        raise click.UsageError(
            "For pKa analysis all 8 output files are required.\n"
            f"Missing gas-phase: "
            f"{', '.join(missing_gas) if missing_gas else 'none'}\n"
            f"Missing solvent SP: "
            f"{', '.join(missing_solv) if missing_solv else 'none'}"
        )

    if reference_pka is None:
        raise click.UsageError(
            "-rp/--reference-pka is required for output-file analysis."
        )


def _detect_program_from_outputs(filepaths):
    detected = {}
    programs = set()
    for fp in filepaths:
        p = get_program_type_from_file(fp)
        detected[fp] = p
        if p != "unknown":
            programs.add(p)

    if not programs:
        raise click.UsageError(
            "Could not detect output-file program type from supplied "
            "files.  Supported: Gaussian and ORCA output files."
        )
    if len(programs) > 1:
        pairs = ", ".join(f"{k}: {v}" for k, v in detected.items())
        raise click.UsageError(
            "Supplied files contain mixed program types.  "
            "Use outputs from a single QC program.\n"
            f"Detected: {pairs}"
        )
    program = next(iter(programs))
    if program not in {"gaussian", "orca"}:
        raise click.UsageError(
            f"Detected unsupported program '{program}'.  "
            "Only Gaussian and ORCA are supported."
        )
    return program


def _auto_discover_pka_files(ha_gas_path, href_gas_path, program=None):
    """Infer companion output paths from HA and HRef gas-phase paths.

    Naming convention produced by the pKa job classes::

        <basename>.<ext>        gas-phase opt+freq
        <basename>_cb.<ext>     conjugate base gas-phase
        <basename>_sp.<ext>     solvent single-point
        <basename>_cb_sp.<ext>  conjugate base solvent SP
    """
    if program is None:
        program = get_program_type_from_file(ha_gas_path)

    extensions = _EXT_MAP.get(program, [".log", ".out"])

    def _find(directory, stem):
        for ext in extensions:
            candidate = os.path.join(directory, stem + ext)
            if os.path.isfile(candidate):
                return candidate
        return None

    def _derive(gas_path, suffix):
        dirpath = os.path.dirname(gas_path) or "."
        stem = os.path.splitext(os.path.basename(gas_path))[0]
        found = _find(dirpath, f"{stem}{suffix}")
        if found is not None:
            return found
        return os.path.join(dirpath, f"{stem}{suffix}{extensions[0]}")

    results = {
        "a": _derive(ha_gas_path, "_cb"),
        "ha_solv": _derive(ha_gas_path, "_sp"),
        "a_solv": _derive(ha_gas_path, "_cb_sp"),
        "ref": _derive(href_gas_path, "_cb"),
        "href_solv": _derive(href_gas_path, "_sp"),
        "ref_solv": _derive(href_gas_path, "_cb_sp"),
    }

    missing = [
        f"  {k}: {v}" for k, v in results.items() if not os.path.isfile(v)
    ]
    if missing:
        raise click.UsageError(
            "Auto-discovery could not find some companion output files.\n"
            "Missing files:\n" + "\n".join(missing) + "\n\n"
            "Provide them explicitly or ensure output files follow:\n"
            "  <basename>_cb.<ext>     (conjugate base)\n"
            "  <basename>_sp.<ext>     (solvent single-point)\n"
            "  <basename>_cb_sp.<ext>  (conjugate base solvent SP)"
        )
    return results


def _resolve_output_cls(output_files):
    """Return the correct pKa output class for the detected backend."""
    program = _detect_program_from_outputs(output_files)
    logger.info(f"Auto-detected output program: {program}")
    if program == "gaussian":
        from chemsmart.io.gaussian.output import Gaussian16pKaOutput

        return Gaussian16pKaOutput
    from chemsmart.io.orca.output import ORCApKaOutput

    return ORCApKaOutput


@click.group(name="pka", cls=MyGroup)
@click.option(
    "-T",
    "--temperature",
    type=float,
    default=298.15,
    show_default=True,
    help="Temperature in Kelvin.",
)
@click.option(
    "-conc",
    "--concentration",
    type=float,
    default=1.0,
    show_default=True,
    help="Concentration in mol/L.",
)
@click.option(
    "-csg",
    "--cutoff-entropy-grimme",
    type=float,
    default=100.0,
    show_default=True,
    help="Cutoff frequency (cm^-1) for entropy (Grimme quasi-RRHO).",
)
@click.option(
    "-ch",
    "--cutoff-enthalpy",
    type=float,
    default=100.0,
    show_default=True,
    help="Cutoff frequency (cm^-1) for enthalpy (Head-Gordon).",
)
@click.pass_context
def pka(
    ctx, temperature, concentration, cutoff_entropy_grimme, cutoff_enthalpy
):
    """Backend-independent pKa output analysis."""
    ctx.ensure_object(dict)
    ctx.obj["pka_shared"] = dict(
        temperature=temperature,
        concentration=concentration,
        cutoff_entropy_grimme=cutoff_entropy_grimme,
        cutoff_enthalpy=cutoff_enthalpy,
    )


@pka.command("analyze", cls=MyCommand)
@click_pka_analyze_options
@click.pass_context
def analyze(
    ctx,
    ha,
    a,
    href,
    ref,
    ha_solv,
    a_solv,
    href_solv,
    ref_solv,
    reference_pka,
    **kwargs,
):
    """Compute pKa from existing output files (auto-detects backend).

    Uses the Dual-level Proton Exchange scheme.

    \b
    Species labels:
        HA   - target acid           HRef - reference acid
        A-   - target conjugate base Ref- - reference conjugate base

    \b
    Reaction:  HA + Ref-  ->  A- + HRef

    \b
    Only -ha and -hr are required.  The remaining six files are
    auto-discovered from the naming convention:
      <basename>_cb.<ext>     conjugate base
      <basename>_sp.<ext>     solvent single-point
      <basename>_cb_sp.<ext>  conjugate base solvent SP
    Override any auto-discovered path with the corresponding flag.

    \b
    Examples:
      chemsmart run pka analyze \\
          -ha acid_opt.log -hr collidine-H_opt.log -rp 6.75

      chemsmart run pka analyze \\
          -ha acid.log -a base.log \\
          -hr ref.log -r ref_base.log \\
          -has acid_sp.log -as base_sp.log \\
          -hrs ref_sp.log -rs ref_base_sp.log \\
          -rp 6.75 -T 298.15
    """
    shared = ctx.obj["pka_shared"]

    # Auto-discover missing companion files when at least -ha and -hr given
    if ha is not None and href is not None:
        optional = {
            "a": a,
            "ha_solv": ha_solv,
            "a_solv": a_solv,
            "ref": ref,
            "href_solv": href_solv,
            "ref_solv": ref_solv,
        }
        if any(v is None for v in optional.values()):
            discovered = _auto_discover_pka_files(ha, href)
            for key in optional:
                if optional[key] is None:
                    optional[key] = discovered[key]
                    logger.info(f"Auto-discovered {key}: {discovered[key]}")
            a = optional["a"]
            ha_solv = optional["ha_solv"]
            a_solv = optional["a_solv"]
            ref = optional["ref"]
            href_solv = optional["href_solv"]
            ref_solv = optional["ref_solv"]

    validate_analyze_files(
        ha, a, href, ref, ha_solv, a_solv, href_solv, ref_solv, reference_pka
    )

    output_files = [ha, a, href, ref, ha_solv, a_solv, href_solv, ref_solv]
    output_cls = _resolve_output_cls(output_files)

    logger.info("Computing pKa (Dual-level Proton Exchange)...")
    output_cls.print_pka_summary(
        ha_gas_file=ha,
        a_gas_file=a,
        href_gas_file=href,
        ref_gas_file=ref,
        ha_solv_file=ha_solv,
        a_solv_file=a_solv,
        href_solv_file=href_solv,
        ref_solv_file=ref_solv,
        pka_reference=reference_pka,
        temperature=shared["temperature"],
        concentration=shared["concentration"],
        cutoff_entropy_grimme=shared["cutoff_entropy_grimme"],
        cutoff_enthalpy=shared["cutoff_enthalpy"],
    )

    # Return None so process_pipeline skips jobrunner execution.
    return None


@pka.command("batch-analyze", cls=MyCommand)
@click.option(
    "-o",
    "--output-table",
    type=click.Path(exists=True),
    required=True,
    help=(
        "Table of precomputed output file paths.  Columns: basename, ha_gas, a_gas, href_gas, ref_gas, ha_sp, a_sp, href_sp, ref_sp, pka_ref."
    ),
)
@click.option(
    "-O",
    "--output-results",
    type=click.Path(),
    default=None,
    help="Path to write results table (.csv/.txt).  Stdout if omitted.",
)
@click.option(
    "-p",
    "--program",
    type=click.Choice(["gaussian", "orca", "auto"]),
    default="auto",
    show_default=True,
    help=(
        "QC program that produced the output files.  'auto' detects from the first file in the table."
    ),
)
@click.pass_context
def batch_analyze(ctx, output_table, output_results, program, **kwargs):
    """Batch pKa computation from a table of precomputed output files.

    Blank reference-acid cells are filled from the previous row.

    \b
    Examples:
      chemsmart run pka batch-analyze -o outputs.csv
      chemsmart run pka batch-analyze -o outputs.csv -O results.csv
    """
    shared = ctx.obj["pka_shared"]

    if program == "auto":
        import csv

        with open(output_table) as fh:
            reader = csv.DictReader(fh)
            first_row = next(reader, None)
        if first_row is None:
            raise click.UsageError("Output table is empty.")
        for col in ["ha_gas", "a_gas", "href_gas", "ref_gas"]:
            val = first_row.get(col)
            if val and os.path.isfile(val):
                program = get_program_type_from_file(val)
                break
        if program == "auto" or program == "unknown":
            raise click.UsageError(
                "Could not auto-detect QC program from output table.  "
                "Use -p gaussian or -p orca."
            )
        logger.info(f"Auto-detected program: {program}")

    run_pka_from_output_table(
        output_table=output_table,
        output_results=output_results,
        temperature=shared["temperature"],
        concentration=shared["concentration"],
        cutoff_entropy_grimme=shared["cutoff_entropy_grimme"],
        cutoff_enthalpy=shared["cutoff_enthalpy"],
        program=program,
    )

    return None


def click_pka_thermo_options(f):
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
        "-hr",
        "--href-file",
        type=click.Path(exists=True),
        default=None,
        help="Path to HRef (reference acid) optimization output file.",
    )
    @click.option(
        "-r",
        "--ref-file",
        type=click.Path(exists=True),
        default=None,
        help="Path to Ref- (reference conjugate base) optimization output file.",
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
        help="Cutoff for entropy (cm^-1), Grimme method (default: 100).",
    )
    @click.option(
        "-ch",
        "--cutoff-enthalpy",
        type=float,
        default=100.0,
        help="Cutoff for enthalpy (cm^-1), Head-Gordon method (default: 100).",
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


def format_thermo_results(results, energy_units):
    lines = [
        "=" * 78,
        "pKa Thermochemistry Extraction",
        "=" * 78,
        f"Temperature: {results['settings']['temperature']} K",
        f"Concentration: {results['settings']['concentration']} mol/L",
        f"Entropy cutoff (Grimme): "
        f"{results['settings']['cutoff_entropy_grimme']} cm⁻¹",
        f"Enthalpy cutoff (Head-Gordon): "
        f"{results['settings']['cutoff_enthalpy']} cm⁻¹",
        f"Energy units: {energy_units}",
        "-" * 78,
        "",
        f"{'Species':<10} {'E':<20} {'qh-G(T)':<20} {'G_corr':<20}",
        "-" * 78,
    ]
    for key in ["HA", "A", "HRef", "Ref"]:
        if key in results:
            sp = results[key]
            E = sp["E"]
            qh_G = sp["qh_G"]
            lines.append(
                f"{sp['name']:<10} {E:<20.10f} "
                f"{qh_G:<20.10f} {qh_G - E:<20.10f}"
            )
    lines.append("=" * 78)
    return "\n".join(lines)
