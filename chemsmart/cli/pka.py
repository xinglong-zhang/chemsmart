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

from chemsmart.cli.thermochemistry.thermochemistry import (
    resolve_entropy_cutoff,
    thermochemistry_cutoff_options,
    thermochemistry_temp_pressure_conc_options,
)
from chemsmart.io.file import PKaCDXFile
from chemsmart.utils.cli import MyCommand, MyGroup
from chemsmart.utils.io import (
    detect_program_type_from_files,
    get_program_output_extensions,
    get_program_type_from_file,
)

logger = logging.getLogger(__name__)


# ── helpers moved from cli/pka_helpers.py ──────────────────────────────


def resolve_pka_entropy_cutoff(cutoff_entropy_grimme, cutoff_entropy_truhlar):
    """Resolve pKa entropy cutoff; default to Grimme 100 cm⁻¹ when unset."""
    s_freq_cutoff, entropy_method = resolve_entropy_cutoff(
        cutoff_entropy_grimme, cutoff_entropy_truhlar
    )
    if s_freq_cutoff is None:
        return 100.0, "grimme"
    return s_freq_cutoff, entropy_method


def click_pka_thermochemistry_options(f):
    """Thermochemistry options reused by pKa submission and analysis."""
    f = thermochemistry_temp_pressure_conc_options(
        f,
        temperature_required=False,
        temperature_default=298.15,
        concentration_default=1.0,
        pressure_default=1.0,
        concentration_short="-conc",
    )
    return thermochemistry_cutoff_options(
        f,
        enthalpy_default=100.0,
    )


def click_pka_shared_options(f):
    f = click_pka_thermochemistry_options(f)

    @click.option(
        "-s",
        "--scheme",
        "scheme",
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
        default=None,
        help="Solvent ID for solution phase SP (default: project setting or water).",
    )
    @functools.wraps(f)
    def wrapper(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper


def click_pka_proton_options(f):
    """Options that identify the proton to remove (-pi / -cc).

    Applied to the ``pka`` *group* so the values are captured once and
    stored in ``ctx.obj``.  The ``submit`` subcommand reads them from
    there rather than re-declaring them.
    """

    @click.option(
        "-pi",
        "--proton-index",
        type=int,
        required=False,
        help=(
            "1-based index of the proton to remove for deprotonation. "
            "Required unless a .cdxml file with a coloured proton is used."
        ),
    )
    @click.option(
        "-cc",
        "--color-code",
        type=int,
        default=None,
        help=(
            "CDXML colour-table index identifying the proton to remove. "
            "If omitted, the uniquely coloured hydrogen is auto-detected."
        ),
    )
    @functools.wraps(f)
    def wrapper(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper


def click_pka_analysis_scheme_options(f):
    """Scheme and proton solvation options for pKa output analysis."""

    @click.option(
        "-s",
        "--scheme",
        "scheme",
        type=click.Choice(["direct", "proton exchange"]),
        default=None,
        help=(
            "Thermodynamic cycle for analysis. 'direct' requires -dG. "
            "Default: proton exchange when omitted."
        ),
    )
    @click.option(
        "-dG",
        "--delta-g-proton",
        "delta_g_proton",
        type=float,
        default=None,
        help=(
            "G_soln(H+) in kcal/mol for the direct cycle (e.g. -265.9 for water). "
            "Used only when both this flag and --scheme direct are specified."
        ),
    )
    @functools.wraps(f)
    def wrapper(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper


def _resolve_pka_analysis_scheme(scheme, delta_g_proton):
    """Validate CLI options and return the active analysis scheme."""
    if scheme is None:
        scheme = "proton exchange"

    if scheme == "direct" and delta_g_proton is None:
        raise click.UsageError(
            "-dG/--delta-g-proton is required when --scheme direct is specified."
        )
    if delta_g_proton is not None and scheme != "direct":
        logger.info(
            "Ignoring -dG/--delta-g-proton because --scheme direct was not "
            "specified; using proton exchange analysis."
        )
    return scheme


def _scheme_display_name(scheme):
    names = {
        "direct": "Direct Dissociation",
        "proton exchange": "Proton Exchange",
    }
    return names.get(scheme, scheme)


def validate_direct_analyze_files(ha, a, ha_solv, a_solv):
    """Validate required files for direct-cycle pKa analysis."""
    required = [
        ("ha", ha),
        ("a", a),
        ("ha-solv", ha_solv),
        ("a-solv", a_solv),
    ]
    missing = [name for name, path in required if path is None]
    if missing:
        raise click.UsageError(
            "For direct-cycle pKa analysis all four output files are required.\n"
            f"Missing: {', '.join(f'--{name}' for name in missing)}"
        )
    for name, path in required:
        if not os.path.isfile(path):
            raise click.UsageError(f"File not found for --{name}: {path}")


def _auto_discover_direct_pka_files(ha_gas_path, program=None):
    """Infer A- and solvent SP paths from the HA gas-phase output path."""
    if program is None:
        program = get_program_type_from_file(ha_gas_path)

    extensions = get_program_output_extensions(program)

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


def resolve_reference_proton(
    reference, reference_proton_index, reference_color_code
):
    """Resolve reference proton index, handling CDXML auto-detection."""
    if reference_proton_index is not None:
        return reference_proton_index
    if not reference or not reference.lower().endswith((".cdx", ".cdxml")):
        return None

    try:
        pka_file = PKaCDXFile(reference)
        # Assuming a single molecule/fragment for the reference
        pka_mol = pka_file.get_pka_molecules(
            color_code=reference_color_code, index="-1"
        )
        return pka_mol.proton_index
    except (ValueError, FileNotFoundError) as exc:
        raise click.UsageError(
            f"Could not auto-detect reference proton from CDXML: {exc}\n"
            "Use -rpi/--reference-proton-index to specify the proton explicitly."
        )


def validate_reference_options(shared):
    reference = shared["reference"]
    if reference is None:
        return

    if shared["scheme"] != "proton exchange":
        raise click.UsageError(
            "Reference acid file can only be used with 'proton exchange' "
            "cycle. Use -s 'proton exchange' or remove the -r option."
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


def _resolve_batch_output_cls(program):
    if program == "gaussian":
        from chemsmart.io.gaussian.output import Gaussian16pKaOutput

        return Gaussian16pKaOutput
    elif program == "orca":
        from chemsmart.io.orca.output import ORCApKaOutput

        return ORCApKaOutput
    else:
        raise ValueError(f"Unsupported program: {program}")


def _first_output_file_from_table(pka_table):
    """Return the first existing output path from a prepared pKa table.

    This supports backend auto-detection without re-reading the raw table file.
    """
    first_entry = next(iter(pka_table), None)
    if first_entry is None:
        raise click.UsageError("Output table is empty.")
    for col in ["ha_gas", "a_gas", "href_gas", "ref_gas"]:
        val = first_entry.get(col)
        if val and os.path.isfile(val):
            return val
    return None


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

    extensions = get_program_output_extensions(program)

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
    program = detect_program_type_from_files(
        output_files,
        allowed_programs={"gaussian", "orca"},
    )
    logger.info(f"Auto-detected output program: {program}")
    if program == "gaussian":
        from chemsmart.io.gaussian.output import Gaussian16pKaOutput

        return Gaussian16pKaOutput
    from chemsmart.io.orca.output import ORCApKaOutput

    return ORCApKaOutput


@click.group(name="pka", cls=MyGroup)
@click_pka_thermochemistry_options
@click_pka_analysis_scheme_options
@click.pass_context
def pka(
    ctx,
    temperature,
    concentration,
    pressure,
    cutoff_entropy_grimme,
    cutoff_entropy_truhlar,
    cutoff_enthalpy,
    scheme,
    delta_g_proton,
):
    """Backend-independent pKa output analysis."""
    s_freq_cutoff, entropy_method = resolve_pka_entropy_cutoff(
        cutoff_entropy_grimme, cutoff_entropy_truhlar
    )

    ctx.ensure_object(dict)
    ctx.obj["pka_shared"] = dict(
        temperature=temperature,
        concentration=concentration,
        pressure=pressure,
        cutoff_entropy_grimme=s_freq_cutoff,
        cutoff_enthalpy=cutoff_enthalpy,
        entropy_method=entropy_method,
        scheme=scheme,
        delta_g_proton=delta_g_proton,
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
    scheme = _resolve_pka_analysis_scheme(
        shared.get("scheme"), shared.get("delta_g_proton")
    )

    if scheme == "direct":
        if ha is None:
            raise click.UsageError(
                "-ha/--ha is required for direct-cycle pKa analysis."
            )
        optional = {"a": a, "ha_solv": ha_solv, "a_solv": a_solv}
        if any(v is None for v in optional.values()):
            discovered = _auto_discover_direct_pka_files(ha)
            for key in optional:
                if optional[key] is None:
                    optional[key] = discovered[key]
                    logger.info(f"Auto-discovered {key}: {discovered[key]}")
            a = optional["a"]
            ha_solv = optional["ha_solv"]
            a_solv = optional["a_solv"]

        validate_direct_analyze_files(ha, a, ha_solv, a_solv)

        output_files = [ha, a, ha_solv, a_solv]
        output_cls = _resolve_output_cls(output_files)

        logger.info("Computing pKa (Direct Dissociation)...")
        output_cls.print_pka_summary(
            ha_gas_file=ha,
            a_gas_file=a,
            ha_solv_file=ha_solv,
            a_solv_file=a_solv,
            scheme="direct",
            delta_G_proton=shared["delta_g_proton"],
            temperature=shared["temperature"],
            concentration=shared["concentration"],
            pressure=shared["pressure"],
            cutoff_entropy_grimme=shared["cutoff_entropy_grimme"],
            cutoff_enthalpy=shared["cutoff_enthalpy"],
            entropy_method=shared["entropy_method"],
        )
        return None

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

    logger.info(f"Computing pKa ({_scheme_display_name(scheme)})...")
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
        scheme=scheme,
        temperature=shared["temperature"],
        concentration=shared["concentration"],
        pressure=shared["pressure"],
        cutoff_entropy_grimme=shared["cutoff_entropy_grimme"],
        cutoff_enthalpy=shared["cutoff_enthalpy"],
        entropy_method=shared["entropy_method"],
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
    scheme = _resolve_pka_analysis_scheme(
        shared.get("scheme"), shared.get("delta_g_proton")
    )

    from chemsmart.utils.utils import PKaOutputTable

    logger.info(f"Reading pKa output table: {output_table}")
    try:
        pka_output_table = PKaOutputTable.from_file(output_table)
        pka_output_table.prepare(check_file_exists=True, scheme=scheme)
    except (FileNotFoundError, ValueError) as e:
        raise click.UsageError(str(e))

    logger.info(f"Found {len(pka_output_table)} entries in output table")

    if program == "auto":
        output_file = _first_output_file_from_table(pka_output_table)
        if output_file is not None:
            program = get_program_type_from_file(output_file)
        if program == "auto" or program == "unknown":
            raise click.UsageError(
                "Could not auto-detect QC program from output table.  "
                "Use -p gaussian or -p orca."
            )
        logger.info(f"Auto-detected program: {program}")

    output_cls = _resolve_batch_output_cls(program)
    logger.info(
        f"Computing pKa ({_scheme_display_name(scheme)}) "
        f"for {len(pka_output_table)} systems "
        f"(T={shared['temperature']}K, program={program})"
    )
    results = pka_output_table.run_pka(
        output_cls=output_cls,
        temperature=shared["temperature"],
        concentration=shared["concentration"],
        pressure=shared["pressure"],
        cutoff_entropy_grimme=shared["cutoff_entropy_grimme"],
        cutoff_enthalpy=shared["cutoff_enthalpy"],
        entropy_method=shared["entropy_method"],
        scheme=scheme,
        delta_G_proton=shared.get("delta_g_proton"),
    )
    output_string = pka_output_table.echo_pka_output_table_results(
        results=results,
        output_results=output_results,
        temperature=shared["temperature"],
        pressure=shared["pressure"],
        scheme=scheme,
    )
    click.echo(output_string)

    return None


def click_pka_thermo_options(f):
    f = click_pka_thermochemistry_options(f)

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
        f"Pressure: {results['settings'].get('pressure', 1.0)} atm",
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
