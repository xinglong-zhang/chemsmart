"""
Backend-independent pKa output analysis.

Registered as ``chemsmart run pka`` — strictly post-processing,
never invokes Gaussian, ORCA, or any other QC backend.

Subcommands
-----------
analyze        Compute pKa from 8 existing output files.
batch-analyze  Batch pKa from a table of output file paths.
"""

import logging
import os

import click

from chemsmart.cli.pka_helpers import (
    click_pka_analyze_options,
    run_pka_from_output_table,
    validate_analyze_files,
)
from chemsmart.utils.cli import MyCommand, MyGroup
from chemsmart.utils.io import get_program_type_from_file

logger = logging.getLogger(__name__)

# Extension preference per backend (used by auto-discovery).
_EXT_MAP = {
    "gaussian": [".log", ".out"],
    "orca": [".out", ".log"],
}


# ── internal helpers (kept in this module, not pka_helpers) ─────────────


def _detect_program_from_outputs(filepaths):
    """Detect and validate a unique QC program across supplied outputs."""
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
    else:
        from chemsmart.io.orca.output import ORCApKaOutput

        return ORCApKaOutput


# ═══════════════════════════════════════════════════════════════════════
# pka group — top-level, backend-independent
# ═══════════════════════════════════════════════════════════════════════


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
    """Backend-independent pKa output analysis.

    \b
    This command is strictly post-processing — it never submits
    or invokes Gaussian / ORCA calculations.
    The QC backend is auto-detected from output-file signatures.

    \b
    Subcommands:
      analyze        Compute pKa from existing output files.
      batch-analyze  Batch pKa from a table of output file paths.

    \b
    For job submission use the backend-specific subcommands:
      chemsmart run gaussian ... pka submit ...
      chemsmart run orca     ... pka submit ...
    """
    ctx.ensure_object(dict)
    ctx.obj["pka_shared"] = dict(
        temperature=temperature,
        concentration=concentration,
        cutoff_entropy_grimme=cutoff_entropy_grimme,
        cutoff_enthalpy=cutoff_enthalpy,
    )


# ═══════════════════════════════════════════════════════════════════════
# pka analyze
# ═══════════════════════════════════════════════════════════════════════


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
    OutputCls = _resolve_output_cls(output_files)

    logger.info("Computing pKa (Dual-level Proton Exchange)...")

    OutputCls.print_pka_summary(
        ha_gas_file=ha,
        a_gas_file=a,
        hb_gas_file=href,
        b_gas_file=ref,
        ha_solv_file=ha_solv,
        a_solv_file=a_solv,
        hb_solv_file=href_solv,
        b_solv_file=ref_solv,
        pka_reference=reference_pka,
        temperature=shared["temperature"],
        concentration=shared["concentration"],
        cutoff_entropy_grimme=shared["cutoff_entropy_grimme"],
        cutoff_enthalpy=shared["cutoff_enthalpy"],
    )

    # Return None so process_pipeline skips jobrunner execution.
    return None


# ═══════════════════════════════════════════════════════════════════════
# pka batch-analyze
# ═══════════════════════════════════════════════════════════════════════


@pka.command("batch-analyze", cls=MyCommand)
@click.option(
    "-o",
    "--output-table",
    type=click.Path(exists=True),
    required=True,
    help=(
        "Table of precomputed output file paths.  Columns: basename, "
        "ha_gas, a_gas, hb_gas, b_gas, ha_sp, a_sp, hb_sp, b_sp, "
        "pka_ref."
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
        "QC program that produced the output files.  "
        "'auto' detects from the first file in the table."
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

    # Resolve 'auto' to a concrete backend
    if program == "auto":
        # Peek at the first data file in the table to detect the program
        import csv

        with open(output_table) as fh:
            reader = csv.DictReader(fh)
            first_row = next(reader, None)
        if first_row is None:
            raise click.UsageError("Output table is empty.")
        # Try the first non-empty file column
        for col in ["ha_gas", "a_gas", "hb_gas", "b_gas"]:
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
