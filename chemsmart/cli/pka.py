import logging
import os

import click

from chemsmart.utils.io import get_program_type_from_file

logger = logging.getLogger(__name__)

# Output extensions tried during auto-discovery, ordered by preference.
_EXT_MAP = {
    "gaussian": [".log", ".out"],
    "orca": [".out", ".log"],
}


def _detect_program_from_outputs(filepaths):
    """Detect and validate a unique QC program across all provided outputs.

    Delegates per-file detection to ``get_program_type_from_file`` so
    output-format parsing remains centralized in ``chemsmart.utils.io``.
    """
    detected = {}
    programs = set()
    for fp in filepaths:
        p = get_program_type_from_file(fp)
        detected[fp] = p
        if p != "unknown":
            programs.add(p)

    if not programs:
        raise click.UsageError(
            "Could not detect output-file program type from supplied files. "
            "Supported for pKa analysis: Gaussian and ORCA output files."
        )

    if len(programs) > 1:
        pairs = ", ".join([f"{k}: {v}" for k, v in detected.items()])
        raise click.UsageError(
            "Supplied files contain mixed program types. "
            "Please use outputs from a single QC program (Gaussian or ORCA).\n"
            f"Detected: {pairs}"
        )

    program = next(iter(programs))
    if program not in {"gaussian", "orca"}:
        raise click.UsageError(
            f"Detected unsupported program '{program}'. "
            "Only Gaussian and ORCA are supported for pKa output analysis."
        )
    return program


def _auto_discover_pka_files(ha_gas_path, hb_gas_path, program=None):
    """Infer the remaining 6 output paths from the HA and HB gas-phase paths.

    Naming convention (produced by the pKa job classes):

    * HA gas  :  ``<acid_basename>.<ext>``
    * A⁻ gas  :  ``<acid_basename>_cb.<ext>``
    * HA solv :  ``<acid_basename>_sp.<ext>``
    * A⁻ solv :  ``<acid_basename>_cb_sp.<ext>``
    * HB gas  :  ``<ref_basename>.<ext>``
    * B⁻ gas  :  ``<ref_basename>_cb.<ext>``
    * HB solv :  ``<ref_basename>_sp.<ext>``
    * B⁻ solv :  ``<ref_basename>_cb_sp.<ext>``

    Parameters
    ----------
    ha_gas_path : str
        Path to the HA gas-phase output (required).
    hb_gas_path : str
        Path to the HB gas-phase output (required).
    program : str or None
        ``"gaussian"`` or ``"orca"``.  If *None* the program is
        auto-detected from *ha_gas_path*.

    Returns
    -------
    dict
        Keys: ``a``, ``b_output``, ``ha_solv``,
        ``a_solv``, ``hb_solv_output``, ``b_solv_output``.
        Each value is an existing file path (str).

    Raises
    ------
    click.UsageError
        If any inferred file does not exist on disk.
    """
    if program is None:
        program = get_program_type_from_file(ha_gas_path)

    extensions = _EXT_MAP.get(program, [".log", ".out"])

    def _find(directory, stem, extensions):
        """Return the first existing path for ``stem`` + one of *extensions*."""
        for ext in extensions:
            candidate = os.path.join(directory, stem + ext)
            if os.path.isfile(candidate):
                return candidate
        return None

    def _derive_companion(gas_path, suffix, extensions):
        """Derive a companion file path from a gas-phase output path.

        ``suffix`` is appended to the stem of *gas_path* (e.g. ``"_cb"``
        or ``"_sp"``).
        """
        dirpath = os.path.dirname(gas_path) or "."
        stem = os.path.splitext(os.path.basename(gas_path))[0]
        companion_stem = f"{stem}{suffix}"
        found = _find(dirpath, companion_stem, extensions)
        if found is not None:
            return found
        # Construct the expected path for the error message
        expected = os.path.join(dirpath, companion_stem + extensions[0])
        return expected  # will be validated later

    results = {
        "a": _derive_companion(ha_gas_path, "_cb", extensions),
        "ha_solv": _derive_companion(ha_gas_path, "_sp", extensions),
        "a_solv": _derive_companion(ha_gas_path, "_cb_sp", extensions),
        "b_output": _derive_companion(hb_gas_path, "_cb", extensions),
        "hb_solv_output": _derive_companion(hb_gas_path, "_sp", extensions),
        "b_solv_output": _derive_companion(hb_gas_path, "_cb_sp", extensions),
    }

    # Validate that all inferred files exist
    missing = [
        f"  {key}: {path}"
        for key, path in results.items()
        if not os.path.isfile(path)
    ]
    if missing:
        raise click.UsageError(
            "Auto-discovery could not find some companion output files.\n"
            "Missing files:\n" + "\n".join(missing) + "\n\n"
            "Either provide them explicitly via CLI options or ensure\n"
            "output files follow the naming convention:\n"
            "  <basename>_cb.<ext>    (conjugate base)\n"
            "  <basename>_sp.<ext>    (solvent single-point)\n"
            "  <basename>_cb_sp.<ext> (conjugate base solvent SP)"
        )

    return results


@click.command(name="pka")
@click.option(
    "-rp",
    "--reference-pka",
    type=float,
    required=True,
    help="Experimental pKa of the reference acid HB used in the proton-exchange scheme.",
)
@click.option(
    "-ha",
    "--ha-output",
    type=click.Path(exists=True),
    required=True,
    help="Path to HA (protonated acid) gas-phase optimization/frequency output file.",
)
@click.option(
    "-a",
    "--a-output",
    type=click.Path(exists=True),
    default=None,
    help=(
        "Path to A⁻ (conjugate base) gas-phase output file.  "
        "Auto-discovered as <HA_basename>_cb.<ext> when omitted."
    ),
)
@click.option(
    "-hb",
    "--hb-output",
    type=click.Path(exists=True),
    required=True,
    help="Path to HB (reference acid) gas-phase optimization/frequency output file.",
)
@click.option(
    "-b",
    "--b-output",
    type=click.Path(exists=True),
    default=None,
    help=(
        "Path to B⁻ (reference conjugate base) gas-phase output file.  "
        "Auto-discovered as <HB_basename>_cb.<ext> when omitted."
    ),
)
@click.option(
    "-has",
    "--ha-solv-output",
    type=click.Path(exists=True),
    default=None,
    help=(
        "Path to HA solvent single-point output file.  "
        "Auto-discovered as <HA_basename>_sp.<ext> when omitted."
    ),
)
@click.option(
    "-as",
    "--a-solv-output",
    type=click.Path(exists=True),
    default=None,
    help=(
        "Path to A⁻ solvent single-point output file.  "
        "Auto-discovered as <HA_basename>_cb_sp.<ext> when omitted."
    ),
)
@click.option(
    "-hbs",
    "--hb-solv-output",
    type=click.Path(exists=True),
    default=None,
    help=(
        "Path to HB solvent single-point output file.  "
        "Auto-discovered as <HB_basename>_sp.<ext> when omitted."
    ),
)
@click.option(
    "-bs",
    "--b-solv-output",
    type=click.Path(exists=True),
    default=None,
    help=(
        "Path to B⁻ solvent single-point output file.  "
        "Auto-discovered as <HB_basename>_cb_sp.<ext> when omitted."
    ),
)
@click.option(
    "-T",
    "--temperature",
    type=float,
    default=298.15,
    show_default=True,
    help="Temperature in Kelvin used in thermochemistry corrections and pKa calculation.",
)
@click.option(
    "-conc",
    "--concentration",
    type=float,
    default=1.0,
    show_default=True,
    help="Concentration in mol/L used in quasi-harmonic thermochemistry.",
)
@click.option(
    "-csg",
    "--cutoff-entropy-grimme",
    type=float,
    default=100.0,
    show_default=True,
    help="Cutoff frequency (cm^-1) for entropy using Grimme's quasi-RRHO method.",
)
@click.option(
    "-ch",
    "--cutoff-enthalpy",
    type=float,
    default=100.0,
    show_default=True,
    help="Cutoff frequency (cm^-1) for enthalpy using Head-Gordon's quasi-RRHO method.",
)
@click.pass_context
def pka(
    ctx,
    reference_pka,
    ha_output,
    a_output,
    hb_output,
    b_output,
    ha_solv_output,
    a_solv_output,
    hb_solv_output,
    b_solv_output,
    temperature,
    concentration,
    cutoff_entropy_grimme,
    cutoff_enthalpy,
):
    """Analyze pKa from existing outputs without selecting Gaussian/ORCA subcommands.

    \b
    This mode is strictly post-processing and does not submit jobs.
    The QC backend is auto-detected from output-file signatures.

    \b
    Only -ha and -hb are required.  The remaining six output files
    are auto-discovered from the naming convention:
      <basename>_cb.<ext>    → conjugate base
      <basename>_sp.<ext>    → solvent single-point
      <basename>_cb_sp.<ext> → conjugate base solvent SP
    Override any auto-discovered path with the corresponding flag.
    """
    # --- auto-discover missing companion files ---
    optional_files = {
        "a": a_output,
        "b_output": b_output,
        "ha_solv": ha_solv_output,
        "a_solv": a_solv_output,
        "hb_solv_output": hb_solv_output,
        "b_solv_output": b_solv_output,
    }

    if any(v is None for v in optional_files.values()):
        discovered = _auto_discover_pka_files(ha_output, hb_output)
        # Only fill in values that were not explicitly provided
        for key in optional_files:
            if optional_files[key] is None:
                optional_files[key] = discovered[key]
                logger.info(f"Auto-discovered {key}: {discovered[key]}")

    a_output = optional_files["a"]
    b_output = optional_files["b_output"]
    ha_solv_output = optional_files["ha_solv"]
    a_solv_output = optional_files["a_solv"]
    hb_solv_output = optional_files["hb_solv_output"]
    b_solv_output = optional_files["b_solv_output"]

    output_files = [
        ha_output,
        a_output,
        hb_output,
        b_output,
        ha_solv_output,
        a_solv_output,
        hb_solv_output,
        b_solv_output,
    ]

    program = _detect_program_from_outputs(output_files)
    logger.info(f"Detected output program: {program}")

    kwargs = dict(
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

    if program == "gaussian":
        from chemsmart.io.gaussian.output import Gaussian16pKaOutput

        Gaussian16pKaOutput.print_pka_summary(**kwargs)
    else:
        from chemsmart.io.orca.output import ORCApKaOutput

        ORCApKaOutput.print_pka_summary(**kwargs)

    # Return None so run.process_pipeline skips jobrunner execution.
    return None
