"""
Shared helpers for Gaussian and ORCA pKa CLI subcommands.

Centralises option decorators, CDXML proton-detection logic,
reference-acid validation, and the batch output-table post-processing
pipeline so that ``cli/gaussian/pka.py`` and ``cli/orca/pka.py`` stay
thin and free of duplication.
"""

import functools
import logging

import click

logger = logging.getLogger(__name__)


# ── shared Click option decorators ──────────────────────────────────────


def click_pka_shared_options(f):
    """Options that apply across *all* pKa execution modes.

    Attached to the ``pka`` group so every subcommand can read them
    from ``ctx.obj["pka_shared"]``.
    """

    @click.option(
        "-t",
        "--thermodynamic-cycle",
        type=click.Choice(["direct", "proton exchange"]),
        default="proton exchange",
        help=(
            "Thermodynamic cycle type. 'proton exchange' uses a reference "
            "acid (default). 'direct' uses absolute free energy of H+ in "
            "water."
        ),
    )
    @click.option(
        "-r",
        "--reference",
        type=click.Path(exists=True),
        default=None,
        help=(
            "Path to geometry file for reference acid (HB) for proton "
            "exchange cycle."
        ),
    )
    @click.option(
        "-rpi",
        "--reference-proton-index",
        type=int,
        default=None,
        help=(
            "1-based index of the proton to remove from reference acid "
            "(HB). Required when --reference is provided."
        ),
    )
    @click.option(
        "-rcc",
        "--reference-color-code",
        type=int,
        default=None,
        help=(
            "CDXML colour-table index identifying the proton in the "
            "reference acid file."
        ),
    )
    @click.option(
        "-rc",
        "--reference-charge",
        type=int,
        default=None,
        help="Charge of the reference acid (HB).",
    )
    @click.option(
        "-rm",
        "--reference-multiplicity",
        type=int,
        default=None,
        help="Multiplicity of the reference acid (HB).",
    )
    @click.option(
        "--reference-conjugate-base-charge",
        type=int,
        default=None,
        help=(
            "Charge of the reference conjugate base (B-). "
            "Defaults to (reference_charge - 1)."
        ),
    )
    @click.option(
        "--reference-conjugate-base-multiplicity",
        type=int,
        default=None,
        help=(
            "Multiplicity of the reference conjugate base (B-). "
            "Defaults to reference_multiplicity."
        ),
    )
    @click.option(
        "-dG",
        "--delta-g-proton",
        type=float,
        default=-265.9,
        help=(
            "Absolute free energy of H+ in water (kcal/mol) for direct "
            "cycle. Default: -265.9 kcal/mol (Tissandier et al., 1998)."
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
            "Multiplicity of the conjugate base (A-). "
            "Defaults to multiplicity."
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
            "Cutoff frequency for entropy (cm^-1) using Grimme's "
            "quasi-RRHO method (default: 100)."
        ),
    )
    @click.option(
        "-ch",
        "--cutoff-enthalpy",
        type=float,
        default=100.0,
        help=(
            "Cutoff frequency for enthalpy (cm^-1) using Head-Gordon's "
            "method (default: 100)."
        ),
    )
    @functools.wraps(f)
    def wrapper(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper


def click_pka_submit_options(f):
    """Options specific to the ``pka submit`` subcommand."""

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
    """Options for the ``pka analyze`` subcommand (output-file parsing).

    \b
    Species labels:
        HA   – target acid (protonated form)
        A-   – target conjugate base
        HRef – reference acid
        Ref- – reference conjugate base

    \b
    Proton-exchange schematic:
        HA + Ref-  →  A- + HRef
    """
    f = click.option(
        "-rp",
        "--reference-pka",
        type=float,
        default=None,
        help=(
            "Experimental pKa of reference acid HRef. "
            "Required for proton exchange analysis."
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


def click_pka_thermo_options(f):
    """Shared options for the ``pka thermo`` subcommand."""

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


# ── CDXML proton-detection helpers ──────────────────────────────────────


def resolve_proton_index(filename, proton_index, color_code):
    """Resolve the proton index for the *target* molecule.

    If *proton_index* is already set, returns it unchanged.  Otherwise
    attempts CDXML colour-based auto-detection.

    Returns
    -------
    proton_index : int
        Resolved 1-based proton index.
    pka_molecules : list or None
        If the CDXML file contains multiple fragments each carrying a
        coloured proton, the list of ``PKaMolecule`` objects is returned
        so the caller can create one job per molecule.  ``None`` when
        a single proton index suffices.

    Raises
    ------
    click.UsageError
        When the proton cannot be determined.
    """
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
                f"Detected {len(pka_mols)} molecules with per-fragment "
                f"proton auto-detection in {filename}."
            )
            return None, pka_mols

        proton_index = pka_mols[0].proton_index
        logger.info(
            f"Detected proton index {proton_index} from CDXML "
            f"colour in {filename}."
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
    """Auto-detect reference acid proton index from CDXML if needed.

    Returns the (possibly updated) *reference_proton_index*.

    Raises
    ------
    click.UsageError
        When auto-detection fails and no explicit index was given.
    """
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
                f"Detected reference proton index "
                f"{reference_proton_index} from CDXML colour "
                f"in {reference}."
            )
            return reference_proton_index
        except ValueError as exc:
            raise click.UsageError(
                f"Could not auto-detect reference proton from "
                f"CDXML colour: {exc}\n"
                "Use -rpi/--reference-proton-index to specify "
                "the proton explicitly."
            )
    elif reference_color_code is not None:
        raise click.UsageError(
            "-rcc/--reference-color-code can only be used when "
            "--reference is a .cdx/.cdxml file."
        )

    return None


def validate_reference_options(shared):
    """Validate that all required reference acid options are present.

    Parameters
    ----------
    shared : dict
        The ``ctx.obj["pka_shared"]`` dictionary.

    Raises
    ------
    click.UsageError
        When required options are missing.
    """
    reference = shared["reference"]
    if reference is None:
        return

    if shared["thermodynamic_cycle"] != "proton exchange":
        raise click.UsageError(
            "Reference acid file can only be used with 'proton exchange' "
            "cycle. Use -t 'proton exchange' or remove the -r option."
        )

    # Attempt CDXML auto-detect for reference proton index
    shared["reference_proton_index"] = resolve_reference_proton(
        reference,
        shared["reference_proton_index"],
        shared["reference_color_code"],
    )

    missing = []
    if shared["reference_proton_index"] is None:
        missing.append(
            "-rpi/--reference-proton-index (or use a .cdxml reference "
            "with a coloured proton)"
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


# ── Batch output-table post-processing ─────────────────────────────────


def run_pka_from_output_table(
    output_table,
    output_results,
    temperature,
    concentration,
    cutoff_entropy_grimme,
    cutoff_enthalpy,
    program="gaussian",
):
    """Compute pKa values from a table of precomputed output files.

    Implements *batch-analyze* mode.  Called from both Gaussian and ORCA
    ``pka batch-analyze`` CLI handlers.

    Parameters
    ----------
    output_table : str
        Path to the output-table file.
    output_results : str or None
        Path to write results table.  If *None*, results print to stdout.
    temperature : float
    concentration : float
    cutoff_entropy_grimme : float
    cutoff_enthalpy : float
    program : ``"gaussian"`` or ``"orca"``
    """
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


# ── Analyze-mode validation helper ──────────────────────────────────────


def validate_analyze_files(
    ha, a, href, ref, ha_solv, a_solv, href_solv, ref_solv, reference_pka
):
    """Validate that all 8 output files and reference pKa are provided.

    Raises
    ------
    click.UsageError
        When any required file or value is missing.
    """
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


# ── Thermo format helper ───────────────────────────────────────────────


def format_thermo_results(results, energy_units):
    """Format thermochemistry results for display."""
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
    for key in ["HA", "A", "HB", "B"]:
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
