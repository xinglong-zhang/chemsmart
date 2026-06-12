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

from chemsmart.analysis.thermochemistry import Thermochemistry
from chemsmart.cli.thermochemistry.thermochemistry import (
    resolve_entropy_cutoff,
    thermochemistry_cutoff_options,
    thermochemistry_temp_pressure_conc_options,
)
from chemsmart.io.file import PKaCDXFile
from chemsmart.utils.cli import MyCommand, MyGroup
from chemsmart.utils.constants import HARTREE_TO_KCAL_MOL, energy_conversion
from chemsmart.utils.io import get_program_type_from_file

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


def _pka_thermochemistry_kwargs(
    temperature,
    concentration,
    pressure,
    cutoff_entropy_grimme,
    cutoff_enthalpy,
    entropy_method,
):
    return {
        "temperature": temperature,
        "concentration": concentration,
        "pressure": pressure,
        "s_freq_cutoff": cutoff_entropy_grimme,
        "h_freq_cutoff": cutoff_enthalpy,
        "entropy_method": entropy_method,
        "energy_units": "hartree",
        "check_imaginary_frequencies": True,
    }


def _require_thermochemistry_property(thermo, filepath, attr, label):
    value = getattr(thermo, attr)
    if value is None:
        raise ValueError(f"Could not extract {label} from file: {filepath}")
    return value


def pka_gas_phase_data(
    filepath,
    temperature=298.15,
    concentration=1.0,
    pressure=1.0,
    cutoff_entropy_grimme=100.0,
    cutoff_enthalpy=100.0,
    entropy_method="grimme",
):
    """Return gas-phase SCF energy and qh-G correction in Hartree."""
    thermo = Thermochemistry.from_filepath(
        filepath,
        **_pka_thermochemistry_kwargs(
            temperature,
            concentration,
            pressure,
            cutoff_entropy_grimme,
            cutoff_enthalpy,
            entropy_method,
        ),
    )
    electronic_energy_j_mol = _require_thermochemistry_property(
        thermo,
        filepath,
        "electronic_energy",
        "SCF energy",
    )
    qh_gibbs_j_mol = _require_thermochemistry_property(
        thermo,
        filepath,
        "qrrho_gibbs_free_energy",
        "quasi-harmonic Gibbs free energy",
    )
    electronic_energy_au = energy_conversion(
        "j/mol", "hartree", electronic_energy_j_mol
    )
    qh_gibbs_au = energy_conversion("j/mol", "hartree", qh_gibbs_j_mol)
    return electronic_energy_au, qh_gibbs_au - electronic_energy_au


def pka_solvent_scf_energy(filepath):
    """Return solvent-phase SCF energy in Hartree."""
    thermo = Thermochemistry.from_filepath(filepath)
    electronic_energy_j_mol = _require_thermochemistry_property(
        thermo,
        filepath,
        "electronic_energy",
        "SCF energy",
    )
    return energy_conversion("j/mol", "hartree", electronic_energy_j_mol)


def compute_pka(
    ha_gas_file,
    a_gas_file,
    href_gas_file=None,
    ref_gas_file=None,
    ha_solv_file=None,
    a_solv_file=None,
    href_solv_file=None,
    ref_solv_file=None,
    pka_reference=None,
    temperature=298.15,
    concentration=1.0,
    pressure=1.0,
    cutoff_entropy_grimme=100.0,
    cutoff_enthalpy=100.0,
    entropy_method="grimme",
    scheme="proton exchange",
    delta_G_proton=None,
):
    """Compute pKa from output files using program-independent thermochemistry."""
    if scheme == "direct":
        if delta_G_proton is None:
            raise ValueError(
                "delta_G_proton is required when scheme='direct'."
            )
    elif pka_reference is None:
        raise ValueError(
            "pka_reference is required when scheme='proton exchange'."
        )
    else:
        missing = [
            name
            for name, value in (
                ("href_gas_file", href_gas_file),
                ("ref_gas_file", ref_gas_file),
                ("ha_solv_file", ha_solv_file),
                ("a_solv_file", a_solv_file),
                ("href_solv_file", href_solv_file),
                ("ref_solv_file", ref_solv_file),
            )
            if value is None
        ]
        if missing:
            raise ValueError(
                "Missing required files for proton exchange scheme: "
                + ", ".join(missing)
            )

    if scheme == "direct" and (ha_solv_file is None or a_solv_file is None):
        raise ValueError(
            "ha_solv_file and a_solv_file are required for scheme='direct'."
        )

    thermo_kwargs = dict(
        temperature=temperature,
        concentration=concentration,
        pressure=pressure,
        cutoff_entropy_grimme=cutoff_entropy_grimme,
        cutoff_enthalpy=cutoff_enthalpy,
        entropy_method=entropy_method,
    )

    E_gas_HA_au, G_corr_HA_au = pka_gas_phase_data(
        ha_gas_file, **thermo_kwargs
    )
    E_gas_A_au, G_corr_A_au = pka_gas_phase_data(a_gas_file, **thermo_kwargs)
    E_solv_HA_au = pka_solvent_scf_energy(ha_solv_file)
    E_solv_A_au = pka_solvent_scf_energy(a_solv_file)
    G_soln_HA_au = E_solv_HA_au + G_corr_HA_au
    G_soln_A_au = E_solv_A_au + G_corr_A_au

    R_kcal = 0.001987204
    ln10 = 2.302585093

    if scheme == "direct":
        G_soln_HA_kcal = G_soln_HA_au * HARTREE_TO_KCAL_MOL
        G_soln_A_kcal = G_soln_A_au * HARTREE_TO_KCAL_MOL
        delta_G_diss_kcal_mol = G_soln_A_kcal + delta_G_proton - G_soln_HA_kcal
        delta_G_diss_au = delta_G_diss_kcal_mol / HARTREE_TO_KCAL_MOL
        pka = delta_G_diss_kcal_mol / (R_kcal * temperature * ln10)
        return {
            "pKa": pka,
            "scheme": "direct",
            "delta_G_proton_kcal_mol": delta_G_proton,
            "delta_G_diss_kcal_mol": delta_G_diss_kcal_mol,
            "delta_G_diss_au": delta_G_diss_au,
            "delta_G_soln_kcal_mol": delta_G_diss_kcal_mol,
            "delta_G_soln_au": delta_G_diss_au,
            "temperature": temperature,
            "G_soln_HA_au": G_soln_HA_au,
            "G_soln_A_au": G_soln_A_au,
            "E_solv_HA_au": E_solv_HA_au,
            "E_solv_A_au": E_solv_A_au,
            "G_corr_HA_au": G_corr_HA_au,
            "G_corr_A_au": G_corr_A_au,
            "E_gas_HA_au": E_gas_HA_au,
            "E_gas_A_au": E_gas_A_au,
        }

    E_gas_HRef_au, G_corr_HRef_au = pka_gas_phase_data(
        href_gas_file, **thermo_kwargs
    )
    E_gas_Ref_au, G_corr_Ref_au = pka_gas_phase_data(
        ref_gas_file, **thermo_kwargs
    )
    E_solv_HRef_au = pka_solvent_scf_energy(href_solv_file)
    E_solv_Ref_au = pka_solvent_scf_energy(ref_solv_file)
    G_soln_HRef_au = E_solv_HRef_au + G_corr_HRef_au
    G_soln_Ref_au = E_solv_Ref_au + G_corr_Ref_au

    delta_G_soln_au = (G_soln_A_au + G_soln_HRef_au) - (
        G_soln_HA_au + G_soln_Ref_au
    )
    delta_G_soln_kcal_mol = delta_G_soln_au * HARTREE_TO_KCAL_MOL
    pka = pka_reference + delta_G_soln_kcal_mol / (R_kcal * temperature * ln10)

    return {
        "pKa": pka,
        "scheme": "proton exchange",
        "pKa_reference": pka_reference,
        "delta_G_soln_kcal_mol": delta_G_soln_kcal_mol,
        "delta_G_soln_au": delta_G_soln_au,
        "temperature": temperature,
        "G_soln_HA_au": G_soln_HA_au,
        "G_soln_A_au": G_soln_A_au,
        "G_soln_HRef_au": G_soln_HRef_au,
        "G_soln_Ref_au": G_soln_Ref_au,
        "E_solv_HA_au": E_solv_HA_au,
        "E_solv_A_au": E_solv_A_au,
        "E_solv_HRef_au": E_solv_HRef_au,
        "E_solv_Ref_au": E_solv_Ref_au,
        "G_corr_HA_au": G_corr_HA_au,
        "G_corr_A_au": G_corr_A_au,
        "G_corr_HRef_au": G_corr_HRef_au,
        "G_corr_Ref_au": G_corr_Ref_au,
        "E_gas_HA_au": E_gas_HA_au,
        "E_gas_A_au": E_gas_A_au,
        "E_gas_HRef_au": E_gas_HRef_au,
        "E_gas_Ref_au": E_gas_Ref_au,
    }


def _thermochemistry_value_in_units(
    thermo, filepath, attr, energy_units, label
):
    value_j_mol = _require_thermochemistry_property(
        thermo, filepath, attr, label
    )
    return energy_conversion("j/mol", energy_units, value_j_mol)


def compute_pka_thermochemistry(
    ha_file=None,
    a_file=None,
    href_file=None,
    ref_file=None,
    temperature=298.15,
    concentration=1.0,
    pressure=1.0,
    cutoff_entropy_grimme=100.0,
    cutoff_enthalpy=100.0,
    energy_units="hartree",
    entropy_method="grimme",
):
    """Extract gas-phase thermochemistry for pKa species from output files."""
    results = {
        "settings": {
            "temperature": temperature,
            "concentration": concentration,
            "pressure": pressure,
            "cutoff_entropy_grimme": cutoff_entropy_grimme,
            "cutoff_enthalpy": cutoff_enthalpy,
            "energy_units": energy_units,
        }
    }
    thermo_kwargs = _pka_thermochemistry_kwargs(
        temperature,
        concentration,
        pressure,
        cutoff_entropy_grimme,
        cutoff_enthalpy,
        entropy_method,
    )

    def get_species_thermo(file_path, name):
        if file_path is None:
            return None
        thermo = Thermochemistry.from_filepath(file_path, **thermo_kwargs)
        return {
            "name": name,
            "E": _thermochemistry_value_in_units(
                thermo,
                file_path,
                "electronic_energy",
                energy_units,
                "SCF energy",
            ),
            "qh_G": _thermochemistry_value_in_units(
                thermo,
                file_path,
                "qrrho_gibbs_free_energy",
                energy_units,
                "quasi-harmonic Gibbs free energy",
            ),
            "ZPE": _thermochemistry_value_in_units(
                thermo,
                file_path,
                "zero_point_energy",
                energy_units,
                "zero-point energy",
            ),
            "H": _thermochemistry_value_in_units(
                thermo, file_path, "enthalpy", energy_units, "enthalpy"
            ),
            "qh_H": _thermochemistry_value_in_units(
                thermo,
                file_path,
                "qrrho_enthalpy",
                energy_units,
                "quasi-harmonic enthalpy",
            ),
            "G": _thermochemistry_value_in_units(
                thermo,
                file_path,
                "gibbs_free_energy",
                energy_units,
                "Gibbs free energy",
            ),
        }

    if ha_file is not None:
        results["HA"] = get_species_thermo(ha_file, "HA")
    if a_file is not None:
        results["A"] = get_species_thermo(a_file, "A-")
    if href_file is not None:
        results["HRef"] = get_species_thermo(href_file, "HRef")
    if ref_file is not None:
        results["Ref"] = get_species_thermo(ref_file, "Ref-")
    return results


def print_pka_summary(
    ha_gas_file,
    a_gas_file,
    href_gas_file=None,
    ref_gas_file=None,
    ha_solv_file=None,
    a_solv_file=None,
    href_solv_file=None,
    ref_solv_file=None,
    pka_reference=None,
    temperature=298.15,
    concentration=1.0,
    pressure=1.0,
    cutoff_entropy_grimme=100.0,
    cutoff_enthalpy=100.0,
    entropy_method="grimme",
    scheme="proton exchange",
    delta_G_proton=None,
):
    """Print a formatted summary of a dual-level pKa calculation."""
    result = compute_pka(
        ha_gas_file=ha_gas_file,
        a_gas_file=a_gas_file,
        href_gas_file=href_gas_file,
        ref_gas_file=ref_gas_file,
        ha_solv_file=ha_solv_file,
        a_solv_file=a_solv_file,
        href_solv_file=href_solv_file,
        ref_solv_file=ref_solv_file,
        pka_reference=pka_reference,
        temperature=temperature,
        concentration=concentration,
        pressure=pressure,
        cutoff_entropy_grimme=cutoff_entropy_grimme,
        cutoff_enthalpy=cutoff_enthalpy,
        entropy_method=entropy_method,
        scheme=scheme,
        delta_G_proton=delta_G_proton,
    )

    if scheme == "direct":
        print("=" * 78)
        print("pKa Calculation - Direct Dissociation Scheme")
        print("=" * 78)
        print("Reaction: HA → A⁻ + H⁺")
        print(f"Temperature: {temperature} K")
        print()
        print("Method:")
        print("  G_corr = qh-G(T) - E_gas  (from gas-phase freq calculation)")
        print("  G_soln = E_solv + G_corr  (solution free energy)")
        print("  ΔG_diss = G_soln(A⁻) + G_soln(H⁺) - G_soln(HA)")
        print("  pKa = ΔG_diss / (2.303 × R × T)")
        print("-" * 78)
        print()
        print("Gas-Phase Electronic Energies (E_gas, au):")
        print(f"  HA:  {result['E_gas_HA_au']:.10f}")
        print(f"  A⁻:  {result['E_gas_A_au']:.10f}")
        print()
        print("Thermal Corrections (G_corr = qh-G - E_gas, au):")
        print(f"  HA:  {result['G_corr_HA_au']:.10f}")
        print(f"  A⁻:  {result['G_corr_A_au']:.10f}")
        print()
        print("Solvent Single-Point Energies (E_solv, au):")
        print(f"  HA:  {result['E_solv_HA_au']:.10f}")
        print(f"  A⁻:  {result['E_solv_A_au']:.10f}")
        print()
        print("Solution Free Energies (G_soln = E_solv + G_corr, au):")
        print(f"  HA:  {result['G_soln_HA_au']:.10f}")
        print(f"  A⁻:  {result['G_soln_A_au']:.10f}")
        print("-" * 78)
        print()
        print("pKa Calculation:")
        print(f"  G_soln(H⁺) = {delta_G_proton:.4f} kcal/mol")
        print(f"  ΔG_diss = {result['delta_G_diss_au']:.10f} au")
        print(f"         = {result['delta_G_diss_kcal_mol']:.4f} kcal/mol")
        print()
        print(f"  *** Computed pKa(HA) = {result['pKa']:.2f} ***")
        print("=" * 78)
        return

    print("=" * 78)
    print("pKa Calculation - Dual-level Proton Exchange Scheme")
    print("=" * 78)
    print("Reaction: HA + Ref⁻ → A⁻ + HRef")
    print(f"Temperature: {temperature} K")
    print()
    print("Method:")
    print("  G_corr = qh-G(T) - E_gas  (from gas-phase freq calculation)")
    print("  G_soln = E_solv + G_corr  (solution free energy)")
    print(
        "  ΔG_soln = [G(A⁻)_soln + G(HRef)_soln] - [G(HA)_soln + G(Ref⁻)_soln]"
    )
    print("  pKa = pKa_ref + ΔG_soln / (RT × ln10)")
    print("-" * 78)
    print()
    print("Gas-Phase Electronic Energies (E_gas, au):")
    print(f"  HA:  {result['E_gas_HA_au']:.10f}")
    print(f"  A⁻:  {result['E_gas_A_au']:.10f}")
    print(f"  HRef:  {result['E_gas_HRef_au']:.10f}")
    print(f"  Ref⁻:  {result['E_gas_Ref_au']:.10f}")
    print()
    print("Thermal Corrections (G_corr = qh-G - E_gas, au):")
    print(f"  HA:  {result['G_corr_HA_au']:.10f}")
    print(f"  A⁻:  {result['G_corr_A_au']:.10f}")
    print(f"  HRef:  {result['G_corr_HRef_au']:.10f}")
    print(f"  Ref⁻:  {result['G_corr_Ref_au']:.10f}")
    print()
    print("Solvent Single-Point Energies (E_solv, au):")
    print(f"  HA:  {result['E_solv_HA_au']:.10f}")
    print(f"  A⁻:  {result['E_solv_A_au']:.10f}")
    print(f"  HRef:  {result['E_solv_HRef_au']:.10f}")
    print(f"  Ref⁻:  {result['E_solv_Ref_au']:.10f}")
    print()
    print("Solution Free Energies (G_soln = E_solv + G_corr, au):")
    print(f"  HA:  {result['G_soln_HA_au']:.10f}")
    print(f"  A⁻:  {result['G_soln_A_au']:.10f}")
    print(f"  HRef:  {result['G_soln_HRef_au']:.10f}")
    print(f"  Ref⁻:  {result['G_soln_Ref_au']:.10f}")
    print("-" * 78)
    print()
    print("pKa Calculation:")
    print(f"  ΔG_soln = {result['delta_G_soln_au']:.10f} au")
    print(f"         = {result['delta_G_soln_kcal_mol']:.4f} kcal/mol")
    print(f"  pKa(HRef)_ref = {pka_reference:.2f}")
    print()
    print(f"  *** Computed pKa(HA) = {result['pKa']:.2f} ***")
    print("=" * 78)


def click_pka_thermochemistry_options(f):
    """Thermochemistry options reused by pKa submission and analysis."""
    f = thermochemistry_temp_pressure_conc_options(
        f,
        temperature_required=False,
        temperature_default=298.15,
        concentration_default=1.0,
        pressure_default=1.0,
        concentration_short="-c",
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


def is_pka_cdxml_input(filename):
    """Return True when *filename* is a ChemDraw CDX/CDXML structure file."""
    return bool(filename) and str(filename).lower().endswith(
        (".cdx", ".cdxml")
    )


def apply_pka_molecule_charge_multiplicity(opt_settings, molecule):
    """Use charge/multiplicity from *molecule* when absent on *opt_settings*."""
    import copy

    updated = copy.copy(opt_settings)
    mol_charge = molecule.charge
    mol_mult = molecule.multiplicity
    if updated.charge is None and mol_charge is not None:
        updated.charge = int(mol_charge)
    if updated.multiplicity is None and mol_mult is not None:
        updated.multiplicity = int(mol_mult)
    return updated


def require_pka_charge_multiplicity(opt_settings, source_hint=""):
    """Raise when charge or multiplicity are still unset after resolution."""
    missing = []
    if opt_settings.charge is None:
        missing.append("-c/--charge")
    if opt_settings.multiplicity is None:
        missing.append("-m/--multiplicity")
    if not missing:
        return
    suffix = f" ({source_hint})" if source_hint else ""
    raise click.UsageError(
        "Charge and multiplicity are required for pKa submission. "
        f"Missing: {', '.join(missing)}. "
        "Provide them on the parent command or use a ChemDraw structure "
        f"from which they can be inferred{suffix}."
    )


def is_pka_batch_invocation(ctx):
    """Return True when the nested ``pka`` command targets ``batch`` mode."""
    if getattr(ctx, "invoked_subcommand", None) != "pka":
        return False

    tokens = []
    current = ctx
    while current is not None:
        if getattr(current, "invoked_subcommand", None) == "pka":
            tokens.extend(str(token) for token in (current.args or []))
        current = current.parent

    if "submit" in tokens:
        return False
    return "batch" in tokens


def resolve_pka_batch_row(filepath, proton_index=None, color_code=None):
    """Resolve proton index and molecule for one pKa submission-table row.

    Each table row maps to a single job. When ``proton_index`` is omitted for a
    single-molecule ``.cdxml`` / ``.cdx`` filepath, the coloured proton is
    auto-detected. An explicit ``proton_index`` always takes precedence. Multi-
    molecule CDXML files are rejected here; pass them directly as ``-f`` with
    ``pka batch`` instead.

    Returns:
        tuple[int, Molecule | PKaMolecule]: Resolved proton index and structure.
    """
    from chemsmart.io.file import PKaCDXFile
    from chemsmart.io.molecules.structure import Molecule

    filepath = str(filepath)
    if proton_index is not None:
        return int(proton_index), Molecule.from_filepath(filepath)

    if not is_pka_cdxml_input(filepath):
        raise ValueError(
            f"Missing proton_index for {filepath}. "
            "Provide proton_index in the table, or use a single-molecule "
            ".cdxml/.cdx file with a coloured proton and leave proton_index blank."
        )

    cdx_file = PKaCDXFile(filepath)
    try:
        pka_molecules = cdx_file.get_pka_molecules(
            color_code=color_code,
            index=":",
            return_list=True,
        )
    except ValueError as exc:
        raise ValueError(
            f"Could not auto-detect proton from CDXML colour for {filepath}: "
            f"{exc}"
        ) from exc

    if len(pka_molecules) != 1:
        raise ValueError(
            f"CDXML file {filepath} contains {len(pka_molecules)} molecules. "
            "Submission-table rows support single-molecule CDXML files only. "
            "Pass a multi-molecule CDXML file directly as -f with pka batch."
        )

    pka_mol = pka_molecules[0]
    return pka_mol.proton_index, pka_mol


def batch_pka_jobs_from_cdxml(
    ctx,
    skip_completed,
    create_jobs_fn,
    invoke_submit_fn,
    **kwargs,
):
    """Create pKa jobs from a CDXML batch input via coloured-proton detection."""
    from chemsmart.io.file import PKaCDXFile

    filename = ctx.obj.get("filename")
    shared = ctx.obj["pka_shared"]
    proton_index, color_code = resolve_pka_submit_proton_options(ctx)
    try:
        proton_index, pka_molecules = PKaCDXFile.resolve_proton_index(
            filename, proton_index, color_code
        )
    except ValueError as exc:
        raise click.UsageError(str(exc)) from exc

    if pka_molecules is not None:
        return create_jobs_fn(
            ctx, pka_molecules, shared, skip_completed, **kwargs
        )

    return invoke_submit_fn(
        ctx,
        skip_completed=skip_completed,
        proton_index=proton_index,
        color_code=color_code,
        **kwargs,
    )


def resolve_pka_submit_proton_options(ctx, proton_index=None, color_code=None):
    """Resolve proton options for ``pka submit`` from multiple Click scopes.

    The ``pka`` group captures ``-pi/--proton-index``, but nested invocation
    via ``chemsmart run`` does not always populate ``ctx.obj`` from group-level
    options.  ``submit`` therefore re-declares the same options and this helper
    merges submit args, parent group params, and ``ctx.obj``.
    """
    parent = ctx.parent
    if proton_index is None and parent is not None:
        proton_index = parent.params.get("proton_index")
    if proton_index is None:
        proton_index = ctx.obj.get("pka_proton_index")

    if color_code is None and parent is not None:
        color_code = parent.params.get("color_code")
    if color_code is None:
        color_code = ctx.obj.get("pka_color_code")

    ctx.obj["pka_proton_index"] = proton_index
    ctx.obj["pka_color_code"] = color_code
    return proton_index, color_code


def click_pka_proton_options(f):
    """Options that identify the proton to remove (-pi / -cc).

    Applied to the ``pka`` group and ``submit`` so values are available whether
    options appear before or after the ``submit`` token (required for per-row
    ``chemsmart run`` scripts generated by ``chemsmart sub`` batch mode).
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
    from chemsmart.utils.utils import (
        PKA_TARGET_SUFFIX_HELP,
        discover_pka_target_companion_outputs,
    )

    results = discover_pka_target_companion_outputs(
        ha_gas_path, program=program
    )

    missing = [
        f"  {k}: {v}" for k, v in results.items() if not os.path.isfile(v)
    ]
    if missing:
        raise click.UsageError(
            "Auto-discovery could not find some companion output files.\n"
            "Missing files:\n" + "\n".join(missing) + "\n\n"
            "Provide them explicitly or ensure output files follow:\n"
            + PKA_TARGET_SUFFIX_HELP
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


def validate_reference_options(shared):
    reference = shared["reference"]
    if reference is None:
        return

    if shared["scheme"] != "proton exchange":
        raise click.UsageError(
            "Reference acid file can only be used with 'proton exchange' "
            "cycle. Use -s 'proton exchange' or remove the -r option."
        )

    shared["reference_proton_index"] = PKaCDXFile.resolve_reference_proton(
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


def _validate_pka_table_program(pka_table, program):
    """Ensure explicit -p matches every output file in the table."""
    output_fields = (
        "ha_gas",
        "a_gas",
        "ha_sp",
        "a_sp",
        "href_gas",
        "ref_gas",
        "href_sp",
        "ref_sp",
    )
    for entry in pka_table.entries:
        for field in output_fields:
            path = entry.get(field)
            if not _is_existing_output_path(path):
                continue
            detected = get_program_type_from_file(str(path))
            if detected != program:
                raise click.UsageError(
                    f"File '{path}' was detected as {detected!r}, but "
                    f"batch-analyze was run with -p {program}."
                )


def _is_existing_output_path(value):
    """Return True when *value* is a non-empty path to an existing file."""
    from chemsmart.utils.utils import normalize_table_cell

    path = normalize_table_cell(value)
    if path is None:
        return False
    return os.path.isfile(str(path))


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

    solv_names = [f"{name}-solv" for name in file_names]
    missing_files = []
    for name, path in zip(file_names, required_gas):
        if path is not None and not os.path.isfile(path):
            missing_files.append(f"  --{name}: {path}")
    for name, path in zip(solv_names, required_solv):
        if path is not None and not os.path.isfile(path):
            missing_files.append(f"  --{name}: {path}")
    if missing_files:
        raise click.UsageError(
            "One or more pKa analysis files do not exist:\n"
            + "\n".join(missing_files)
        )


def _auto_discover_pka_files(ha_gas_path, href_gas_path, program=None):
    """Infer companion output paths from HA and HRef gas-phase paths."""
    from chemsmart.utils.utils import (
        PKA_REFERENCE_SUFFIX_HELP,
        PKA_TARGET_SUFFIX_HELP,
        discover_pka_reference_companion_outputs,
        discover_pka_target_companion_outputs,
    )

    if program is None:
        program = get_program_type_from_file(ha_gas_path)

    results = discover_pka_target_companion_outputs(
        ha_gas_path, program=program
    )
    results.update(discover_pka_reference_companion_outputs(href_gas_path))

    missing = [
        f"  {k}: {v}" for k, v in results.items() if not os.path.isfile(v)
    ]
    if missing:
        raise click.UsageError(
            "Auto-discovery could not find some companion output files.\n"
            "Missing files:\n" + "\n".join(missing) + "\n\n"
            "Provide them explicitly or ensure output files follow:\n"
            + PKA_TARGET_SUFFIX_HELP
            + "\n"
            + PKA_REFERENCE_SUFFIX_HELP
        )
    return results


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
    auto-discovered from the same naming convention as batch-analyze:
      <basename>_pka_A_opt.<ext>   conjugate base
      <basename>_pka_HA_sp.<ext>   HA solvent single-point
      <basename>_pka_A_sp.<ext>    conjugate base solvent SP
      (and the corresponding _pka_Ref_* / _pka_HRef_sp files for HRef)
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

        logger.info("Computing pKa (Direct Dissociation)...")
        print_pka_summary(
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

    logger.info(f"Computing pKa ({_scheme_display_name(scheme)})...")
    print_pka_summary(
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
    help="Path to write the formatted results report.  Stdout if omitted.",
)
@click.option(
    "-p",
    "--program",
    type=click.Choice(["gaussian", "orca", "auto"]),
    default="auto",
    show_default=True,
    help=(
        "Require every populated output path to match this backend.  "
        "'auto' (default) parses each file independently and supports "
        "mixed Gaussian/ORCA tables."
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

    if program != "auto":
        _validate_pka_table_program(pka_output_table, program)

    program_label = "auto" if program == "auto" else program

    logger.info(
        f"Computing pKa ({_scheme_display_name(scheme)}) "
        f"for {len(pka_output_table)} systems "
        f"(T={shared['temperature']}K, program={program_label})"
    )
    results = pka_output_table.run_pka(
        output_cls=compute_pka,
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
