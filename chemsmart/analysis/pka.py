"""Dual-level pKa analysis from output files.

Post-processes finished HA / A− (and optional HRef / Ref−) calculations.
Job submission lives under ``chemsmart sub … gaussian … pka``,
``chemsmart sub … orca … pka``, or ``chemsmart sub … chain … pka``;
analysis is exposed as ``chemsmart run pka``.
"""

from chemsmart.analysis.thermochemistry import (
    Thermochemistry,
    _extract_thermochemistry_property,
    _thermochemistry_kwargs,
    gas_phase_data,
    solvent_scf_energy,
)
from chemsmart.utils.constants import HARTREE_TO_KCAL_MOL, energy_conversion

pka_gas_phase_data = gas_phase_data
pka_solvent_scf_energy = solvent_scf_energy


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

    E_gas_HA_au, G_corr_HA_au = gas_phase_data(ha_gas_file, **thermo_kwargs)
    E_gas_A_au, G_corr_A_au = gas_phase_data(a_gas_file, **thermo_kwargs)
    E_solv_HA_au = solvent_scf_energy(ha_solv_file)
    E_solv_A_au = solvent_scf_energy(a_solv_file)
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

    E_gas_HRef_au, G_corr_HRef_au = gas_phase_data(
        href_gas_file, **thermo_kwargs
    )
    E_gas_Ref_au, G_corr_Ref_au = gas_phase_data(ref_gas_file, **thermo_kwargs)
    E_solv_HRef_au = solvent_scf_energy(href_solv_file)
    E_solv_Ref_au = solvent_scf_energy(ref_solv_file)
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
    value_j_mol = _extract_thermochemistry_property(
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
    thermo_kwargs = _thermochemistry_kwargs(
        temperature,
        concentration,
        pressure,
        cutoff_entropy_grimme,
        cutoff_enthalpy,
        entropy_method,
        energy_units=energy_units,
    )

    def get_species_thermo(filepath, name):
        if filepath is None:
            return None
        thermo = Thermochemistry(filename=filepath, **thermo_kwargs)
        return {
            "name": name,
            "E": _thermochemistry_value_in_units(
                thermo,
                filepath,
                "electronic_energy",
                energy_units,
                "SCF energy",
            ),
            "qh_G": _thermochemistry_value_in_units(
                thermo,
                filepath,
                "qrrho_gibbs_free_energy",
                energy_units,
                "quasi-harmonic Gibbs free energy",
            ),
            "ZPE": _thermochemistry_value_in_units(
                thermo,
                filepath,
                "zero_point_energy",
                energy_units,
                "zero-point energy",
            ),
            "H": _thermochemistry_value_in_units(
                thermo, filepath, "enthalpy", energy_units, "enthalpy"
            ),
            "qh_H": _thermochemistry_value_in_units(
                thermo,
                filepath,
                "qrrho_enthalpy",
                energy_units,
                "quasi-harmonic enthalpy",
            ),
            "G": _thermochemistry_value_in_units(
                thermo,
                filepath,
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
