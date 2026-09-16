"""Dual-level exchange redox analysis from output files.

Post-processes finished oxidized / reduced target and reference calculations.
Job submission lives under ``chemsmart sub … gaussian … redox``,
``chemsmart sub … orca … redox``, or ``chemsmart sub … chain … redox``;
analysis is exposed as ``chemsmart run redox analyze``.
"""

from __future__ import annotations

from dataclasses import dataclass

from chemsmart.analysis.thermochemistry import (
    gas_phase_data,
    solvent_scf_energy,
)
from chemsmart.utils.constants import FARADAY, energy_conversion

_REGISTRY: dict[str, RedoxReference] = {}


def require_positive_n_electrons(n_electrons):
    """Return ``n_electrons`` after requiring a positive integer."""
    if n_electrons < 1:
        raise ValueError("n_electrons must be a positive integer.")
    return n_electrons


def resolve_n_electrons(n_electrons, reference=None, default=None):
    """Resolve electron count from an explicit value, couple, or default."""
    if n_electrons is None:
        if reference is not None:
            n_electrons = reference.n_electrons
        elif default is not None:
            n_electrons = default
        else:
            raise ValueError("n_electrons is required.")
    elif reference is not None and n_electrons != reference.n_electrons:
        raise ValueError(
            f"n_electrons ({n_electrons}) must match the reference "
            f"couple {reference.name!r} (n={reference.n_electrons})."
        )
    return require_positive_n_electrons(n_electrons)


@dataclass(frozen=True)
class RedoxReference:
    """One experimental or computational redox reference couple.

    Attributes:
        name: Registry key, e.g. ``fc_fc+``.
        E_ref_V: Reference reduction potential in volts.
        n_electrons: Electrons transferred in the couple.
        scale: Reporting scale, e.g. ``Fc/Fc+`` or ``SHE``.
        couple_label: Chemical couple, e.g. ``Fc/Fc+``.
        ox_file: Optional path to the oxidized geometry.
        red_file: Optional path to the reduced geometry.
        ox_charge: Optional charge of the oxidized species.
        ox_multiplicity: Optional multiplicity of the oxidized species.
        red_charge: Optional charge of the reduced species.
        red_multiplicity: Optional multiplicity of the reduced species.
        solvent: Optional solvent label associated with ``E_ref_V``.
    """

    name: str
    E_ref_V: float
    n_electrons: int
    scale: str
    couple_label: str
    ox_file: str | None = None
    red_file: str | None = None
    ox_charge: int | None = None
    ox_multiplicity: int | None = None
    red_charge: int | None = None
    red_multiplicity: int | None = None
    solvent: str | None = None

    def __post_init__(self):
        if not self.name:
            raise ValueError("RedoxReference.name is required.")
        require_positive_n_electrons(self.n_electrons)


def register_redox_reference(reference):
    """Register a :class:`RedoxReference` under ``reference.name``.

    Args:
        reference (RedoxReference): Couple to add to the registry.

    Raises:
        TypeError: If ``reference`` is not a :class:`RedoxReference`.
        ValueError: If ``reference.name`` is already registered.
    """
    if not isinstance(reference, RedoxReference):
        raise TypeError(
            "reference must be a RedoxReference, "
            f"got {type(reference).__name__}."
        )
    if reference.name in _REGISTRY:
        raise ValueError(
            f"Redox reference {reference.name!r} is already registered."
        )
    _REGISTRY[reference.name] = reference


def get_redox_reference(name):
    """Return the registered :class:`RedoxReference` for ``name``.

    Args:
        name (str): Registry key.

    Returns:
        RedoxReference: The registered couple.

    Raises:
        ValueError: If ``name`` is not registered.
    """
    reference = _REGISTRY.get(name)
    if reference is None:
        available = ", ".join(sorted(_REGISTRY)) or "(none)"
        raise ValueError(
            f"Unknown redox reference {name!r}. Available: {available}."
        )
    return reference


def list_redox_references():
    """Return registered couples sorted by name."""
    return tuple(_REGISTRY[name] for name in sorted(_REGISTRY))


def resolve_redox_reference(reference=None):
    """Return a :class:`RedoxReference` from an object, name, or default.

    Args:
        reference: A :class:`RedoxReference`, registry name, or ``None``
            (defaults to ``fc_fc+``).

    Returns:
        RedoxReference: Resolved couple.

    Raises:
        TypeError: If ``reference`` is not a supported type.
    """
    if reference is None:
        return get_redox_reference("fc_fc+")
    if isinstance(reference, RedoxReference):
        return reference
    if isinstance(reference, str):
        return get_redox_reference(reference)
    raise TypeError(
        "reference must be a RedoxReference, registry name, or None, "
        f"got {type(reference).__name__}."
    )


def _species_terms(gas_file, solv_file, thermo_kwargs):
    """Return gas, correction, solvent, and solution G terms in Hartree."""
    e_gas_au, g_corr_au = gas_phase_data(gas_file, **thermo_kwargs)
    e_solv_au = solvent_scf_energy(solv_file)
    return {
        "E_gas_au": e_gas_au,
        "G_corr_au": g_corr_au,
        "E_solv_au": e_solv_au,
        "G_soln_au": e_solv_au + g_corr_au,
    }


def _require_files(**files):
    missing = [name for name, path in files.items() if not path]
    if missing:
        raise ValueError(
            "Missing required files for redox exchange analysis: "
            + ", ".join(missing)
        )


def compute_redox_potential(
    ox_gas_file,
    red_gas_file,
    ref_ox_gas_file,
    ref_red_gas_file,
    ox_solv_file,
    red_solv_file,
    ref_ox_solv_file,
    ref_red_solv_file,
    reference=None,
    e_ref=None,
    n_electrons=None,
    temperature=298.15,
    concentration=1.0,
    pressure=1.0,
    cutoff_entropy_grimme=100.0,
    cutoff_enthalpy=100.0,
    entropy_method="grimme",
):
    """Compute the exchange redox potential from output files.

    Dual-level solution free energies follow the pKa helpers:
    ``G_soln = E_solv + (qh-G − E_gas)``.

    Reaction: ``Ox + Ref_red → Red + Ref_ox``

    ``ΔG_exchange = G(Red) + G(Ref_ox) − G(Ox) − G(Ref_red)``

    ``E_target = E_ref − ΔG_exchange / (n F)``

    Args:
        ox_gas_file: Gas-phase opt+freq output for the oxidized target.
        red_gas_file: Gas-phase opt+freq output for the reduced target.
        ref_ox_gas_file: Gas-phase opt+freq output for the oxidized reference.
        ref_red_gas_file: Gas-phase opt+freq output for the reduced reference.
        ox_solv_file: Solution-phase SP output for the oxidized target.
        red_solv_file: Solution-phase SP output for the reduced target.
        ref_ox_solv_file: Solution-phase SP output for the oxidized reference.
        ref_red_solv_file: Solution-phase SP output for the reduced reference.
        reference: :class:`RedoxReference` or registry name. Used for
            ``E_ref`` when ``e_ref`` is omitted, and for ``n`` / labels.
        e_ref: Reference reduction potential in volts. Required when
            ``reference`` is omitted.
        n_electrons: Electrons transferred. Defaults to the reference
            couple, or 1 when only ``e_ref`` is given. Must match the
            couple when both are given.
        temperature: Temperature in Kelvin for quasi-harmonic G.
        concentration: Concentration in mol/L.
        pressure: Pressure in atm.
        cutoff_entropy_grimme: Entropy cutoff in cm⁻¹ (Grimme qRRHO).
        cutoff_enthalpy: Enthalpy cutoff in cm⁻¹.
        entropy_method: Quasi-RRHO entropy method (``grimme`` or ``truhlar``).

    Returns:
        dict: Exchange ΔG (au, kcal/mol, eV, J/mol), ``E_ref_V``,
        ``E_target_V``, reference metadata, and per-species G terms.
    """
    _require_files(
        ox_gas_file=ox_gas_file,
        red_gas_file=red_gas_file,
        ref_ox_gas_file=ref_ox_gas_file,
        ref_red_gas_file=ref_red_gas_file,
        ox_solv_file=ox_solv_file,
        red_solv_file=red_solv_file,
        ref_ox_solv_file=ref_ox_solv_file,
        ref_red_solv_file=ref_red_solv_file,
    )
    if e_ref is None and reference is None:
        raise ValueError(
            "--e-ref is required for analyze (reference reduction "
            "potential in volts), or pass a registry couple."
        )
    resolved = (
        resolve_redox_reference(reference) if reference is not None else None
    )
    e_ref_v = resolved.E_ref_V if e_ref is None else float(e_ref)
    n_electrons = resolve_n_electrons(
        n_electrons, reference=resolved, default=1
    )

    thermo_kwargs = dict(
        temperature=temperature,
        concentration=concentration,
        pressure=pressure,
        cutoff_entropy_grimme=cutoff_entropy_grimme,
        cutoff_enthalpy=cutoff_enthalpy,
        entropy_method=entropy_method,
    )
    ox = _species_terms(ox_gas_file, ox_solv_file, thermo_kwargs)
    red = _species_terms(red_gas_file, red_solv_file, thermo_kwargs)
    ref_ox = _species_terms(ref_ox_gas_file, ref_ox_solv_file, thermo_kwargs)
    ref_red = _species_terms(
        ref_red_gas_file, ref_red_solv_file, thermo_kwargs
    )

    delta_g_exchange_au = (
        red["G_soln_au"]
        + ref_ox["G_soln_au"]
        - ox["G_soln_au"]
        - ref_red["G_soln_au"]
    )
    delta_g_exchange_j_mol = energy_conversion(
        "hartree", "j/mol", delta_g_exchange_au
    )
    delta_g_exchange_kcal_mol = energy_conversion(
        "hartree", "kcal/mol", delta_g_exchange_au
    )
    delta_g_exchange_ev = energy_conversion(
        "hartree", "ev", delta_g_exchange_au
    )
    e_target_v = e_ref_v - (delta_g_exchange_j_mol / (n_electrons * FARADAY))
    return {
        "delta_G_exchange_au": delta_g_exchange_au,
        "delta_G_exchange_kcal_mol": delta_g_exchange_kcal_mol,
        "delta_G_exchange_eV": delta_g_exchange_ev,
        "delta_G_exchange_J_mol": delta_g_exchange_j_mol,
        "E_ref_V": e_ref_v,
        "E_target_V": e_target_v,
        "scale": resolved.scale if resolved is not None else "",
        "couple_label": resolved.couple_label if resolved is not None else "",
        "reference_name": resolved.name if resolved is not None else None,
        "n_electrons": n_electrons,
        "temperature": temperature,
        "G_soln_ox_au": ox["G_soln_au"],
        "G_soln_red_au": red["G_soln_au"],
        "G_soln_ref_ox_au": ref_ox["G_soln_au"],
        "G_soln_ref_red_au": ref_red["G_soln_au"],
        "E_solv_ox_au": ox["E_solv_au"],
        "E_solv_red_au": red["E_solv_au"],
        "E_solv_ref_ox_au": ref_ox["E_solv_au"],
        "E_solv_ref_red_au": ref_red["E_solv_au"],
        "G_corr_ox_au": ox["G_corr_au"],
        "G_corr_red_au": red["G_corr_au"],
        "G_corr_ref_ox_au": ref_ox["G_corr_au"],
        "G_corr_ref_red_au": ref_red["G_corr_au"],
        "E_gas_ox_au": ox["E_gas_au"],
        "E_gas_red_au": red["E_gas_au"],
        "E_gas_ref_ox_au": ref_ox["E_gas_au"],
        "E_gas_ref_red_au": ref_red["E_gas_au"],
    }


def format_redox_summary(result):
    """Return a formatted summary of an exchange redox calculation."""
    if result.get("reference_name"):
        reference_line = (
            f"Reference: {result['reference_name']} "
            f"({result['couple_label']}, {result['scale']})"
        )
    else:
        reference_line = f"Reference: E_ref = {result['E_ref_V']:.4f} V"
    scale = result.get("scale") or ""
    scale_suffix = f" ({scale})" if scale else ""
    lines = [
        "=" * 78,
        "Redox Potential - Dual-level Exchange Scheme",
        "=" * 78,
        "Reaction: Ox + Ref_red → Red + Ref_ox",
        reference_line,
        f"n = {result['n_electrons']}",
        f"Temperature: {result['temperature']} K",
        "",
        "Method:",
        "  G_corr = qh-G(T) - E_gas  (from gas-phase freq calculation)",
        "  G_soln = E_solv + G_corr  (solution free energy)",
        "  ΔG_exchange = G(Red) + G(Ref_ox) − G(Ox) − G(Ref_red)",
        "  E_target = E_ref − ΔG_exchange / (n F)",
        "-" * 78,
        "",
        "Gas-Phase Electronic Energies (E_gas, au):",
        f"  Ox:      {result['E_gas_ox_au']:.10f}",
        f"  Red:     {result['E_gas_red_au']:.10f}",
        f"  Ref_ox:  {result['E_gas_ref_ox_au']:.10f}",
        f"  Ref_red: {result['E_gas_ref_red_au']:.10f}",
        "",
        "Thermal Corrections (G_corr = qh-G - E_gas, au):",
        f"  Ox:      {result['G_corr_ox_au']:.10f}",
        f"  Red:     {result['G_corr_red_au']:.10f}",
        f"  Ref_ox:  {result['G_corr_ref_ox_au']:.10f}",
        f"  Ref_red: {result['G_corr_ref_red_au']:.10f}",
        "",
        "Solvent Single-Point Energies (E_solv, au):",
        f"  Ox:      {result['E_solv_ox_au']:.10f}",
        f"  Red:     {result['E_solv_red_au']:.10f}",
        f"  Ref_ox:  {result['E_solv_ref_ox_au']:.10f}",
        f"  Ref_red: {result['E_solv_ref_red_au']:.10f}",
        "",
        "Solution Free Energies (G_soln = E_solv + G_corr, au):",
        f"  Ox:      {result['G_soln_ox_au']:.10f}",
        f"  Red:     {result['G_soln_red_au']:.10f}",
        f"  Ref_ox:  {result['G_soln_ref_ox_au']:.10f}",
        f"  Ref_red: {result['G_soln_ref_red_au']:.10f}",
        "-" * 78,
        "",
        "Redox Potential:",
        f"  ΔG_exchange = {result['delta_G_exchange_au']:.10f} au",
        f"              = {result['delta_G_exchange_kcal_mol']:.4f} kcal/mol",
        f"              = {result['delta_G_exchange_eV']:.4f} eV",
        f"              = {result['delta_G_exchange_J_mol']:.4f} J/mol",
        f"  E_ref = {result['E_ref_V']:.4f} V{scale_suffix}",
        "",
        f"  *** E_target = {result['E_target_V']:.4f} V{scale_suffix} ***",
        "=" * 78,
    ]
    return "\n".join(lines)


register_redox_reference(
    RedoxReference(
        name="fc_fc+",
        E_ref_V=0.0,
        n_electrons=1,
        scale="Fc/Fc+",
        couple_label="Fc/Fc+",
    )
)
