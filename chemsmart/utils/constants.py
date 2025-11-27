"""
Physical constants and energy unit conversions.

Convenience constants and a small, explicit energy unit conversion helper
that complement ASE's units. Includes common pressure/length/mass factors
and robust conversions among Hartree, eV, kcal/mol, kJ/mol, and J/mol.

Highlights:
- `atm_to_pa`, `R`, `bohr_to_meter`, `amu_to_kg`, `hartree_to_joules`.
- `energy_conversion(from_unit, to_unit, value)` to convert between energy
  units (case-insensitive).
"""

import logging

from ase import units

logger = logging.getLogger(__name__)


# Pressure and fundamental constants
atm_to_pa = 101325  # 1 atm = 101325 Pa
R = units._k * units._Nav  # Ideal gas constant
bohr_to_meter = (
    1 * units.Bohr / units.m
)  # 1 Bohr = 0.52917721067 Å = 0.52917721067e-10 m
amu_to_kg = 1 * units._amu  # 1 amu = 1.66053906660e-27 kg
hartree_to_joules = 4.35974434e-18  # 1 Hartree = 4.35974434 × 10^-18 Joules
cal_to_joules = 4.184  # 1 Calorie = 4.184 Joules

# Conversion factors for energy units
joule_per_mol_to_eV = 1.0364269574711572e-05  # J/mol to eV
joule_per_mol_to_kcal_per_mol = 1 / 4184  # J/mol to kcal/mol
joule_per_mol_to_kJ_per_mol = 0.001  # J/mol to kJ/mol
joule_per_mol_to_hartree = 1 / (
    hartree_to_joules * units._Nav
)  # J/mol to Hartree
# joule_per_mol_to_hartree = 3.8087991196914175e-07


def energy_conversion(from_unit, to_unit, value=1.0):
    """
    Convert an energy value between supported units (case-insensitive).

    Supports: Hartree, eV, kcal/mol, kJ/mol, and J/mol. Unit strings are
    matched case-insensitively; canonical lower-case keys are
    "hartree", "ev", "kcal/mol", "kj/mol", and "j/mol".

    Args:
        from_unit (str): Unit to convert from. Options: 'hartree', 'eV',
            'kcal/mol', 'kJ/mol', 'J/mol' (case-insensitive).
        to_unit (str): Unit to convert to. Options: 'hartree', 'eV',
            'kcal/mol', 'kJ/mol', 'J/mol' (case-insensitive).
        value (float, optional): Energy value to convert. Defaults to 1.0
            (i.e., returns the conversion factor).

    Returns:
        float or None: Converted value, or None if `value` is None.

    Raises:
        ValueError: If either unit is unsupported.
    """
    if value is None:
        return None
    valid_units = ["hartree", "ev", "kcal/mol", "kj/mol", "j/mol"]
    if from_unit.lower() not in valid_units:
        raise ValueError(
            f"Unsupported from_unit: {from_unit}. Choose from {valid_units}"
        )
    if to_unit.lower() not in valid_units:
        raise ValueError(
            f"Unsupported to_unit: {to_unit}. Choose from {valid_units}"
        )

    # Conversion factors to J/mol (SI units used internally in Thermochemistry)
    to_j_per_mol = {
        "hartree": hartree_to_joules * units._Nav,  # Hartree to J/mol
        "ev": 1 / joule_per_mol_to_eV,  # eV to J/mol
        "kcal/mol": 4184.0,  # kcal/mol to J/mol (1 kcal = 4184 J)
        "kj/mol": 1000.0,  # kJ/mol to J/mol
        "j/mol": 1.0,  # J/mol to J/mol
    }

    # Convert input value to J/mol
    value_in_j_per_mol = value * to_j_per_mol[from_unit.lower()]

    # Convert from J/mol to target unit
    conversion_factor = to_j_per_mol[to_unit.lower()]
    converted_value = value_in_j_per_mol / conversion_factor

    logger.debug(
        f"Converted {value} {from_unit} to {converted_value:.6f} {to_unit}"
    )
    return converted_value
