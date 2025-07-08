"""Constants not found in ase units."""

from ase import units

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
