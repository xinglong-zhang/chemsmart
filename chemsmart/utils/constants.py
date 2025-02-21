"""Constants not found in ase units."""

from ase import units

atm_to_pa = 101325  # 1 atm = 101325 Pa
R = units._k * units._Nav  # Ideal gas constant
bohr_to_meter = (
    1 * units.Bohr / units.m
)  # 1 Bohr = 0.52917721067 Å = 0.52917721067e-10 m
amu_to_kg = 1 * units._amu  # 1 amu = 1.66053906660e-27 kg
