"""Conversions for units."""
from ase.units import create_units as ase_create_units

__codata_version__ = '2022'  # here we use 2022 version of CODATA


# create a dictionary of units of <class 'ase.units.Units'>
ase_units = ase_create_units(__codata_version__)

# add additional units and constants
ase_units['atm_to_pa'] = 101325.0  # Atmosphere to Pascal
ase_units['R'] = 8.314462618  # Gas constant, J/(mol K)

# explicitly define fundamental constants, in SI units
ase_units['kB'] = 1.380649e-23  # Boltzmann constant, J/K
ase_units['e'] = 1.602176634e-19  # Elementary charge, C
ase_units['_hplanck'] = 6.62607015e-34  # Planck constant, J s
ase_units['me'] = 9.10938356e-31  # Electron mass, kg
ase_units['mp'] = 1.6726219e-27  # Proton mass, kg
ase_units['Nav'] = 6.02214076e23  # Avogadro number, mol^-1
ase_units['amu'] = 1.66053906660e-27  # Atomic mass unit, kg
ase_units['c'] = 299792458.0  # Speed of light, m/s
ase_units['mu0'] = 4e-7 * 3.141592653589793  # Permeability of vacuum = 1.25663706127(20)×10−6 N⋅A−2
