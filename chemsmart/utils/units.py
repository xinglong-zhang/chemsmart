"""Conversions for units."""
from ase.units import create_units

__codata_version__ = '2022'  # here we use 2022 version of CODATA
ase_units = create_units(__codata_version__)

class CustomUnits(ase_units):
    """Subclass of ase.units.Units with additional units and
    constants that are not present in ase."""

    # Generate type hints dynamically for ASE unit attributes
    #'_c', '_mu0', '_Grav', '_hplanck', '_e', '_me', '_mp', '_Nav', '_k', '_amu'


    __annotations__ = {key: float for key in ase_units.keys()}  # Add type hints

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # explicitly declare parent class attributes
        for k, v in ase_units.items():
            setattr(self, k, v)
        #
        # self._c = ase_units['_c']
        # self._mu0 = ase_units['_mu0']
        # self._Grav = ase_units['_Grav']
        # self._hplanck = ase_units['_hplanck']
        # self._e = ase_units['_e']
        # self._me = ase_units['_me']
        # self._mp = ase_units['_mp']
        # self._Nav = ase_units['_Nav']
        # self._k = ase_units['_k']
        # self._amu = ase_units['_amu']

        # units derived from the CODATA values
        self._eps0 = ase_units['_eps0']
        self._hbar = ase_units['_hbar']
        self.Ang = ase_units['Ang']
        self.Angstrom = ase_units['Angstrom']
        self.nm = ase_units['nm']
        self.Bohr = ase_units['Bohr']
        self.eV = ase_units['eV']
        self.Hartree = ase_units['Hartree']
        self.kJ = ase_units['kJ']
        self.kcal = ase_units['kcal']
        self.mol = ase_units['mol']
        self.Rydberg = ase_units['Rydberg']
        self.Ry = ase_units['Ry']
        self.Ha = ase_units['Ha']
        self.second = ase_units['second']
        self.fs = ase_units['fs']
        self.kB = ase_units['kB']

        #     u['_eps0'] = (1 / u['_mu0'] / u['_c']**2)  # permittivity of vacuum
        #     u['_hbar'] = u['_hplanck'] / (2 * pi)  # Planck constant / 2pi, J s
        #
        #     u['Ang'] = u['Angstrom'] = 1.0
        #     u['nm'] = 10.0
        #     u['Bohr'] = (4e10 * pi * u['_eps0'] * u['_hbar']**2 /
        #                  u['_me'] / u['_e']**2)  # Bohr radius
        #
        #     u['eV'] = 1.0
        #     u['Hartree'] = (u['_me'] * u['_e']**3 / 16 / pi**2 /
        #                     u['_eps0']**2 / u['_hbar']**2)
        #     u['kJ'] = 1000.0 / u['_e']
        #     u['kcal'] = 4.184 * u['kJ']
        #     u['mol'] = u['_Nav']
        #     u['Rydberg'] = 0.5 * u['Hartree']
        #     u['Ry'] = u['Rydberg']
        #     u['Ha'] = u['Hartree']
        #
        #     u['second'] = 1e10 * sqrt(u['_e'] / u['_amu'])
        #     u['fs'] = 1e-15 * u['second']
        #
        #     u['kB'] = u['_k'] / u['_e']  # Boltzmann constant, eV/K
        #
        #     u['Pascal'] = (1 / u['_e']) / 1e30  # J/m^3
        #     u['GPa'] = 1e9 * u['Pascal']
        #     u['bar'] = 1e5 * u['Pascal']
        #
        #     u['Debye'] = 1.0 / 1e11 / u['_e'] / u['_c']
        #     u['alpha'] = (u['_e']**2 / (4 * pi * u['_eps0']) /
        #                   u['_hbar'] / u['_c'])  # fine structure constant
        #     u['invcm'] = (100 * u['_c'] * u['_hplanck'] /
        #                   u['_e'])  # cm^-1 energy unit
        #
        #     # Derived atomic units that have no assigned name:
        #     # atomic unit of time, s:
        #     u['_aut'] = u['_hbar'] / (u['alpha']**2 * u['_me'] * u['_c']**2)
        #     # atomic unit of velocity, m/s:
        #     u['_auv'] = u['_e']**2 / u['_hbar'] / (4 * pi * u['_eps0'])
        #     # atomic unit of force, N:
        #     u['_auf'] = u['alpha']**3 * u['_me']**2 * u['_c']**3 / u['_hbar']
        #     # atomic unit of pressure, Pa:
        #     u['_aup'] = u['alpha']**5 * u['_me']**4 * u['_c']**5 / u['_hbar']**3
        #
        #     u['AUT'] = u['second'] * u['_aut']
        #
        #     # SI units
        #     u['m'] = 1e10 * u['Ang']  # metre
        #     u['kg'] = 1. / u['_amu']  # kilogram
        #     u['s'] = u['second']  # second
        #     u['A'] = 1.0 / u['_e'] / u['s']  # ampere
        #     # derived
        #     u['J'] = u['kJ'] / 1000  # Joule = kg * m**2 / s**2
        #     u['C'] = 1.0 / u['_e']  # Co


        # self.data = ase_units.CODATA[__codata_version__]  # Store CODATA constants separately
        # Additional unit conversions and constants
        self.atm_to_pa = 101325.0  # Atmosphere to Pascal
        self.R = 8.314462618  # Gas constant, J/(mol K)

        # Overriding or explicitly defining fundamental constants
        self.kB = 1.380649e-23  # Boltzmann constant, J/K
        self.e = 1.602176634e-19  # Elementary charge, C
        self.hplanck = 6.62607015e-34  # Planck constant, J s
        self.me = 9.10938356e-31  # Electron mass, kg
        self.mp = 1.6726219e-27  # Proton mass, kg
        self.Nav = 6.02214076e23  # Avogadro number, mol^-1
        self.amu = 1.66053906660e-27  # Atomic mass unit, kg
        self.c = 299792458.0  # Speed of light, m/s
        self.mu0 = 4e-7 * 3.141592653589793  # Permeability of vacuum

    # def __getattr__(self, name):
    #     """Fallback to ASE's parent class if the attribute isn't found."""
    #     try:
    #         data = object.__getattribute__(self, "data")  # Avoid recursion
    #         if name in data:
    #             return data[name]  # Fetch from CODATA dictionary
    #     except AttributeError:
    #         pass
    #
    #     # Check if the attribute exists in the ASE base class
    #     if hasattr(super(), name):
    #         return getattr(super(), name)  # Inherit from ase.units.Units
    #
    #     raise AttributeError(f"'{self.__class__.__name__}' object has no attribute '{name}'")

    # def __getattr__(self, name):
    #     """Ensure ASE attributes are accessible and recognized by IDEs."""
    #     if name in self.data:
    #         return self.data[name]  # Fetch from ASE parent class dictionary
    #     raise AttributeError(f"'{self.__class__.__name__}' object has no attribute '{name}'")
    #

units = CustomUnits(ase_units.CODATA[ase_units.__codata_version__])
print(ase_units.CODATA[ase_units.__codata_version__])
print(units._hplanck)
print(ase_units.CODATA[__codata_version__].keys())
