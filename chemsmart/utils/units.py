"""Conversions for units."""
from ase.units import create_units as ase_create_units

__codata_version__ = '2022'  # here we use 2022 version of CODATA


class CustomUnits(ase_units.Units):
    """Subclass of ase.units.Units with additional units and
    constants that are not present in ase."""

    # Generate type hints dynamically for ASE unit attributes
    #'_c', '_mu0', '_Grav', '_hplanck', '_e', '_me', '_mp', '_Nav', '_k', '_amu'
    _hplanck: float

    __annotations__ = {key: float for key in ase_units.CODATA[__codata_version__].keys()}  # Add type hints

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
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
