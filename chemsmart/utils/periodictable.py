"""
Periodic table utilities for element properties and conversions.

Provides a comprehensive interface for accessing element data including
atomic numbers, masses, radii, and isotope information. Integrates with
ASE data sources and ChemSmart isotope data for complete element handling.
"""

from ase.data import chemical_symbols as elements
from ase.data import covalent_radii
from ase.data.vdw import vdw_radii

from chemsmart.utils.isotopes_data import isotopes


class PeriodicTable:
    """
    Periodic table interface for element data and conversions.

    Provides methods for converting between element symbols, atomic numbers,
    and various atomic properties. Includes access to atomic masses,
    radii, and isotope data for all elements.
    """

    PERIODIC_TABLE = [str(element) for element in elements]

    @property
    def atomic_masses(self):
        """
        Get atomic masses for all elements.

        Returns:
            numpy.ndarray: Array of atomic masses indexed by atomic number.
        """
        from ase.data import atomic_masses_iupac2016 as atomic_masses

        return atomic_masses

    def to_element(self, element_str):
        """
        Normalize an element symbol's capitalization.

        Returns a symbol with the first character uppercased. For multi-
        letter symbols, the remainder is preserved as provided (i.e., this
        does not force lowercase on the second letter).

        Args:
            element_str (str): Element symbol in any case format.

        Returns:
            str: Properly capitalized element symbol.
        """
        # if element_str.upper() == "TV":
        #     pass
        return (
            element_str.upper()
            if len(element_str) == 1
            else f"{element_str[0].upper()}{element_str[1:]}"
        )

    def sorted_periodic_table_list(self, list_of_elements):
        """
        Sort elements by atomic number order.

        Arranges a list of element symbols according to their position
        in the periodic table (atomic number order).

        Args:
            list_of_elements (list[str]): Element symbols to sort (case-sensitive;
                ensure canonical capitalization as in `PERIODIC_TABLE`).

        Returns:
            list[str]: Element symbols sorted by atomic number.
        """
        return sorted(
            list_of_elements, key=lambda x: self.PERIODIC_TABLE.index(x)
        )

    def to_atomic_number(self, symbol):
        """
        Convert element symbol to atomic number.

        Args:
            symbol (str): Element symbol (e.g., 'H', 'He', 'Li'). Case-sensitive;
                use `to_element` beforehand to normalize if needed.

        Returns:
            int: Atomic number of the element.
        """
        # if symbol.upper() == 'TV':
        #     pass
        return self.PERIODIC_TABLE.index(symbol)

    def to_symbol(self, atomic_number):
        """
        Convert atomic number to element symbol.

        Args:
            atomic_number (int): Atomic number of the element.

        Returns:
            str: Properly formatted element symbol.
        """
        element_str = self.PERIODIC_TABLE[atomic_number]
        return self.to_element(element_str)

    def to_atomic_mass(self, symbol):
        """
        Get the standard atomic mass of an element.

        Args:
            symbol (str): Element symbol.

        Returns:
            float: Standard atomic mass in atomic mass units.
        """
        return self.atomic_masses[self.to_atomic_number(symbol)]

    def to_weighted_atomic_mass_by_abundance(self, symbol):
        """
        Get abundance-weighted atomic mass of an element.

        Calculates atomic mass based on natural isotope abundances.

        Args:
            symbol (str): Element symbol.

        Returns:
            float: Abundance-weighted atomic mass in atomic mass units.
        """
        return isotopes[self.to_atomic_number(symbol)]["weighted_atomic_mass"]

    def to_most_abundant_atomic_mass(self, symbol):
        """
        Get atomic mass of the most abundant isotope.

        Args:
            symbol (str): Element symbol.

        Returns:
            float: Atomic mass of the most abundant isotope.
        """
        return isotopes[self.to_atomic_number(symbol)]["most_abundant"]["mass"]

    def vdw_radius(self, symbol):
        """
        Get van der Waals radius of an element.

        Args:
            symbol (str): Element symbol.

        Returns:
            float: Van der Waals radius in Angstroms.
        """
        # obtain the van der waals radius of the element
        return vdw_radii[self.to_atomic_number(symbol)]

    def covalent_radius(self, symbol):
        """
        Get covalent radius of an element.

        Args:
            symbol (str): Element symbol.

        Returns:
            float: Covalent radius in Angstroms.
        """
        # obtain the covalent radius of the element
        return covalent_radii[self.to_atomic_number(symbol)]
