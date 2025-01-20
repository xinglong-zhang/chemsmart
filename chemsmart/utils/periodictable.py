from ase.data import chemical_symbols as elements
from ase.data.vdw import vdw_radii
from ase.data import covalent_radii


class PeriodicTable:
    PERIODIC_TABLE = [str(element) for element in elements]

    def to_element(self, element_str):
        """Function to convert lower case element to element with proper caps in periodic table."""
        # if element_str.upper() == "TV":
        #     pass
        return (
            element_str.upper()
            if len(element_str) == 1
            else f"{element_str[0].upper()}{element_str[1:]}"
        )

    def sorted_periodic_table_list(self, list_of_elements):
        """Sort the list of elements according to atomic number."""
        return sorted(
            list_of_elements, key=lambda x: self.PERIODIC_TABLE.index(x)
        )

    def to_atomic_number(self, symbol):
        # if symbol.upper() == 'TV':
        #     pass
        return self.PERIODIC_TABLE.index(symbol)

    def to_symbol(self, atomic_number):
        element_str = self.PERIODIC_TABLE[atomic_number]
        return self.to_element(element_str)

    def vdw_radius(self, symbol):
        # obtain the van der waals radius of the element
        return vdw_radii[self.to_atomic_number(symbol)]

    def covalent_radius(self, symbol):
        # obtain the covalent radius of the element
        return covalent_radii[self.to_atomic_number(symbol)]
