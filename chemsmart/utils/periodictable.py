from ase.data import chemical_symbols as elements


class PeriodicTable:
    PERIODIC_TABLE = [str(element) for element in elements]

    def to_element(self, element_str):
        """Function to convert lower case element to element with proper caps in periodic table."""
        return element_str.upper() if len(element_str) == 1 else f'{element_str[0].upper()}{element_str[1:]}'

    def sorted_periodic_table_list(self, list_of_elements):
        return sorted(list_of_elements, key=lambda x: self.PERIODIC_TABLE.index(x))

    def to_atomic_number(self, symbol):
        return self.PERIODIC_TABLE.index(symbol)

    def to_symbol(self, atomic_number):
        element_str = self.PERIODIC_TABLE[atomic_number]
        return self.to_element(element_str)
