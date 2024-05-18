from pymatgen.io.gaussian import GaussianOutput
from chemsmart.utils.mixins import FileMixin
from functools import cached_property


class Gaussian16Output(FileMixin):
    def __init__(self, filename):
        self.filename = filename
        self.gaussian_object = GaussianOutput(
            filename=filename
        )  # Gaussian output object from pymatgen

    @cached_property
    def tddft_transitions(self):
        """
        Read a excitation energies after a TD-DFT calculation.

        Returns:
            A list: A list of tuple for each transition such as
                    [(energie (eV), lambda (nm), oscillatory strength), ... ]
        """
        return self.gaussian_object.read_excitation_energies()

    @cached_property
    def excitation_energies_eV(self):
        """
        Read TDDFT transitions and return the transition energies in eV as a list
        """
        excitation_energies_eV = []
        for i in self.tddft_transitions:
            excitation_energies_eV.append(i[0])
        return excitation_energies_eV
