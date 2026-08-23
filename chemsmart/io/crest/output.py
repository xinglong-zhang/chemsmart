import logging
from functools import cached_property

from chemsmart.io.crest.file import CRESTEnergiesFile, CRESTMainOut
from chemsmart.io.crest.folder import CRESTFolder
from chemsmart.io.molecules.structure import Molecule
from chemsmart.utils.constants import kcal_per_mol_to_hartree
from chemsmart.utils.mixins import FolderOutputMixin

logger = logging.getLogger(__name__)


class CRESTOutput(FolderOutputMixin):
    """
    Integrates and exposes all relevant data from a CREST calculation directory.

    This class discovers and parses all relevant CREST output files in a given
    folder, then provides convenient access to conformer ensembles, energies,
    and run metadata.

    Args:
        folder (str or CRESTFolder): Path to CREST calculation folder or
            CRESTFolder instance.
    """

    FILE_PARSERS = (
        "main_out",
        "energies_file",
    )

    def __init__(self, folder):
        self.folder = (
            folder if isinstance(folder, CRESTFolder) else CRESTFolder(folder)
        )

    @cached_property
    def main_out(self):
        """Main CREST output file parser (CRESTMainOut)."""
        path = self.folder.crest_out_filepath
        return CRESTMainOut(path) if path else None

    @cached_property
    def energies_file(self):
        """Energies file parser (CRESTEnergiesFile)."""
        path = self.folder.energies_filepath
        return CRESTEnergiesFile(path) if path else None

    @property
    def normal_termination(self):
        """Whether CREST terminated normally."""
        if self.main_out is not None:
            return self.main_out.normal_termination
        return False

    @cached_property
    def conformers(self):
        """List of conformer Molecule objects from crest_conformers.xyz.

        Returns:
            list[Molecule]: All conformers sorted by energy, or empty list.
        """
        path = self.folder.conformers_xyz_filepath
        if path is None:
            return []
        try:
            molecules = Molecule.from_filepath(
                filepath=path, index=":", return_list=True
            )
            for molecule in molecules:
                molecule.charge = self.charge
                molecule.multiplicity = self.multiplicity
            return molecules
        except Exception as exc:
            logger.error(f"Error reading conformers from {path}: {exc}")
            return []

    @cached_property
    def rotamers(self):
        """List of rotamer Molecule objects from crest_rotamers.xyz.

        Rotamers are the larger pre-CREGEN ensemble sampled during the CREST
        run. This set is typically much larger than the final conformer list.

        Returns:
            list[Molecule]: All rotamers with energies from the xyz comment
                lines, or empty list.
        """
        path = self.folder.rotamers_xyz_filepath
        if path is None:
            return []
        try:
            molecules = Molecule.from_filepath(
                filepath=path, index=":", return_list=True
            )
            for molecule in molecules:
                molecule.charge = self.charge
                molecule.multiplicity = self.multiplicity
            return molecules
        except Exception as exc:
            logger.error(f"Error reading rotamers from {path}: {exc}")
            return []

    @cached_property
    def best_conformer(self):
        """Lowest-energy conformer Molecule from crest_best.xyz.

        Returns:
            Molecule or None: The best conformer, or None if unavailable.
        """
        path = self.folder.best_xyz_filepath
        if path is None:
            if self.conformers:
                return self.conformers[0]
            return None
        try:
            molecule = Molecule.from_filepath(filepath=path, index="-1")
            molecule.charge = self.charge
            molecule.multiplicity = self.multiplicity
            return molecule
        except Exception as exc:
            logger.error(f"Error reading best conformer from {path}: {exc}")
            return None

    @cached_property
    def energies(self):
        """Absolute conformer energies in Hartree.

        Returns:
            list[float]: Absolute conformer energies in Hartree, or empty list.
        """
        conformer_energies = [
            mol.energy for mol in self.conformers if mol.energy is not None
        ]
        if conformer_energies and len(conformer_energies) == len(
            self.conformers
        ):
            return conformer_energies

        if (
            self.energies_file is not None
            and self.main_out is not None
            and self.main_out.lowest_energy is not None
        ):
            return [
                self.main_out.lowest_energy
                + rel_energy * kcal_per_mol_to_hartree
                for rel_energy in self.energies_file.relative_energies
            ]
        return []

    @cached_property
    def num_conformers(self):
        """Number of unique conformers found.

        Returns:
            int: Number of conformers.
        """
        if (
            self.main_out is not None
            and self.main_out.num_conformers is not None
        ):
            return self.main_out.num_conformers
        if self.energies_file is not None:
            return self.energies_file.num_conformers
        return len(self.conformers)

    @cached_property
    def num_rotamers(self):
        """Number of unique rotamers sampled before CREGEN filtering.

        Returns:
            int: Number of rotamers.
        """
        if (
            self.main_out is not None
            and self.main_out.num_rotamers is not None
        ):
            return self.main_out.num_rotamers
        return len(self.rotamers)
