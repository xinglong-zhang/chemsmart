import re
from functools import cached_property

from chemsmart.utils.mixins import FileMixin
from chemsmart.utils.utils import string2index_1based


class XYZFile(FileMixin):
    """xyz file object."""

    def __init__(self, filename):
        self.filename = filename

    def __repr__(self):
        return f"XYZFile({self.filename})"

    def __str__(self):
        return f"XYZFile object with filename: {self.filename}"

    @cached_property
    def num_atoms(self):
        return int(self.contents[0])

    @cached_property
    def molecule(self):
        return self.get_molecules(index="-1")

    @cached_property
    def comments(self):
        return self.get_comments(index="-1")

    def _get_molecules_and_comments(self, index=":", return_list=False):
        """Return a molecule object or a list of molecule objects from an xyz file.
        The xzy file can either contain a single molecule, as conventionally, or a list
        of molecules, such as those in crest_conformers.xyz file."""
        from chemsmart.io.molecules.structure import Molecule

        all_molecules = []
        comments = []
        i = 0
        while i < len(self.contents):
            # Read number of atoms
            num_atoms = int(self.contents[i].strip())
            if num_atoms == 0:
                raise ValueError("Number of atoms in the xyz file is zero!")
            i += 1
            # Read comment line
            comment = self.contents[i].strip()
            comments.append(comment)
            i += 1
            # Read the coordinate block
            coordinate_block = self.contents[i : i + num_atoms]
            i += num_atoms
            molecule = Molecule.from_coordinate_block_text(coordinate_block)

            # Store the molecule data
            all_molecules.append(molecule)

        molecules = all_molecules[string2index_1based(index)]
        comments = comments[string2index_1based(index)]
        if return_list and not isinstance(molecules, list):
            return [molecules], [comments]
        else:
            return molecules, comments

    def get_molecules(self, index=":", return_list=False):
        # Ensure that when return_list=False, molecules is always treated as a list before iteration:
        molecules, comments = self._get_molecules_and_comments(
            index=index, return_list=True
        )

        # Ensures energy is assigned before returning a single molecule:
        if len(comments) != 0:
            for i, comment in enumerate(comments):
                energy = None
                try:
                    energy = float(comment)
                except ValueError:
                    # example case:
                    # "Empirical formula: C191H241Cu2N59O96P14    Energy(Hartree): -25900.214629"
                    match = re.search(r"(-?\d+\.\d+)", comment)
                    if match:
                        energy = float(match.group(1))
                    else:
                        # No energy found, skip
                        continue
                if energy:
                    molecules[i].energy = energy  # Assign energy

        if return_list:
            return molecules
        else:
            return (
                molecules[0] if molecules else None
            )  # Return a single molecule if list has one item

    def get_comments(self, index=":", return_list=False):
        _, comments = self._get_molecules_and_comments(
            index=index, return_list=return_list
        )
        return comments
