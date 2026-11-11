import re
from functools import cached_property

from chemsmart.utils.mixins import FileMixin
from chemsmart.utils.repattern import raw_energy_value_pattern
from chemsmart.utils.utils import string2index_1based


class XYZFile(FileMixin):
    """
    Parser for XYZ coordinate files.

    This class handles XYZ files containing molecular coordinates, supporting
    both single-molecule files and multi-molecule files (e.g., conformer
    ensembles). It can extract molecular geometries, comments, and energies
    from XYZ file formats.
    """

    def __init__(self, filename):
        """
        Initialize XYZ file parser.

        Args:
            filename (str): Path to the XYZ file to parse
        """
        self.filename = filename

    def __repr__(self):
        return f"XYZFile({self.filename})"

    def __str__(self):
        return f"XYZFile object with filename: {self.filename}"

    @cached_property
    def num_atoms(self):
        """
        Get number of atoms from the first line of XYZ file.

        Returns:
            int: Number of atoms in the molecular structure
        """
        return int(self.contents[0])

    @cached_property
    def molecule(self):
        """
        Get the last molecule from the XYZ file.

        Returns:
            Molecule: Last molecular structure in the file
        """
        return self.get_molecules(index="-1")

    @cached_property
    def comments(self):
        """
        Get the last comment line from the XYZ file.

        Returns:
            str: Last comment line in the file
        """
        return self.get_comments(index="-1")

    def iter_molecule_results(self, index=":"):
        """Read frames sequentially, yielding original one-based positions.

        The selection syntax matches ``get_molecules``. Negative positions
        require a counting pass, but neither pass stores all frames.
        """
        from chemsmart.io.molecules.structure import Molecule
        from chemsmart.utils.utils import parse_index_specification

        selector = parse_index_specification(str(index))
        needs_count = (
            isinstance(selector, int)
            and selector < 0
            or isinstance(selector, list)
            and any(i < 0 for i in selector)
            or isinstance(selector, slice)
            and (
                selector.step is not None
                and selector.step < 0
                or selector.start is not None
                and selector.start < 0
                or selector.stop is not None
                and selector.stop < 0
            )
        )
        count = (
            sum(1 for _ in self._iter_frame_blocks()) if needs_count else None
        )
        if count is not None:
            selected = parse_index_specification(str(index), total_count=count)
        else:
            selected = selector
        if isinstance(selected, int):
            wanted = {selected}
        elif isinstance(selected, list):
            wanted = set(selected)
        else:
            wanted = None

        matched = False
        for offset, comment, block in self._iter_frame_blocks():
            if wanted is not None:
                include = offset in wanted
            elif selected.step is not None and selected.step < 0:
                include = offset in range(*selected.indices(count))
            else:
                start = selected.start or 0
                stop = selected.stop
                step = selected.step or 1
                include = (
                    offset >= start
                    and (stop is None or offset < stop)
                    and (offset - start) % step == 0
                )
            if not include:
                continue
            matched = True
            try:
                molecule = Molecule.from_coordinate_block_text(
                    block, structure_index_in_file=offset + 1
                )
                match = re.findall(raw_energy_value_pattern, comment)
                if match:
                    molecule.energy = float(match[0])
                yield offset + 1, molecule, None
            except (ValueError, TypeError, IndexError) as error:
                yield offset + 1, None, error
        if not matched:
            raise IndexError(
                f"Index {index} selects no structures in {self.filename}."
            )

    def _iter_frame_blocks(self):
        """Yield one XYZ frame at a time without loading the full file."""
        with open(self.filename) as stream:
            offset = 0
            while header := stream.readline():
                atoms = int(header.strip())
                if atoms <= 0:
                    raise ValueError(
                        "Number of atoms in the xyz file is zero!"
                    )
                comment = stream.readline()
                if not comment:
                    raise ValueError(f"Incomplete XYZ frame {offset + 1}.")
                block = [stream.readline() for _ in range(atoms)]
                if any(not line for line in block):
                    raise ValueError(f"Incomplete XYZ frame {offset + 1}.")
                yield offset, comment.strip(), block
                offset += 1

    def _get_molecules_and_comments(self, index=":", return_list=False):
        """
        Extract molecules and comments from XYZ file.

        The XYZ file can contain a single molecule (conventional format) or
        multiple molecules (e.g., conformer ensembles from CREST).

        Args:
            index (str): Index specification for molecule selection
            return_list (bool): Whether to always return lists

        Returns:
            tuple: (molecules, comments) - molecules and their comment lines

        Raises:
            ValueError: If number of atoms is zero
        """
        from chemsmart.io.molecules.structure import Molecule

        all_molecules = []
        comments = []
        i = 0
        frame_idx = 0  # frame counter across the entire file (assigned to molecules as 1-based)
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

            # Assign frame index (1-based) where the structure appears in file
            frame_idx += 1
            molecule = Molecule.from_coordinate_block_text(
                coordinate_block, structure_index_in_file=frame_idx
            )

            # Store the molecule data
            all_molecules.append(molecule)

        molecules = all_molecules[string2index_1based(index)]
        comments = comments[string2index_1based(index)]
        if return_list and not isinstance(molecules, list):
            return [molecules], [comments]
        else:
            return molecules, comments

    def get_molecules(self, index=":", return_list=False):
        """
        Extract molecular structures from XYZ file with energy assignment.

        Parses comment lines to extract energy values and assigns them to
        the corresponding Molecule objects.

        Args:
            index (str): Index specification for molecule selection
            return_list (bool): Whether to return list format

        Returns:
            Molecule or list: Single molecule
            or list of molecules with energies
        """
        # Ensure that when return_list=False, molecules
        # is always treated as a list before iteration:
        molecules, comments = self._get_molecules_and_comments(
            index=index, return_list=True
        )

        # Ensures energy is assigned before returning a single molecule:
        if len(comments) != 0:
            for i, comment in enumerate(comments):
                # will extract the first float number in the line.
                # example case 1: "Empirical formula:
                # C191H241Cu2N59O96P14 Energy(Hartree): -25900.214629"
                # energy will be -25900.214629.
                match = re.findall(raw_energy_value_pattern, comment)
                if match:
                    molecules[i].energy = float(match[0])
                    # Assign energy to the only or
                    # the first negative float number
                else:
                    # No energy found, skip
                    continue
        if return_list:
            return molecules
        else:
            return (
                molecules[0] if molecules else None
            )  # Return a single molecule if list has one item

    def get_comments(self, index=":", return_list=False):
        """
        Extract comment lines from XYZ file.

        Args:
            index (str): Index specification for comment selection
            return_list (bool): Whether to return list format

        Returns:
            str or list: Single comment or list of comments
        """
        _, comments = self._get_molecules_and_comments(
            index=index, return_list=return_list
        )
        return comments
