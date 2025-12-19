import logging
import re

import numpy as np

from chemsmart.io.molecules.structure import Molecule
from chemsmart.utils.mixins import FileMixin

logger = logging.getLogger(__name__)


class SDFFile(FileMixin):
    """
    SDF file object.
    """

    def __init__(self, filename):
        self.filename = filename

    @property
    def molecule(self):
        return self.get_molecule()

    def get_molecule(self):
        list_of_symbols = []
        cart_coords = []
        # sdf line pattern containing coordinates and element type
        from chemsmart.utils.repattern import sdf_pattern

        for line in self.contents:
            match = re.match(sdf_pattern, line)
            if match:
                x = float(match.group(1))
                y = float(match.group(2))
                z = float(match.group(3))
                atom_type = str(match.group(4))
                list_of_symbols.append(atom_type)
                cart_coords.append((x, y, z))

        cart_coords = np.array(cart_coords)

        if len(list_of_symbols) == 0 or len(cart_coords) == 0:
            raise ValueError("No coordinates found in the SDF file!")

        return Molecule.from_symbols_and_positions_and_pbc_conditions(
            list_of_symbols=list_of_symbols, positions=cart_coords
        )


class CDXFile(FileMixin):
    """
    ChemDraw file object for reading .cdx and .cdxml files.

    Supports both binary (.cdx) and XML-based (.cdxml) ChemDraw formats.
    Uses RDKit for parsing and generates 3D coordinates using EmbedMolecule.

    Args:
        filename (str or pathlib.Path): Path to the ChemDraw file.
    """

    def __init__(self, filename):
        from pathlib import Path

        # Accept both str and pathlib.Path
        self.filename = (
            str(filename) if isinstance(filename, Path) else filename
        )

    @property
    def molecules(self):
        """
        Return all molecules from the ChemDraw file.
        """
        return self._parse_chemdraw_file()

    def _parse_chemdraw_file(self):
        """
        Parse the ChemDraw file and return a list of Molecule objects.

        Uses RDKit to parse the file and generate 3D coordinates.
        Falls back to Open Babel (via obtain_mols_from_cdx_via_obabel) for .cdx
        if RDKit cannot read it.

        Returns:
            list[Molecule]: List of Molecule objects with 3D coordinates.

        Raises:
            ValueError: If no molecules can be read or no valid 3D structures
            can be generated from the ChemDraw file.
        """
        from pathlib import Path

        from rdkit import Chem
        from rdkit.Chem import AllChem

        suffix = Path(self.filename).suffix.lower()

        rdkit_mols = []
        try:
            # NOTE: RDKit's MolsFromCDXMLFile always supports CDXML.
            # CDX files are only supported if RDKit was built with
            # ChemDraw CDX support.
            rdkit_mols = list(
                Chem.MolsFromCDXMLFile(
                    self.filename, sanitize=False, removeHs=False
                )
            )
        except Exception as e:
            logger.debug(
                f"RDKit MolsFromCDXMLFile failed for {self.filename}: {e}"
            )

        # Fallback for .cdx: use Open Babel helper if RDKit gave nothing.
        if not rdkit_mols and suffix == ".cdx":
            logger.debug(
                f"RDKit did not return any molecules for {self.filename}; "
                "falling back to Open Babel CDX reader."
            )
            try:
                from chemsmart.utils.io import obtain_mols_from_cdx_via_obabel

                rdkit_mols = obtain_mols_from_cdx_via_obabel(self.filename)
            except Exception as e:
                logger.debug(
                    f"Open Babel CDX fallback failed for {self.filename}: {e}"
                )

        if not rdkit_mols:
            raise ValueError(
                f"No molecules could be read from ChemDraw file: {self.filename}"
            )

        molecules = []
        for rdkit_mol in rdkit_mols:
            if rdkit_mol is None:
                continue

            # Handle organometallic complexes: normalize metal bonds and sanitize
            try:
                from chemsmart.utils.io import (
                    normalize_metal_bonds,
                    safe_sanitize,
                )

                rdkit_mol = normalize_metal_bonds(rdkit_mol)
                rdkit_mol = safe_sanitize(rdkit_mol)
            except Exception as e:
                logger.warning(
                    f"Error normalizing metal bonds or sanitizing molecule "
                    f"in {self.filename}: {e}. Continuing with original molecule."
                )

            # Add explicit hydrogens for proper structure
            rdkit_mol = Chem.AddHs(rdkit_mol)

            # Generate 3D coordinates
            try:
                # Try to embed the molecule to get 3D coordinates
                result = AllChem.EmbedMolecule(rdkit_mol, randomSeed=42)
                if result == -1:
                    # Embedding failed, try with random coordinates
                    result = AllChem.EmbedMolecule(
                        rdkit_mol,
                        useRandomCoords=True,
                        randomSeed=42,
                    )
                    if result == -1:
                        logger.warning(
                            f"Could not generate 3D coordinates for a molecule "
                            f"in {self.filename}. Skipping this molecule."
                        )
                        continue

                # Optimize the geometry (may fail for exotic atom types)
                try:
                    AllChem.MMFFOptimizeMolecule(rdkit_mol)
                except Exception as e:
                    logger.debug(
                        f"MMFF optimization failed for a molecule in {self.filename}: {e}"
                    )

            except Exception as e:
                logger.warning(
                    f"Error generating 3D coordinates for molecule in {self.filename}: {e}"
                )
                continue

            # Convert RDKit Mol to Molecule
            mol = Molecule.from_rdkit_mol(rdkit_mol)
            molecules.append(mol)

        if not molecules:
            raise ValueError(
                "No valid molecules with 3D coordinates could be generated "
                f"from ChemDraw file: {self.filename}"
            )

        return molecules

    def get_molecules(self, index="-1", return_list=False):
        """
        Get molecule(s) from the ChemDraw file.

        Args:
            index (str or int): Index specification:
                - "-1": Return the last molecule (default)
                - ":": Return all molecules
                - "1": Return the first molecule (1-based indexing)
                - "1:3": Return molecules 1 to 3 (1-based indexing)
            return_list (bool): If True, always return a list.

        Returns:
            Molecule or list[Molecule]: Single Molecule or list of Molecules.
        """
        molecules = self.molecules

        # Handle index specification
        if isinstance(index, int):
            index = str(index)

        if index == ":":
            # Return all molecules
            return molecules if return_list else molecules

        if index == "-1":
            # Return last molecule
            mol = molecules[-1]
            return [mol] if return_list else mol

        # Handle 1-based indexing
        if index.isdigit() or (index.startswith("-") and index[1:].isdigit()):
            idx = int(index)
            if idx > 0:
                # 1-based positive indexing
                mol = molecules[idx - 1]
            else:
                # Negative indexing
                mol = molecules[idx]
            return [mol] if return_list else mol

        # Handle slice notation (e.g., "1:3")
        if ":" in index:
            parts = index.split(":")
            start = int(parts[0]) - 1 if parts[0] else 0
            end = int(parts[1]) if parts[1] else len(molecules)
            return molecules[start:end]

        raise ValueError(f"Invalid index specification: {index}")
