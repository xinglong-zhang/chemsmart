"""
RDKit isomorphism-based molecular grouping.

Groups molecules by graph isomorphism using RDKit.
"""

import logging
import os
from collections import defaultdict
from typing import Iterable, List, Optional, Tuple

import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolHash

from chemsmart.io.molecules.structure import Molecule

from .runner import MoleculeGrouper

logger = logging.getLogger(__name__)


class RDKitIsomorphismGrouper(MoleculeGrouper):
    """
    Group molecules by RDKit graph isomorphism.

    Uses RDKit's molecular hashing and isomorphism detection to group
    molecules with identical connectivity patterns. Efficient for large
    datasets due to hash-based pre-filtering.

    Attributes:
        molecules (Iterable[Molecule]): Inherited; collection of molecules to
            group.
        num_procs (int): Inherited; number of worker processes.
    """

    def __init__(
        self,
        molecules: Iterable[Molecule],
        num_procs: int = 1,
        ignore_hydrogens: bool = False,
        label: str = None,
        conformer_ids: List[str] = None,
        **kwargs,
    ):
        """
        Initialize RDKit isomorphism-based molecular grouper.

        Args:
            molecules (Iterable[Molecule]): Collection of molecules to group.
            num_procs (int): Number of processes for parallel computation.
            ignore_hydrogens (bool): Whether to remove hydrogens before
                isomorphism comparison. Defaults to False.
            label (str): Label/name for output files. Defaults to None.
            conformer_ids (list[str]): Custom IDs for each molecule.
        """
        super().__init__(
            molecules, num_procs, label=label, conformer_ids=conformer_ids
        )
        self.ignore_hydrogens = ignore_hydrogens

    def _mol_to_rdkit(self, mol: Molecule):
        """Convert Molecule to RDKit mol object."""
        try:
            # Generate XYZ block string
            lines = [str(mol.num_atoms), ""]
            for symbol, pos in zip(mol.chemical_symbols, mol.positions):
                lines.append(
                    f"{symbol} {pos[0]:.6f} {pos[1]:.6f} {pos[2]:.6f}"
                )
            xyz_string = "\n".join(lines)

            rdkit_mol = Chem.MolFromXYZBlock(xyz_string)
            if rdkit_mol is not None:
                # Determine connectivity
                try:
                    Chem.rdDetermineBonds.DetermineConnectivity(rdkit_mol)
                except AttributeError:
                    # Fallback for older RDKit versions
                    from rdkit.Chem import rdDetermineBonds

                    rdDetermineBonds.DetermineConnectivity(rdkit_mol)

                # Remove hydrogens if requested
                if self.ignore_hydrogens:
                    rdkit_mol = Chem.RemoveHs(rdkit_mol)

            return rdkit_mol
        except Exception as e:
            logger.warning(f"Failed to convert molecule to RDKit: {e}")
            return None

    def _get_mol_hash(self, rdkit_mol) -> Optional[str]:
        """Get canonical hash for RDKit molecule."""
        if rdkit_mol is None:
            return None
        try:
            return rdMolHash.MolHash(
                rdkit_mol, rdMolHash.HashFunction.CanonicalSmiles
            )
        except Exception:
            return None

    def group(self) -> Tuple[List[List[Molecule]], List[List[int]]]:
        """
        Group molecules by RDKit graph isomorphism.

        Converts molecules to RDKit format, computes canonical hashes,
        and groups molecules with identical hashes. Also saves results
        to an Excel file.

        Returns:
            Tuple[List[List[Molecule]], List[List[int]]]: Tuple containing:
                - List of molecule groups (each group is isomorphic)
                - List of index groups (corresponding indices for each group)
        """
        import time

        grouping_start_time = time.time()

        # Convert all molecules to RDKit and compute hashes
        rdkit_mols = []
        mol_hashes = []
        for mol in self.molecules:
            rdkit_mol = self._mol_to_rdkit(mol)
            rdkit_mols.append(rdkit_mol)
            mol_hashes.append(self._get_mol_hash(rdkit_mol))

        # Group by hash
        hash_groups = defaultdict(list)
        hash_indices = defaultdict(list)

        molecules_list = list(self.molecules)
        for i, (mol, mol_hash) in enumerate(zip(molecules_list, mol_hashes)):
            if mol_hash is not None:
                hash_groups[mol_hash].append(mol)
                hash_indices[mol_hash].append(i)
            else:
                # Invalid molecules get their own group
                hash_groups[f"invalid_{i}"].append(mol)
                hash_indices[f"invalid_{i}"].append(i)

        groups = list(hash_groups.values())
        index_groups = list(hash_indices.values())
        hashes = list(hash_groups.keys())

        grouping_time = time.time() - grouping_start_time

        # Save results to Excel
        self._save_isomorphism_results(
            hashes, groups, index_groups, grouping_time
        )

        # Cache results
        self._cached_groups = groups
        self._cached_group_indices = index_groups

        logger.info(f"Found {len(groups)} isomorphism groups")

        return groups, index_groups

    def _save_isomorphism_results(
        self,
        hashes: List[str],
        groups: List[List[Molecule]],
        index_groups: List[List[int]],
        grouping_time: float = None,
    ):
        """Save isomorphism grouping results to Excel file."""
        n = sum(len(g) for g in groups)

        # Create output directory
        if self.label:
            output_dir = f"{self.label}_group_result"
        else:
            output_dir = "group_result"
        os.makedirs(output_dir, exist_ok=True)

        # Create filename
        label_prefix = f"{self.label}_" if self.label else ""
        filename = os.path.join(
            output_dir,
            f"{label_prefix}{self.__class__.__name__}.xlsx",
        )

        with pd.ExcelWriter(filename, engine="openpyxl") as writer:
            # Sheet 1: Summary with header info
            summary_data = []
            for i, (mol_hash, group) in enumerate(zip(hashes, groups)):
                display_hash = (
                    mol_hash
                    if not mol_hash.startswith("invalid_")
                    else "Invalid"
                )
                summary_data.append(
                    {
                        "Group": i + 1,
                        "Hash/SMILES": (
                            display_hash[:50] + "..."
                            if len(str(display_hash)) > 50
                            else display_hash
                        ),
                        "Count": len(group),
                    }
                )

            summary_df = pd.DataFrame(summary_data)
            summary_df.to_excel(
                writer,
                sheet_name="Isomorphism_Groups",
                startrow=6,
                index=False,
            )

            worksheet = writer.sheets["Isomorphism_Groups"]
            row = 1
            worksheet[f"A{row}"] = (
                f"Isomorphism Grouping Results - {self.__class__.__name__}"
            )
            row += 1
            worksheet[f"A{row}"] = f"Total Molecules: {n}"
            row += 1
            worksheet[f"A{row}"] = f"Unique Isomorphism Classes: {len(groups)}"
            row += 1
            if grouping_time is not None:
                worksheet[f"A{row}"] = (
                    f"Grouping Time: {grouping_time:.2f} seconds"
                )
                row += 1

            # Sheet 2: Groups detail
            groups_data = []
            for i, indices in enumerate(index_groups):
                if self.conformer_ids is not None:
                    member_labels = [
                        self.conformer_ids[idx] for idx in indices
                    ]
                else:
                    member_labels = [str(idx + 1) for idx in indices]

                groups_data.append(
                    {
                        "Group": i + 1,
                        "Members": ", ".join(member_labels),
                    }
                )

            groups_df = pd.DataFrame(groups_data)
            groups_df.to_excel(writer, sheet_name="Groups", index=False)

        logger.info(f"Isomorphism grouping results saved to {filename}")

    def __repr__(self):
        return f"{self.__class__.__name__}(num_procs={self.num_procs})"
