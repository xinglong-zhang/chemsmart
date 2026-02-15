"""
RDKit isomorphism-based molecular grouping.

Groups molecules by graph isomorphism using RDKit.
"""

import logging
from collections import defaultdict
from typing import Iterable, List, Optional, Tuple

import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolHash

from chemsmart.io.molecules.structure import Molecule

from .base import MoleculeGrouper

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
        output_format: str = "xlsx",
        **kwargs,
    ):
        """
        Initialize RDKit isomorphism-based molecular grouper.

        Args:
            molecules (Iterable[Molecule]): Collection of molecules to group.
            num_procs (int): Number of processes for parallel computation.
            ignore_hydrogens (bool): Whether to remove hydrogens before
                isomorphism comparison. Defaults to False. Warning: For some
                molecules, removing hydrogens may cause kekulization
                issues. If errors occur, try setting this to False.
            label (str): Label/name for output files. Defaults to None.
            conformer_ids (list[str]): Custom IDs for each molecule.
            output_format (str): Output format ('xlsx', 'csv', 'txt'). Defaults to 'xlsx'.
        """
        super().__init__(
            molecules,
            num_procs,
            label=label,
            conformer_ids=conformer_ids,
            output_format=output_format,
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
                    try:
                        rdkit_mol = Chem.RemoveHs(rdkit_mol, sanitize=False)
                        # Re-sanitize without kekulization to avoid aromatic ring issues
                        Chem.SanitizeMol(
                            rdkit_mol,
                            sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL
                            ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE,
                        )
                    except Exception as e:
                        logger.warning(
                            f"Failed to remove hydrogens: {e}. Using original molecule."
                        )
                        # Re-convert without removing hydrogens
                        rdkit_mol = Chem.MolFromXYZBlock(xyz_string)
                        if rdkit_mol is not None:
                            Chem.rdDetermineBonds.DetermineConnectivity(
                                rdkit_mol
                            )

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
        """Save isomorphism grouping results to file using ResultsRecorder."""
        n = sum(len(g) for g in groups)

        # Build header info
        header_info = [
            ("", f"Isomorphism Grouping Results - {self.__class__.__name__}"),
            ("Total Molecules", n),
            ("Unique Isomorphism Classes", len(groups)),
            ("Ignore Hydrogens", self.ignore_hydrogens),
            ("Num Procs", self.num_procs),
        ]

        if grouping_time is not None:
            header_info.append(
                ("Grouping Time", f"{grouping_time:.2f} seconds")
            )

        # Use ResultsRecorder to save
        recorder = self._get_results_recorder()

        # Build summary data
        summary_data = []
        for i, (mol_hash, group) in enumerate(zip(hashes, groups)):
            display_hash = (
                mol_hash if not mol_hash.startswith("invalid_") else "Invalid"
            )
            summary_data.append(
                {
                    "Group": i + 1,
                    "Hash/SMILES": display_hash,
                    "Count": len(group),
                }
            )
        summary_df = pd.DataFrame(summary_data)

        # Build sheets data using recorder's method
        sheets_data = {
            "Isomorphism_Groups": summary_df,
            "Groups": recorder.build_groups_dataframe(index_groups, n),
        }

        recorder.record_results(
            grouper_name=self.__class__.__name__,
            header_info=header_info,
            sheets_data=sheets_data,
            matrix_data=None,
            suffix=None,
            startrow=7,
        )

    def __repr__(self):
        return f"{self.__class__.__name__}(num_procs={self.num_procs})"
