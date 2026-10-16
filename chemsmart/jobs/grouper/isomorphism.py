"""
RDKit isomorphism-based molecular grouping.

Groups molecules by graph isomorphism using RDKit.
"""

import logging
from collections import defaultdict
from typing import Iterable, List, Tuple

import pandas as pd

from chemsmart.io.molecules.structure import Molecule

from .base import MoleculeGrouper

logger = logging.getLogger(__name__)


class RDKitIsomorphismGrouper(MoleculeGrouper):
    """
    Group molecules by RDKit graph isomorphism.

    Uses RDKit bond-order/aromaticity perception and molecular hashing to group
    molecules with identical chemical graphs. Efficient for large
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
        matrix_format: str = "xlsx",
        energy_type: str = "E",
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
            matrix_format (str): Output format ('xlsx', 'csv', 'txt'). Defaults to 'xlsx'.
        """
        super().__init__(
            molecules,
            num_procs,
            label=label,
            conformer_ids=conformer_ids,
            matrix_format=matrix_format,
            energy_type=energy_type,
            **kwargs,
        )
        self.ignore_hydrogens = ignore_hydrogens

    def group(self) -> Tuple[List[List[Molecule]], List[List[int]]]:
        """
        Group molecules by RDKit graph isomorphism.

        Builds RDKit molecules from coordinates, lets RDKit perceive bond orders
        and aromaticity, computes canonical hashes, and groups molecules with
        identical chemical graphs. Also saves results
        to an Excel file.

        Returns:
            Tuple[List[List[Molecule]], List[List[int]]]: Tuple containing:
                - List of molecule groups (each group is isomorphic)
                - List of index groups (corresponding indices for each group)
        """
        import time

        grouping_start_time = time.time()

        # Convert all molecules to RDKit and compute hashes
        mol_hashes = []
        for mol in self.molecules:
            try:
                mol_hash = mol.get_rdkit_hash(
                    ignore_hydrogens=self.ignore_hydrogens
                )
            except Exception as e:
                logger.warning(f"RDKit molecular hashing failed: {e}")
                mol_hash = None
            mol_hashes.append(mol_hash)

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

        index_groups = list(hash_indices.values())
        groups, index_groups = self._apply_lowest_representative_ordering(
            index_groups
        )
        hashes = list(hash_groups.keys())

        grouping_time = time.time() - grouping_start_time

        # Save results through unified record entrypoint
        self.record(
            hashes=hashes,
            groups=groups,
            index_groups=index_groups,
            grouping_time=grouping_time,
        )

        # Cache results
        self._cached_groups = groups
        self._cached_group_indices = index_groups

        logger.info(f"Found {len(groups)} isomorphism groups")

        return groups, index_groups

    def _record_results(
        self,
        hashes: List[str],
        groups: List[List[Molecule]],
        index_groups: List[List[int]],
        grouping_time: float = None,
    ):
        """Strategy-specific result writer for isomorphism grouping."""
        n = sum(len(g) for g in groups)

        # Build header info
        header_info = [
            ("", f"Isomorphism Grouping Results - {self.__class__.__name__}"),
            ("Total Molecules", n),
            ("Unique Isomorphism Classes", len(groups)),
            ("Energy Type", self.energy_type),
            ("Ignore Hydrogens", self.ignore_hydrogens),
            ("Num Procs", self.num_procs),
        ]

        if grouping_time is not None:
            header_info.append(
                ("Grouping Time", f"{grouping_time:.2f} seconds")
            )

        self._append_input_usage_header(header_info)
        self._append_thermo_header(header_info)

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
            startrow=len(header_info) + 2,
        )

    def __repr__(self):
        return f"{self.__class__.__name__}(num_procs={self.num_procs})"
