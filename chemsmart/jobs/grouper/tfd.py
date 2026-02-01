"""
Torsion Fingerprint Deviation (TFD) based molecular grouping.

Groups molecular conformers based on torsion angle similarity using TFD.
"""

import logging
import os
from typing import Iterable, List, Tuple

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import TorsionFingerprints

from chemsmart.io.molecules.structure import Molecule

from .runner import MoleculeGrouper

logger = logging.getLogger(__name__)


class TorsionFingerprintGrouper(MoleculeGrouper):
    """
    Groups conformers based on Torsion Fingerprint Deviation (TFD).

    TFD is a measure of the similarity of the torsion angles of rotatable
    bonds between conformers. Lower values indicate higher similarity.

    Reference: Schulz-Gasch et al., JCIM, 1499-1512 (2012)

    Attributes:
        threshold (float): TFD threshold for grouping (lower values = more similar).
        num_groups (int): Number of groups to create (alternative to threshold).
        num_procs (int): Number of processes for parallel computation.
        use_weights (bool): Whether to use torsion weights in TFD calculation.
        max_dev (str): Normalization method ('equal' or 'spec').
        symm_radius (int): Radius for calculating atom invariants.
        ignore_colinear_bonds (bool): Whether to ignore single bonds adjacent to triple bonds.
    """

    def __init__(
        self,
        molecules: Iterable[Molecule],
        threshold: float = None,
        num_groups: int = None,
        num_procs: int = 1,
        use_weights: bool = True,
        max_dev: str = "equal",
        symm_radius: int = 2,
        ignore_colinear_bonds: bool = True,
        ignore_hydrogens: bool = False,
        label: str = None,
        conformer_ids: List[str] = None,
        **kwargs,
    ):
        """
        Initialize TFD-based conformer grouper.

        Args:
            molecules (Iterable[Molecule]): Collection of molecule conformers to group.
            threshold (float): TFD threshold for grouping. Lower values indicate more
                similar torsion patterns. Defaults to 0.1. Ignored if num_groups is specified.
            num_groups (int): Number of groups to create. When specified,
                automatically determines threshold to create this many groups.
            num_procs (int): Number of processes for parallel computation.
            use_weights (bool): Whether to use torsion weights in TFD calculation. Defaults to True.
            max_dev (str): Normalization method:
                - 'equal': all torsions normalized using 180.0 (default)
                - 'spec': each torsion normalized using specific maximal deviation
            symm_radius (int): Radius for calculating atom invariants. Defaults to 2.
            ignore_colinear_bonds (bool): If True, single bonds adjacent to triple bonds
                are ignored. Defaults to True.
            ignore_hydrogens (bool): Whether to exclude hydrogen atoms from TFD calculation.
                Defaults to False.
            label (str): Label/name for output files. Defaults to None.
            conformer_ids (list[str]): Custom IDs for each molecule (e.g., ['c1', 'c2']).
        """
        super().__init__(
            molecules, num_procs, label=label, conformer_ids=conformer_ids
        )

        # Validate that threshold and num_groups are mutually exclusive
        if threshold is not None and num_groups is not None:
            raise ValueError(
                "Cannot specify both threshold (-t) and num_groups (-N). Please use only one."
            )

        if threshold is None and num_groups is None:
            threshold = 0.1
        self.threshold = threshold
        self.num_groups = num_groups
        self.use_weights = use_weights
        self.max_dev = max_dev
        self.symm_radius = symm_radius
        self.ignore_colinear_bonds = ignore_colinear_bonds
        self.ignore_hydrogens = ignore_hydrogens

        # Convert to list for indexing
        self.molecules = list(molecules)
        self._prepare_conformer_molecule()

    def _prepare_conformer_molecule(self):
        """
        Prepare a single RDKit molecule with multiple conformers from input molecules.

        This assumes all input molecules are conformers of the same chemical structure.
        """
        if not self.molecules:
            self.rdkit_mol = None
            self.valid_conformer_ids = []
            return

        # Use the first molecule as the base structure
        base_mol = self.molecules[0]
        self.rdkit_mol = base_mol.to_rdkit()

        if self.rdkit_mol is None:
            self.valid_conformer_ids = []
            return

        # Remove hydrogens if ignore_hydrogens is True
        if self.ignore_hydrogens:
            self.rdkit_mol = Chem.RemoveHs(self.rdkit_mol)
            logger.info("Removed hydrogen atoms for TFD calculation")

        # Clear existing conformers and add all input molecules as conformers
        self.rdkit_mol.RemoveAllConformers()
        self.valid_conformer_ids = []

        # Get heavy atom indices for filtering positions if ignore_hydrogens
        if self.ignore_hydrogens:
            heavy_atom_indices = [
                i
                for i, sym in enumerate(base_mol.chemical_symbols)
                if sym != "H"
            ]
        else:
            heavy_atom_indices = None

        for i, mol in enumerate(self.molecules):
            try:
                # Get positions, filtering out hydrogens if needed
                if self.ignore_hydrogens and heavy_atom_indices is not None:
                    positions = mol.positions[heavy_atom_indices]
                    num_atoms = len(heavy_atom_indices)
                else:
                    positions = mol.positions
                    num_atoms = mol.num_atoms

                # Convert molecule positions to RDKit conformer
                conf = Chem.Conformer(num_atoms)
                for atom_idx, pos in enumerate(positions):
                    conf.SetAtomPosition(atom_idx, pos)

                # Add conformer to the molecule
                conf_id = self.rdkit_mol.AddConformer(conf, assignId=True)
                self.valid_conformer_ids.append(conf_id)

            except Exception as e:
                logger.warning(f"Failed to add conformer {i}: {str(e)}")
                continue

        logger.info(
            f"Prepared molecule with {len(self.valid_conformer_ids)} valid conformers"
            + (" (hydrogens ignored)" if self.ignore_hydrogens else "")
        )

    def _calculate_tfd(self, conf_pair: Tuple[int, int]) -> float:
        """
        Calculate TFD between two conformers.

        Args:
            conf_pair (Tuple[int, int]): Indices of conformers to compare.

        Returns:
            float: TFD value between the conformers.
        """
        if self.rdkit_mol is None or len(self.valid_conformer_ids) == 0:
            return float("inf")

        i, j = conf_pair

        try:
            conf_id1 = self.valid_conformer_ids[i]
            conf_id2 = self.valid_conformer_ids[j]

            tfd_values = TorsionFingerprints.GetTFDBetweenConformers(
                self.rdkit_mol,
                confIds1=[conf_id1],
                confIds2=[conf_id2],
                useWeights=self.use_weights,
                maxDev=self.max_dev,
                symmRadius=self.symm_radius,
                ignoreColinearBonds=self.ignore_colinear_bonds,
            )

            return tfd_values[0] if tfd_values else float("inf")

        except Exception as e:
            logger.warning(
                f"TFD calculation failed for conformers {i}, {j}: {str(e)}"
            )
            return float("inf")

    def group(self) -> Tuple[List[List[Molecule]], List[List[int]]]:
        """
        Group conformers based on TFD similarity.

        Returns:
            Tuple[List[List[Molecule]], List[List[int]]]: Tuple containing:
                - List of molecule groups (each group is a list of molecules)
                - List of index groups (corresponding indices for each group)
        """
        import time

        grouping_start_time = time.time()

        n = len(self.molecules)

        if n == 0:
            return [], []

        if n == 1:
            return [self.molecules], [[0]]

        if self.rdkit_mol is None or len(self.valid_conformer_ids) == 0:
            logger.warning(
                "No valid conformers found, each molecule becomes its own group"
            )
            return [[mol] for mol in self.molecules], [[i] for i in range(n)]

        # Generate conformer pairs for TFD calculation
        indices = [(i, j) for i in range(n) for j in range(i + 1, n)]
        total_pairs = len(indices)

        print(
            f"[{self.__class__.__name__}] Starting TFD calculation for {n} conformers ({total_pairs} pairs)"
        )
        print(f"  - TFD threshold: {self.threshold}")
        print(f"  - Use weights: {self.use_weights}")
        print(f"  - Max deviation: {self.max_dev}")
        print(f"  - Symmetry radius: {self.symm_radius}")
        print(f"  - Ignore colinear bonds: {self.ignore_colinear_bonds}")

        # Calculate TFD values with real-time output
        tfd_values = []
        for idx, (i, j) in enumerate(indices):
            tfd = self._calculate_tfd((i, j))
            tfd_values.append(tfd)
            print(
                f"The {idx+1}/{total_pairs} pair (conformer{i+1}, conformer{j+1}) calculation finished, TFD= {tfd:.7f}"
            )

        # Build full TFD matrix
        tfd_matrix = np.zeros((n, n))
        for (i, j), tfd in zip(indices, tfd_values):
            tfd_matrix[i, j] = tfd_matrix[j, i] = tfd

        # Choose grouping strategy
        if self.num_groups is not None:
            groups, index_groups = self._group_by_num_groups(
                tfd_values, indices, n
            )
        else:
            groups, index_groups = self._group_by_threshold(
                tfd_values, indices, n
            )

        grouping_time = time.time() - grouping_start_time

        # Save TFD matrix
        if self.label:
            output_dir = f"{self.label}_group_result"
        else:
            output_dir = "group_result"
        os.makedirs(output_dir, exist_ok=True)

        label_prefix = f"{self.label}_" if self.label else ""
        if self.num_groups is not None:
            matrix_filename = os.path.join(
                output_dir,
                f"{label_prefix}{self.__class__.__name__}_N{self.num_groups}.xlsx",
            )
        else:
            matrix_filename = os.path.join(
                output_dir,
                f"{label_prefix}{self.__class__.__name__}_T{self.threshold}.xlsx",
            )

        self._save_tfd_matrix(
            tfd_matrix, matrix_filename, grouping_time, groups, index_groups
        )

        # Cache results
        self._cached_groups = groups
        self._cached_group_indices = index_groups

        print(
            f"[{self.__class__.__name__}] Found {len(groups)} groups using TFD"
        )

        return groups, index_groups

    def _group_by_threshold(self, tfd_values, indices, n):
        """Threshold-based grouping for TFD using complete linkage."""
        adj_matrix = np.zeros((n, n), dtype=bool)
        for (i, j), tfd in zip(indices, tfd_values):
            if tfd <= self.threshold:
                adj_matrix[i, j] = adj_matrix[j, i] = True

        groups, index_groups = self._complete_linkage_grouping(adj_matrix, n)
        return groups, index_groups

    def _complete_linkage_grouping(self, adj_matrix, n):
        """Perform complete linkage grouping."""
        assigned = [False] * n
        groups = []
        index_groups = []

        for i in range(n - 1, -1, -1):
            if assigned[i]:
                continue
            current_group = [i]
            assigned[i] = True

            for j in range(i - 1, -1, -1):
                if assigned[j]:
                    continue
                can_join = all(
                    adj_matrix[j, member] for member in current_group
                )
                if can_join:
                    current_group.append(j)
                    assigned[j] = True

            current_group.sort()
            groups.append([self.molecules[idx] for idx in current_group])
            index_groups.append(current_group)

        groups.reverse()
        index_groups.reverse()
        return groups, index_groups

    def _group_by_num_groups(self, tfd_values, indices, n):
        """Automatic grouping to create specified number of groups."""
        if self.num_groups >= n:
            print(
                f"[{self.__class__.__name__}] Requested {self.num_groups} groups but only {n} molecules. Creating {n} groups."
            )
            groups = [[mol] for mol in self.molecules]
            index_groups = [[i] for i in range(n)]
            return groups, index_groups

        threshold = self._find_optimal_tfd_threshold(tfd_values, indices, n)
        self._auto_threshold = threshold

        print(
            f"[{self.__class__.__name__}] Auto-determined threshold: {threshold:.7f} to create {self.num_groups} groups"
        )

        adj_matrix = np.zeros((n, n), dtype=bool)
        for (i, j), tfd in zip(indices, tfd_values):
            if tfd <= threshold:
                adj_matrix[i, j] = adj_matrix[j, i] = True

        groups, index_groups = self._complete_linkage_grouping(adj_matrix, n)
        actual_groups = len(groups)

        print(
            f"[{self.__class__.__name__}] Created {actual_groups} groups (requested: {self.num_groups})"
        )

        if actual_groups > self.num_groups:
            groups, index_groups = self._merge_groups_to_target(
                groups, index_groups
            )

        return groups, index_groups

    def _find_optimal_tfd_threshold(self, tfd_values, indices, n):
        """Find threshold using binary search."""
        sorted_tfd = sorted([tfd for tfd in tfd_values if not np.isinf(tfd)])

        if not sorted_tfd:
            return 0.0

        low, high = 0, len(sorted_tfd) - 1
        best_threshold = sorted_tfd[-1]

        while low <= high:
            mid = (low + high) // 2
            threshold = sorted_tfd[mid]

            adj_matrix = np.zeros((n, n), dtype=bool)
            for (idx_i, idx_j), tfd in zip(indices, tfd_values):
                if tfd <= threshold:
                    adj_matrix[idx_i, idx_j] = adj_matrix[idx_j, idx_i] = True

            num_groups_found = self._count_groups(adj_matrix, n)

            if num_groups_found == self.num_groups:
                return threshold
            elif num_groups_found > self.num_groups:
                low = mid + 1
            else:
                best_threshold = threshold
                high = mid - 1

        return best_threshold

    def _count_groups(self, adj_matrix, n):
        """Count number of groups using complete linkage."""
        assigned = [False] * n
        num_groups = 0
        for i in range(n):
            if assigned[i]:
                continue
            current_group = [i]
            assigned[i] = True
            for j in range(i + 1, n):
                if assigned[j]:
                    continue
                can_join = all(adj_matrix[j, m] for m in current_group)
                if can_join:
                    current_group.append(j)
                    assigned[j] = True
            num_groups += 1
        return num_groups

    def _merge_groups_to_target(self, groups, index_groups):
        """Merge groups to reach target number."""
        while len(groups) > self.num_groups:
            min_idx = min(range(len(groups)), key=lambda i: len(groups[i]))
            largest_idx = max(
                range(len(groups)),
                key=lambda i: len(groups[i]) if i != min_idx else -1,
            )
            groups[largest_idx].extend(groups[min_idx])
            index_groups[largest_idx].extend(index_groups[min_idx])
            groups.pop(min_idx)
            index_groups.pop(min_idx)
        return groups, index_groups

    def _save_tfd_matrix(
        self,
        tfd_matrix: np.ndarray,
        filename: str,
        grouping_time: float = None,
        groups: List[List[Molecule]] = None,
        index_groups: List[List[int]] = None,
    ):
        """Save TFD matrix to Excel file."""
        n = tfd_matrix.shape[0]

        if filename.endswith(".txt"):
            filename = filename[:-4] + ".xlsx"
        elif not filename.endswith(".xlsx"):
            filename = filename + ".xlsx"

        if self.conformer_ids is not None and len(self.conformer_ids) == n:
            row_labels = self.conformer_ids
            col_labels = self.conformer_ids
        else:
            row_labels = [str(i + 1) for i in range(n)]
            col_labels = [str(j + 1) for j in range(n)]

        matrix_display = np.where(np.isinf(tfd_matrix), np.nan, tfd_matrix)
        df = pd.DataFrame(matrix_display, index=row_labels, columns=col_labels)

        with pd.ExcelWriter(filename, engine="openpyxl") as writer:
            df.to_excel(
                writer,
                sheet_name="TFD_Matrix",
                startrow=12,
                float_format="%.7f",
            )

            worksheet = writer.sheets["TFD_Matrix"]

            row = 1
            worksheet[f"A{row}"] = (
                f"Full TFD Matrix ({n}x{n}) - {self.__class__.__name__}"
            )
            row += 1
            worksheet[f"A{row}"] = (
                "Based on Schulz-Gasch et al., JCIM, 1499-1512 (2012)"
            )
            row += 1

            if hasattr(self, "num_groups") and self.num_groups is not None:
                worksheet[f"A{row}"] = (
                    f"Requested Groups (-N): {self.num_groups}"
                )
                row += 1
                if (
                    hasattr(self, "_auto_threshold")
                    and self._auto_threshold is not None
                ):
                    worksheet[f"A{row}"] = (
                        f"Auto-determined Threshold: {self._auto_threshold:.7f}"
                    )
                    row += 1
            else:
                worksheet[f"A{row}"] = f"Threshold: {self.threshold}"
                row += 1

            worksheet[f"A{row}"] = f"Use Weights: {self.use_weights}"
            row += 1
            worksheet[f"A{row}"] = f"Max Deviation: {self.max_dev}"
            row += 1
            worksheet[f"A{row}"] = f"Symmetry Radius: {self.symm_radius}"
            row += 1
            worksheet[f"A{row}"] = (
                f"Ignore Colinear Bonds: {self.ignore_colinear_bonds}"
            )
            row += 1

            if grouping_time is not None:
                worksheet[f"A{row}"] = (
                    f"Grouping Time: {grouping_time:.2f} seconds"
                )
                row += 1

            worksheet[f"A{row}"] = (
                "Lower values indicate higher torsional similarity. Empty cells indicate calculation failures."
            )
            row += 1

            # Auto-adjust column widths
            for column in worksheet.columns:
                max_length = 0
                column_letter = column[0].column_letter
                for cell in column:
                    try:
                        if cell.value:
                            max_length = max(max_length, len(str(cell.value)))
                    except (TypeError, AttributeError):
                        pass
                adjusted_width = min(max_length + 2, 18)
                worksheet.column_dimensions[column_letter].width = (
                    adjusted_width
                )

            # Add Groups sheet if groups are provided
            if groups is not None and index_groups is not None:
                self._write_groups_sheet(writer, groups, index_groups)

        logger.info(f"TFD matrix saved to {filename}")

    def _write_groups_sheet(self, writer, groups, index_groups):
        """Write groups information to a separate sheet."""
        groups_data = []
        for i, (group, indices) in enumerate(zip(groups, index_groups)):
            if self.conformer_ids is not None:
                member_labels = [self.conformer_ids[idx] for idx in indices]
            else:
                member_labels = [str(idx + 1) for idx in indices]

            groups_data.append(
                {
                    "Group": i + 1,
                    "Members": len(group),
                    "Indices": str(member_labels),
                    "Representative": (
                        member_labels[0] if member_labels else "N/A"
                    ),
                }
            )

        groups_df = pd.DataFrame(groups_data)
        groups_df.to_excel(writer, sheet_name="Groups", index=False)

    def __repr__(self):
        if self.num_groups is not None:
            return (
                f"{self.__class__.__name__}(num_groups={self.num_groups}, "
                f"num_procs={self.num_procs}, use_weights={self.use_weights}, "
                f"max_dev='{self.max_dev}', symm_radius={self.symm_radius})"
            )
        else:
            return (
                f"{self.__class__.__name__}(threshold={self.threshold}, "
                f"num_procs={self.num_procs}, use_weights={self.use_weights}, "
                f"max_dev='{self.max_dev}', symm_radius={self.symm_radius})"
            )
