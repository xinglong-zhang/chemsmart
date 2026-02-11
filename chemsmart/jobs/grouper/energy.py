"""
Energy-based molecular grouping algorithm.

This module contains the EnergyGrouper implementation which groups molecules
based on energy differences using a threshold or target number of groups.
"""

import logging
import os
import time
from typing import Iterable, List, Tuple

import numpy as np
import pandas as pd

from chemsmart.io.molecules.structure import Molecule

from .runner import MoleculeGrouper

logger = logging.getLogger(__name__)

# Conversion factor: 1 Hartree = 627.509474 kcal/mol
HARTREE_TO_KCAL = 627.509474
KCAL_TO_HARTREE = 1.0 / HARTREE_TO_KCAL


class EnergyGrouper(MoleculeGrouper):
    """
    Energy-based molecular grouping.

    Groups molecules based on energy differences. Molecules with energy
    differences below the threshold are grouped together using complete
    linkage clustering.

    Attributes:
        molecules (Iterable[Molecule]): Collection of molecules to group.
        num_procs (int): Number of worker processes/threads.
        threshold (float): Energy difference threshold for grouping (in kcal/mol).
        threshold_hartree (float): Threshold converted to Hartree for internal use.
        num_groups (int): Alternative to threshold - target number of groups.
        label (str): Label for output files.
        conformer_ids (list[str]): Custom conformer IDs for labeling.
    """

    def __init__(
        self,
        molecules: Iterable[Molecule],
        threshold: float = None,
        num_groups: int = None,
        num_procs: int = 1,
        label: str = None,
        conformer_ids: List[str] = None,
        **kwargs,
    ):
        """
        Initialize Energy-based molecular grouper.

        Args:
            molecules (Iterable[Molecule]): Collection of molecules to group.
            threshold (float): Energy difference threshold for grouping in kcal/mol.
                Defaults to 1.0 kcal/mol.
                Ignored if num_groups is specified.
            num_groups (int): Number of groups to create. When specified,
                automatically determines threshold to create this many groups.
            num_procs (int): Number of processes (not used for energy calculation
                but kept for API consistency). Defaults to 1.
            label (str): Label/name for output files. Defaults to None.
            conformer_ids (list[str]): Custom IDs for each molecule (e.g., ['c1', 'c2']).

        Note:
            Uses complete linkage clustering: a structure joins a group only if
            its energy difference to ALL existing members is below the threshold.
        """
        super().__init__(
            molecules, num_procs, label=label, conformer_ids=conformer_ids
        )

        # Validate that threshold and num_groups are mutually exclusive
        if threshold is not None and num_groups is not None:
            raise ValueError(
                "Cannot specify both threshold (-T) and num_groups (-N). Please use only one."
            )

        if threshold is None and num_groups is None:
            threshold = 1.0  # Default: 1 kcal/mol

        # Store threshold in kcal/mol (user-friendly unit)
        self.threshold = threshold
        # Convert to Hartree for internal calculations (only if threshold is set)
        self.threshold_hartree = (
            threshold * KCAL_TO_HARTREE if threshold is not None else None
        )
        self.num_groups = num_groups
        self._auto_threshold = (
            None  # Will store auto-determined threshold in kcal/mol
        )

        # Validate that all molecules have energy information
        self._validate_energies()

    def _validate_energies(self) -> None:
        """
        Validate that all molecules have energy information.

        Raises:
            ValueError: If any molecule is missing energy information.
        """
        missing_energy = []
        for i, mol in enumerate(self.molecules):
            if mol.energy is None:
                missing_energy.append(i + 1)

        if missing_energy:
            if len(missing_energy) <= 5:
                raise ValueError(
                    f"Molecules at indices {missing_energy} are missing energy information. "
                    "Energy grouping requires all molecules to have energy values."
                )
            else:
                raise ValueError(
                    f"Found {len(missing_energy)} molecules without energy information. "
                    "Energy grouping requires all molecules to have energy values."
                )

    def _calculate_energy_diff(
        self, idx_pair: Tuple[int, int]
    ) -> Tuple[float, float]:
        """
        Calculate energy difference between two molecules.

        Args:
            idx_pair (Tuple[int, int]): Tuple of molecule indices (i, j).

        Returns:
            Tuple[float, float]: (relative energy difference, absolute energy difference) in Hartree.
                - Relative: E_j - E_i (positive means j has higher energy than i)
                - Absolute: |E_j - E_i| (for threshold comparison)
        """
        i, j = idx_pair
        energy_i = self.molecules[i].energy
        energy_j = self.molecules[j].energy
        relative_diff = energy_j - energy_i  # Positive if j has higher energy
        return relative_diff, abs(relative_diff)

    def group(self) -> Tuple[List[List[Molecule]], List[List[int]]]:
        """
        Group molecules by energy similarity using complete linkage clustering.

        Computes pairwise energy differences between all molecules and groups
        those within the specified threshold using complete linkage clustering,
        or automatically determines threshold to create the specified number
        of groups.

        Returns:
            Tuple[List[List[Molecule]], List[List[int]]]: Tuple containing:
                - List of molecule groups (each group is a list of molecules)
                - List of index groups (corresponding indices for each group)
        """
        grouping_start_time = time.time()

        n = len(self.molecules)
        indices = [(i, j) for i in range(n) for j in range(i + 1, n)]
        total_pairs = len(indices)

        logger.info(
            f"[{self.__class__.__name__}] Starting calculation for {n} molecules ({total_pairs} pairs)"
        )

        # Calculate energy differences (both relative and absolute)
        energy_diff_relative = []  # For output matrix (with sign)
        energy_diff_absolute = []  # For threshold comparison
        for idx, (i, j) in enumerate(indices):
            rel_diff, abs_diff = self._calculate_energy_diff((i, j))
            energy_diff_relative.append(rel_diff)
            energy_diff_absolute.append(abs_diff)
            logger.info(
                f"The {idx+1}/{total_pairs} pair (conformer{i+1}, conformer{j+1}) calculation finished, "
                f"Energy Diff= {rel_diff:+.10f} Hartree ({rel_diff * HARTREE_TO_KCAL:+.4f} kcal/mol)"
            )

        # Build full energy difference matrix for output (with sign, relative to smaller index)
        # matrix[i,j] = E_j - E_i (positive means j has higher energy)
        # matrix[j,i] = E_i - E_j (negative of the above)
        energy_matrix = np.zeros((n, n))
        for (i, j), rel_diff in zip(indices, energy_diff_relative):
            energy_matrix[i, j] = rel_diff  # E_j - E_i
            energy_matrix[j, i] = -rel_diff  # E_i - E_j

        # Create output directory
        if self.label:
            output_dir = f"{self.label}_group_result"
        else:
            output_dir = "group_result"
        os.makedirs(output_dir, exist_ok=True)

        # Create filename based on label, grouper type and threshold/num_groups
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

        # Choose grouping strategy based on parameters
        if self.num_groups is not None:
            groups, index_groups = self._group_by_num_groups(
                energy_matrix, energy_diff_absolute, indices
            )
        else:
            groups, index_groups = self._group_by_threshold(
                energy_diff_absolute, indices
            )

        # Calculate total grouping time
        grouping_end_time = time.time()
        grouping_time = grouping_end_time - grouping_start_time

        # Cache the results BEFORE saving (so Groups sheet can be populated)
        self._cached_groups = groups
        self._cached_group_indices = index_groups

        # Save energy difference matrix
        self._save_energy_matrix(energy_matrix, matrix_filename, grouping_time)

        return groups, index_groups

    def _group_by_threshold(
        self, energy_diff_values: List[float], indices: List[Tuple[int, int]]
    ) -> Tuple[List[List[Molecule]], List[List[int]]]:
        """
        Threshold-based grouping using complete linkage.

        A structure joins a group only if its energy difference to ALL
        existing members is below the threshold.
        """
        n = len(self.molecules)

        # Build adjacency matrix (use threshold_hartree for comparison since energy_diff is in Hartree)
        adj_matrix = np.zeros((n, n), dtype=bool)
        for (i, j), diff in zip(indices, energy_diff_values):
            if diff < self.threshold_hartree:
                adj_matrix[i, j] = adj_matrix[j, i] = True

        # Complete linkage grouping
        groups, index_groups = self._complete_linkage_grouping(adj_matrix, n)

        return groups, index_groups

    def _complete_linkage_grouping(
        self, adj_matrix: np.ndarray, n: int
    ) -> Tuple[List[List[Molecule]], List[List[int]]]:
        """
        Perform complete linkage grouping (from last to first).

        A structure joins a group only if its energy difference to ALL
        existing members of that group is below the threshold.

        Args:
            adj_matrix: Boolean adjacency matrix where adj_matrix[i,j] = True
                       if energy diff between i and j is below threshold
            n: Number of molecules

        Returns:
            Tuple of (groups, index_groups)
        """
        assigned = [False] * n
        groups = []
        index_groups = []

        # Iterate from last to first (higher energy structures first as seeds)
        for i in range(n - 1, -1, -1):
            if assigned[i]:
                continue

            # Start a new group with molecule i
            current_group = [i]
            assigned[i] = True

            # Try to add unassigned molecules with lower indices (lower energy)
            for j in range(i - 1, -1, -1):
                if assigned[j]:
                    continue

                # Check if j is connected to ALL members in current_group
                can_join = all(
                    adj_matrix[j, member] for member in current_group
                )

                if can_join:
                    current_group.append(j)
                    assigned[j] = True

            # Sort indices within group (lowest index = lowest energy first)
            current_group.sort()
            groups.append([self.molecules[idx] for idx in current_group])
            index_groups.append(current_group)

        # Reverse groups so that groups containing lower-energy structures come first
        groups.reverse()
        index_groups.reverse()

        return groups, index_groups

    def _group_by_num_groups(
        self,
        energy_matrix: np.ndarray,
        energy_diff_values: List[float],
        indices: List[Tuple[int, int]],
    ) -> Tuple[List[List[Molecule]], List[List[int]]]:
        """Automatic grouping to create specified number of groups."""
        n = len(self.molecules)

        if self.num_groups >= n:
            logger.info(
                f"[{self.__class__.__name__}] Requested {self.num_groups} groups but only {n} molecules. Creating {n} groups."
            )
            groups = [[mol] for mol in self.molecules]
            index_groups = [[i] for i in range(n)]
            return groups, index_groups

        # Find appropriate threshold to create desired number of groups (returns in Hartree)
        threshold_hartree = self._find_optimal_threshold(
            energy_diff_values, indices, n
        )
        # Store in kcal/mol for reporting
        self._auto_threshold = threshold_hartree * HARTREE_TO_KCAL

        logger.info(
            f"[{self.__class__.__name__}] Auto-determined threshold: {self._auto_threshold:.4f} kcal/mol "
            f"({threshold_hartree:.10f} Hartree) to create {self.num_groups} groups"
        )

        # Build adjacency matrix with the determined threshold
        adj_matrix = np.zeros((n, n), dtype=bool)
        for (i, j), diff in zip(indices, energy_diff_values):
            if diff < threshold_hartree:
                adj_matrix[i, j] = adj_matrix[j, i] = True

        # Use complete linkage grouping
        groups, index_groups = self._complete_linkage_grouping(adj_matrix, n)
        actual_groups = len(groups)

        logger.info(
            f"[{self.__class__.__name__}] Created {actual_groups} groups (requested: {self.num_groups})"
        )

        if actual_groups > self.num_groups:
            groups, index_groups = self._merge_groups_to_target(
                groups, index_groups, adj_matrix
            )

        return groups, index_groups

    def _find_optimal_threshold(
        self,
        energy_diff_values: List[float],
        indices: List[Tuple[int, int]],
        n: int,
    ) -> float:
        """Find threshold that creates approximately the desired number of groups using binary search."""
        sorted_diffs = sorted(
            [diff for diff in energy_diff_values if not np.isinf(diff)]
        )

        if not sorted_diffs:
            return 0.0

        # Binary search for optimal threshold
        low, high = 0, len(sorted_diffs) - 1
        best_threshold = sorted_diffs[-1]

        while low <= high:
            mid = (low + high) // 2
            threshold = sorted_diffs[mid]

            # Build adjacency matrix with this threshold
            adj_matrix = np.zeros((n, n), dtype=bool)
            for (idx_i, idx_j), diff in zip(indices, energy_diff_values):
                if diff < threshold:
                    adj_matrix[idx_i, idx_j] = adj_matrix[idx_j, idx_i] = True

            # Use complete linkage to count groups
            groups, _ = self._complete_linkage_grouping(adj_matrix, n)
            num_groups_found = len(groups)

            if num_groups_found == self.num_groups:
                return threshold
            elif num_groups_found > self.num_groups:
                # Too many groups, need higher threshold (more permissive)
                low = mid + 1
            else:
                # Too few groups, need lower threshold (more restrictive)
                best_threshold = threshold
                high = mid - 1

        return best_threshold

    def _merge_groups_to_target(
        self,
        groups: List[List[Molecule]],
        index_groups: List[List[int]],
        adj_matrix: np.ndarray,
    ) -> Tuple[List[List[Molecule]], List[List[int]]]:
        """
        Merge groups to reach target number when using complete linkage.

        Args:
            groups: List of molecule groups
            index_groups: List of index groups
            adj_matrix: Boolean adjacency matrix

        Returns:
            Tuple of (merged_groups, merged_index_groups)
        """
        while len(groups) > self.num_groups:
            # Find the smallest group
            min_size = float("inf")
            min_idx = 0
            for i, g in enumerate(groups):
                if len(g) < min_size:
                    min_size = len(g)
                    min_idx = i

            # Find the best group to merge with (most connections)
            best_merge_idx = -1
            best_connection_count = -1

            for i, target_indices in enumerate(index_groups):
                if i == min_idx:
                    continue

                # Count connections between source group and target group
                connection_count = 0
                for src_idx in index_groups[min_idx]:
                    for tgt_idx in target_indices:
                        if adj_matrix[src_idx, tgt_idx]:
                            connection_count += 1

                if connection_count > best_connection_count:
                    best_connection_count = connection_count
                    best_merge_idx = i

            # Merge the smallest group into the best target
            if best_merge_idx >= 0:
                groups[best_merge_idx].extend(groups[min_idx])
                index_groups[best_merge_idx].extend(index_groups[min_idx])
            else:
                # No connections found, merge into the largest group
                largest_idx = max(
                    range(len(groups)),
                    key=lambda i: len(groups[i]) if i != min_idx else -1,
                )
                groups[largest_idx].extend(groups[min_idx])
                index_groups[largest_idx].extend(index_groups[min_idx])

            # Remove the merged group
            groups.pop(min_idx)
            index_groups.pop(min_idx)

        return groups, index_groups

    def _save_energy_matrix(
        self,
        energy_matrix: np.ndarray,
        filename: str,
        grouping_time: float = None,
    ):
        """Save energy difference matrix to Excel file (in kcal/mol units)."""
        n = energy_matrix.shape[0]

        # Ensure .xlsx extension
        if filename.endswith(".txt"):
            filename = filename[:-4] + ".xlsx"
        elif not filename.endswith(".xlsx"):
            filename = filename + ".xlsx"

        # Create DataFrame with labels - use conformer_ids if available
        if self.conformer_ids is not None and len(self.conformer_ids) == n:
            row_labels = self.conformer_ids
            col_labels = self.conformer_ids
        else:
            row_labels = [str(i + 1) for i in range(n)]
            col_labels = [str(j + 1) for j in range(n)]

        # Convert energy matrix from Hartree to kcal/mol for output
        energy_matrix_kcal = energy_matrix * HARTREE_TO_KCAL
        df = pd.DataFrame(
            energy_matrix_kcal, index=row_labels, columns=col_labels
        )

        # Create Excel writer
        with pd.ExcelWriter(filename, engine="openpyxl") as writer:
            # Write matrix to 'Energy_Matrix' sheet starting from row 8
            df.to_excel(
                writer,
                sheet_name="Energy_Matrix",
                startrow=8,
                float_format="%.4f",
            )

            # Get the worksheet to add header information
            worksheet = writer.sheets["Energy_Matrix"]

            # Add header information
            row = 1
            worksheet[f"A{row}"] = (
                f"Relative Energy Difference Matrix ({n}x{n}) - {self.__class__.__name__}"
            )
            row += 1
            worksheet[f"A{row}"] = (
                "Values are relative energy differences in kcal/mol: matrix[i,j] = E_j - E_i"
            )
            row += 1
            worksheet[f"A{row}"] = (
                "Positive value means column molecule has higher energy than row molecule"
            )
            row += 1
            worksheet[f"A{row}"] = (
                f"Conversion: 1 Hartree = {HARTREE_TO_KCAL} kcal/mol"
            )
            row += 1

            # Threshold or num_groups
            if self.num_groups is not None:
                worksheet[f"A{row}"] = (
                    f"Requested Groups (-N): {self.num_groups}"
                )
                row += 1
                if self._auto_threshold is not None:
                    # _auto_threshold is already in kcal/mol
                    worksheet[f"A{row}"] = (
                        f"Auto-determined Threshold: {self._auto_threshold:.4f} kcal/mol "
                        f"({self._auto_threshold * KCAL_TO_HARTREE:.10f} Hartree)"
                    )
                    row += 1
            else:
                # self.threshold is in kcal/mol (user input)
                worksheet[f"A{row}"] = (
                    f"Threshold: {self.threshold:.4f} kcal/mol "
                    f"({self.threshold_hartree:.10f} Hartree)"
                )
                row += 1

            # Number of processors
            worksheet[f"A{row}"] = f"Num Procs: {self.num_procs}"
            row += 1

            if grouping_time is not None:
                worksheet[f"A{row}"] = (
                    f"Grouping Time: {grouping_time:.2f} seconds"
                )
                row += 1

            # Auto-adjust column widths
            for column in worksheet.columns:
                max_length = 0
                column_letter = column[0].column_letter
                for cell in column:
                    try:
                        if len(str(cell.value)) > max_length:
                            max_length = len(str(cell.value))
                    except Exception:
                        # Skip cells that cannot be converted to string for width calculation
                        pass
                adjusted_width = min(max_length + 2, 20)
                worksheet.column_dimensions[column_letter].width = (
                    adjusted_width
                )

            # Now add Groups sheet
            groups = self._cached_groups
            index_groups = self._cached_group_indices

            if groups is not None and index_groups is not None:
                groups_data = []
                for i, indices in enumerate(index_groups):
                    # Get conformer IDs if available
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

                logger.info(
                    f"Groups sheet added with {len(groups)} groups to {filename}"
                )

        logger.info(f"Energy matrix saved to {filename}")

    def __repr__(self):
        if self.num_groups is not None:
            return (
                f"{self.__class__.__name__}(num_groups={self.num_groups}, "
                f"num_procs={self.num_procs})"
            )
        else:
            return (
                f"{self.__class__.__name__}(threshold={self.threshold} kcal/mol, "
                f"num_procs={self.num_procs})"
            )
