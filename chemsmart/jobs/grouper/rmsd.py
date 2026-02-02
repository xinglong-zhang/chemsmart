"""
RMSD-based molecular grouping algorithms.

This module contains all RMSD-based grouper implementations:
- RMSDGrouper: Abstract base class for RMSD groupers
- BasicRMSDGrouper: Standard Kabsch RMSD
- HungarianRMSDGrouper: Hungarian algorithm for atom assignment
- SpyRMSDGrouper: Using graph isomorphism for symmetry correction
- IRMSDGrouper: Invariant RMSD with APSP algorithm
- PymolRMSDGrouper: PyMOL-based alignment
- RMSDGrouperSharedMemory: Shared memory optimization
"""

import logging
import multiprocessing
import os
from abc import abstractmethod
from collections import defaultdict
from multiprocessing import RawArray
from typing import Iterable, List, Optional, Tuple

import networkx as nx
import numpy as np
import pandas as pd
from networkx.algorithms import isomorphism
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist

from chemsmart.io.molecules.structure import Molecule
from chemsmart.utils.periodictable import PeriodicTable
from chemsmart.utils.utils import kabsch_align

from .runner import MoleculeGrouper

logger = logging.getLogger(__name__)


class RMSDGrouper(MoleculeGrouper):
    """
    Abstract base class for RMSD-based molecular grouping.

    Groups molecules based on geometric similarity of atomic positions using
    various RMSD calculation methods. This base class provides common
    functionality for all RMSD-based groupers while allowing subclasses
    to implement specific RMSD calculation algorithms.

    Follows the same design pattern as JobRunner - provides the core
    grouping logic while allowing subclasses to implement specific
    RMSD calculation methods via _calculate_rmsd_core().

    Attributes:
        molecules (Iterable[Molecule]): Inherited; collection of molecules to
            group.
        num_procs (int): Inherited; number of worker processes/threads.
        threshold (float): RMSD threshold for grouping molecules.
        align_molecules (bool): Whether to align molecules before RMSD calculation.
        ignore_hydrogens (bool): Whether to exclude hydrogen atoms from RMSD.
    """

    def __init__(
        self,
        molecules: Iterable[Molecule],
        threshold=None,  # RMSD threshold for grouping
        num_groups=None,  # Number of groups to create (alternative to threshold)
        num_procs: int = 1,
        align_molecules: bool = True,
        ignore_hydrogens: bool = False,
        label: str = None,  # Label for output files
        conformer_ids: List[str] = None,  # Custom conformer IDs for labeling
        **kwargs,  # Option to ignore H atoms for grouping
    ):
        """
        Initialize RMSD-based molecular grouper.

        Args:
            molecules (Iterable[Molecule]): Collection of molecules to group.
            threshold (float): RMSD threshold for grouping. Defaults to 0.5.
                Ignored if num_groups is specified.
            num_groups (int): Number of groups to create. When specified,
                automatically determines threshold to create this many groups.
            num_procs (int): Number of processes for parallel computation.
            align_molecules (bool): Whether to align molecules using Kabsch
                algorithm before RMSD calculation. Defaults to True.
            ignore_hydrogens (bool): Whether to exclude hydrogen atoms from
                RMSD calculation. Defaults to False.
            label (str): Label/name for output files. Defaults to None.
            conformer_ids (list[str]): Custom IDs for each molecule (e.g., ['c1', 'c2']).

        Note:
            Uses complete linkage clustering: a structure joins a group only if
            its RMSD to ALL existing members is below the threshold. This prevents
            the chaining effect where dissimilar structures end up in the same
            group through intermediate "bridge" structures.
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
            threshold = 0.5
        self.threshold = threshold  # RMSD threshold for grouping
        self.num_groups = num_groups  # Number of groups to create
        self.align_molecules = align_molecules
        self.ignore_hydrogens = ignore_hydrogens
        # Cache sorted chemical symbols as sets for faster comparison
        self._chemical_symbol_sets = [
            set(mol.chemical_symbols) for mol in molecules
        ]

    def _get_heavy_atoms(self, mol: Molecule) -> Tuple[np.ndarray, List[str]]:
        """
        Extract heavy atoms (non-hydrogen) if ignore_hydrogens is enabled.

        Args:
            mol (Molecule): Molecule to process.

        Returns:
            Tuple[np.ndarray, List[str]]: Tuple containing positions array
                and chemical symbols list (filtered or full based on settings).
        """
        if self.ignore_hydrogens:
            non_h_indices = [
                i for i, sym in enumerate(mol.chemical_symbols) if sym != "H"
            ]
            return mol.positions[non_h_indices], [
                mol.chemical_symbols[i] for i in non_h_indices
            ]
        return (
            mol.positions,
            mol.chemical_symbols,
        )  # Use all atoms if flag is False

    @abstractmethod
    def _calculate_rmsd(self, idx_pair: Tuple[int, int]) -> float:
        """Calculate RMSD between two molecules. Must be implemented by subclasses."""
        pass

    def group(self) -> Tuple[List[List[Molecule]], List[List[int]]]:
        """
        Group molecules by geometric similarity using RMSD clustering.

        Computes pairwise RMSD values between all molecules and groups
        those within the specified threshold using connected components
        clustering, or automatically determines threshold to create
        the specified number of groups. Automatically saves RMSD matrix
        to group_result folder.

        Returns:
            Tuple[List[List[Molecule]], List[List[int]]]: Tuple containing:
                - List of molecule groups (each group is a list of molecules)
                - List of index groups (corresponding indices for each group)
        """
        import time

        # Record start time for grouping process
        grouping_start_time = time.time()

        n = len(self.molecules)
        indices = [(i, j) for i in range(n) for j in range(i + 1, n)]
        total_pairs = len(indices)

        print(
            f"[{self.__class__.__name__}] Starting calculation for {n} molecules ({total_pairs} pairs)"
        )

        # For real-time output, calculate one by one instead of using multiprocessing
        rmsd_values = []
        for idx, (i, j) in enumerate(indices):
            rmsd = self._calculate_rmsd((i, j))
            rmsd_values.append(rmsd)
            print(
                f"The {idx+1}/{total_pairs} pair (conformer{i+1}, conformer{j+1}) calculation finished, RMSD= {rmsd:.7f}"
            )

        # Build full RMSD matrix for output
        rmsd_matrix = np.zeros((n, n))
        for (i, j), rmsd in zip(indices, rmsd_values):
            rmsd_matrix[i, j] = rmsd_matrix[j, i] = rmsd

        # Save RMSD matrix to group_result folder
        # Use label for folder name if provided
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

        # Choose grouping strategy based on parameters (do this first to set _auto_threshold)
        if self.num_groups is not None:
            groups, index_groups = self._group_by_num_groups(
                rmsd_matrix, rmsd_values, indices
            )
        else:
            groups, index_groups = self._group_by_threshold(
                rmsd_values, indices
            )

        # Calculate total grouping time
        grouping_end_time = time.time()
        grouping_time = grouping_end_time - grouping_start_time

        # Cache the results BEFORE saving (so Groups sheet can be populated)
        self._cached_groups = groups
        self._cached_group_indices = index_groups

        # Save full matrix (after grouping to include auto-determined threshold)
        self._save_rmsd_matrix(rmsd_matrix, matrix_filename, grouping_time)

        return groups, index_groups

    def _group_by_threshold(self, rmsd_values, indices):
        """
        Threshold-based grouping using complete linkage.

        A structure joins a group only if its RMSD to ALL existing members
        is below the threshold. This prevents the chaining effect.
        """
        n = len(self.molecules)

        # Build adjacency matrix for clustering
        adj_matrix = np.zeros((n, n), dtype=bool)
        for (i, j), rmsd in zip(indices, rmsd_values):
            if rmsd < self.threshold:
                adj_matrix[i, j] = adj_matrix[j, i] = True

        # Complete linkage grouping
        groups, index_groups = self._complete_linkage_grouping(adj_matrix, n)

        return groups, index_groups

    def _complete_linkage_grouping(self, adj_matrix, n):
        """
        Perform complete linkage grouping (from last to first).

        A structure joins a group only if its RMSD to ALL existing members
        of that group is below the threshold (i.e., all edges exist in adj_matrix).
        This prevents the chaining effect where dissimilar structures end up
        in the same group through intermediate "bridge" structures.

        Iterates from last to first structure so that higher-energy structures
        (typically at the end) are more likely to form singleton groups,
        preserving lower-energy structures in larger groups.

        Args:
            adj_matrix: Boolean adjacency matrix where adj_matrix[i,j] = True
                       if RMSD between i and j is below threshold
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
            # Add the completed group
            groups.append([self.molecules[idx] for idx in current_group])
            index_groups.append(current_group)

        # Reverse groups so that groups containing lower-energy structures come first
        groups.reverse()
        index_groups.reverse()

        return groups, index_groups

    def _group_by_num_groups(self, rmsd_matrix, rmsd_values, indices):
        """Automatic grouping to create specified number of groups."""
        n = len(self.molecules)

        if self.num_groups >= n:
            # If requesting more groups than molecules, each molecule is its own group
            print(
                f"[{self.__class__.__name__}] Requested {self.num_groups} groups but only {n} molecules. Creating {n} groups."
            )
            groups = [[mol] for mol in self.molecules]
            index_groups = [[i] for i in range(n)]
            return groups, index_groups

        # Find appropriate threshold to create desired number of groups
        threshold = self._find_optimal_threshold(rmsd_values, indices, n)

        # Store the auto-determined threshold for summary reporting
        self._auto_threshold = threshold

        print(
            f"[{self.__class__.__name__}] Auto-determined threshold: {threshold:.7f} to create {self.num_groups} groups"
        )

        # Build adjacency matrix with the determined threshold
        adj_matrix = np.zeros((n, n), dtype=bool)
        for (i, j), rmsd in zip(indices, rmsd_values):
            if rmsd < threshold:
                adj_matrix[i, j] = adj_matrix[j, i] = True

        # Use complete linkage grouping
        groups, index_groups = self._complete_linkage_grouping(adj_matrix, n)
        actual_groups = len(groups)

        print(
            f"[{self.__class__.__name__}] Created {actual_groups} groups (requested: {self.num_groups})"
        )

        if actual_groups > self.num_groups:
            # If we have too many groups, merge the smallest ones
            groups, index_groups = self._merge_groups_to_target(
                groups, index_groups, adj_matrix
            )

        return groups, index_groups

    def _find_optimal_threshold(self, rmsd_values, indices, n):
        """Find threshold that creates approximately the desired number of groups using binary search."""
        # Sort RMSD values to try different thresholds
        sorted_rmsd = sorted(
            [rmsd for rmsd in rmsd_values if not np.isinf(rmsd)]
        )

        if not sorted_rmsd:
            return 0.0

        # Binary search for optimal threshold
        low, high = 0, len(sorted_rmsd) - 1
        best_threshold = sorted_rmsd[-1]

        while low <= high:
            mid = (low + high) // 2
            threshold = sorted_rmsd[mid]

            # Build adjacency matrix with this threshold
            adj_matrix = np.zeros((n, n), dtype=bool)
            for (idx_i, idx_j), rmsd in zip(indices, rmsd_values):
                if rmsd < threshold:
                    adj_matrix[idx_i, idx_j] = adj_matrix[idx_j, idx_i] = True

            # Use complete linkage to count groups
            groups, _ = self._complete_linkage_grouping(adj_matrix, n)
            num_groups_found = len(groups)

            if num_groups_found == self.num_groups:
                # Found exact match
                return threshold
            elif num_groups_found > self.num_groups:
                # Too many groups, need higher threshold (more permissive)
                low = mid + 1
            else:
                # Too few groups, need lower threshold (more restrictive)
                best_threshold = threshold
                high = mid - 1

        return best_threshold

    def _merge_groups_to_target(self, groups, index_groups, adj_matrix):
        """
        Merge groups to reach target number when using complete linkage.

        Merges smallest groups into the most compatible larger groups,
        where compatibility is determined by the number of edges in adj_matrix.

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

    def __repr__(self):
        if self.num_groups is not None:
            return (
                f"{self.__class__.__name__}(num_groups={self.num_groups}, "
                f"num_procs={self.num_procs}, align_molecules={self.align_molecules}, "
                f"ignore_hydrogens={self.ignore_hydrogens})"
            )
        else:
            return (
                f"{self.__class__.__name__}(threshold={self.threshold}, "
                f"num_procs={self.num_procs}, align_molecules={self.align_molecules}, "
                f"ignore_hydrogens={self.ignore_hydrogens})"
            )

    def calculate_full_rmsd_matrix(
        self, output_file: Optional[str] = None
    ) -> np.ndarray:
        """
        Calculate the full RMSD matrix for all molecule pairs.

        Args:
            output_file (str, optional): Path to save RMSD matrix as text file

        Returns:
            np.ndarray: Symmetric RMSD matrix (n x n)
        """
        n = len(self.molecules)
        rmsd_matrix = np.zeros((n, n))

        logger.info(f"Calculating full RMSD matrix for {n} molecules")

        # Calculate upper triangular matrix (symmetric)
        indices = [(i, j) for i in range(n) for j in range(i + 1, n)]

        # Use multiprocessing for efficiency
        with multiprocessing.Pool(self.num_procs) as pool:
            rmsd_values = pool.map(self._calculate_rmsd, indices)

        # Fill the matrix
        for (i, j), rmsd in zip(indices, rmsd_values):
            rmsd_matrix[i, j] = rmsd_matrix[j, i] = rmsd

        # Save to file if requested
        if output_file:
            os.makedirs(
                (
                    os.path.dirname(output_file)
                    if os.path.dirname(output_file)
                    else "."
                ),
                exist_ok=True,
            )
            self._save_rmsd_matrix(rmsd_matrix, output_file)

        return rmsd_matrix

    def _save_rmsd_matrix(
        self,
        rmsd_matrix: np.ndarray,
        filename: str,
        grouping_time: float = None,
    ):
        """Save RMSD matrix to Excel file with 7 decimal precision."""
        n = rmsd_matrix.shape[0]

        # Change extension to .xlsx if it's .txt
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

        # Replace inf with string "∞" for display
        matrix_display = np.where(np.isinf(rmsd_matrix), np.nan, rmsd_matrix)
        df = pd.DataFrame(matrix_display, index=row_labels, columns=col_labels)

        # Create Excel writer with openpyxl engine
        with pd.ExcelWriter(filename, engine="openpyxl") as writer:
            # Write matrix to 'RMSD_Matrix' sheet starting from row 8 to leave room for header info
            df.to_excel(
                writer,
                sheet_name="RMSD_Matrix",
                startrow=8,
                float_format="%.7f",
            )

            # Get the worksheet to add header information
            worksheet = writer.sheets["RMSD_Matrix"]

            # Add header information
            row = 1
            worksheet[f"A{row}"] = (
                f"Full RMSD Matrix ({n}x{n}) - {self.__class__.__name__}"
            )
            row += 1

            # Threshold or num_groups
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
                        f"Auto-determined Threshold: {self._auto_threshold:.7f} Å"
                    )
                    row += 1
            else:
                worksheet[f"A{row}"] = f"Threshold: {self.threshold} Å"
                row += 1

            # Parameters specific to different groupers
            worksheet[f"A{row}"] = (
                f"Align Molecules: {getattr(self, 'align_molecules', 'N/A')}"
            )
            row += 1
            worksheet[f"A{row}"] = (
                f"Ignore Hydrogens: {getattr(self, 'ignore_hydrogens', False)}"
            )
            row += 1

            # IRMSDGrouper specific parameters
            if hasattr(self, "check_stereo"):
                worksheet[f"A{row}"] = f"Check Stereo: {self.check_stereo}"
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
                        if cell.value:
                            max_length = max(max_length, len(str(cell.value)))
                    except (TypeError, AttributeError):
                        pass
                adjusted_width = min(max_length + 2, 18)
                worksheet.column_dimensions[column_letter].width = (
                    adjusted_width
                )

            # Add Groups sheet
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

        logger.info(f"RMSD matrix saved to {filename}")


class BasicRMSDGrouper(RMSDGrouper):
    """
    Basic RMSD grouper using standard Euclidean distance calculation.

    Implements the most straightforward RMSD calculation using the standard
    formula: sqrt(mean(sum((pos1 - pos2)^2))). This is the classic RMSD
    implementation that compares atomic positions directly.
    """

    def _calculate_rmsd(self, idx_pair: Tuple[int, int]) -> float:
        """Calculate RMSD between two molecules without atom reordering."""
        i, j = idx_pair
        mol1, mol2 = self.molecules[i], self.molecules[j]

        if self.ignore_hydrogens:
            pos1, symbols1 = self._get_heavy_atoms(mol1)
            pos2, symbols2 = self._get_heavy_atoms(mol2)
        else:
            pos1, symbols1 = mol1.positions, list(mol1.chemical_symbols)
            pos2, symbols2 = mol2.positions, list(mol2.chemical_symbols)

        # Quick check for compatibility
        if len(symbols1) != len(symbols2) or sorted(symbols1) != sorted(
            symbols2
        ):
            return np.inf

        if self.align_molecules:
            logger.debug("Aligning molecules using Kabsch algorithm.")
            pos1, pos2, _, _, _ = kabsch_align(pos1, pos2)
        rmsd = np.sqrt(np.mean(np.sum((pos1 - pos2) ** 2, axis=1)))
        return rmsd


class HungarianRMSDGrouper(RMSDGrouper):
    """
    Hungarian RMSD grouper for optimal atom assignment.

    Uses the Hungarian algorithm (Kuhn-Munkres algorithm) to find the optimal
    assignment of atoms between two molecules that minimizes RMSD. This approach
    handles cases where atom ordering might differ between molecules of the same
    chemical structure.
    """

    def __init__(
        self,
        molecules,
        threshold=None,
        num_groups=None,
        num_procs: int = 1,
        align_molecules: bool = True,
        ignore_hydrogens: bool = False,
        label: str = None,
        **kwargs,
    ):
        super().__init__(
            molecules,
            threshold,
            num_groups,
            num_procs,
            align_molecules,
            ignore_hydrogens,
            label=label,
        )

    def _calculate_rmsd(self, idx_pair):
        i, j = idx_pair
        mol1, mol2 = self.molecules[i], self.molecules[j]

        if self.ignore_hydrogens:
            pos1, symbols1 = self._get_heavy_atoms(mol1)
            pos2, symbols2 = self._get_heavy_atoms(mol2)
        else:
            pos1, symbols1 = mol1.positions, list(mol1.chemical_symbols)
            pos2, symbols2 = mol2.positions, list(mol2.chemical_symbols)

        if len(symbols1) != len(symbols2) or sorted(symbols1) != sorted(
            symbols2
        ):
            return np.inf

        # Use Hungarian algorithm for optimal atom matching
        elements = sorted(set(symbols1))
        matched_idx1 = []
        matched_idx2 = []
        for elem in elements:
            idxs1 = [k for k, s in enumerate(symbols1) if s == elem]
            idxs2 = [k for k, s in enumerate(symbols2) if s == elem]

            if len(idxs1) == 1 and len(idxs2) == 1:
                matched_idx1.extend(idxs1)
                matched_idx2.extend(idxs2)
            else:
                pos1_elem = pos1[idxs1]
                pos2_elem = pos2[idxs2]
                dist_matrix = cdist(pos1_elem, pos2_elem, metric="sqeuclidean")
                row_ind, col_ind = linear_sum_assignment(dist_matrix)
                matched_idx1.extend([idxs1[r] for r in row_ind])
                matched_idx2.extend([idxs2[c] for c in col_ind])

        pos1_matched = pos1[matched_idx1]
        pos2_matched = pos2[matched_idx2]

        if self.align_molecules:
            logger.debug("Aligning molecules using Kabsch algorithm.")
            pos1_matched, pos2_matched, _, _, rmsd = kabsch_align(
                pos1_matched, pos2_matched
            )

        rmsd = np.sqrt(
            np.mean(np.sum((pos1_matched - pos2_matched) ** 2, axis=1))
        )
        return rmsd


class SpyRMSDGrouper(RMSDGrouper):
    """
    SpyRMSD grouper with graph isomorphism symmetry correction.

    Uses graph isomorphism for advanced RMSD calculations with symmetry correction.
    This approach handles molecular symmetries and atom permutations more
    comprehensively than basic Hungarian assignment.
    """

    def __init__(
        self,
        molecules,
        threshold=None,
        num_groups=None,
        num_procs: int = 1,
        align_molecules: bool = True,
        ignore_hydrogens: bool = False,
        cache: bool = True,
        label: str = None,
        **kwargs,
    ):
        super().__init__(
            molecules,
            threshold,
            num_groups,
            num_procs,
            align_molecules,
            ignore_hydrogens,
            label=label,
        )
        self.cache = cache
        self.periodic_table = PeriodicTable()
        self.best_isomorphisms = {}

    @staticmethod
    def _center_coordinates(pos: np.ndarray) -> np.ndarray:
        """Calculate the center of geometry (centroid) of atomic coordinates."""
        return np.mean(pos, axis=0)

    def _symbol_to_atomicnum(self, symbol: str) -> int:
        """Convert element symbol to atomic number using PeriodicTable."""
        try:
            return self.periodic_table.to_atomic_number(symbol)
        except (ValueError, IndexError):
            logger.warning(f"Unknown element symbol: {symbol}")
            return 0

    def _symmrmsd(
        self,
        pos1: np.ndarray,
        pos2: np.ndarray,
        symbols1: list,
        symbols2: list,
        adj_matrix1: np.ndarray,
        adj_matrix2: np.ndarray,
        mol_idx_pair: Tuple[int, int] = None,
    ) -> float:
        """Calculate symmetry-corrected RMSD using graph isomorphism."""
        if len(pos1) != len(pos2):
            if mol_idx_pair:
                self.best_isomorphisms[mol_idx_pair] = None
            return np.inf

        if sorted(symbols1) != sorted(symbols2):
            if mol_idx_pair:
                self.best_isomorphisms[mol_idx_pair] = None
            return np.inf

        G1 = nx.from_numpy_array(adj_matrix1)
        G2 = nx.from_numpy_array(adj_matrix2)

        for i, symbol in enumerate(symbols1):
            G1.nodes[i]["element"] = self._symbol_to_atomicnum(symbol)
        for i, symbol in enumerate(symbols2):
            G2.nodes[i]["element"] = self._symbol_to_atomicnum(symbol)

        node_match = isomorphism.categorical_node_match("element", 0)
        GM = isomorphism.GraphMatcher(G1, G2, node_match=node_match)

        if not GM.is_isomorphic():
            if mol_idx_pair:
                self.best_isomorphisms[mol_idx_pair] = None
            return np.inf

        all_isomorphisms = list(GM.isomorphisms_iter())

        min_rmsd = np.inf
        best_isomorphism = None

        for mapping in all_isomorphisms:
            idx1 = list(range(len(symbols1)))
            idx2 = [mapping[i] for i in range(len(symbols2))]
            reordered_pos2 = np.array(
                [pos2[idx2[i]] for i in range(len(pos2))]
            )

            if self.align_molecules:
                _, _, _, _, rmsd = kabsch_align(pos1, reordered_pos2)
            else:
                diff = pos1 - reordered_pos2
                rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))

            if rmsd < min_rmsd:
                min_rmsd = rmsd
                best_isomorphism = (idx1, idx2)

        if mol_idx_pair:
            self.best_isomorphisms[mol_idx_pair] = best_isomorphism

        return min_rmsd

    def _calculate_rmsd(self, idx_pair):
        i, j = idx_pair
        mol1, mol2 = self.molecules[i], self.molecules[j]

        if self.ignore_hydrogens:
            pos1, symbols1 = self._get_heavy_atoms(mol1)
            pos2, symbols2 = self._get_heavy_atoms(mol2)
        else:
            pos1, symbols1 = mol1.positions, list(mol1.chemical_symbols)
            pos2, symbols2 = mol2.positions, list(mol2.chemical_symbols)

        if len(symbols1) != len(symbols2) or sorted(symbols1) != sorted(
            symbols2
        ):
            self.best_isomorphisms[(i, j)] = None
            return np.inf

        adj_matrix1 = (
            mol1.adjacency_matrix
            if hasattr(mol1, "adjacency_matrix")
            else None
        )
        adj_matrix2 = (
            mol2.adjacency_matrix
            if hasattr(mol2, "adjacency_matrix")
            else None
        )

        if adj_matrix1 is None or adj_matrix2 is None:
            try:
                if adj_matrix1 is None:
                    graph1 = mol1.to_graph()
                    adj_matrix1 = nx.adjacency_matrix(graph1).toarray()
                if adj_matrix2 is None:
                    graph2 = mol2.to_graph()
                    adj_matrix2 = nx.adjacency_matrix(graph2).toarray()
            except Exception as e:
                logger.warning(f"Failed to generate adjacency matrices: {e}")

        if adj_matrix1 is not None and adj_matrix2 is not None:
            try:
                rmsd = self._symmrmsd(
                    pos1=pos1,
                    pos2=pos2,
                    symbols1=symbols1,
                    symbols2=symbols2,
                    adj_matrix1=adj_matrix1,
                    adj_matrix2=adj_matrix2,
                    mol_idx_pair=(i, j),
                )
                return rmsd
            except Exception as e:
                logger.warning(f"Symmetry-corrected RMSD failed: {e}")
                self.best_isomorphisms[(i, j)] = None
                return np.inf

        self.best_isomorphisms[(i, j)] = None
        return np.inf

    def get_best_isomorphism(
        self, mol_idx1: int, mol_idx2: int
    ) -> Optional[Tuple[List[int], List[int]]]:
        """Get the best isomorphism mapping between two molecules."""
        if (mol_idx1, mol_idx2) in self.best_isomorphisms:
            return self.best_isomorphisms[(mol_idx1, mol_idx2)]
        elif (mol_idx2, mol_idx1) in self.best_isomorphisms:
            mapping = self.best_isomorphisms[(mol_idx2, mol_idx1)]
            if mapping is not None:
                return mapping[1], mapping[0]
            return None
        else:
            return None


class IRMSDGrouper(RMSDGrouper):
    """
    Invariant RMSD (iRMSD) Grouper.

    This grouper computes the permutation-invariant RMSD between molecular
    structures using canonical atom identification based on molecular graph
    topology (All-Pair-Shortest-Path algorithm, similar to Morgan algorithm).

    The implementation follows the iRMSD algorithm from:
    J. Chem. Inf. Model. 2025, 65, 4501-4511
    """

    # Standard rotation matrices for symmetry operations
    Rx180 = np.array([[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]])
    Ry180 = np.array([[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]])
    Rz180 = np.array([[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]])
    Rx90 = np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, -1.0, 0.0]])
    Ry90 = np.array([[0.0, 0.0, -1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]])
    Rz90 = np.array([[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    Rx90T = Rx90.T
    Ry90T = Ry90.T
    Rz90T = Rz90.T

    # D3 covalent radii (in Angstroms) for connectivity detection
    _COVALENT_RADII_D3 = {
        1: 0.32,
        2: 0.46,
        3: 1.20,
        4: 0.94,
        5: 0.77,
        6: 0.75,
        7: 0.71,
        8: 0.63,
        9: 0.64,
        10: 0.67,
        11: 1.40,
        12: 1.25,
        13: 1.13,
        14: 1.04,
        15: 1.10,
        16: 1.02,
        17: 0.99,
        18: 0.96,
        19: 1.76,
        20: 1.54,
        21: 1.33,
        22: 1.22,
        23: 1.21,
        24: 1.10,
        25: 1.07,
        26: 1.04,
        27: 1.00,
        28: 0.99,
        29: 1.01,
        30: 1.09,
        31: 1.12,
        32: 1.09,
        33: 1.15,
        34: 1.10,
        35: 1.14,
        36: 1.17,
        37: 1.89,
        38: 1.67,
        39: 1.47,
        40: 1.39,
        41: 1.32,
        42: 1.24,
        43: 1.15,
        44: 1.13,
        45: 1.13,
        46: 1.08,
        47: 1.15,
        48: 1.23,
        49: 1.28,
        50: 1.26,
        51: 1.26,
        52: 1.23,
        53: 1.32,
        54: 1.31,
        55: 2.09,
        56: 1.76,
        57: 1.62,
        58: 1.47,
        59: 1.58,
        60: 1.57,
        61: 1.56,
        62: 1.55,
        63: 1.51,
        64: 1.52,
        65: 1.51,
        66: 1.50,
        67: 1.49,
        68: 1.49,
        69: 1.48,
        70: 1.53,
        71: 1.46,
        72: 1.37,
        73: 1.31,
        74: 1.23,
        75: 1.18,
        76: 1.16,
        77: 1.11,
        78: 1.12,
        79: 1.13,
        80: 1.32,
        81: 1.30,
        82: 1.30,
        83: 1.36,
        84: 1.31,
        85: 1.38,
        86: 1.42,
        87: 2.01,
        88: 1.81,
        89: 1.67,
        90: 1.58,
        91: 1.52,
        92: 1.53,
        93: 1.54,
        94: 1.55,
    }

    _CN_K1 = 16.0
    _CN_K2 = 4.0 / 3.0

    def __init__(
        self,
        molecules: Iterable[Molecule],
        threshold=None,
        num_groups=None,
        num_procs: int = 1,
        align_molecules: bool = False,
        ignore_hydrogens: bool = False,
        check_stereo: str = "auto",
        label: str = None,
        **kwargs,
    ):
        if threshold is None and num_groups is None:
            threshold = 0.125
        super().__init__(
            molecules,
            threshold,
            num_groups,
            num_procs,
            align_molecules,
            ignore_hydrogens,
            label=label,
        )
        self._pt = PeriodicTable()
        self.check_stereo = (
            check_stereo.lower() if isinstance(check_stereo, str) else "auto"
        )
        self._check_stereo_enabled = self.check_stereo != "off"
        logger.info(
            "Using native Python iRMSD implementation with APSP algorithm"
        )

    def _get_covalent_radius(self, atomic_number: int) -> float:
        return self._COVALENT_RADII_D3.get(atomic_number, 1.5)

    def _compute_cn_d3(
        self, atomic_numbers: np.ndarray, positions: np.ndarray
    ) -> np.ndarray:
        n_atoms = len(atomic_numbers)
        cn = np.zeros(n_atoms, dtype=np.float64)
        cov_radii = np.array(
            [self._get_covalent_radius(z) for z in atomic_numbers]
        )
        for i in range(n_atoms):
            for j in range(i + 1, n_atoms):
                rij = np.linalg.norm(positions[i] - positions[j])
                if rij < 1e-6:
                    continue
                rcov_sum = cov_radii[i] + cov_radii[j]
                arg = -self._CN_K1 * (self._CN_K2 * rcov_sum / rij - 1.0)
                cn_contrib = 1.0 / (1.0 + np.exp(arg))
                cn[i] += cn_contrib
                cn[j] += cn_contrib
        return cn

    def _build_connectivity_matrix(
        self,
        atomic_numbers: np.ndarray,
        positions: np.ndarray,
        cn_threshold: float = 0.5,
    ) -> np.ndarray:
        n_atoms = len(atomic_numbers)
        connectivity = np.zeros((n_atoms, n_atoms), dtype=np.int32)
        cov_radii = np.array(
            [self._get_covalent_radius(z) for z in atomic_numbers]
        )
        for i in range(n_atoms):
            for j in range(i + 1, n_atoms):
                rij = np.linalg.norm(positions[i] - positions[j])
                if rij < 1e-6:
                    continue
                rcov_sum = cov_radii[i] + cov_radii[j]
                arg = -self._CN_K1 * (self._CN_K2 * rcov_sum / rij - 1.0)
                cn_contrib = 1.0 / (1.0 + np.exp(arg))
                if cn_contrib > cn_threshold:
                    connectivity[i, j] = 1
                    connectivity[j, i] = 1
        return connectivity

    def _floyd_warshall_apsp(self, connectivity: np.ndarray) -> np.ndarray:
        n = connectivity.shape[0]
        INF = n + 1
        dist = np.full((n, n), INF, dtype=np.int32)
        np.fill_diagonal(dist, 0)
        dist[connectivity == 1] = 1
        for k in range(n):
            for i in range(n):
                for j in range(n):
                    if dist[i, k] + dist[k, j] < dist[i, j]:
                        dist[i, j] = dist[i, k] + dist[k, j]
        return dist

    def _compute_apsp_invariants(
        self, atomic_numbers: np.ndarray, dist_matrix: np.ndarray
    ) -> np.ndarray:
        n = len(atomic_numbers)
        invariants = np.zeros(n, dtype=np.int64)
        P = 31
        M = 10**18 + 9
        for i in range(n):
            inv_list = [
                (int(dist_matrix[i, j]), int(atomic_numbers[j]))
                for j in range(n)
                if i != j
            ]
            inv_list.sort()
            h = int(atomic_numbers[i])
            p_pow = P
            for d, z in inv_list:
                h = (h + (d * 1000 + z) * p_pow) % M
                p_pow = (p_pow * P) % M
            invariants[i] = h
        return invariants

    def _compute_canonical_ranks(
        self, atomic_numbers: np.ndarray, positions: np.ndarray
    ) -> np.ndarray:
        connectivity = self._build_connectivity_matrix(
            atomic_numbers, positions
        )
        dist_matrix = self._floyd_warshall_apsp(connectivity)
        invariants = self._compute_apsp_invariants(atomic_numbers, dist_matrix)
        unique_invs = sorted(set(invariants))
        inv_to_rank = {inv: rank + 1 for rank, inv in enumerate(unique_invs)}
        return np.array(
            [inv_to_rank[inv] for inv in invariants], dtype=np.int32
        )

    def _get_atomic_masses(self, symbols):
        return np.array([self._pt.to_atomic_mass(s) for s in symbols])

    def _compute_cma_and_shift(self, coords, masses):
        total_mass = np.sum(masses) + 1e-20
        cma = np.sum(coords * masses[:, np.newaxis], axis=0) / total_mass
        return coords - cma, cma

    def _compute_inertia_tensor(self, coords, masses):
        t = np.array([i * 1e-10 for i in range(1, 7)])
        for m, r in zip(masses, coords):
            x, y, z = r
            t[0] += m * (y**2 + z**2)
            t[1] -= m * x * y
            t[2] += m * (z**2 + x**2)
            t[3] -= m * z * x
            t[4] -= m * y * z
            t[5] += m * (x**2 + y**2)
        return np.array(
            [[t[0], t[1], t[3]], [t[1], t[2], t[4]], [t[3], t[4], t[5]]]
        )

    def _compute_rotational_constants(self, inertia):
        icm2MHz = 2.9979245e4
        Aamu2icm = 16.8576522
        eig, evec = np.linalg.eigh(inertia)
        evec[np.abs(evec) < 1e-9] = 0.0
        rot = np.zeros(3)
        for i in range(3):
            rot[i] = 0.0 if eig[i] < 3e-4 else icm2MHz * Aamu2icm / eig[i]
        return rot, evec

    def _check_unique_axes(self, rot, thr=0.01):
        unique = np.array([False, False, False])
        if rot[0] < 1e-10 or rot[1] < 1e-10:
            return unique, 3
        diff = np.array(
            [
                abs(rot[1] / rot[0] - 1.0),
                abs(rot[2] / rot[0] - 1.0),
                abs(rot[2] / rot[1] - 1.0),
            ]
        )
        if diff[0] > thr and diff[1] > thr:
            unique[0] = True
        if diff[0] > thr and diff[2] > thr:
            unique[1] = True
        if diff[1] > thr and diff[2] > thr:
            unique[2] = True
        n_unique = np.sum(unique)
        if n_unique == 3:
            return unique, 0
        elif n_unique == 1:
            return unique, 1 if unique[0] else (2 if unique[2] else 3)
        return unique, 3

    def _align_to_principal_axes(self, coords, masses):
        shifted_coords, _ = self._compute_cma_and_shift(coords, masses)
        inertia = self._compute_inertia_tensor(shifted_coords, masses)
        rot, evec = self._compute_rotational_constants(inertia)
        if np.linalg.det(evec) < 0:
            evec[:, 0] = -evec[:, 0]
        return shifted_coords @ evec, rot, evec

    def _rmsd_quaternion(self, ref_coords, mol_coords):
        n_atoms = len(ref_coords)
        x = ref_coords.T.copy()
        y = mol_coords.T.copy()
        x = x - np.mean(x, axis=1, keepdims=True)
        y = y - np.mean(y, axis=1, keepdims=True)
        x_norm = np.sum(x**2)
        y_norm = np.sum(y**2)
        R = x @ y.T
        S = np.array(
            [
                [
                    R[0, 0] + R[1, 1] + R[2, 2],
                    R[1, 2] - R[2, 1],
                    R[2, 0] - R[0, 2],
                    R[0, 1] - R[1, 0],
                ],
                [
                    R[1, 2] - R[2, 1],
                    R[0, 0] - R[1, 1] - R[2, 2],
                    R[0, 1] + R[1, 0],
                    R[0, 2] + R[2, 0],
                ],
                [
                    R[2, 0] - R[0, 2],
                    R[0, 1] + R[1, 0],
                    -R[0, 0] + R[1, 1] - R[2, 2],
                    R[1, 2] + R[2, 1],
                ],
                [
                    R[0, 1] - R[1, 0],
                    R[0, 2] + R[2, 0],
                    R[1, 2] + R[2, 1],
                    -R[0, 0] - R[1, 1] + R[2, 2],
                ],
            ]
        )
        eigenvalues, _ = np.linalg.eigh(S)
        lambda_max = eigenvalues[-1]
        rmsd_sq = max(0.0, (x_norm + y_norm - 2.0 * lambda_max)) / n_atoms
        return np.sqrt(rmsd_sq)

    def _apsp_permutation(self, ref_coords, mol_coords, atomic_numbers, ranks):
        n_atoms = len(atomic_numbers)
        perm = np.arange(n_atoms)
        total_cost = 0.0
        rank_groups = defaultdict(list)
        for i, rank in enumerate(ranks):
            rank_groups[rank].append(i)
        for rank, indices in rank_groups.items():
            if len(indices) <= 1:
                continue
            indices = np.array(indices)
            cost_matrix = cdist(
                ref_coords[indices], mol_coords[indices], metric="sqeuclidean"
            )
            row_ind, col_ind = linear_sum_assignment(cost_matrix)
            perm[indices[row_ind]] = indices[col_ind]
            total_cost += cost_matrix[row_ind, col_ind].sum()
        return perm, total_cost

    def _test_orientations_apsp(
        self, ref_coords, mol_coords, atomic_numbers, ranks, uniqueness_case
    ):
        best_rmsd = np.inf
        best_perm = None
        n_iterations = {0: 1, 1: 2, 2: 2}.get(uniqueness_case, 4)
        mol_working = mol_coords.copy()
        for ii in range(n_iterations):
            for rot_matrix in [None, self.Rx180, self.Ry180, self.Rx180]:
                if rot_matrix is not None:
                    mol_working = mol_working @ rot_matrix.T
                perm, _ = self._apsp_permutation(
                    ref_coords, mol_working, atomic_numbers, ranks
                )
                rmsd = self._rmsd_quaternion(ref_coords, mol_working[perm])
                if rmsd < best_rmsd:
                    best_rmsd = rmsd
                    best_perm = perm.copy()
                if best_rmsd < 1e-10:
                    return best_rmsd, best_perm
            mol_working = mol_working @ self.Ry180.T
            if uniqueness_case == 0:
                break
            elif uniqueness_case == 1 and ii == 0:
                mol_working = mol_working @ self.Rx90.T
            elif uniqueness_case == 2 and ii == 0:
                mol_working = mol_working @ self.Rz90.T
            elif uniqueness_case == 3:
                if ii == 0:
                    mol_working = mol_working @ self.Rz90.T
                elif ii == 1:
                    mol_working = mol_working @ self.Rz90T.T @ self.Ry90.T
                elif ii == 2:
                    mol_working = mol_working @ self.Ry90T.T @ self.Rx90.T
        return best_rmsd, best_perm

    def _find_initial_correspondence(
        self, atomic_numbers1, ranks1, atomic_numbers2, ranks2
    ):
        n_atoms = len(atomic_numbers1)
        perm = np.zeros(n_atoms, dtype=np.int32)
        sig1 = [(ranks1[i], atomic_numbers1[i]) for i in range(n_atoms)]
        sig2 = [(ranks2[i], atomic_numbers2[i]) for i in range(n_atoms)]
        sig_to_indices1 = defaultdict(list)
        sig_to_indices2 = defaultdict(list)
        for i, sig in enumerate(sig1):
            sig_to_indices1[sig].append(i)
        for i, sig in enumerate(sig2):
            sig_to_indices2[sig].append(i)
        if set(sig_to_indices1.keys()) != set(sig_to_indices2.keys()):
            return None
        for sig in sig_to_indices1:
            if len(sig_to_indices1[sig]) != len(sig_to_indices2[sig]):
                return None
        for sig in sig_to_indices1:
            for idx1, idx2 in zip(sig_to_indices1[sig], sig_to_indices2[sig]):
                perm[idx1] = idx2
        return perm

    def _irmsd_core(self, pos1, pos2, atomic_numbers1, atomic_numbers2):
        symbols1 = [self._pt.to_symbol(int(z)) for z in atomic_numbers1]
        symbols2 = [self._pt.to_symbol(int(z)) for z in atomic_numbers2]
        masses1 = np.array([self._pt.to_atomic_mass(s) for s in symbols1])
        masses2 = np.array([self._pt.to_atomic_mass(s) for s in symbols2])
        ranks1 = self._compute_canonical_ranks(atomic_numbers1, pos1)
        ranks2 = self._compute_canonical_ranks(atomic_numbers2, pos2)
        initial_perm = self._find_initial_correspondence(
            atomic_numbers1, ranks1, atomic_numbers2, ranks2
        )
        if initial_perm is None:
            return np.inf
        pos2_reordered = pos2[initial_perm]
        masses2_reordered = masses2[initial_perm]
        ranks = ranks1
        ref_aligned, rot_ref, evec_ref = self._align_to_principal_axes(
            pos1, masses1
        )
        mol_aligned, rot_mol, evec_mol = self._align_to_principal_axes(
            pos2_reordered, masses2_reordered
        )
        _, uniqueness_case = self._check_unique_axes(rot_ref)
        best_rmsd, _ = self._test_orientations_apsp(
            ref_aligned, mol_aligned, atomic_numbers1, ranks, uniqueness_case
        )
        if self._check_stereo_enabled:
            mol_inverted = -mol_aligned
            rmsd_inv, _ = self._test_orientations_apsp(
                ref_aligned,
                mol_inverted,
                atomic_numbers1,
                ranks,
                uniqueness_case,
            )
            if rmsd_inv < best_rmsd:
                best_rmsd = rmsd_inv
        return best_rmsd

    def _calculate_rmsd(self, mol_idx_pair: Tuple[int, int]) -> float:
        i, j = mol_idx_pair
        mol1, mol2 = self.molecules[i], self.molecules[j]
        pos1, pos2 = mol1.positions, mol2.positions
        symbols1, symbols2 = mol1.chemical_symbols, mol2.chemical_symbols
        if sorted(symbols1) != sorted(symbols2):
            logger.warning(f"Molecules {i} and {j} have different atom types")
            return np.inf
        if pos1.shape != pos2.shape:
            logger.warning(f"Molecules {i} and {j} have different shapes")
            return np.inf
        if self.ignore_hydrogens:
            non_h_mask1 = np.array([s != "H" for s in symbols1])
            non_h_mask2 = np.array([s != "H" for s in symbols2])
            pos1, pos2 = pos1[non_h_mask1], pos2[non_h_mask2]
            symbols1 = [s for s in symbols1 if s != "H"]
            symbols2 = [s for s in symbols2 if s != "H"]
        atomic_numbers1 = np.array(
            [self._pt.to_atomic_number(s) for s in symbols1]
        )
        atomic_numbers2 = np.array(
            [self._pt.to_atomic_number(s) for s in symbols2]
        )
        return self._irmsd_core(pos1, pos2, atomic_numbers1, atomic_numbers2)


class PymolRMSDGrouper(RMSDGrouper):
    """Group molecules using PyMOL's align command for RMSD calculation."""

    def __init__(
        self,
        molecules: Iterable[Molecule],
        threshold: float = 0.5,
        num_groups=None,
        num_procs: int = 1,
        ignore_hydrogens: bool = False,
        label: str = None,
        **kwargs,
    ):
        super().__init__(
            molecules=molecules,
            threshold=threshold,
            num_groups=num_groups,
            num_procs=1,
            align_molecules=True,
            ignore_hydrogens=ignore_hydrogens,
            label=label,
            **kwargs,
        )
        self._temp_dir = None
        self._xyz_files = []
        self._mol_names = []
        self._alignment_cache = {}
        self._init_pymol()
        self._prepare_molecules()

    def _init_pymol(self):
        try:
            import pymol
            from pymol import cmd

            pymol.finish_launching(["pymol", "-qc"])
            self.cmd = cmd
            self.cmd.reinitialize()
        except ImportError as e:
            raise ImportError(
                f"PyMOL not available: {e}\nPlease install PyMOL: conda install -c conda-forge pymol-open-source"
            )

    def _prepare_molecules(self):
        import tempfile

        self._temp_dir = tempfile.mkdtemp(prefix="pymol_rmsd_")
        for i, mol in enumerate(self.molecules):
            mol_name = f"mol_{i}"
            xyz_path = os.path.join(self._temp_dir, f"{mol_name}.xyz")
            self._write_xyz(mol, xyz_path)
            self.cmd.load(xyz_path, mol_name)
            self._xyz_files.append(xyz_path)
            self._mol_names.append(mol_name)

    def _write_xyz(self, mol: Molecule, filepath: str):
        positions = mol.positions
        symbols = mol.chemical_symbols
        if self.ignore_hydrogens:
            non_h_indices = [i for i, s in enumerate(symbols) if s != "H"]
            positions = positions[non_h_indices]
            symbols = [symbols[i] for i in non_h_indices]
        with open(filepath, "w") as f:
            f.write(f"{len(symbols)}\n")
            f.write("Generated by PymolRMSDGrouper\n")
            for pos, sym in zip(positions, symbols):
                f.write(f"{sym} {pos[0]:.6f} {pos[1]:.6f} {pos[2]:.6f}\n")

    def _calculate_rmsd(self, idx_pair: Tuple[int, int]) -> float:
        i, j = idx_pair
        cache_key = tuple(sorted([i, j]))
        if cache_key in self._alignment_cache:
            return self._alignment_cache[cache_key]
        mol_i, mol_j = self._mol_names[i], self._mol_names[j]
        try:
            result = self.cmd.align(mol_i, mol_j)
            rmsd = (
                float("inf")
                if result is None or not isinstance(result, (list, tuple))
                else result[0]
            )
            if rmsd < 0 or not np.isfinite(rmsd):
                rmsd = float("inf")
        except Exception as e:
            logger.warning(f"PyMOL align failed for ({i}, {j}): {e}")
            rmsd = float("inf")
        self._alignment_cache[cache_key] = rmsd
        return rmsd

    def __del__(self):
        import shutil

        try:
            if hasattr(self, "cmd"):
                self.cmd.quit()
        except Exception:
            pass
        try:
            if hasattr(self, "_temp_dir") and self._temp_dir:
                shutil.rmtree(self._temp_dir, ignore_errors=True)
        except Exception:
            pass

    def __repr__(self):
        return (
            f"{self.__class__.__name__}(num_groups={self.num_groups})"
            if self.num_groups is not None
            else f"{self.__class__.__name__}(threshold={self.threshold})"
        )


class RMSDGrouperSharedMemory(MoleculeGrouper):
    """Group molecules based on RMSD using shared memory optimization."""

    def __init__(
        self,
        molecules: Iterable[Molecule],
        threshold: float = 0.5,
        num_procs: int = 1,
        align_molecules: bool = True,
        ignore_hydrogens: bool = False,
        label: str = None,
        conformer_ids: List[str] = None,
        **kwargs,
    ):
        super().__init__(
            molecules, num_procs, label=label, conformer_ids=conformer_ids
        )
        self.threshold = threshold
        self.align_molecules = align_molecules
        self.ignore_hydrogens = ignore_hydrogens

    def group(self) -> Tuple[List[List[Molecule]], List[List[int]]]:
        n = len(self.molecules)
        num_atoms = self.molecules[0].positions.shape[0]
        shared_pos = RawArray("d", n * num_atoms * 3)
        pos_np = np.frombuffer(shared_pos, dtype=np.float64).reshape(
            n, num_atoms, 3
        )
        for i, mol in enumerate(self.molecules):
            pos_np[i] = mol.positions
        indices = [(i, j) for i in range(n) for j in range(i + 1, n)]
        with multiprocessing.Pool(
            self.num_procs,
            initializer=self._init_worker,
            initargs=(shared_pos, (n, num_atoms, 3)),
        ) as pool:
            rmsd_values = pool.map(self._calculate_rmsd, indices)
        adj_matrix = np.zeros((n, n), dtype=bool)
        for (i, j), rmsd in zip(indices, rmsd_values):
            if rmsd < self.threshold:
                adj_matrix[i, j] = adj_matrix[j, i] = True
        groups, index_groups = self._complete_linkage_grouping(adj_matrix, n)
        self._cached_groups = groups
        self._cached_group_indices = index_groups
        return groups, index_groups

    def _complete_linkage_grouping(self, adj_matrix, n):
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
                if all(adj_matrix[j, member] for member in current_group):
                    current_group.append(j)
                    assigned[j] = True
            current_group.sort()
            groups.append([self.molecules[idx] for idx in current_group])
            index_groups.append(current_group)
        groups.reverse()
        index_groups.reverse()
        return groups, index_groups

    @staticmethod
    def _init_worker(shared_pos, pos_shape):
        global shared_positions
        shared_positions = np.frombuffer(shared_pos, dtype=np.float64).reshape(
            pos_shape
        )

    def _calculate_rmsd(self, idx_pair: Tuple[int, int]) -> float:
        i, j = idx_pair
        pos1 = np.array(shared_positions[i])
        pos2 = np.array(shared_positions[j])
        if pos1.shape != pos2.shape:
            return np.inf
        if self.align_molecules:
            pos1, pos2, _, _, _ = kabsch_align(pos1, pos2)
        return np.sqrt(np.mean(np.sum((pos1 - pos2) ** 2, axis=1)))


__all__ = [
    "RMSDGrouper",
    "BasicRMSDGrouper",
    "HungarianRMSDGrouper",
    "SpyRMSDGrouper",
    "IRMSDGrouper",
    "PymolRMSDGrouper",
    "RMSDGrouperSharedMemory",
]
