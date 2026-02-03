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
from multiprocessing import RawArray
from typing import Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist

from chemsmart.io.molecules.structure import Molecule
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
                "Cannot specify both threshold (-T) and num_groups (-N). Please use only one."
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
            actual_inversion = getattr(self, "_actual_inversion", None)
            if actual_inversion is not None:
                worksheet[f"A{row}"] = f"Inversion: {actual_inversion}"
                row += 1

            # SpyRMSDGrouper specific parameters
            cache = getattr(self, "cache", None)
            if cache is not None:
                worksheet[f"A{row}"] = f"Cache: {cache}"
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
    SpyRMSD grouper using the spyrmsd package for symmetry-corrected RMSD.

    Uses the spyrmsd library for advanced RMSD calculations with symmetry
    correction via graph isomorphism. This approach handles molecular
    symmetries and atom permutations comprehensively.

    Attributes:
        minimize (bool): Whether to minimize RMSD by centering and rotating.
        cache (bool): Whether to cache graph isomorphisms for performance.
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
        conformer_ids: List[str] = None,
        **kwargs,
    ):
        """
        Initialize SpyRMSD grouper.

        Args:
            molecules: Collection of molecules to group.
            threshold: RMSD threshold for grouping.
            num_groups: Number of groups to create (alternative to threshold).
            num_procs: Number of processes for parallel computation.
            align_molecules: Whether to minimize RMSD (center and rotate).
            ignore_hydrogens: Whether to exclude hydrogen atoms.
            cache: Whether to cache graph isomorphisms.
            label: Label for output files.
            conformer_ids: Custom IDs for each molecule.
        """
        super().__init__(
            molecules,
            threshold,
            num_groups,
            num_procs,
            align_molecules,
            ignore_hydrogens,
            label=label,
            conformer_ids=conformer_ids,
        )
        self.cache = cache
        self.minimize = align_molecules  # spyrmsd uses 'minimize' parameter
        self.best_isomorphisms = {}

        # Pre-compute molecule data for spyrmsd
        self._prepare_spyrmsd_data()

    def _prepare_spyrmsd_data(self):
        """Pre-compute coordinates, atomic numbers, and adjacency matrices.

        Uses spyrmsd's native io.loadmol (via OpenBabel) to ensure
        adjacency matrices are generated exactly the same way as the
        original spyrmsd package.
        """
        import os
        import tempfile

        from spyrmsd import io as spy_io

        self._coords_list = []
        self._anum_list = []
        self._adj_list = []

        for mol in self.molecules:
            if self.ignore_hydrogens:
                pos, symbols = self._get_heavy_atoms(mol)
            else:
                pos = mol.positions
                symbols = list(mol.chemical_symbols)

            # Write molecule to temp xyz file and load with spyrmsd
            # This ensures adjacency matrix is generated exactly like original spyrmsd
            with tempfile.NamedTemporaryFile(
                mode="w", suffix=".xyz", delete=False
            ) as f:
                f.write(f"{len(symbols)}\n")
                f.write("temp\n")
                for sym, p in zip(symbols, pos):
                    f.write(f"{sym} {p[0]:.10f} {p[1]:.10f} {p[2]:.10f}\n")
                tmp_file = f.name

            try:
                # Load with spyrmsd's io (uses OpenBabel internally)
                spy_mol = spy_io.loadmol(tmp_file)
                coords = spy_mol.coordinates
                anum = spy_mol.atomicnums
                adj = spy_mol.adjacency_matrix
            finally:
                # Clean up temp file
                os.unlink(tmp_file)

            self._coords_list.append(coords)
            self._anum_list.append(anum)
            self._adj_list.append(adj)

    def _calculate_rmsd(self, idx_pair):
        """Calculate symmetry-corrected RMSD using spyrmsd package."""
        from spyrmsd import rmsd as spyrmsd_rmsd

        i, j = idx_pair

        coords_ref = self._coords_list[i]
        coords_other = self._coords_list[j]
        anum_ref = self._anum_list[i]
        anum_other = self._anum_list[j]
        adj_ref = self._adj_list[i]
        adj_other = self._adj_list[j]

        # Check compatibility
        if len(coords_ref) != len(coords_other):
            self.best_isomorphisms[(i, j)] = None
            return np.inf

        if not np.array_equal(np.sort(anum_ref), np.sort(anum_other)):
            self.best_isomorphisms[(i, j)] = None
            return np.inf

        try:
            # Use spyrmsd.symmrmsd with single coordinate (not list)
            # symmrmsd expects coords as single array or list of arrays
            result = spyrmsd_rmsd.symmrmsd(
                coords_ref,
                coords_other,  # Single coordinate array
                anum_ref,
                anum_other,
                adj_ref,
                adj_other,
                minimize=self.minimize,
                cache=self.cache,
                return_best_isomorphism=True,
            )

            # Result is (rmsd, isomorphism) when return_best_isomorphism=True
            if isinstance(result, tuple):
                rmsd_value, best_iso = result
                self.best_isomorphisms[(i, j)] = best_iso
            else:
                rmsd_value = result
                self.best_isomorphisms[(i, j)] = None

            return float(rmsd_value)

        except Exception as e:
            logger.warning(
                f"spyrmsd calculation failed for pair ({i}, {j}): {e}"
            )
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
    Invariant RMSD (iRMSD) Grouper using external irmsd package.

    This grouper computes the permutation-invariant RMSD between molecular
    structures by calling the external 'irmsd' command-line tool via subprocess.

    The iRMSD algorithm:
    - Assigns canonical atom identities independent of input atom order
    - Performs symmetry-aware alignment using principal axes
    - Solves the linear sum assignment problem (LSAP, Hungarian algorithm)
    - Handles false enantiomers via z-mirror checking

    Requirements:
        The 'irmsd' package must be installed in a separate conda environment
        with numpy>=2.0. Set environment variable to specify the irmsd location:

        Option 1: Set IRMSD_CONDA_ENV to the conda environment name
            export IRMSD_CONDA_ENV=irmsd_env

        Option 2: Set IRMSD_PATH to the full path of the irmsd executable
            export IRMSD_PATH=/path/to/conda/envs/irmsd_env/bin/irmsd

        Setup example:
            conda create -n irmsd_env python=3.10 numpy>=2.0
            conda activate irmsd_env
            pip install irmsd
            conda deactivate
            export IRMSD_CONDA_ENV=irmsd_env

    Reference: J. Chem. Inf. Model. 2025, 65, 4501-4511
    """

    def __init__(
        self,
        molecules: Iterable[Molecule],
        threshold=None,
        num_groups=None,
        num_procs: int = 1,
        align_molecules: bool = False,
        ignore_hydrogens: bool = False,
        inversion: str = "auto",
        label: str = None,
        conformer_ids: List[str] = None,
        **kwargs,
    ):
        """
        Initialize IRMSDGrouper.

        Args:
            molecules: Collection of molecules to group.
            threshold: RMSD threshold for grouping (default: 0.125).
            num_groups: Number of groups to create (alternative to threshold).
            num_procs: Number of processes (note: external calls are sequential).
            align_molecules: Not used (irmsd handles alignment internally).
            ignore_hydrogens: Whether to use only heavy atoms (--heavy flag).
            inversion: Inversion check mode: 'on', 'off', or 'auto' (default).
            label: Label for output files.
            conformer_ids: Custom IDs for each molecule.
        """
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
            conformer_ids=conformer_ids,
        )
        self.inversion = (
            inversion.lower() if isinstance(inversion, str) else "auto"
        )
        self._irmsd_cmd = self._find_irmsd_command()
        if self._irmsd_cmd:
            logger.info(f"Using irmsd command: {self._irmsd_cmd}")
        else:
            raise RuntimeError(
                "irmsd command not found. IRMSDGrouper requires the irmsd package.\n\n"
                "Quick setup (one-line command):\n"
                "  conda create -n irmsd_env python=3.10 'numpy>=2.0' -y && "
                "conda run -n irmsd_env pip install irmsd\n\n"
                "Then set environment variable:\n"
                "  export IRMSD_CONDA_ENV=irmsd_env\n\n"
                "To make it permanent, add to ~/.zshrc or ~/.bashrc,\n"
                "then run: source ~/.zshrc or source ~/.bashrc"
            )

    def _find_irmsd_command(self) -> Optional[str]:
        """
        Find the irmsd command, checking multiple sources.

        Priority:
        1. IRMSD_PATH environment variable (full path to irmsd executable)
        2. IRMSD_CONDA_ENV environment variable (conda environment name)
        3. irmsd in current PATH

        Returns:
            Full path to irmsd command, or None if not found.
        """
        import os
        import shutil

        # Option 1: IRMSD_PATH environment variable
        irmsd_path = os.environ.get("IRMSD_PATH")
        if (
            irmsd_path
            and os.path.isfile(irmsd_path)
            and os.access(irmsd_path, os.X_OK)
        ):
            return irmsd_path

        # Option 2: IRMSD_CONDA_ENV environment variable
        conda_env = os.environ.get("IRMSD_CONDA_ENV")
        if conda_env:
            # Try to find conda base path
            conda_base = os.environ.get("CONDA_PREFIX_1") or os.environ.get(
                "CONDA_PREFIX"
            )
            if conda_base:
                # Go up to envs directory
                if "envs" in conda_base:
                    envs_dir = conda_base.rsplit("envs", 1)[0] + "envs"
                else:
                    envs_dir = os.path.join(
                        os.path.dirname(conda_base), "envs"
                    )

                irmsd_path = os.path.join(envs_dir, conda_env, "bin", "irmsd")
                if os.path.isfile(irmsd_path) and os.access(
                    irmsd_path, os.X_OK
                ):
                    return irmsd_path

            # Try common conda locations
            for base in [
                os.path.expanduser("~/anaconda3"),
                os.path.expanduser("~/miniconda3"),
                "/opt/anaconda3",
                "/opt/miniconda3",
                os.path.expanduser("~/.conda"),
            ]:
                irmsd_path = os.path.join(
                    base, "envs", conda_env, "bin", "irmsd"
                )
                if os.path.isfile(irmsd_path) and os.access(
                    irmsd_path, os.X_OK
                ):
                    return irmsd_path

        # Option 3: irmsd in current PATH
        irmsd_in_path = shutil.which("irmsd")
        if irmsd_in_path:
            return irmsd_in_path

        return None

    def _write_two_molecules_xyz(
        self, mol1: Molecule, mol2: Molecule, filepath: str
    ) -> None:
        """Write two molecules to a single XYZ file for irmsd compare."""
        with open(filepath, "w") as f:
            for mol in [mol1, mol2]:
                n_atoms = len(mol.chemical_symbols)
                energy = mol.energy if mol.energy is not None else 0.0
                f.write(f"{n_atoms}\n")
                f.write(f"Energy = {energy}\n")
                for symbol, pos in zip(mol.chemical_symbols, mol.positions):
                    f.write(
                        f"{symbol:2s} {pos[0]:15.8f} {pos[1]:15.8f} {pos[2]:15.8f}\n"
                    )

    def _parse_irmsd_output(
        self, output: str, parse_inversion: bool = True
    ) -> Tuple[float, Optional[str]]:
        """Parse iRMSD value and inversion setting from irmsd compare output.

        Args:
            output: The stdout from irmsd compare command.
            parse_inversion: Whether to parse inversion setting (skip if already obtained).

        Returns:
            Tuple[float, Optional[str]]: RMSD value and actual inversion setting used.
        """
        rmsd_value = np.inf
        actual_inversion = None

        for line in output.split("\n"):
            if "iRMSD:" in line:
                parts = line.split()
                try:
                    rmsd_value = float(parts[1])
                except (IndexError, ValueError):
                    pass
            elif parse_inversion and "Inversion check:" in line:
                # Parse "Inversion check: on/off/auto"
                parts = line.split(":")
                if len(parts) >= 2:
                    actual_inversion = parts[1].strip()

        return rmsd_value, actual_inversion

    def _calculate_rmsd(self, mol_idx_pair: Tuple[int, int]) -> float:
        """Calculate iRMSD between two molecules using external irmsd tool."""
        import subprocess
        import tempfile

        i, j = mol_idx_pair
        mol1, mol2 = self.molecules[i], self.molecules[j]

        # Check compatibility
        symbols1, symbols2 = mol1.chemical_symbols, mol2.chemical_symbols
        if sorted(symbols1) != sorted(symbols2):
            logger.warning(f"Molecules {i} and {j} have different atom types")
            return np.inf

        # Create temporary XYZ file with both molecules
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".xyz", delete=False
        ) as tmp_file:
            tmp_path = tmp_file.name
            self._write_two_molecules_xyz(mol1, mol2, tmp_path)

        try:
            # Build irmsd compare command
            cmd = [self._irmsd_cmd, "compare", tmp_path, "--ref-idx", "0"]

            # Add inversion option
            if self.inversion in ["on", "off", "auto"]:
                cmd.extend(["--inversion", self.inversion])

            # Add heavy atom option if ignoring hydrogens
            if self.ignore_hydrogens:
                cmd.append("--heavy")

            # Run irmsd compare
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=60,  # 60 second timeout per pair
            )

            if result.returncode != 0:
                logger.warning(
                    f"irmsd failed for pair ({i}, {j}): {result.stderr}"
                )
                return np.inf

            # Parse output (only parse inversion on first call)
            need_inversion = not hasattr(self, "_actual_inversion")
            rmsd_value, actual_inversion = self._parse_irmsd_output(
                result.stdout, parse_inversion=need_inversion
            )

            # Store the actual inversion setting (only need to do this once)
            if actual_inversion:
                self._actual_inversion = actual_inversion

            return rmsd_value

        except subprocess.TimeoutExpired:
            logger.warning(f"irmsd timed out for pair ({i}, {j})")
            return np.inf
        except Exception as e:
            logger.warning(f"irmsd error for pair ({i}, {j}): {e}")
            return np.inf
        finally:
            # Clean up temporary file
            import os

            try:
                os.unlink(tmp_path)
            except OSError:
                pass


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
        # PyMOL only supports single-threaded operation
        if num_procs > 1:
            raise ValueError(
                f"PymolRMSDGrouper only supports single-threaded operation (num_procs=1), "
                f"got num_procs={num_procs}. PyMOL cannot run in parallel."
            )

        super().__init__(
            molecules=molecules,
            threshold=threshold,
            num_groups=num_groups,
            num_procs=1,  # Always force to 1 for PyMOL
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
