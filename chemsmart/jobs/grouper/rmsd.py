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
    Invariant RMSD (iRMSD) Grouper.

    This grouper computes the permutation-invariant RMSD between molecular
    structures. The implementation exactly follows the Fortran CREST code:
    https://github.com/crest-lab/crest (src/sorting/irmsd_module.f90)

    The iRMSD algorithm:
    - Assigns canonical atom identities independent of input atom order
    - Performs symmetry-aware alignment using principal axes
    - Solves the linear sum assignment problem (LSAP, Hungarian algorithm)
    - Handles false enantiomers via z-mirror checking

    Reference: J. Chem. Inf. Model. 2025, 65, 4501-4511
    """

    # Rotation matrices - exactly as defined in Fortran CREST
    # Fortran uses column-major order, these are transposed for Python's row-major
    # In Fortran: reshape([1,0,0, 0,-1,0, 0,0,-1], [3,3]) gives column-major
    # To apply: mol_xyz = matmul(R, mol_xyz) in Fortran = mol_xyz @ R.T in Python (for row vectors)
    # But Fortran stores xyz as (3, nat), so matmul(R, xyz) works directly
    # In Python with xyz as (nat, 3), we do xyz @ R.T
    Rx180 = np.array([[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]])
    Ry180 = np.array([[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]])
    Rz180 = np.array([[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]])
    Rx90 = np.array([[1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]])
    Ry90 = np.array([[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]])
    Rz90 = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
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

    def _compute_apsp_invariants_fortran(
        self, atomic_numbers: np.ndarray, dist_matrix: np.ndarray
    ) -> np.ndarray:
        """
        Compute APSP invariants following the Fortran CREST implementation.

        The Fortran algorithm:
        1. Start with rinv = 1 for all atoms
        2. For each distance d from 1 to maxdist:
           - For each atom j, sum up rinv[k] for all k where dist[k,j] == d
           - Add this to rinv[j]
        3. Rank atoms by their final rinv values (higher = lower rank number)
        """
        n = len(atomic_numbers)
        max_dist = int(
            np.max(dist_matrix[dist_matrix < n + 1])
        )  # Exclude INF values

        rinv = np.ones(n, dtype=np.float64)

        for d in range(1, max_dist + 1):
            tmp_rinv = np.zeros(n, dtype=np.float64)
            for j in range(n):
                for k in range(n):
                    if dist_matrix[k, j] == d:
                        tmp_rinv[j] += rinv[k]
            rinv += tmp_rinv

        # Convert to invariants (higher rinv = lower rank number)
        # We'll use negative rinv so that sorting gives correct order
        return rinv

    def _compute_canonical_ranks(
        self, atomic_numbers: np.ndarray, positions: np.ndarray
    ) -> np.ndarray:
        """
        Compute canonical atom ranks following Fortran CREST implementation.

        Uses APSP (All-Pairs-Shortest-Path) based invariants with:
        - APSP distance-based invariant
        - Atomic number
        - Number of hydrogen neighbors

        This matches Fortran: inv = inv*1000 + ati*10 + hneigh
        """
        connectivity = self._build_connectivity_matrix(
            atomic_numbers, positions
        )
        dist_matrix = self._floyd_warshall_apsp(connectivity)

        # Get APSP-based invariants
        rinv = self._compute_apsp_invariants_fortran(
            atomic_numbers, dist_matrix
        )

        # Count hydrogen neighbors for each atom
        h_neighbors = np.zeros(len(atomic_numbers), dtype=np.int32)
        for i in range(len(atomic_numbers)):
            for j in range(len(atomic_numbers)):
                if (
                    connectivity[i, j] == 1 and atomic_numbers[j] == 1
                ):  # H is atomic number 1
                    h_neighbors[i] += 1

        # Update invariants with atomic number and H neighbor count (like Fortran)
        # inv = inv * 1000 + ati * 10 + hneigh
        invariants = rinv * 1000 + atomic_numbers * 10 + h_neighbors

        # Rank by invariant value (higher invariant = lower rank number)
        sorted_indices = np.argsort(-invariants)
        ranks = np.zeros(len(atomic_numbers), dtype=np.int32)

        current_rank = 1
        prev_inv = None
        for i, idx in enumerate(sorted_indices):
            if prev_inv is not None and not np.isclose(
                invariants[idx], prev_inv, rtol=1e-9
            ):
                current_rank = i + 1
            ranks[idx] = current_rank
            prev_inv = invariants[idx]

        return ranks

    def _compute_canonical_ranks_both(
        self,
        atomic_numbers1: np.ndarray,
        positions1: np.ndarray,
        atomic_numbers2: np.ndarray,
        positions2: np.ndarray,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute canonical ranks for both molecules.

        Returns:
            ranks1, ranks2: canonical rank arrays for both molecules
        """
        ranks1 = self._compute_canonical_ranks(atomic_numbers1, positions1)
        ranks2 = self._compute_canonical_ranks(atomic_numbers2, positions2)
        return ranks1, ranks2

    def _fallback_ranks(self, atomic_numbers: np.ndarray) -> np.ndarray:
        """
        Fallback ranking using only atom types (like Fortran fallbackranks).

        This is used when canonical ranking is not available or fails.
        """
        unique_types = sorted(set(atomic_numbers))
        type_to_rank = {t: i + 1 for i, t in enumerate(unique_types)}
        return np.array(
            [type_to_rank[z] for z in atomic_numbers], dtype=np.int32
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
        """Check which rotational axes are unique.

        Exactly matches Fortran uniqueax subroutine:

        diff(1) = abs(rot(2)/rot(1) - 1.0)  ! B/A
        diff(2) = abs(rot(3)/rot(1) - 1.0)  ! C/A
        diff(3) = abs(rot(3)/rot(2) - 1.0)  ! C/B

        unique(1) = diff(1) > thr AND diff(2) > thr  (A unique)
        unique(2) = diff(1) > thr AND diff(3) > thr  (B unique)
        unique(3) = diff(2) > thr AND diff(3) > thr  (C unique)

        In Python (0-indexed):
        rot[0] = A, rot[1] = B, rot[2] = C
        """
        unique = np.array([False, False, False])

        # Avoid division by zero
        if rot[0] < 1e-10 or rot[1] < 1e-10:
            return unique, 3

        # Match Fortran exactly (note: Fortran is 1-indexed)
        # diff(1) = rot(2)/rot(1) = rot[1]/rot[0] = B/A
        # diff(2) = rot(3)/rot(1) = rot[2]/rot[0] = C/A
        # diff(3) = rot(3)/rot(2) = rot[2]/rot[1] = C/B
        diff = np.array(
            [
                abs(rot[1] / rot[0] - 1.0),  # B/A
                abs(rot[2] / rot[0] - 1.0),  # C/A
                abs(rot[2] / rot[1] - 1.0),  # C/B
            ]
        )

        # unique(1) = diff(1) > thr AND diff(2) > thr
        if diff[0] > thr and diff[1] > thr:
            unique[0] = True  # A is unique
        # unique(2) = diff(1) > thr AND diff(3) > thr
        if diff[0] > thr and diff[2] > thr:
            unique[1] = True  # B is unique
        # unique(3) = diff(2) > thr AND diff(3) > thr
        if diff[1] > thr and diff[2] > thr:
            unique[2] = True  # C is unique

        n_unique = np.sum(unique)

        # Fortran logic:
        # select case(nunique)
        # case ( 3 ) -> uniquenesscase = 0
        # case ( 1 ) -> if(unique(1)) uniquenesscase = 1; if(unique(3)) uniquenesscase = 2
        # case ( 0 ) -> uniquenesscase = 3
        # end select
        # Note: B unique (nunique==1, unique(2)==True) falls through with uniquenesscase = 0

        if n_unique == 3:
            return unique, 0
        elif n_unique == 1:
            if unique[0]:  # A unique (Fortran unique(1))
                return unique, 1
            elif unique[2]:  # C unique (Fortran unique(3))
                return unique, 2
            else:  # B unique - Fortran leaves uniquenesscase = 0
                return unique, 0
        elif n_unique == 0:
            return unique, 3
        else:  # n_unique == 2
            return unique, 0  # Fortran default

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

    def _fallback_ranks_both(self, atomic_numbers1, atomic_numbers2):
        """
        Compute fallback ranks for both molecules based on atom types.

        This exactly matches Fortran fallbackranks subroutine.
        Creates a unified typemap from both molecules, then assigns ranks.

        Returns:
            ranks1, ranks2: rank arrays for both molecules
        """

        # Build typemap: collect unique atom types from both molecules
        typemap = []
        for at in atomic_numbers1:
            if at not in typemap:
                typemap.append(at)
        for at in atomic_numbers2:
            if at not in typemap:
                typemap.append(at)

        # Build reverse typemap: atom_type -> rank
        rtypemap = {at: i + 1 for i, at in enumerate(typemap)}

        # Assign ranks
        ranks1 = np.array(
            [rtypemap[at] for at in atomic_numbers1], dtype=np.int32
        )
        ranks2 = np.array(
            [rtypemap[at] for at in atomic_numbers2], dtype=np.int32
        )

        return ranks1, ranks2

    def _rank_to_order(self, ranks):
        """
        Convert rank array to order array.

        Exactly matches Fortran rank_2_order subroutine.
        """
        n = len(ranks)
        order = np.zeros(n, dtype=np.int32)
        max_rank = np.max(ranks) if len(ranks) > 0 else 0
        k = 0
        for rank in range(1, max_rank + 1):
            for j in range(n):
                if ranks[j] == rank:
                    k += 1
                    order[j] = k
        return order

    def _mol_atom_sort(self, coords, masses, current_order, target_order):
        """
        Sort molecule atoms from current_order to target_order.

        Exactly matches Fortran molatomsort subroutine.
        """
        n = len(current_order)
        coords_out = coords.copy()
        masses_out = masses.copy()

        # Create index map: target_order[i] tells us what position atom i should go to
        index_map = np.zeros(n, dtype=np.int32)
        for i in range(n):
            index_map[current_order[i] - 1] = (
                i  # -1 because Fortran is 1-indexed
            )

        # Reorder atoms
        for i in range(n):
            correct_atom = target_order[i] - 1  # -1 for 0-indexing
            current_position = index_map[correct_atom]

            if i != current_position:
                # Swap
                coords_out[[i, current_position]] = coords_out[
                    [current_position, i]
                ]
                masses_out[[i, current_position]] = masses_out[
                    [current_position, i]
                ]

                # Update index_map
                index_map[current_order[i] - 1] = current_position
                index_map[current_order[current_position] - 1] = i

                # Update current_order
                current_order[i], current_order[current_position] = (
                    current_order[current_position],
                    current_order[i],
                )

        return coords_out, masses_out

    def _compute_lsap_for_rank_groups(
        self, ref_coords, mol_coords, ranks, ngroup, nranks
    ):
        """
        Iterate through rank groups and solve LSAP for each.

        Exactly matches Fortran min_rmsd_iterate_through_groups subroutine.

        Returns:
            iwork: permutation array
            total_cost: total LSAP cost
        """
        n_atoms = len(ranks)
        iwork = np.arange(n_atoms, dtype=np.int32)  # Initialize to identity
        total_cost = 0.0

        for rr in range(1, nranks + 1):
            if ngroup[rr] <= 1:
                # Skip ranks with only one atom (already assigned)
                continue

            # Get indices of atoms with this rank
            ref_indices = np.where(ranks == rr)[0]
            mol_indices = (
                ref_indices.copy()
            )  # Same indices since ranks are equal

            rnknat = len(ref_indices)

            # Build cost matrix (squared distances)
            cost_matrix = np.zeros((rnknat, rnknat), dtype=np.float32)
            for ii, ref_idx in enumerate(ref_indices):
                for jj, mol_idx in enumerate(mol_indices):
                    diff = ref_coords[ref_idx] - mol_coords[mol_idx]
                    cost_matrix[ii, jj] = np.sum(diff**2)

            # Solve LSAP using scipy
            row_ind, col_ind = linear_sum_assignment(cost_matrix)

            # Update iwork with the mapping
            for ii, jj in zip(row_ind, col_ind):
                iwork[ref_indices[ii]] = mol_indices[jj]
                total_cost += cost_matrix[ii, jj]

        return iwork, total_cost

    def _min_rmsd_rotcheck_permute(
        self,
        ref_coords,
        mol_coords,
        ranks,
        ngroup,
        nranks,
        step,
        uniquenesscase,
    ):
        """
        Test different orientations and compute LSAP costs.

        Exactly matches Fortran min_rmsd_rotcheck_permute subroutine.

        Returns:
            vals: array of 16 LSAP costs
            order_backup: list of 16 permutation arrays
            mol_coords_backup: list of 16 mol coordinate arrays (for reconstructing later)
        """
        INF = np.inf
        vals = np.full(16, INF)
        order_backup = [None] * 16
        mol_coords_backup = [
            None
        ] * 16  # Store mol coordinates for each orientation

        mol = mol_coords.copy()

        # ALIGNLOOP: do ii=1,4
        for ii in range(4):
            # Test 4 orientations with 180-degree rotations

            # Orientation 1: no rotation (from current state)
            iwork, cost = self._compute_lsap_for_rank_groups(
                ref_coords, mol, ranks, ngroup, nranks
            )
            vals[0 + 4 * ii] = cost
            order_backup[0 + 4 * ii] = iwork.copy()
            mol_coords_backup[0 + 4 * ii] = mol.copy()

            # Orientation 2: Rx180
            mol = mol @ self.Rx180.T
            iwork, cost = self._compute_lsap_for_rank_groups(
                ref_coords, mol, ranks, ngroup, nranks
            )
            vals[1 + 4 * ii] = cost
            order_backup[1 + 4 * ii] = iwork.copy()
            mol_coords_backup[1 + 4 * ii] = mol.copy()

            # Orientation 3: Ry180 (from Rx180 state)
            mol = mol @ self.Ry180.T
            iwork, cost = self._compute_lsap_for_rank_groups(
                ref_coords, mol, ranks, ngroup, nranks
            )
            vals[2 + 4 * ii] = cost
            order_backup[2 + 4 * ii] = iwork.copy()
            mol_coords_backup[2 + 4 * ii] = mol.copy()

            # Orientation 4: Rx180 (from Rx180+Ry180 state)
            mol = mol @ self.Rx180.T
            iwork, cost = self._compute_lsap_for_rank_groups(
                ref_coords, mol, ranks, ngroup, nranks
            )
            vals[3 + 4 * ii] = cost
            order_backup[3 + 4 * ii] = iwork.copy()
            mol_coords_backup[3 + 4 * ii] = mol.copy()

            # Restore: Ry180
            mol = mol @ self.Ry180.T

            # Handle uniqueness cases - exactly as Fortran
            if uniquenesscase == 0:
                # 3 unique principal axes - exit after first iteration
                break
            elif uniquenesscase == 1:
                # Only A unique
                if ii == 1:
                    mol = mol @ self.Rx90T.T
                    break
                mol = mol @ self.Rx90.T
            elif uniquenesscase == 2:
                # Only C unique
                if ii == 1:
                    mol = mol @ self.Rz90T.T
                    break
                mol = mol @ self.Rz90.T
            elif uniquenesscase == 3:
                # Rotationally ambiguous
                if ii == 0:
                    mol = mol @ self.Rz90.T
                elif ii == 1:
                    mol = mol @ self.Rz90T.T
                    mol = mol @ self.Ry90.T
                elif ii == 2:
                    mol = mol @ self.Ry90T.T
                    mol = mol @ self.Rx90.T
                else:
                    mol = mol @ self.Rx90T.T
                    break

        return vals, order_backup, mol_coords_backup

    def _irmsd_core(self, pos1, pos2, atomic_numbers1, atomic_numbers2):
        """
        Core iRMSD calculation - exactly follows Fortran min_rmsd subroutine.

        Fortran flow:
        1. Get/compute ranks (we use canonical ranks from ref)
        2. If ranks differ, sort mol to match ref (skipped since we use same ranks)
        3. Count ngroup for each rank
        4. Call axis() on mol ONLY to align to principal axes
        5. Call min_rmsd_rotcheck_permute for step 1 (no z-mirror)
        6. If stereocheck, mirror z, call axis() again, call min_rmsd_rotcheck_permute for step 2
        7. Select orientation with minimum LSAP cost (minloc of tmprmsd_sym)
        8. Apply permutation to reorder atoms
        9. Compute final RMSD with rmsd() function
        """
        symbols2 = [self._pt.to_symbol(int(z)) for z in atomic_numbers2]
        masses2 = np.array([self._pt.to_atomic_mass(s) for s in symbols2])

        # Step 1: Compute canonical ranks from ref molecule only
        ranks = self._compute_canonical_ranks(atomic_numbers1, pos1)

        # Step 2: Skipped - we use same ranks for both molecules

        # Step 3: Count ngroup for each rank
        nranks = int(np.max(ranks))
        ngroup = np.zeros(nranks + 1, dtype=np.int32)
        for rank in ranks:
            if rank > 0:
                ngroup[rank] += 1

        # Working copy of mol coordinates (ref stays unchanged!)
        mol_xyz = pos2.copy()

        # Step 4: Initial alignment of mol to principal axes
        mol_aligned, rotconst, _ = self._align_to_principal_axes(
            mol_xyz, masses2
        )
        _, uniquenesscase = self._check_unique_axes(rotconst)

        # Initialize storage for 32 orientations (16 per step × 2 steps)
        INF = np.inf
        tmprmsd_sym = np.full(32, INF)
        order_backup = [None] * 32
        mol_coords_backup = [None] * 32

        # Step 5: Run min_rmsd_rotcheck_permute for step 1 (no z-mirror)
        vals1, orders1, mols1 = self._min_rmsd_rotcheck_permute(
            pos1, mol_aligned.copy(), ranks, ngroup, nranks, 1, uniquenesscase
        )
        tmprmsd_sym[0:16] = vals1
        for i in range(16):
            order_backup[i] = orders1[i]
            mol_coords_backup[i] = mols1[i]

        # Step 6: If stereocheck, mirror z and re-run
        if self._check_stereo_enabled:
            mol_inverted = mol_aligned.copy()
            mol_inverted[:, 2] = -mol_inverted[:, 2]  # Mirror z
            mol_inverted, _, _ = self._align_to_principal_axes(
                mol_inverted, masses2
            )

            vals2, orders2, mols2 = self._min_rmsd_rotcheck_permute(
                pos1,
                mol_inverted.copy(),
                ranks,
                ngroup,
                nranks,
                2,
                uniquenesscase,
            )
            tmprmsd_sym[16:32] = vals2
            for i in range(16):
                order_backup[16 + i] = orders2[i]
                mol_coords_backup[16 + i] = mols2[i]

        # Step 7: Select orientation with minimum LSAP cost
        best_idx = np.argmin(tmprmsd_sym)
        best_perm = order_backup[best_idx]
        best_mol_coords = mol_coords_backup[best_idx]

        # Step 8: Apply permutation to reorder atoms
        mol_permuted = best_mol_coords[best_perm]

        # Step 9: Compute final RMSD
        final_rmsd = self._rmsd_quaternion(pos1, mol_permuted)

        return final_rmsd

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
