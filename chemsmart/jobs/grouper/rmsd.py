"""
RMSD-based molecular grouping algorithms.

This module contains all RMSD-based grouper implementations:
- RMSDGrouper: Abstract base class for RMSD groupers
- BasicRMSDGrouper: Standard Kabsch RMSD
- HungarianRMSDGrouper: Hungarian algorithm for atom assignment
- SpyRMSDGrouper: Using graph isomorphism for symmetry correction
- IRMSDGrouper: Invariant RMSD with APSP algorithm
- PymolRMSDGrouper: PyMOL-based alignment
"""

import importlib
import logging
import multiprocessing
import os
from abc import abstractmethod
from typing import Iterable, List, Optional, Tuple

import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist

from chemsmart.io.molecules.structure import Molecule
from chemsmart.utils.utils import find_irmsd_command, kabsch_align

from .base import MatrixGrouper

logger = logging.getLogger(__name__)


class RMSDGrouper(MatrixGrouper):
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

    supports_multiprocessing = True

    def __init__(
        self,
        molecules: Iterable[Molecule],
        threshold=None,
        num_groups=None,
        num_procs: int = 1,
        align_molecules: bool = True,
        ignore_hydrogens: bool = False,
        label: str = None,
        conformer_ids: List[str] = None,
        matrix_format: str = "xlsx",
        energy_type: str = "E",
        **kwargs,
    ):
        """
        Initialize RMSD-based molecular grouper.

        Args:
            molecules (Iterable[Molecule]): Collection of molecules to group.
            threshold (float): RMSD threshold for grouping. Defaults to 0.5.
                Ignored if num_groups is specified. In threshold mode,
                clusters are formed by hierarchical complete-linkage
                clustering with linkage distance <= threshold.
            num_groups (int): Requested number of groups. When specified,
                complete-linkage merges are applied level by level. If the
                next full complete-linkage distance level would reduce the
                number of groups below the requested value, the previous
                level is retained.
            num_procs (int): Number of processes for parallel computation.
            align_molecules (bool): Whether to align molecules using Kabsch
                algorithm before RMSD calculation. Defaults to True.
            ignore_hydrogens (bool): Whether to exclude hydrogen atoms from
                RMSD calculation. Defaults to False.
            label (str): Label/name for output files. Defaults to None.
            conformer_ids (list[str]): Custom IDs for each molecule (e.g., ['c1', 'c2']).
            matrix_format (str): Output format ('xlsx', 'csv', 'txt'). Defaults to 'xlsx'.
            energy_type (str): Energy type for output files. Defaults to 'E'.

        Note:
            Uses standard hierarchical agglomerative complete-linkage
            clustering based on the full pairwise RMSD matrix.
        """
        super().__init__(
            molecules,
            threshold=threshold,
            num_groups=num_groups,
            num_procs=num_procs,
            label=label,
            conformer_ids=conformer_ids,
            matrix_format=matrix_format,
            energy_type=energy_type,
            **kwargs,
        )
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

    def calculate_rmsd_pair(self, i: int, j: int) -> float:
        """Public API: calculate RMSD for one molecule pair by index.

        This wrapper keeps ``_calculate_rmsd`` as an internal subclass contract
        while exposing a stable external method with shared validation and
        consistent error messages.
        """
        i_idx, j_idx = self._validate_pair_indices(i, j)
        if i_idx == j_idx:
            logger.debug(
                f"[{self.__class__.__name__}] RMSD requested for identical indices "
                f"({i_idx}, {j_idx}); returning 0.0"
            )
            return 0.0

        return float(self._calculate_rmsd((i_idx, j_idx)))

    def _calculate_pair_payload(self, idx_pair: Tuple[int, int]):
        """Worker payload for one pair; subclasses can return extra metadata."""
        return self._calculate_rmsd(idx_pair)

    def _calculate_pair_chunk(
        self, pairs: List[Tuple[int, int]]
    ) -> List[object]:
        """Worker payload for a chunk of pairs, preserving per-pair payload types."""
        return [self._calculate_pair_payload(pair) for pair in pairs]

    def _consume_pair_payload(
        self, idx_pair: Tuple[int, int], payload
    ) -> float:
        """Convert worker payload to RMSD and merge any parent-side metadata."""
        return float(payload)

    def _calculate_pairwise_rmsd_values(
        self, indices: List[Tuple[int, int]]
    ) -> List[float]:
        """Run pairwise RMSD jobs (serial or multiprocessing) with shared progress logging."""
        total_pairs = len(indices)
        rmsd_values: List[float] = []
        next_progress = 10

        if self.num_procs > 1 and total_pairs > 0:
            chunk_size = max(1, total_pairs // (self.num_procs * 8))
            pair_chunks = [
                indices[start : start + chunk_size]
                for start in range(0, total_pairs, chunk_size)
            ]

            completed = 0
            with multiprocessing.Pool(processes=self.num_procs) as pool:
                for pair_chunk, payload_chunk in zip(
                    pair_chunks,
                    pool.imap(
                        self._calculate_pair_chunk,
                        pair_chunks,
                        chunksize=1,
                    ),
                ):
                    for pair, payload in zip(pair_chunk, payload_chunk):
                        rmsd_values.append(
                            self._consume_pair_payload(pair, payload)
                        )

                    completed += len(pair_chunk)
                    next_progress = self._log_progress(
                        completed, total_pairs, next_progress
                    )
        else:
            for idx, pair in enumerate(indices):
                payload = self._calculate_pair_payload(pair)
                rmsd_values.append(self._consume_pair_payload(pair, payload))
                next_progress = self._log_progress(
                    idx + 1, total_pairs, next_progress
                )

        return rmsd_values

    def group(self) -> Tuple[List[List[Molecule]], List[List[int]]]:
        """
        Group molecules by geometric similarity using hierarchical
        complete-linkage RMSD clustering.

        Computes pairwise RMSD values between all molecules and groups
        them either by a linkage-distance threshold or by a requested
        number of groups. In num_groups mode, if the next full
        complete-linkage distance level would reduce the number of groups
        below the requested value, the previous level is retained.
        Automatically saves the full RMSD matrix to the group_result
        folder.

        Returns:
            Tuple[List[List[Molecule]], List[List[int]]]: Tuple containing:
                - List of molecule groups (each group is a list of molecules)
                - List of index groups (corresponding indices for each group)
        """
        import time

        # Record start time for grouping process
        grouping_start_time = time.time()
        self._auto_threshold = None

        n = len(self.molecules)
        indices = [(i, j) for i in range(n) for j in range(i + 1, n)]
        total_pairs = len(indices)

        logger.info(
            f"[{self.__class__.__name__}] Starting calculation for {n} molecules ({total_pairs} pairs)"
        )

        rmsd_values = self._calculate_pairwise_rmsd_values(indices)

        # Build full RMSD matrix for output
        rmsd_matrix = np.zeros((n, n))
        for (i, j), rmsd in zip(indices, rmsd_values):
            rmsd_matrix[i, j] = rmsd_matrix[j, i] = rmsd

        # Choose grouping strategy based on parameters (do this first to set _auto_threshold)
        if self.num_groups is not None:
            groups, index_groups = self.group_by_num_groups(rmsd_matrix)
        else:
            groups, index_groups = self.group_by_threshold(rmsd_matrix)

        # Calculate total grouping time
        grouping_end_time = time.time()
        grouping_time = grouping_end_time - grouping_start_time

        # Cache the results BEFORE saving (so Groups sheet can be populated)
        self._cached_groups = groups
        self._cached_group_indices = index_groups

        # Save full matrix using unified record entrypoint
        self.record(rmsd_matrix=rmsd_matrix, grouping_time=grouping_time)

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

        rmsd_values = self._calculate_pairwise_rmsd_values(indices)

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
            self.record(rmsd_matrix=rmsd_matrix)

        return rmsd_matrix

    def _record_results(
        self,
        rmsd_matrix: np.ndarray,
        grouping_time: float = None,
    ):
        """Strategy-specific result writer for RMSD grouping."""
        n = rmsd_matrix.shape[0]

        # Use ResultsRecorder to save
        recorder = self._get_results_recorder()
        labels = recorder.get_labels(n)

        # Build header info
        header_info = [
            ("", f"Full RMSD Matrix ({n}x{n}) - {self.__class__.__name__}"),
        ]

        if self.num_groups is not None:
            header_info.append(("Requested Groups (-N)", self.num_groups))
            header_info.append(
                (
                    "Actual Groups",
                    (
                        len(self._cached_group_indices)
                        if self._cached_group_indices is not None
                        else 0
                    ),
                )
            )
            if self._auto_threshold is not None:
                header_info.append(
                    (
                        "Auto-determined Threshold",
                        f"{self._auto_threshold:.7f} Å",
                    )
                )
        else:
            header_info.append(("Threshold", f"{self.threshold} Å"))

        header_info.append(("Align Molecules", self.align_molecules))
        header_info.append(("Energy Type", self.energy_type))
        header_info.append(("Ignore Hydrogens", self.ignore_hydrogens))

        # IRMSDGrouper specific parameters
        if isinstance(self, IRMSDGrouper):
            header_info.append(("Inversion", self._actual_inversion))

        # SpyRMSDGrouper specific parameters
        if isinstance(self, SpyRMSDGrouper):
            header_info.append(("Cache", self.cache))

        header_info.append(("Num Procs", self.num_procs))

        if grouping_time is not None:
            header_info.append(
                ("Grouping Time", f"{grouping_time:.2f} seconds")
            )

        self._append_input_usage_header(header_info)
        self._append_thermo_header(header_info)

        # Build sheets data
        sheets_data = {}
        index_groups = self._cached_group_indices
        if index_groups is not None:
            sheets_data["Groups"] = recorder.build_groups_dataframe(
                index_groups, n
            )

        # Determine suffix
        if self.num_groups is not None:
            suffix = f"N{self.num_groups}"
        else:
            suffix = f"T{self.threshold}"

        recorder.record_results(
            grouper_name=self.__class__.__name__,
            header_info=header_info,
            sheets_data=sheets_data,
            matrix_data=("RMSD_Matrix", rmsd_matrix, labels),
            suffix=suffix,
            startrow=len(header_info) + 2,
        )


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
        conformer_ids: List[str] = None,
        matrix_format: str = "xlsx",
        energy_type: str = "E",
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
            conformer_ids=conformer_ids,
            matrix_format=matrix_format,
            energy_type=energy_type,
            **kwargs,
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
        matrix_format: str = "xlsx",
        energy_type: str = "E",
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
            matrix_format: Output format ('xlsx', 'csv', 'txt').
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
            matrix_format=matrix_format,
            energy_type=energy_type,
            **kwargs,
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

    def _calculate_rmsd_with_isomorphism(
        self, idx_pair
    ) -> Tuple[float, Optional[Tuple[List[int], List[int]]]]:
        """Calculate symmetry-corrected RMSD and return best isomorphism."""
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
            return np.inf, None

        if not np.array_equal(np.sort(anum_ref), np.sort(anum_other)):
            return np.inf, None

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
            else:
                rmsd_value = result
                best_iso = None

            return float(rmsd_value), best_iso

        except Exception as e:
            logger.warning(
                f"spyrmsd calculation failed for pair ({i}, {j}): {e}"
            )
            return np.inf, None

    def _calculate_rmsd(self, idx_pair):
        """Calculate symmetry-corrected RMSD and store best isomorphism."""
        rmsd_value, best_iso = self._calculate_rmsd_with_isomorphism(idx_pair)
        self.best_isomorphisms[idx_pair] = best_iso
        return rmsd_value

    def _calculate_pair_payload(self, idx_pair):
        """Return SpyRMSD payload used by multiprocessing workers."""
        return self._calculate_rmsd_with_isomorphism(idx_pair)

    def _consume_pair_payload(self, idx_pair, payload) -> float:
        """Merge best isomorphism metadata in parent process and return RMSD."""
        rmsd_value, best_iso = payload
        self.best_isomorphisms[idx_pair] = best_iso
        return float(rmsd_value)

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
    Invariant RMSD (iRMSD) Grouper using the irmsd package.

    When the irmsd Python API is importable in the current environment, this
    grouper calls its compiled backend directly and avoids per-pair subprocess
    and temporary-file overhead. The existing external CLI implementation is
    retained as a compatibility fallback.

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
        matrix_format: str = "xlsx",
        energy_type: str = "E",
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
            matrix_format: Output format ('xlsx', 'csv', 'txt').
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
            matrix_format=matrix_format,
            energy_type=energy_type,
            **kwargs,
        )
        self.inversion = (
            inversion.lower() if isinstance(inversion, str) else "auto"
        )
        self._actual_inversion = None
        self._use_direct_api = False
        self._irmsd_cmd = None

        # Prefer the in-process Python/Fortran API. This avoids launching one
        # external process and writing one temporary XYZ file for every pair.
        try:
            irmsd_module = importlib.import_module("irmsd.api.irmsd_exposed")
            if not callable(irmsd_module.get_irmsd):
                raise AttributeError("irmsd API does not provide get_irmsd")
            self._use_direct_api = True
            logger.info("Using direct irmsd Python/Fortran API")
        except (ImportError, OSError, AttributeError) as exc:
            logger.info(
                f"Direct irmsd Python API is unavailable ({exc}); "
                "falling back to the external irmsd CLI"
            )
            self._irmsd_cmd = find_irmsd_command()

        if not self._use_direct_api and not self._irmsd_cmd:
            raise RuntimeError(
                "irmsd is not available through either the Python API or the "
                "external command. IRMSDGrouper requires the irmsd package.\n\n"
                "For the fastest implementation, install irmsd in the same "
                "environment as CHEMSMART so that "
                "irmsd.api.irmsd_exposed.get_irmsd can be imported.\n\n"
                "Alternatively, configure the existing CLI fallback with:\n"
                "  export IRMSD_CONDA_ENV=irmsd_env\n"
                "or:\n"
                "  export IRMSD_PATH=/path/to/irmsd"
            )

    @staticmethod
    def _inversion_code(inversion: str) -> int:
        """Map the public inversion mode to the irmsd low-level API code."""
        return {"auto": 0, "on": 1, "off": 2}[inversion]

    def _prepare_irmsd_arrays(
        self, mol: Molecule
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Return atomic numbers and coordinates for direct iRMSD calls."""
        from ase.data import atomic_numbers

        symbols = list(mol.chemical_symbols)
        positions = np.asarray(mol.positions, dtype=np.float64)

        if self.ignore_hydrogens:
            keep = [idx for idx, symbol in enumerate(symbols) if symbol != "H"]
            symbols = [symbols[idx] for idx in keep]
            positions = positions[keep]

        numbers = np.asarray(
            [atomic_numbers[symbol] for symbol in symbols], dtype=np.int32
        )
        return numbers, positions

    def _calculate_rmsd_direct(self, mol_idx_pair: Tuple[int, int]) -> float:
        """Calculate one iRMSD pair through the in-process compiled API."""
        irmsd_module = importlib.import_module("irmsd.api.irmsd_exposed")
        get_irmsd = irmsd_module.get_irmsd

        i, j = mol_idx_pair
        mol1, mol2 = self.molecules[i], self.molecules[j]

        symbols1 = list(mol1.chemical_symbols)
        symbols2 = list(mol2.chemical_symbols)
        if self.ignore_hydrogens:
            symbols1 = [symbol for symbol in symbols1 if symbol != "H"]
            symbols2 = [symbol for symbol in symbols2 if symbol != "H"]

        if sorted(symbols1) != sorted(symbols2):
            logger.warning(f"Molecules {i} and {j} have different atom types")
            return np.inf

        try:
            numbers1, positions1 = self._prepare_irmsd_arrays(mol1)
            numbers2, positions2 = self._prepare_irmsd_arrays(mol2)
            rmsd_value, _, _, _, _ = get_irmsd(
                numbers1,
                positions1,
                numbers2,
                positions2,
                iinversion=self._inversion_code(self.inversion),
            )

            # The low-level API receives the requested inversion mode directly;
            # unlike the CLI it does not return a separate text field to parse.
            if self._actual_inversion is None:
                self._actual_inversion = self.inversion

            rmsd_value = float(rmsd_value)
            return rmsd_value if np.isfinite(rmsd_value) else np.inf
        except Exception as exc:
            logger.warning(
                f"irmsd API calculation failed for pair ({i}, {j}): {exc}"
            )
            return np.inf

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
                    # Failed to parse RMSD value from output line, continue searching
                    pass
            elif parse_inversion and "Inversion check:" in line:
                # Parse "Inversion check: on/off/auto"
                parts = line.split(":")
                if len(parts) >= 2:
                    actual_inversion = parts[1].strip()

        return rmsd_value, actual_inversion

    def _calculate_rmsd_cli(self, mol_idx_pair: Tuple[int, int]) -> float:
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
            need_inversion = self._actual_inversion is None
            rmsd_value, actual_inversion = self._parse_irmsd_output(
                result.stdout, parse_inversion=need_inversion
            )

            # Store the actual inversion setting (only need to do this once)
            if actual_inversion:
                self._actual_inversion = actual_inversion

            return float(rmsd_value) if rmsd_value is not None else np.inf

        except subprocess.TimeoutExpired:
            logger.warning(f"irmsd timed out for pair ({i}, {j})")
            return np.inf
        except Exception as e:
            logger.warning(f"irmsd error for pair ({i}, {j}): {e}")
            return np.inf
        finally:
            # Clean up temporary file
            try:
                os.unlink(tmp_path)
            except OSError:
                # Ignore cleanup errors
                pass

    def _calculate_rmsd(self, mol_idx_pair: Tuple[int, int]) -> float:
        """Calculate iRMSD using the direct API when available, else the CLI."""
        if self._use_direct_api:
            return self._calculate_rmsd_direct(mol_idx_pair)
        return self._calculate_rmsd_cli(mol_idx_pair)

    def _calculate_pair_payload(self, idx_pair):
        """Return iRMSD value together with resolved inversion metadata."""
        rmsd_value = self._calculate_rmsd(idx_pair)
        return rmsd_value, self._actual_inversion

    def _consume_pair_payload(self, idx_pair, payload) -> float:
        """Merge iRMSD inversion metadata in the parent process."""
        rmsd_value, actual_inversion = payload
        if self._actual_inversion is None and actual_inversion is not None:
            self._actual_inversion = actual_inversion
        return float(rmsd_value)


class PymolRMSDGrouper(RMSDGrouper):
    """Group molecules using PyMOL's align command for RMSD calculation."""

    supports_multiprocessing = False

    def __init__(
        self,
        molecules: Iterable[Molecule],
        threshold=None,
        num_groups=None,
        num_procs: int = 1,
        ignore_hydrogens: bool = False,
        label: str = None,
        conformer_ids: List[str] = None,
        matrix_format: str = "xlsx",
        energy_type: str = "E",
        **kwargs,
    ):
        super().__init__(
            molecules=molecules,
            threshold=threshold,
            num_groups=num_groups,
            num_procs=num_procs,
            align_molecules=True,
            ignore_hydrogens=ignore_hydrogens,
            label=label,
            conformer_ids=conformer_ids,
            matrix_format=matrix_format,
            energy_type=energy_type,
            **kwargs,
        )
        self._temp_dir = None
        self._xyz_files = []
        self._mol_names = []
        self._alignment_cache = {}
        self.cmd = None  # Will be set in _init_pymol()
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

        # Clean up PyMOL session
        try:
            if self.cmd is not None:
                self.cmd.quit()
        except Exception:
            logger.debug(
                "Error while quitting PyMOL cmd in __del__", exc_info=True
            )

        # Clean up temporary directory
        try:
            if self._temp_dir is not None:
                shutil.rmtree(self._temp_dir, ignore_errors=True)
        except Exception:
            logger.debug(
                "Error while removing temporary directory in __del__",
                exc_info=True,
            )

    def __repr__(self):
        return (
            f"{self.__class__.__name__}(num_groups={self.num_groups})"
            if self.num_groups is not None
            else f"{self.__class__.__name__}(threshold={self.threshold})"
        )


__all__ = [
    "RMSDGrouper",
    "BasicRMSDGrouper",
    "HungarianRMSDGrouper",
    "SpyRMSDGrouper",
    "IRMSDGrouper",
    "PymolRMSDGrouper",
]
