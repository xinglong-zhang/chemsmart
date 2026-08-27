"""
Energy-based molecular grouping algorithm.

This module contains the EnergyGrouper implementation which groups molecules
based on energy differences using a threshold or target number of groups.
"""

import logging
import time
from typing import Iterable, List, Tuple

import numpy as np

from chemsmart.io.molecules.structure import Molecule

from .base import MatrixGrouper

logger = logging.getLogger(__name__)

# Conversion factor: 1 Hartree = 627.509474 kcal/mol
HARTREE_TO_KCAL = 627.509474
KCAL_TO_HARTREE = 1.0 / HARTREE_TO_KCAL


class EnergyGrouper(MatrixGrouper):
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
        matrix_format: str = "xlsx",
        energy_type: str = "E",
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
            matrix_format (str): Output format ('xlsx', 'csv', 'txt'). Defaults to 'xlsx'.

        Note:
            Uses hierarchical complete-linkage clustering based on the full
            pairwise energy-difference matrix.
        """
        if threshold is None and num_groups is None:
            threshold = 1.0  # Default: 1 kcal/mol

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

        # Store threshold in kcal/mol (user-friendly unit)
        self.threshold = threshold
        # Convert to Hartree for internal calculations (only if threshold is set)
        self.threshold_hartree = (
            threshold * KCAL_TO_HARTREE if threshold is not None else None
        )
        self._auto_threshold = (
            None  # Will store auto-determined threshold in kcal/mol
        )

        # Validate that all molecules have energy information
        self._validate_energies()

    def _validate_energies(self) -> None:
        """
        Validate that all molecules have usable energy information.

        Raises:
            ValueError: If any molecule is missing or has invalid energy information.
        """
        missing_energy = []
        for i, mol in enumerate(self.molecules):
            energy = mol.energy
            if energy is None:
                missing_energy.append(i + 1)
                continue

            try:
                energy_value = float(energy)
            except (TypeError, ValueError):
                missing_energy.append(i + 1)
                continue

            # Reject NaN/inf so downstream matrix/grouping math is always valid.
            if not np.isfinite(energy_value):
                missing_energy.append(i + 1)

        if missing_energy:
            if len(missing_energy) <= 5:
                raise ValueError(
                    f"Molecules at indices {missing_energy} are missing energy information or have invalid values. "
                    "Energy grouping requires all molecules to have finite energy values."
                )
            else:
                raise ValueError(
                    f"Found {len(missing_energy)} molecules with missing energy information or invalid values. "
                    "Energy grouping requires all molecules to have finite energy values."
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
        Group molecules by energy similarity using hierarchical
        complete-linkage clustering.

        Computes pairwise energy differences between all molecules and groups
        those within the specified threshold, or automatically selects the
        deepest complete-linkage distance level that does not reduce the
        number of groups below the requested value.

        Returns:
            Tuple[List[List[Molecule]], List[List[int]]]: Tuple containing:
                - List of molecule groups (each group is a list of molecules)
                - List of index groups (corresponding indices for each group)
        """
        grouping_start_time = time.time()
        self._auto_threshold = None

        n = len(self.molecules)
        indices = [(i, j) for i in range(n) for j in range(i + 1, n)]
        total_pairs = len(indices)

        logger.info(
            f"[{self.__class__.__name__}] Starting calculation for {n} molecules ({total_pairs} pairs)"
        )

        # Calculate energy differences (both relative and absolute)
        energy_diff_relative = []  # For output matrix (with sign)
        energy_diff_absolute = []  # For threshold comparison
        next_progress = 10
        for idx, (i, j) in enumerate(indices):
            rel_diff, abs_diff = self._calculate_energy_diff((i, j))
            energy_diff_relative.append(rel_diff)
            energy_diff_absolute.append(abs_diff)
            next_progress = self._log_progress(
                idx + 1, total_pairs, next_progress
            )

        # Build full energy difference matrix for output (with sign, relative to smaller index)
        # matrix[i,j] = E_j - E_i (positive means j has higher energy)
        # matrix[j,i] = E_i - E_j (negative of the above)
        energy_matrix = np.zeros((n, n))
        energy_distance_matrix = np.zeros((n, n))
        for (i, j), rel_diff in zip(indices, energy_diff_relative):
            energy_matrix[i, j] = rel_diff  # E_j - E_i
            energy_matrix[j, i] = -rel_diff  # E_i - E_j
        for (i, j), abs_diff in zip(indices, energy_diff_absolute):
            energy_distance_matrix[i, j] = abs_diff
            energy_distance_matrix[j, i] = abs_diff

        # Choose grouping strategy based on parameters
        if self.num_groups is not None:
            groups, index_groups = self._group_by_num_groups(
                energy_distance_matrix
            )
        else:
            groups, index_groups = self._group_by_threshold(
                energy_distance_matrix
            )

        # Calculate total grouping time
        grouping_end_time = time.time()
        grouping_time = grouping_end_time - grouping_start_time

        # Cache the results BEFORE saving (so Groups sheet can be populated)
        self._cached_groups = groups
        self._cached_group_indices = index_groups

        # Save energy difference matrix through unified record entrypoint
        self.record(energy_matrix=energy_matrix, grouping_time=grouping_time)

        return groups, index_groups

    def _group_by_threshold(
        self, distance_matrix: np.ndarray
    ) -> Tuple[List[List[Molecule]], List[List[int]]]:
        original_threshold = self.threshold
        result: Tuple[List[List[Molecule]], List[List[int]]] | None = None
        try:
            self.threshold = self.threshold_hartree
            result = super()._group_by_threshold(distance_matrix)
        finally:
            self.threshold = original_threshold

        assert result is not None, "Energy grouping failed to produce groups."
        return result

    def _build_complete_linkage_tree(
        self, distance_matrix: np.ndarray
    ) -> np.ndarray:
        try:
            return super()._build_complete_linkage_tree(distance_matrix)
        except ValueError as exc:
            if "non-finite" in str(exc):
                raise ValueError(
                    "Pairwise energy-difference matrix contains non-finite values. "
                    "Energy grouping requires all molecules to have finite energy values."
                ) from exc
            raise

    def _group_by_num_groups(
        self, distance_matrix: np.ndarray
    ) -> Tuple[List[List[Molecule]], List[List[int]]]:
        original_threshold = self.threshold
        result: Tuple[List[List[Molecule]], List[List[int]]] | None = None
        try:
            self.threshold = self.threshold_hartree
            result = super()._group_by_num_groups(distance_matrix)
        finally:
            self.threshold = original_threshold

        assert (
            result is not None
        ), "Energy grouping failed to produce groups in num_groups mode."

        groups, index_groups = result

        if self._auto_threshold is not None:
            self._auto_threshold *= HARTREE_TO_KCAL

        return groups, index_groups

    def _record_results(
        self,
        energy_matrix: np.ndarray,
        grouping_time: float = None,
    ):
        """Strategy-specific result writer for energy grouping."""
        n = energy_matrix.shape[0]

        # Use ResultsRecorder to save
        recorder = self._get_results_recorder()
        labels = recorder.get_labels(n)

        # Convert energy matrix from Hartree to kcal/mol for output
        energy_matrix_kcal = energy_matrix * HARTREE_TO_KCAL

        # Build header info
        header_info = [
            (
                "",
                f"Relative Energy Difference Matrix ({n}x{n}) - {self.__class__.__name__}",
            ),
            (
                "",
                "Values are relative energy differences in kcal/mol: matrix[i,j] = E_j - E_i",
            ),
            (
                "",
                "Positive value means column molecule has higher energy than row molecule",
            ),
            ("Conversion", f"1 Hartree = {HARTREE_TO_KCAL} kcal/mol"),
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
                        f"{self._auto_threshold:.4f} kcal/mol ({self._auto_threshold * KCAL_TO_HARTREE:.10f} Hartree)",
                    )
                )
        else:
            header_info.append(
                (
                    "Threshold",
                    f"{self.threshold:.4f} kcal/mol ({self.threshold_hartree:.10f} Hartree)",
                )
            )
        header_info.append(("Energy Type", self.energy_type))
        header_info.append(("Num Procs", self.num_procs))

        if grouping_time is not None:
            header_info.append(
                ("Grouping Time", f"{grouping_time:.2f} seconds")
            )

        self._append_input_usage_header(header_info)
        self._append_thermo_header(header_info)

        # Build sheets data using recorder's method
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
            matrix_data=("Energy_Matrix", energy_matrix_kcal, labels),
            suffix=suffix,
            startrow=len(header_info) + 2,
            float_format="%.4f",
        )

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
