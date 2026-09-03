"""
Base classes for molecular structure grouping algorithms.

This module contains the core abstract base class and utilities:
- MoleculeGrouper: Abstract base class for all grouper implementations
- MatrixGrouper: Base class for pairwise-matrix-based grouping algorithms
- ResultsRecorder: Unified results output with multiple format support
- StructureGrouperConfig: Configuration container for structure matching
"""

import logging
import os
from abc import ABC, abstractmethod
from typing import Any, Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

from chemsmart.io.molecules.structure import Molecule

logger = logging.getLogger(__name__)


class ResultsRecorder:
    """
    Unified results recorder for grouper output with multiple format support.

    Handles all output file generation (xlsx/csv/txt), including:
    - Matrix data output
    - Groups DataFrame generation
    - Header information formatting
    - Column width auto-adjustment

    Attributes:
        output_dir (str): Directory for output files.
        label (str): Label prefix for output filenames.
        matrix_format (str): Output format ('xlsx', 'csv', 'txt').
        conformer_ids (list[str]): Custom IDs for labeling molecules.
    """

    SUPPORTED_FORMATS = {"xlsx", "csv", "txt"}

    def __init__(
        self,
        output_dir: str,
        label: str = None,
        matrix_format: str = "xlsx",
        conformer_ids: List[str] = None,
    ):
        """
        Initialize the results recorder.

        Args:
            output_dir (str): Directory for output files.
            label (str): Label prefix for output filenames.
            matrix_format (str): Output format. Defaults to 'xlsx'.
                Supported: 'xlsx', 'csv', 'txt'.
            conformer_ids (list[str]): Custom IDs for labeling molecules.
        """
        self.output_dir = output_dir
        self.label = label
        self.matrix_format = matrix_format.lower().lstrip(".")
        self.conformer_ids = conformer_ids

        if self.matrix_format not in self.SUPPORTED_FORMATS:
            raise ValueError(
                f"Unsupported output format: '{matrix_format}'. "
                f"Supported formats: {self.SUPPORTED_FORMATS}"
            )

        os.makedirs(output_dir, exist_ok=True)

    def get_labels(self, n: int) -> List[str]:
        """
        Get labels for matrix rows/columns.

        Args:
            n (int): Number of molecules.

        Returns:
            List[str]: Labels from conformer_ids if available, otherwise "1", "2", ...
        """
        if self.conformer_ids is not None and len(self.conformer_ids) == n:
            return self.conformer_ids
        else:
            return [str(i + 1) for i in range(n)]

    def build_groups_dataframe(
        self,
        index_groups: List[List[int]],
        n_molecules: int,
        extra_columns: Dict[str, List[Any]] = None,
    ) -> pd.DataFrame:
        """
        Build a DataFrame for the Groups sheet.

        Args:
            index_groups: List of index groups from grouping.
            n_molecules: Total number of molecules (for label generation).
            extra_columns: Optional dict of additional columns {column_name: values_list}.

        Returns:
            pd.DataFrame with Group, Members columns (and any extra columns).
        """
        groups_data = []
        labels = self.get_labels(n_molecules)

        for i, indices in enumerate(index_groups):
            member_labels = [labels[idx] for idx in indices]
            row = {
                "Group": i + 1,
                "Members": ", ".join(member_labels),
            }
            # Add extra columns if provided
            if extra_columns:
                for col_name, values in extra_columns.items():
                    if i < len(values):
                        row[col_name] = values[i]
            groups_data.append(row)

        return pd.DataFrame(groups_data)

    def get_filename(self, grouper_name: str, suffix: str = None) -> str:
        """
        Generate output filename based on label, grouper name, and format.

        Args:
            grouper_name (str): Name of the grouper class.
            suffix (str): Optional suffix (e.g., 'T0.5' for threshold).

        Returns:
            str: Full path to output file.
        """
        label_prefix = f"{self.label}_" if self.label else ""
        suffix_part = f"_{suffix}" if suffix else ""
        filename = (
            f"{label_prefix}{grouper_name}{suffix_part}.{self.matrix_format}"
        )
        return os.path.join(self.output_dir, filename)

    def record_results(
        self,
        grouper_name: str,
        header_info: List[Tuple[str, Any]],
        sheets_data: Dict[str, pd.DataFrame],
        matrix_data: Optional[Tuple[str, np.ndarray, List[str]]] = None,
        suffix: str = None,
        startrow: int = None,
        float_format: str = "%.7f",
    ) -> str:
        """
        Record grouping results to file.

        Args:
            grouper_name (str): Name of the grouper class.
            header_info (list): List of (key, value) tuples for header information.
            sheets_data (dict): Dict of {sheet_name: DataFrame} for tabular data.
            matrix_data (tuple): Optional (sheet_name, matrix, labels) for matrix output.
            suffix (str): Optional suffix for filename.
            startrow (int): Row to start writing matrix data. If None, calculated
                from header_info length + 2.
            float_format (str): Format string for floating point numbers. Default "%.7f".

        Returns:
            str: Path to the output file.
        """
        filename = self.get_filename(grouper_name, suffix)

        if self.matrix_format == "xlsx":
            self._write_xlsx(
                filename,
                header_info,
                sheets_data,
                matrix_data,
                startrow,
                float_format,
            )
        elif self.matrix_format == "csv":
            self._write_csv(filename, header_info, sheets_data, matrix_data)
        elif self.matrix_format == "txt":
            self._write_txt(filename, header_info, sheets_data, matrix_data)
        else:
            raise ValueError(
                f"Unsupported output format: '{self.matrix_format}'"
            )

        logger.info(f"Results saved to {filename}")
        return filename

    def _write_xlsx(
        self,
        filename: str,
        header_info: List[Tuple[str, Any]],
        sheets_data: Dict[str, pd.DataFrame],
        matrix_data: Optional[Tuple[str, np.ndarray, List[str]]] = None,
        startrow: int = None,
        float_format: str = "%.7f",
    ):
        """Write results to Excel format with multiple sheets."""
        with pd.ExcelWriter(filename, engine="openpyxl") as writer:
            # Write matrix sheet first if provided
            if matrix_data is not None:
                sheet_name, matrix, labels = matrix_data
                # Replace inf with NaN for display
                matrix_display = np.where(np.isinf(matrix), np.nan, matrix)
                df = pd.DataFrame(matrix_display, index=labels, columns=labels)

                # Calculate startrow if not provided
                if startrow is None:
                    startrow = len(header_info) + 2

                df.to_excel(
                    writer,
                    sheet_name=sheet_name,
                    startrow=startrow,
                    float_format=float_format,
                )

                # Add header information
                worksheet = writer.sheets[sheet_name]
                for row, (key, value) in enumerate(header_info, start=1):
                    # Format: "Key: Value" or just "Value" if key is empty/title
                    # Note: value can be False (boolean), so check for None explicitly
                    if key and value is not None and value != "":
                        worksheet[f"A{row}"] = f"{key}: {value}"
                    elif value is not None and value != "":
                        worksheet[f"A{row}"] = str(value)
                    elif key:
                        worksheet[f"A{row}"] = str(key)
                    # Skip if both key and value are empty/None

                # Auto-adjust column widths for matrix sheet
                self._auto_adjust_columns(worksheet)

            # Write other sheets
            first_sheet_with_header = (
                matrix_data is None
            )  # First sheet gets header if no matrix
            for sheet_name, df in sheets_data.items():
                if matrix_data is not None and sheet_name == matrix_data[0]:
                    continue  # Skip if already written as matrix

                # First sheet (when no matrix) gets header info
                if first_sheet_with_header and sheet_name != "Groups":
                    # Calculate startrow if not provided
                    actual_startrow = (
                        startrow if startrow else len(header_info) + 2
                    )
                    df.to_excel(
                        writer,
                        sheet_name=sheet_name,
                        startrow=actual_startrow,
                        index=False,
                    )

                    # Add header information
                    worksheet = writer.sheets[sheet_name]
                    for row, (key, value) in enumerate(header_info, start=1):
                        # Note: value can be False (boolean), so check for None explicitly
                        if key and value is not None and value != "":
                            worksheet[f"A{row}"] = f"{key}: {value}"
                        elif value is not None and value != "":
                            worksheet[f"A{row}"] = str(value)
                        elif key:
                            worksheet[f"A{row}"] = str(key)
                        # Skip if both key and value are empty/None

                    first_sheet_with_header = (
                        False  # Only first non-Groups sheet
                    )
                else:
                    # Groups sheet and other sheets don't need header info
                    df.to_excel(writer, sheet_name=sheet_name, index=False)

                # Auto-adjust column widths
                self._auto_adjust_columns(writer.sheets[sheet_name])

    def _auto_adjust_columns(self, worksheet, max_width: int = 50):
        """Auto-adjust column widths for a worksheet."""
        for column_cells in worksheet.columns:
            max_length = 0
            column_letter = column_cells[0].column_letter
            for cell in column_cells:
                try:
                    cell_length = len(str(cell.value)) if cell.value else 0
                    max_length = max(max_length, cell_length)
                except (TypeError, AttributeError):
                    pass
            worksheet.column_dimensions[column_letter].width = min(
                max_length + 2, max_width
            )

    def _write_csv(
        self,
        filename: str,
        header_info: List[Tuple[str, Any]],
        sheets_data: Dict[str, pd.DataFrame],
        matrix_data: Optional[Tuple[str, np.ndarray, List[str]]] = None,
    ):
        """Write results to CSV format (one file per sheet)."""
        base_name = filename[:-4]  # Remove .csv extension

        # Write header info to main file
        with open(filename, "w") as f:
            for key, value in header_info:
                # Handle empty key (title lines)
                if key and value is not None and value != "":
                    f.write(f"# {key}: {value}\n")
                elif value is not None and value != "":
                    f.write(f"# {value}\n")
                elif key:
                    f.write(f"# {key}\n")
            f.write("\n")

        # Write matrix if provided
        if matrix_data is not None:
            sheet_name, matrix, labels = matrix_data
            matrix_display = np.where(np.isinf(matrix), np.nan, matrix)
            df = pd.DataFrame(matrix_display, index=labels, columns=labels)
            df.to_csv(filename, mode="a")

        # Write other sheets to separate files
        for sheet_name, df in sheets_data.items():
            if matrix_data is not None and sheet_name == matrix_data[0]:
                continue
            sheet_filename = f"{base_name}_{sheet_name}.csv"
            df.to_csv(sheet_filename, index=False)
            logger.info(f"Sheet '{sheet_name}' saved to {sheet_filename}")

    def _write_txt(
        self,
        filename: str,
        header_info: List[Tuple[str, Any]],
        sheets_data: Dict[str, pd.DataFrame],
        matrix_data: Optional[Tuple[str, np.ndarray, List[str]]] = None,
    ):
        """Write results to plain text format."""
        with open(filename, "w") as f:
            # Write header info
            f.write("=" * 60 + "\n")
            for key, value in header_info:
                # Handle empty key (title lines)
                if key and value is not None and value != "":
                    f.write(f"{key}: {value}\n")
                elif value is not None and value != "":
                    f.write(f"{value}\n")
                elif key:
                    f.write(f"{key}\n")
            f.write("=" * 60 + "\n\n")

            # Write matrix if provided
            if matrix_data is not None:
                sheet_name, matrix, labels = matrix_data
                f.write(f"--- {sheet_name} ---\n")
                # Write column headers
                f.write(" " * 12)
                for label in labels:
                    f.write(f"{label:>12}")
                f.write("\n")
                # Write rows
                for i, row_label in enumerate(labels):
                    f.write(f"{row_label:>12}")
                    for j in range(len(labels)):
                        val = matrix[i, j]
                        if np.isinf(val):
                            f.write(f"{'inf':>12}")
                        else:
                            f.write(f"{val:>12.7f}")
                    f.write("\n")
                f.write("\n")

            # Write other sheets
            for sheet_name, df in sheets_data.items():
                if matrix_data is not None and sheet_name == matrix_data[0]:
                    continue
                f.write(f"--- {sheet_name} ---\n")
                f.write(df.to_string(index=False))
                f.write("\n\n")


class StructureGrouperConfig:
    """
    Configuration container for StructureMatcher parameters.

    Stores tolerance parameters for structure matching algorithms.
    Default values are optimized for heterogeneous molecular systems
    and may need adjustment for specific molecular types.

    Attributes:
        ltol (float): Length tolerance for structure matching.
        stol (float): Site tolerance for atomic position matching.
        angle_tol (float): Angle tolerance in degrees for structure matching.
    """

    def __init__(self, ltol=0.1, stol=0.18, angle_tol=1):
        """
        Initialize structure grouper configuration.

        Args:
            ltol (float): Length tolerance. Defaults to 0.1.
            stol (float): Site tolerance. Defaults to 0.18.
            angle_tol (float): Angle tolerance in degrees. Defaults to 1.
        """
        self.ltol = ltol
        self.stol = stol
        self.angle_tol = angle_tol


class MoleculeGrouper(ABC):
    """
    Abstract base class for molecular structure grouping algorithms.

    Defines the common interface that all molecular grouping strategies
    must implement. Cannot be directly instantiated and designed to
    ensure consistent behavior across different grouping methods.

    Attributes:
        molecules (Iterable[Molecule]): Collection of molecules to group.
        num_procs (int): Number of processes for parallel computation.
        label (str): Label/name for this grouping task (used in output filenames).
        conformer_ids (list[str]): Optional custom IDs for each molecule (e.g., ['c1', 'c2']).
        matrix_format (str): Output format for results ('xlsx', 'csv', 'txt').
        output_dir (str): Base directory for output files.
    """

    # Concrete subclasses opt in when they are safe to run pair calculations
    # in multiple processes.
    supports_multiprocessing = False
    supports_matrix_representative_selection = False
    SUPPORTED_REPRESENTATIVE_STRATEGIES = {"lowest", "center", "top3"}

    def __init__(
        self,
        molecules: Iterable[Molecule],
        num_procs: int = 1,
        label: str = None,
        conformer_ids: List[str] = None,
        skipped_ids: List[str] = None,
        matrix_format: str = "xlsx",
        output_dir: str = None,
        energy_type: str = "E",
        thermo_parameters: str = None,
        representative_strategy: str = "lowest",
    ):
        """
        Initialize the molecular grouper.

        Args:
            molecules (Iterable[Molecule]): Collection of molecules to group.
            num_procs (int): Number of processes for parallel computation.
                Defaults to 1.
            label (str): Label/name for this grouping task. Used in output folder
                and file names. Defaults to None.
            conformer_ids (list[str]): Optional custom IDs for each molecule (e.g., ['c1', 'c2']).
                If provided, these are used as labels in matrix output instead of numeric indices.
            skipped_ids (list[str]): Optional list of IDs that were skipped before grouping.
            matrix_format (str): Output format for results ('xlsx', 'csv', 'txt').
                Defaults to 'xlsx'.
            output_dir (str): Base directory for output files. If None, uses current
                working directory. Defaults to None.
            energy_type (str): Energy type label propagated from CLI/job.
            thermo_parameters (str): Thermochemistry parameters from CLI.
        """
        self.molecules = list(molecules)
        requested_num_procs = int(max(1, num_procs))
        if requested_num_procs > 1 and not self.supports_multiprocessing:
            logger.warning(
                f"{self.__class__.__name__} does not support multiprocessing; "
                "using num_procs=1."
            )
            requested_num_procs = 1
        self.num_procs = requested_num_procs
        self.label = label
        self.conformer_ids = conformer_ids
        self.skipped_ids = skipped_ids if skipped_ids is not None else []
        self.matrix_format = matrix_format
        self.output_dir = output_dir
        self.energy_type = energy_type
        self.thermo_parameters = thermo_parameters
        self.representative_strategy = str(
            representative_strategy or "lowest"
        ).lower()

        # Cache for avoiding repeated grouping calculations
        self._cached_groups = None
        self._cached_group_indices = None
        self._matrix_skipped_ids = []
        self._matrix_skipped_indices = []

        self._validate_inputs()
        self._validate_representative_strategy()

    def _validate_representative_strategy(self) -> None:
        """Validate representative strategy and strategy compatibility."""
        if (
            self.representative_strategy
            not in self.SUPPORTED_REPRESENTATIVE_STRATEGIES
        ):
            raise ValueError(
                f"Unsupported representative strategy '{self.representative_strategy}'. "
                "Supported values: lowest, center, top3."
            )

        if (
            self.representative_strategy in {"center", "top3"}
            and not self.supports_matrix_representative_selection
        ):
            raise ValueError(
                f"Representative strategy '{self.representative_strategy}' requires a pairwise distance matrix. "
                f"Use '-r lowest' for {self.__class__.__name__}."
            )

    def _energy_sort_key(self, original_index: int) -> Tuple[int, float, int]:
        """Deterministic key: finite-energy first, then energy, then index."""
        energy = self.molecules[original_index].energy
        if energy is None:
            return (1, float("inf"), original_index)

        try:
            value = float(energy)
        except (TypeError, ValueError):
            return (1, float("inf"), original_index)

        if not np.isfinite(value):
            return (1, float("inf"), original_index)

        return (0, value, original_index)

    def _order_index_group_by_energy(
        self, index_group: List[int]
    ) -> List[int]:
        """Order one index group by energy with deterministic tie-breaking."""
        return sorted(index_group, key=self._energy_sort_key)

    def _rebuild_groups_from_index_groups(
        self, index_groups: List[List[int]]
    ) -> List[List[Molecule]]:
        """Rebuild molecule groups from index groups preserving order."""
        return [
            [self.molecules[idx] for idx in group] for group in index_groups
        ]

    def _apply_lowest_representative_ordering(
        self, index_groups: List[List[int]]
    ) -> Tuple[List[List[Molecule]], List[List[int]]]:
        """Apply lowest-energy representative ordering to each group."""
        ordered_index_groups = [
            self._order_index_group_by_energy(group) for group in index_groups
        ]
        return (
            self._rebuild_groups_from_index_groups(ordered_index_groups),
            ordered_index_groups,
        )

    def _get_output_dir(self) -> str:
        """Get the output directory path for this grouper."""
        if self.label:
            folder_name = f"{self.label}_group_result"
        else:
            folder_name = "group_result"

        if self.output_dir:
            return os.path.join(self.output_dir, folder_name)
        return folder_name

    def _get_results_recorder(self) -> ResultsRecorder:
        """Create a ResultsRecorder instance for this grouper."""
        return ResultsRecorder(
            output_dir=self._get_output_dir(),
            label=self.label,
            matrix_format=self.matrix_format,
            conformer_ids=self.conformer_ids,
        )

    def _validate_inputs(self) -> None:
        """
        Validate input molecules for grouping.

        Ensures that the input is an iterable collection and all items
        are valid Molecule instances.

        Raises:
            TypeError: If molecules is not iterable or contains non-Molecule items.
        """
        if not isinstance(self.molecules, Iterable):
            raise TypeError("Molecules must be an iterable collection")
        if not all(isinstance(m, Molecule) for m in self.molecules):
            raise TypeError("All items must be Molecule instances")

    def _validate_pair_indices(self, i: int, j: int) -> Tuple[int, int]:
        """Validate and normalize two molecule indices used for pairwise APIs."""
        try:
            i_idx = int(i)
            j_idx = int(j)
        except (TypeError, ValueError) as exc:
            raise TypeError("Pair indices must be integers.") from exc

        n = len(self.molecules)
        if not (0 <= i_idx < n) or not (0 <= j_idx < n):
            raise IndexError(
                f"Pair indices out of range for {n} molecules: ({i_idx}, {j_idx})."
            )

        return i_idx, j_idx

    @abstractmethod
    def group(self) -> Tuple[List[List[Molecule]], List[List[int]]]:
        """
        Main grouping method to return grouped molecules and their indices.

        Must be implemented by subclasses to define specific grouping logic.

        Returns:
            Tuple[List[List[Molecule]], List[List[int]]]: Tuple containing:
                - List of molecule groups (each group is a list of molecules)
                - List of index groups (corresponding indices for each group)
        """
        pass

    def record(self, **kwargs):
        """Public output entrypoint that delegates to strategy-specific recorder."""
        return self._record_results(**kwargs)

    def _record_results(self, **kwargs):
        """Strategy-specific result writer; implemented by concrete groupers."""
        raise NotImplementedError(
            f"{self.__class__.__name__} must implement _record_results()"
        )

    def _append_thermo_header(
        self, header_info: List[Tuple[str, Any]]
    ) -> None:
        """Helper to append thermochemistry parameters to header info."""
        thermo_energy_types = {"QHH", "QHG", "SP_QHG"}
        if (
            self.thermo_parameters
            and self.energy_type
            and self.energy_type.upper() in thermo_energy_types
        ):
            header_info.append(
                ("Thermochemistry Parameters", self.thermo_parameters)
            )

    def _append_input_usage_header(
        self, header_info: List[Tuple[str, Any]]
    ) -> None:
        """Append molecule usage statistics (used/skipped) to header info."""
        combined_skipped_ids = []
        for skipped_id in [*self.skipped_ids, *self._matrix_skipped_ids]:
            if skipped_id not in combined_skipped_ids:
                combined_skipped_ids.append(skipped_id)

        total_molecules = len(self.molecules) + len(self.skipped_ids)
        used_molecules = len(self.molecules) - len(
            self._matrix_skipped_indices
        )
        skipped_ids_text = (
            ", ".join(combined_skipped_ids) if combined_skipped_ids else "None"
        )
        header_info.extend(
            [
                ("Representative Strategy", self.representative_strategy),
                ("Total Molecules", total_molecules),
                ("Used Molecules", used_molecules),
                ("Skipped Molecules", len(combined_skipped_ids)),
                ("Skipped IDs", skipped_ids_text),
            ]
        )

    def unique(
        self, output_dir: str = None, prefix: str = "group"
    ) -> List[Molecule]:
        """
        Get unique representative molecules from each group.

        Returns the selected representative molecule from each group
        (always the first member in the finalized group ordering).
        Also generates XYZ files for each group in that same ordering
        in a dedicated subfolder.

        Args:
            output_dir (str): Base directory for output. If None, uses self.output_dir
                or current directory. Default is None.
            prefix (str): Prefix for output XYZ files. Default is "group".

        Returns:
            List[Molecule]: List of unique representative molecules (first member from each group).
        """

        # Use cached results if available, otherwise compute and cache
        if (
            self._cached_groups is not None
            and self._cached_group_indices is not None
        ):
            logger.info(
                f"[{self.__class__.__name__}] Using cached grouping results"
            )
            groups, group_indices = (
                self._cached_groups,
                self._cached_group_indices,
            )
        else:
            logger.info(
                f"[{self.__class__.__name__}] Computing groups for unique method"
            )
            groups, group_indices = self.group()
            # Cache the results
            self._cached_groups = groups
            self._cached_group_indices = group_indices

        unique_molecules = []

        # Determine base output directory
        base_dir = (
            output_dir if output_dir is not None else (self.output_dir or ".")
        )

        # Create dedicated subfolder for XYZ files (include label if provided)
        if self.label:
            result_folder = f"{self.label}_group_result"
        else:
            result_folder = "group_result"
        full_output_path = os.path.join(base_dir, result_folder)
        os.makedirs(full_output_path, exist_ok=True)

        logger.info(f"Creating XYZ files in folder: {full_output_path}")

        # Determine the file prefix (include label if provided)
        if self.label:
            file_prefix = f"{self.label}_{prefix}"
        else:
            file_prefix = prefix

        for i, (group, indices) in enumerate(zip(groups, group_indices)):
            # Group ordering is finalized by the representative strategy.
            ordered_pairs = list(zip(group, indices))

            # Write group XYZ file preserving representative-defined ordering
            group_filename = os.path.join(
                full_output_path, f"{file_prefix}_{i+1}.xyz"
            )
            with open(group_filename, "w") as f:
                for j, (mol, original_idx) in enumerate(ordered_pairs):
                    # Write the molecule coordinates
                    f.write(f"{mol.num_atoms}\n")

                    # Determine original index label (use conformer_id if available)
                    if self.conformer_ids is not None:
                        original_label = self.conformer_ids[original_idx]
                    else:
                        original_label = str(original_idx + 1)

                    # Create comment line with energy info and original molecule index
                    if mol.energy is not None:
                        comment = f"Group {i+1} Member {j+1} Original_Index: {original_label} {self.energy_type}(Hartree): {mol.energy:.8f}"
                    else:
                        comment = f"Group {i+1} Member {j+1} Original_Index: {original_label} {self.energy_type}: N/A"

                    f.write(f"{comment}\n")

                    # Write coordinates
                    for symbol, position in zip(
                        mol.chemical_symbols, mol.positions
                    ):
                        f.write(
                            f"{symbol:2s} {position[0]:15.10f} {position[1]:15.10f} {position[2]:15.10f}\n"
                        )

            logger.info(
                f"Written group {i+1} with {len(ordered_pairs)} molecules to {group_filename}"
            )

            # Representative is always the first member by convention.
            if ordered_pairs:
                unique_molecules.append(ordered_pairs[0][0])

        logger.info(
            f"Generated {len(groups)} group XYZ files in {full_output_path}"
        )

        return unique_molecules


class MatrixGrouper(MoleculeGrouper):
    """
    Base class for groupers that operate on a pairwise distance matrix.

    Provides shared matrix validation, hierarchical complete-linkage
    clustering, threshold-based grouping, requested-num-groups grouping,
    and representative-aware group ordering.
    """

    supports_matrix_representative_selection = True

    def __init__(
        self,
        molecules: Iterable[Molecule],
        threshold=None,
        num_groups=None,
        num_procs: int = 1,
        label: str = None,
        conformer_ids: List[str] = None,
        skipped_ids: List[str] = None,
        matrix_format: str = "xlsx",
        output_dir: str = None,
        energy_type: str = "E",
        thermo_parameters: str = None,
        representative_strategy: str = "lowest",
    ):
        super().__init__(
            molecules,
            num_procs=num_procs,
            label=label,
            conformer_ids=conformer_ids,
            skipped_ids=skipped_ids,
            matrix_format=matrix_format,
            output_dir=output_dir,
            energy_type=energy_type,
            thermo_parameters=thermo_parameters,
            representative_strategy=representative_strategy,
        )

        if threshold is not None and num_groups is not None:
            raise ValueError(
                "Cannot specify both threshold (-T) and num_groups (-N). Please use only one."
            )

        if threshold is not None and threshold < 0:
            raise ValueError("threshold must be non-negative.")

        if num_groups is not None and num_groups < 1:
            raise ValueError("num_groups must be at least 1.")

        if threshold is None and num_groups is None:
            threshold = 0.5

        self.threshold = threshold
        self.num_groups = num_groups
        self._auto_threshold = None

    def _log_progress(
        self, completed: int, total: int, next_progress: int
    ) -> int:
        """
        Log matrix calculation progress at 10% intervals.
        """
        if total <= 0:
            return next_progress

        current_percent = (completed * 100) // total
        while current_percent >= next_progress and next_progress <= 100:
            logger.info(
                f"[{self.__class__.__name__}] Matrix calculation progress: {next_progress}% ({completed}/{total})"
            )
            next_progress += 10

        return next_progress

    def _build_singleton_groups(
        self,
    ) -> Tuple[List[List[Molecule]], List[List[int]]]:
        """Build singleton groups for the current molecule collection."""
        index_groups = [[i] for i in range(len(self.molecules))]
        groups = [[self.molecules[i]] for i in range(len(self.molecules))]
        return groups, index_groups

    def _get_original_index_label(self, original_index: int) -> str:
        """Return the display label for an original molecule index."""
        if self.conformer_ids is not None and 0 <= original_index < len(
            self.conformer_ids
        ):
            return self.conformer_ids[original_index]
        return str(original_index + 1)

    def _set_matrix_skipped_indices(self, skipped_indices: List[int]) -> None:
        """Store matrix-based skipped structure indices and display IDs."""
        self._matrix_skipped_indices = list(skipped_indices)
        self._matrix_skipped_ids = [
            self._get_original_index_label(idx) for idx in skipped_indices
        ]

    def _validate_distance_matrix(
        self, distance_matrix: np.ndarray
    ) -> np.ndarray:
        """Validate a pairwise distance matrix while allowing +inf placeholders."""
        clustering_matrix = np.asarray(distance_matrix, dtype=float)

        if clustering_matrix.ndim != 2:
            raise ValueError("Distance matrix must be 2-dimensional.")

        n_rows, n_cols = clustering_matrix.shape
        if n_rows != n_cols:
            raise ValueError("Distance matrix must be square.")

        if np.any(np.isnan(clustering_matrix)):
            raise ValueError(
                "Pairwise distance matrix contains non-finite values. "
                "Check that all pairwise calculations completed successfully."
            )

        if np.any(np.isneginf(clustering_matrix)):
            raise ValueError(
                "Distance matrix cannot contain negative distances."
            )

        if np.any(clustering_matrix < 0):
            raise ValueError(
                "Distance matrix cannot contain negative distances."
            )

        if not np.allclose(clustering_matrix, clustering_matrix.T):
            raise ValueError("Distance matrix must be symmetric.")

        if not np.allclose(np.diag(clustering_matrix), 0.0):
            raise ValueError("Distance matrix diagonal must be zero.")

        return clustering_matrix

    def _prepare_clustering_submatrix(
        self,
        distance_matrix: np.ndarray,
        original_indices: Optional[List[int]] = None,
    ) -> Tuple[np.ndarray, List[int]]:
        """Remove structures responsible for +inf pairs and return a finite submatrix."""
        clustering_matrix = self._validate_distance_matrix(distance_matrix)
        n = clustering_matrix.shape[0]

        if original_indices is None:
            original_indices = list(range(n))
        elif len(original_indices) != n:
            raise ValueError(
                "original_indices must match the size of the distance matrix."
            )

        remaining_local_indices = list(range(n))
        skipped_local_indices = []

        while remaining_local_indices:
            current_matrix = clustering_matrix[
                np.ix_(remaining_local_indices, remaining_local_indices)
            ]
            inf_mask = np.isposinf(current_matrix).copy()
            np.fill_diagonal(inf_mask, False)

            if not np.any(inf_mask):
                break

            inf_counts = np.sum(inf_mask, axis=1)
            worst_current_index = int(np.argmax(inf_counts))
            skipped_local_indices.append(
                remaining_local_indices.pop(worst_current_index)
            )

        skipped_original_indices = [
            original_indices[idx] for idx in skipped_local_indices
        ]
        self._set_matrix_skipped_indices(skipped_original_indices)

        valid_local_indices = [
            idx for idx in range(n) if idx not in skipped_local_indices
        ]
        valid_original_indices = [
            original_indices[idx] for idx in valid_local_indices
        ]

        if not valid_local_indices:
            return np.zeros((0, 0), dtype=float), valid_original_indices

        clean_submatrix = clustering_matrix[
            np.ix_(valid_local_indices, valid_local_indices)
        ]
        return clean_submatrix, valid_original_indices

    def _build_complete_linkage_tree(
        self, distance_matrix: np.ndarray
    ) -> Optional[np.ndarray]:
        """Build a SciPy hierarchical complete-linkage tree from a distance matrix."""
        clustering_matrix = self._validate_distance_matrix(distance_matrix)

        if clustering_matrix.shape[0] < 2:
            return None

        if np.any(np.isposinf(clustering_matrix)):
            raise ValueError(
                "Distance matrix for linkage must be finite after removing skipped structures."
            )

        condensed_distances = squareform(clustering_matrix, checks=False)
        return linkage(condensed_distances, method="complete")

    def _build_groups_from_cluster_labels(
        self,
        cluster_labels: np.ndarray,
        original_indices: Optional[List[int]] = None,
    ) -> Tuple[List[List[Molecule]], List[List[int]]]:
        """Convert cluster labels to deterministically ordered groups."""
        clusters: dict[int, List[int]] = {}
        cluster_ids = np.asarray(cluster_labels, dtype=int).reshape(-1)
        for idx, cluster_id in enumerate(cluster_ids.tolist()):
            if cluster_id not in clusters:
                clusters[cluster_id] = []
            clusters[cluster_id].append(idx)

        if original_indices is None:
            original_indices = list(range(len(cluster_ids)))

        index_groups = sorted(
            (
                sorted(original_indices[idx] for idx in indices)
                for indices in clusters.values()
            ),
            key=lambda indices: indices[0],
        )
        groups = [
            [self.molecules[idx] for idx in index_group]
            for index_group in index_groups
        ]
        return groups, index_groups

    def _matrix_score_key(
        self,
        original_index: int,
        score: float,
    ) -> Tuple[float, int, float, int]:
        """Deterministic score key for center/top3 representative selection."""
        energy_missing, energy_value, _ = self._energy_sort_key(original_index)
        return (score, energy_missing, energy_value, original_index)

    def _mean_distance_to_group(
        self,
        candidate_original_index: int,
        group_original_indices: List[int],
        distance_matrix: np.ndarray,
        original_to_matrix_index: Dict[int, int],
    ) -> float:
        """Mean pairwise distance from one candidate to the rest of a group."""
        if len(group_original_indices) <= 1:
            return 0.0

        candidate_matrix_index = original_to_matrix_index[
            candidate_original_index
        ]
        distances = []
        for other_original_index in group_original_indices:
            if other_original_index == candidate_original_index:
                continue
            other_matrix_index = original_to_matrix_index[other_original_index]
            distances.append(
                float(
                    distance_matrix[candidate_matrix_index, other_matrix_index]
                )
            )

        return float(np.mean(distances)) if distances else 0.0

    def _order_index_group_by_center(
        self,
        index_group: List[int],
        distance_matrix: np.ndarray,
        original_to_matrix_index: Dict[int, int],
    ) -> List[int]:
        """Order a group by centrality score ascending."""
        if len(index_group) <= 1:
            return list(index_group)

        scored = []
        for original_index in index_group:
            score = self._mean_distance_to_group(
                original_index,
                index_group,
                distance_matrix,
                original_to_matrix_index,
            )
            scored.append((original_index, score))

        scored.sort(key=lambda item: self._matrix_score_key(item[0], item[1]))
        return [original_index for original_index, _ in scored]

    def _order_index_group_by_top3(
        self,
        index_group: List[int],
        distance_matrix: np.ndarray,
        original_to_matrix_index: Dict[int, int],
    ) -> List[int]:
        """Order group using top3 representative selection with lowest fallback."""
        if len(index_group) < 3:
            return self._order_index_group_by_energy(index_group)

        energy_ordered = self._order_index_group_by_energy(index_group)
        candidate_indices = energy_ordered[:3]

        best_candidate = min(
            candidate_indices,
            key=lambda original_index: self._matrix_score_key(
                original_index,
                self._mean_distance_to_group(
                    original_index,
                    index_group,
                    distance_matrix,
                    original_to_matrix_index,
                ),
            ),
        )

        remaining = [
            original_index
            for original_index in energy_ordered
            if original_index != best_candidate
        ]
        return [best_candidate, *remaining]

    def _order_groups_by_representative_strategy(
        self,
        index_groups: List[List[int]],
        distance_matrix: np.ndarray,
        matrix_original_indices: Optional[List[int]] = None,
    ) -> Tuple[List[List[Molecule]], List[List[int]]]:
        """Order each group so that index_group[0] is the selected representative."""
        if matrix_original_indices is None:
            matrix_original_indices = list(range(distance_matrix.shape[0]))

        original_to_matrix_index = {
            original_index: matrix_index
            for matrix_index, original_index in enumerate(
                matrix_original_indices
            )
        }

        ordered_index_groups: List[List[int]] = []
        for index_group in index_groups:
            if self.representative_strategy == "lowest":
                ordered_group = self._order_index_group_by_energy(index_group)
            elif self.representative_strategy == "center":
                ordered_group = self._order_index_group_by_center(
                    index_group,
                    distance_matrix,
                    original_to_matrix_index,
                )
            else:  # top3
                ordered_group = self._order_index_group_by_top3(
                    index_group,
                    distance_matrix,
                    original_to_matrix_index,
                )

            ordered_index_groups.append(ordered_group)

        return (
            self._rebuild_groups_from_index_groups(ordered_index_groups),
            ordered_index_groups,
        )

    def _cut_tree_by_threshold(
        self, linkage_matrix: np.ndarray, threshold: float
    ) -> np.ndarray:
        """Cut a hierarchical linkage tree at distance <= threshold."""
        return fcluster(linkage_matrix, t=threshold, criterion="distance")

    def _resolve_threshold(self, threshold: Optional[float]) -> float:
        """Resolve and validate a configured or per-call threshold."""
        effective_threshold = (
            self.threshold if threshold is None else threshold
        )
        if effective_threshold is None:
            raise ValueError(
                "A threshold must be provided either when constructing the "
                "grouper or when calling group_by_threshold()."
            )

        try:
            effective_threshold = float(effective_threshold)
        except (TypeError, ValueError) as exc:
            raise ValueError("threshold must be a finite number.") from exc

        if not np.isfinite(effective_threshold):
            raise ValueError("threshold must be a finite number.")
        if effective_threshold < 0:
            raise ValueError("threshold must be non-negative.")
        return effective_threshold

    def _resolve_num_groups(self, num_groups: Optional[int]) -> int:
        """Resolve and validate a configured or per-call num_groups."""
        effective_num_groups = (
            self.num_groups if num_groups is None else num_groups
        )
        if effective_num_groups is None:
            raise ValueError(
                "num_groups must be provided either when constructing the "
                "grouper or when calling group_by_num_groups()."
            )

        try:
            effective_num_groups = int(effective_num_groups)
        except (TypeError, ValueError) as exc:
            raise ValueError("num_groups must be an integer >= 1.") from exc

        if effective_num_groups < 1:
            raise ValueError("num_groups must be at least 1.")
        return effective_num_groups

    def group_by_threshold(
        self,
        matrix: np.ndarray,
        threshold: Optional[float] = None,
    ) -> Tuple[List[List[Molecule]], List[List[int]]]:
        """Group a precomputed output distance matrix at a chosen threshold.

        The base implementation accepts a symmetric, non-negative distance
        matrix whose rows and columns follow ``self.molecules``. Subclasses
        whose recorded matrices use different
        semantics, such as similarity or signed energy differences, normalize
        their matrices before delegating to the private clustering helper.

        The most recent result replaces the grouping cache used by ``unique()``;
        the input matrix and the grouper's configured threshold are not modified.
        """
        (
            distance_matrix,
            distance_threshold,
            original_indices,
            output_skipped_indices,
        ) = self._normalize_output_matrix_for_threshold(matrix, threshold)

        self._auto_threshold = None

        groups, index_groups = self._cluster_distance_matrix_by_threshold(
            distance_matrix,
            distance_threshold,
            original_indices=original_indices,
        )
        self._set_matrix_skipped_indices(
            sorted(
                set(output_skipped_indices) | set(self._matrix_skipped_indices)
            )
        )
        self._cached_groups = groups
        self._cached_group_indices = index_groups
        return groups, index_groups

    def _normalize_output_matrix_for_threshold(
        self,
        matrix: np.ndarray,
        threshold: Optional[float],
    ) -> Tuple[np.ndarray, float, Optional[List[int]], List[int]]:
        """Normalize a recorded output matrix for distance clustering."""
        matrix_array = np.asarray(matrix)
        if matrix_array.ndim == 2 and matrix_array.shape[0] != len(
            self.molecules
        ):
            raise ValueError(
                "Matrix dimensions must match the number of molecules."
            )

        return (
            matrix_array,
            self._resolve_threshold(threshold),
            None,
            [],
        )

    def group_by_num_groups(
        self,
        matrix: np.ndarray,
        num_groups: Optional[int] = None,
    ) -> Tuple[List[List[Molecule]], List[List[int]]]:
        """Group a precomputed output matrix by requested group count."""
        (
            distance_matrix,
            effective_num_groups,
            original_indices,
            output_skipped_indices,
        ) = self._normalize_output_matrix_for_num_groups(matrix, num_groups)
        groups, index_groups = self._cluster_distance_matrix_by_num_groups(
            distance_matrix,
            effective_num_groups,
            original_indices=original_indices,
        )
        self._set_matrix_skipped_indices(
            sorted(
                set(output_skipped_indices) | set(self._matrix_skipped_indices)
            )
        )
        self._cached_groups = groups
        self._cached_group_indices = index_groups
        self._postprocess_auto_threshold_for_num_groups()
        return groups, index_groups

    def _normalize_output_matrix_for_num_groups(
        self,
        matrix: np.ndarray,
        num_groups: Optional[int],
    ) -> Tuple[np.ndarray, int, Optional[List[int]], List[int]]:
        """Normalize a recorded output matrix for num-groups clustering."""
        matrix_array = np.asarray(matrix)
        if matrix_array.ndim == 2 and matrix_array.shape[0] != len(
            self.molecules
        ):
            raise ValueError(
                "Matrix dimensions must match the number of molecules."
            )

        return (
            matrix_array,
            self._resolve_num_groups(num_groups),
            None,
            [],
        )

    def _postprocess_auto_threshold_for_num_groups(self) -> None:
        """Hook for subclasses to transform auto threshold display units."""
        return None

    def _cluster_distance_matrix_by_threshold(
        self,
        distance_matrix: np.ndarray,
        distance_threshold: float,
        original_indices: Optional[List[int]] = None,
    ) -> Tuple[List[List[Molecule]], List[List[int]]]:
        """Cluster a normalized distance matrix by complete linkage."""
        clean_submatrix, valid_original_indices = (
            self._prepare_clustering_submatrix(
                distance_matrix, original_indices=original_indices
            )
        )

        if len(valid_original_indices) == 0:
            return [], []

        if len(valid_original_indices) == 1:
            original_index = valid_original_indices[0]
            return [[self.molecules[original_index]]], [[original_index]]

        linkage_matrix = self._build_complete_linkage_tree(clean_submatrix)
        cluster_labels = self._cut_tree_by_threshold(
            linkage_matrix, distance_threshold
        )
        groups, index_groups = self._build_groups_from_cluster_labels(
            cluster_labels,
            original_indices=valid_original_indices,
        )
        return self._order_groups_by_representative_strategy(
            index_groups,
            clean_submatrix,
            matrix_original_indices=valid_original_indices,
        )

    def _cluster_distance_matrix_by_num_groups(
        self,
        distance_matrix: np.ndarray,
        effective_num_groups: int,
        original_indices: Optional[List[int]] = None,
    ) -> Tuple[List[List[Molecule]], List[List[int]]]:
        """Group by requested number of groups using complete-linkage levels."""
        clean_submatrix, valid_original_indices = (
            self._prepare_clustering_submatrix(
                distance_matrix, original_indices=original_indices
            )
        )
        n = len(valid_original_indices)

        if n == 0:
            self._auto_threshold = None
            return [], []

        if n == 1:
            self._auto_threshold = None
            original_index = valid_original_indices[0]
            return [[self.molecules[original_index]]], [[original_index]]

        if effective_num_groups >= n:
            logger.info(
                f"[{self.__class__.__name__}] Requested {effective_num_groups} groups but only {n} molecules are available. Creating {n} groups."
            )
            self._auto_threshold = None
            groups = [[self.molecules[idx]] for idx in valid_original_indices]
            index_groups = [[idx] for idx in valid_original_indices]
            return groups, index_groups

        linkage_matrix = self._build_complete_linkage_tree(clean_submatrix)
        unique_distances = np.unique(linkage_matrix[:, 2])

        best_labels = np.arange(1, n + 1)
        best_distance = None

        for distance in unique_distances:
            cluster_labels = self._cut_tree_by_threshold(
                linkage_matrix, distance
            )
            num_groups = len(np.unique(cluster_labels))

            if num_groups < effective_num_groups:
                break

            best_labels = cluster_labels
            best_distance = float(distance)

            if num_groups == effective_num_groups:
                break

        self._auto_threshold = best_distance

        logger.info(
            f"[{self.__class__.__name__}] Selected complete-linkage distance level: "
            f"{best_distance:.7f} to keep {len(np.unique(best_labels))} groups for requested={effective_num_groups}"
            if best_distance is not None
            else f"[{self.__class__.__name__}] No linkage merges applied; retaining {n} groups for requested={effective_num_groups}"
        )

        groups, index_groups = self._build_groups_from_cluster_labels(
            best_labels,
            original_indices=valid_original_indices,
        )
        actual_groups = len(groups)

        logger.info(
            f"[{self.__class__.__name__}] Created {actual_groups} groups (requested: {effective_num_groups})"
        )

        return self._order_groups_by_representative_strategy(
            index_groups,
            clean_submatrix,
            matrix_original_indices=valid_original_indices,
        )
