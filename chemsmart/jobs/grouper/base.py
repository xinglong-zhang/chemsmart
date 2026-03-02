"""
Base classes for molecular structure grouping algorithms.

This module contains the core abstract base class and utilities:
- MoleculeGrouper: Abstract base class for all grouper implementations
- ResultsRecorder: Unified results output with multiple format support
- StructureGrouperConfig: Configuration container for structure matching
"""

import logging
import os
from abc import ABC, abstractmethod
from typing import Any, Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd

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
        output_format (str): Output format ('xlsx', 'csv', 'txt').
        conformer_ids (list[str]): Custom IDs for labeling molecules.
    """

    SUPPORTED_FORMATS = {"xlsx", "csv", "txt"}

    def __init__(
        self,
        output_dir: str,
        label: str = None,
        output_format: str = "xlsx",
        conformer_ids: List[str] = None,
    ):
        """
        Initialize the results recorder.

        Args:
            output_dir (str): Directory for output files.
            label (str): Label prefix for output filenames.
            output_format (str): Output format. Defaults to 'xlsx'.
                Supported: 'xlsx', 'csv', 'txt'.
            conformer_ids (list[str]): Custom IDs for labeling molecules.
        """
        self.output_dir = output_dir
        self.label = label
        self.output_format = output_format.lower().lstrip(".")
        self.conformer_ids = conformer_ids

        if self.output_format not in self.SUPPORTED_FORMATS:
            raise ValueError(
                f"Unsupported output format: '{output_format}'. "
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
            f"{label_prefix}{grouper_name}{suffix_part}.{self.output_format}"
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

        if self.output_format == "xlsx":
            self._write_xlsx(
                filename,
                header_info,
                sheets_data,
                matrix_data,
                startrow,
                float_format,
            )
        elif self.output_format == "csv":
            self._write_csv(filename, header_info, sheets_data, matrix_data)
        elif self.output_format == "txt":
            self._write_txt(filename, header_info, sheets_data, matrix_data)
        else:
            raise ValueError(
                f"Unsupported output format: '{self.output_format}'"
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
        output_format (str): Output format for results ('xlsx', 'csv', 'txt').
        output_dir (str): Base directory for output files.
    """

    def __init__(
        self,
        molecules: Iterable[Molecule],
        num_procs: int = 1,
        label: str = None,
        conformer_ids: List[str] = None,
        output_format: str = "xlsx",
        output_dir: str = None,
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
            output_format (str): Output format for results ('xlsx', 'csv', 'txt').
                Defaults to 'xlsx'.
            output_dir (str): Base directory for output files. If None, uses current
                working directory. Defaults to None.
        """
        self.molecules = molecules
        self.num_procs = int(max(1, num_procs))
        self.label = label
        self.conformer_ids = conformer_ids
        self.output_format = output_format
        self.output_dir = output_dir

        # Cache for avoiding repeated grouping calculations
        self._cached_groups = None
        self._cached_group_indices = None

        self._validate_inputs()

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
            output_format=self.output_format,
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

    def unique(
        self, output_dir: str = None, prefix: str = "group"
    ) -> List[Molecule]:
        """
        Get unique representative molecules from each group.

        Returns the lowest energy molecule from each group as a representative
        of that structural family. Also generates XYZ files for each group,
        sorted by energy, in a dedicated subfolder.

        Args:
            output_dir (str): Base directory for output. If None, uses self.output_dir
                or current directory. Default is None.
            prefix (str): Prefix for output XYZ files. Default is "group".

        Returns:
            List[Molecule]: List of unique representative molecules (lowest energy from each group).
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
            # Create tuples of (molecule, original_index) for tracking
            mol_index_pairs = list(zip(group, indices))

            # Filter molecules that have energy information and sort by energy
            molecules_with_energy = [
                (mol, idx)
                for mol, idx in mol_index_pairs
                if mol.energy is not None
            ]
            molecules_without_energy = [
                (mol, idx)
                for mol, idx in mol_index_pairs
                if mol.energy is None
            ]

            # Sort molecules with energy by energy (ascending - lowest first)
            if molecules_with_energy:
                sorted_pairs = sorted(
                    molecules_with_energy, key=lambda pair: pair[0].energy
                )
                # Add molecules without energy at the end
                sorted_pairs.extend(molecules_without_energy)
            else:
                # If no molecules have energy, use original group order
                sorted_pairs = mol_index_pairs

            # Write group XYZ file with all molecules sorted by energy
            group_filename = os.path.join(
                full_output_path, f"{file_prefix}_{i+1}.xyz"
            )
            with open(group_filename, "w") as f:
                for j, (mol, original_idx) in enumerate(sorted_pairs):
                    # Write the molecule coordinates
                    f.write(f"{mol.num_atoms}\n")

                    # Determine original index label (use conformer_id if available)
                    if self.conformer_ids is not None:
                        original_label = self.conformer_ids[original_idx]
                    else:
                        original_label = str(original_idx + 1)

                    # Create comment line with energy info and original molecule index
                    if mol.energy is not None:
                        comment = f"Group {i+1} Member {j+1} Original_Index: {original_label} Energy(Hartree): {mol.energy:.8f}"
                    else:
                        comment = f"Group {i+1} Member {j+1} Original_Index: {original_label} Energy: N/A"

                    f.write(f"{comment}\n")

                    # Write coordinates
                    for symbol, position in zip(
                        mol.chemical_symbols, mol.positions
                    ):
                        f.write(
                            f"{symbol:2s} {position[0]:15.10f} {position[1]:15.10f} {position[2]:15.10f}\n"
                        )

            logger.info(
                f"Written group {i+1} with {len(sorted_pairs)} molecules to {group_filename}"
            )

            # Add the lowest energy molecule (first in sorted pairs) as representative
            unique_molecules.append(sorted_pairs[0][0])

        logger.info(
            f"Generated {len(groups)} group XYZ files in {full_output_path}"
        )

        return unique_molecules
