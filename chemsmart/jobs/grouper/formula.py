"""
Formula-based molecular grouping.

Groups molecules by their molecular formula.
"""

import logging
from collections import defaultdict
from typing import Iterable, List, Tuple

import pandas as pd

from chemsmart.io.molecules.structure import Molecule

from .base import MoleculeGrouper

logger = logging.getLogger(__name__)


class FormulaGrouper(MoleculeGrouper):
    """
    Group molecules by molecular formula.

    Groups molecules that share identical molecular formulas. Uses the
    `Molecule.formula` property to determine chemical composition.

    Attributes:
        molecules (Iterable[Molecule]): Inherited; collection of molecules to
            group.
        num_procs (int): Inherited; number of worker processes (not used
            since formula comparison is fast).
    """

    def __init__(
        self,
        molecules: Iterable[Molecule],
        num_procs: int = 1,
        label: str = None,
        conformer_ids: List[str] = None,
        output_format: str = "xlsx",
        **kwargs,
    ):
        """
        Initialize formula-based molecular grouper.

        Args:
            molecules (Iterable[Molecule]): Collection of molecules to group.
            num_procs (int): Number of processes (unused but kept for API
                consistency).
            label (str): Label/name for output files. Defaults to None.
            conformer_ids (list[str]): Custom IDs for each molecule.
            output_format (str): Output format ('xlsx', 'csv', 'txt'). Defaults to 'xlsx'.
        """
        super().__init__(
            molecules,
            num_procs,
            label=label,
            conformer_ids=conformer_ids,
            output_format=output_format,
        )

    def group(self) -> Tuple[List[List[Molecule]], List[List[int]]]:
        """
        Group molecules by molecular formula.

        Analyzes molecular formulas and groups molecules with identical
        compositions. Also saves results to an Excel file.

        Returns:
            Tuple[List[List[Molecule]], List[List[int]]]: Tuple containing:
                - List of molecule groups (each group has same formula)
                - List of index groups (corresponding indices for each group)
        """
        import time

        grouping_start_time = time.time()

        formula_groups = defaultdict(list)
        formula_indices = defaultdict(list)

        for i, mol in enumerate(self.molecules):
            formula = mol.chemical_formula
            formula_groups[formula].append(mol)
            formula_indices[formula].append(i)

        # Convert to lists
        groups = list(formula_groups.values())
        index_groups = list(formula_indices.values())
        formulas = list(formula_groups.keys())

        grouping_time = time.time() - grouping_start_time

        # Save results to Excel
        self._save_formula_results(
            formulas, groups, index_groups, grouping_time
        )

        # Cache results
        self._cached_groups = groups
        self._cached_group_indices = index_groups

        logger.info(f"Found {len(groups)} unique formulas")
        for formula, group in zip(formulas, groups):
            logger.info(f"  {formula}: {len(group)} molecules")

        return groups, index_groups

    def _save_formula_results(
        self,
        formulas: List[str],
        groups: List[List[Molecule]],
        index_groups: List[List[int]],
        grouping_time: float = None,
    ):
        """Save formula grouping results to file using ResultsRecorder."""
        n = len(list(self.molecules))

        # Build header info
        header_info = [
            ("", f"Formula Grouping Results - {self.__class__.__name__}"),
            ("Total Molecules", n),
            ("Unique Formulas", len(formulas)),
            ("Num Procs", self.num_procs),
        ]

        if grouping_time is not None:
            header_info.append(
                ("Grouping Time", f"{grouping_time:.2f} seconds")
            )

        # Use ResultsRecorder to save
        recorder = self._get_results_recorder()

        # Build formula summary data
        formula_data = []
        for i, formula in enumerate(formulas):
            formula_data.append(
                {
                    "Formula": formula,
                    "Count": len(groups[i]),
                }
            )
        formula_df = pd.DataFrame(formula_data)

        # Build groups dataframe with Formula column using recorder's method
        groups_df = recorder.build_groups_dataframe(
            index_groups, n, extra_columns={"Formula": formulas}
        )
        # Reorder columns to put Formula before Members
        groups_df = groups_df[["Group", "Formula", "Members"]]

        # Build sheets data
        sheets_data = {
            "Formulas": formula_df,
            "Groups": groups_df,
        }

        recorder.record_results(
            grouper_name=self.__class__.__name__,
            header_info=header_info,
            sheets_data=sheets_data,
            matrix_data=None,
            suffix=None,
            startrow=6,
        )

    def __repr__(self):
        return f"{self.__class__.__name__}(num_procs={self.num_procs})"
