"""
Formula-based molecular grouping.

Groups molecules by their molecular formula.
"""

import logging
import os
from collections import defaultdict
from typing import Iterable, List, Tuple

import pandas as pd

from chemsmart.io.molecules.structure import Molecule

from .runner import MoleculeGrouper

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
        """
        super().__init__(
            molecules, num_procs, label=label, conformer_ids=conformer_ids
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
        """Save formula grouping results to Excel file."""
        n = len(list(self.molecules))

        # Create output directory
        if self.label:
            output_dir = f"{self.label}_group_result"
        else:
            output_dir = "group_result"
        os.makedirs(output_dir, exist_ok=True)

        # Create filename
        label_prefix = f"{self.label}_" if self.label else ""
        filename = os.path.join(
            output_dir,
            f"{label_prefix}{self.__class__.__name__}.xlsx",
        )

        with pd.ExcelWriter(filename, engine="openpyxl") as writer:
            # Sheet 1: Formula summary
            formula_data = []
            for i, formula in enumerate(formulas):
                formula_data.append(
                    {
                        "Formula": formula,
                        "Count": len(groups[i]),
                    }
                )

            formula_df = pd.DataFrame(formula_data)

            # Write with header info
            formula_df.to_excel(
                writer,
                sheet_name="Formulas",
                startrow=6,
                index=False,
            )

            worksheet = writer.sheets["Formulas"]
            row = 1
            worksheet[f"A{row}"] = (
                f"Formula Grouping Results - {self.__class__.__name__}"
            )
            row += 1
            worksheet[f"A{row}"] = f"Total Molecules: {n}"
            row += 1
            worksheet[f"A{row}"] = f"Unique Formulas: {len(formulas)}"
            row += 1
            if grouping_time is not None:
                worksheet[f"A{row}"] = (
                    f"Grouping Time: {grouping_time:.2f} seconds"
                )
                row += 1

            # Sheet 2: Groups detail
            groups_data = []
            for i, (formula, indices) in enumerate(
                zip(formulas, index_groups)
            ):
                if self.conformer_ids is not None:
                    member_labels = [
                        self.conformer_ids[idx] for idx in indices
                    ]
                else:
                    member_labels = [str(idx + 1) for idx in indices]

                groups_data.append(
                    {
                        "Group": i + 1,
                        "Formula": formula,
                        "Members": len(indices),
                        "Indices": str(member_labels),
                        "Representative": (
                            member_labels[0] if member_labels else "N/A"
                        ),
                    }
                )

            groups_df = pd.DataFrame(groups_data)
            groups_df.to_excel(writer, sheet_name="Groups", index=False)

        logger.info(f"Formula grouping results saved to {filename}")

    def __repr__(self):
        return f"{self.__class__.__name__}(num_procs={self.num_procs})"
