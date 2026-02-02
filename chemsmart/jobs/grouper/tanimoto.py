"""
Tanimoto fingerprint similarity-based molecular grouping.

Groups molecules based on fingerprint similarity using Tanimoto coefficient.
"""

import logging
import os
from multiprocessing.pool import ThreadPool
from typing import Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.rdFingerprintGenerator import GetRDKitFPGenerator

from chemsmart.io.molecules.structure import Molecule

from .runner import MoleculeGrouper

logger = logging.getLogger(__name__)


class TanimotoSimilarityGrouper(MoleculeGrouper):
    """
    Groups molecules based on fingerprint similarity using Tanimoto coefficient.

    This class supports different fingerprint types and uses connected components
    clustering to group structurally similar molecules.

    Supported fingerprint types:
    - "rdkit": RDKit topological fingerprint (default)
    - "rdk": Legacy RDKit fingerprint
    - "morgan": Morgan (circular) fingerprint (radius=2)
    - "maccs": MACCS keys (166 bits)
    - "atompair": Atom pair fingerprint
    - "torsion": Topological torsion fingerprint
    - "usr": Ultrafast Shape Recognition (3D descriptor)
    - "usrcat": USR with CREDO Atom Types (3D descriptor)
    """

    def __init__(
        self,
        molecules: Iterable[Molecule],
        threshold=None,  # Tanimoto similarity threshold
        num_groups=None,  # Number of groups to create (alternative to threshold)
        num_procs: int = 1,
        fingerprint_type: str = "rdkit",
        use_rdkit_fp: bool = None,  # Legacy support
        label: str = None,  # Label for output files
        ignore_hydrogens: bool = False,
        conformer_ids: List[str] = None,
        **kwargs,
    ):
        """
        Initialize Tanimoto similarity-based molecular grouper.

        Args:
            molecules (Iterable[Molecule]): Collection of molecules to group.
            threshold (float): Tanimoto similarity threshold. Defaults to 0.9.
                Ignored if num_groups is specified.
            num_groups (int): Number of groups to create. When specified,
                automatically determines threshold to create this many groups.
            num_procs (int): Number of processes for parallel computation.
            fingerprint_type (str): Type of fingerprint to use.
                Options: "rdkit", "rdk", "morgan", "maccs", "atompair",
                "torsion", "usr", "usrcat". Defaults to "rdkit".
            use_rdkit_fp (bool): Legacy parameter. If True, sets fingerprint_type="rdkit".
                If False, sets fingerprint_type="rdk".
            label (str): Label/name for output files. Defaults to None.
            ignore_hydrogens (bool): Whether to remove hydrogens before fingerprint
                calculation. Defaults to False.
            conformer_ids (list[str]): Custom IDs for each molecule (e.g., ['c1', 'c2']).
        """
        super().__init__(
            molecules, num_procs, label=label, conformer_ids=conformer_ids
        )

        self.ignore_hydrogens = ignore_hydrogens

        # Validate that threshold and num_groups are mutually exclusive
        if threshold is not None and num_groups is not None:
            raise ValueError(
                "Cannot specify both threshold (-t) and num_groups (-N). Please use only one."
            )

        if threshold is None and num_groups is None:
            threshold = 0.9
        self.threshold = threshold
        self.num_groups = num_groups

        if use_rdkit_fp is not None:
            self.fingerprint_type = "rdkit" if use_rdkit_fp else "rdk"
        else:
            self.fingerprint_type = fingerprint_type.lower()

        # Convert molecules to list for indexing
        self.molecules = list(molecules)

        # Convert valid molecules to RDKit format
        self.rdkit_molecules = []
        self.valid_molecules = []
        for mol in self.molecules:
            rdkit_mol = mol.to_rdkit()
            if rdkit_mol is not None:
                # Remove hydrogens if requested
                if self.ignore_hydrogens:
                    rdkit_mol = Chem.RemoveHs(rdkit_mol)
                self.rdkit_molecules.append(rdkit_mol)
                self.valid_molecules.append(mol)

    def _get_fingerprint(self, rdkit_mol: Chem.Mol) -> Optional[object]:
        """
        Generate a fingerprint for a molecule.

        Args:
            rdkit_mol (Chem.Mol): RDKit molecule object.

        Returns:
            Optional[object]: Molecular fingerprint (BitVect or np.ndarray)
                or None if generation fails.
        """
        try:
            if self.fingerprint_type == "rdkit":
                return GetRDKitFPGenerator().GetFingerprint(rdkit_mol)
            elif self.fingerprint_type == "rdk":
                return Chem.RDKFingerprint(rdkit_mol)
            elif self.fingerprint_type == "morgan":
                return rdMolDescriptors.GetMorganFingerprintAsBitVect(
                    rdkit_mol, 2
                )
            elif self.fingerprint_type == "maccs":
                return rdMolDescriptors.GetMACCSKeysFingerprint(rdkit_mol)
            elif self.fingerprint_type == "atompair":
                return rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(
                    rdkit_mol
                )
            elif self.fingerprint_type == "torsion":
                return rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect(
                    rdkit_mol
                )
            elif self.fingerprint_type == "usr":
                conf = rdkit_mol.GetConformer()
                return np.array(
                    rdMolDescriptors.GetUSR(rdkit_mol, confId=conf.GetId())
                )
            elif self.fingerprint_type == "usrcat":
                conf = rdkit_mol.GetConformer()
                return np.array(
                    rdMolDescriptors.GetUSRCAT(rdkit_mol, confId=conf.GetId())
                )
            else:
                logger.warning(
                    f"Unknown fingerprint type: {self.fingerprint_type}, using RDKit default."
                )
                return GetRDKitFPGenerator().GetFingerprint(rdkit_mol)
        except Exception as e:
            logger.warning(f"Fingerprint generation failed: {str(e)}")
            return None

    def group(self) -> Tuple[List[List[Molecule]], List[List[int]]]:
        """
        Groups molecules based on Tanimoto similarity clustering.

        Computes fingerprints for all molecules, calculates pairwise
        Tanimoto similarities, and groups molecules using connected
        components clustering.

        Returns:
            Tuple[List[List[Molecule]], List[List[int]]]: Tuple containing:
                - List of molecule groups (each group is a list of molecules)
                - List of index groups (corresponding indices for each group)
        """
        import time

        # Record start time for grouping process
        grouping_start_time = time.time()

        n = len(self.molecules)
        print(
            f"[{self.__class__.__name__}] Starting fingerprint calculation for {n} molecules using {self.fingerprint_type} fingerprints"
        )

        # Compute fingerprints in parallel
        with ThreadPool(self.num_procs) as pool:
            fingerprints = pool.map(
                self._get_fingerprint, self.rdkit_molecules
            )

        # Filter valid fingerprints
        valid_indices = [
            i for i, fp in enumerate(fingerprints) if fp is not None
        ]
        valid_fps = [fingerprints[i] for i in valid_indices]
        num_valid = len(valid_indices)

        if num_valid == 0:
            return [], []  # No valid molecules

        print(
            f"[{self.__class__.__name__}] Computing Tanimoto similarities for {num_valid} valid molecules"
        )

        # Compute similarity matrix
        similarity_matrix = np.zeros((num_valid, num_valid), dtype=np.float32)

        # Check if we are using numpy arrays (USR/USRCAT)
        if valid_fps and isinstance(valid_fps[0], np.ndarray):
            # Calculate Tanimoto for continuous variables (vectors)
            fps_array = np.array(valid_fps)
            dot_products = np.dot(fps_array, fps_array.T)
            norms_sq = np.diag(dot_products)
            denominator = norms_sq[:, None] + norms_sq[None, :] - dot_products
            denominator[denominator == 0] = 1e-9
            similarity_matrix = dot_products / denominator
        else:
            # Use RDKit DataStructs for BitVects
            pairs = [
                (i, j)
                for i in range(num_valid)
                for j in range(i + 1, num_valid)
            ]

            with ThreadPool(self.num_procs) as pool:
                similarities = pool.starmap(
                    DataStructs.FingerprintSimilarity,
                    [(valid_fps[i], valid_fps[j]) for i, j in pairs],
                )

            # Fill similarity matrix
            for (i, j), sim in zip(pairs, similarities):
                similarity_matrix[i, j] = similarity_matrix[j, i] = sim

        # Choose grouping strategy based on parameters
        if self.num_groups is not None:
            groups, index_groups = self._group_by_num_groups(
                similarity_matrix, valid_indices
            )
        else:
            groups, index_groups = self._group_by_threshold(
                similarity_matrix, valid_indices
            )

        # Calculate total grouping time
        grouping_end_time = time.time()
        grouping_time = grouping_end_time - grouping_start_time

        # Save Tanimoto matrix to group_result folder
        if self.label:
            output_dir = f"{self.label}_group_result"
        else:
            output_dir = "group_result"
        os.makedirs(output_dir, exist_ok=True)

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

        # Save full matrix
        self._save_tanimoto_matrix(
            similarity_matrix,
            matrix_filename,
            valid_indices,
            grouping_time,
            groups,
            index_groups,
        )

        # Cache results
        self._cached_groups = groups
        self._cached_group_indices = index_groups

        return groups, index_groups

    def _group_by_threshold(self, similarity_matrix, valid_indices):
        """Threshold-based grouping for Tanimoto similarity using complete linkage."""
        adj_matrix = similarity_matrix >= self.threshold
        n = len(valid_indices)
        mol_groups, idx_groups = self._complete_linkage_grouping(
            adj_matrix, n, valid_indices
        )
        return mol_groups, idx_groups

    def _complete_linkage_grouping(self, adj_matrix, n, valid_indices):
        """Perform complete linkage grouping for Tanimoto similarity."""
        assigned = [False] * n
        mol_groups = []
        idx_groups = []

        for i in range(n - 1, -1, -1):
            if assigned[i]:
                continue
            current_group = [i]
            assigned[i] = True

            for j in range(i - 1, -1, -1):
                if assigned[j]:
                    continue
                can_join = all(
                    adj_matrix[j, member] for member in current_group
                )
                if can_join:
                    current_group.append(j)
                    assigned[j] = True

            current_group.sort()
            mol_groups.append(
                [self.valid_molecules[idx] for idx in current_group]
            )
            idx_groups.append([valid_indices[idx] for idx in current_group])

        mol_groups.reverse()
        idx_groups.reverse()
        return mol_groups, idx_groups

    def _group_by_num_groups(self, similarity_matrix, valid_indices):
        """Automatic grouping to create specified number of groups."""
        n = len(valid_indices)

        if self.num_groups >= n:
            print(
                f"[{self.__class__.__name__}] Requested {self.num_groups} groups but only {n} molecules. Creating {n} groups."
            )
            groups = [[self.valid_molecules[i]] for i in range(n)]
            index_groups = [[valid_indices[i]] for i in range(n)]
            return groups, index_groups

        # Extract similarity values for threshold finding
        similarity_values = []
        for i in range(n):
            for j in range(i + 1, n):
                similarity_values.append(similarity_matrix[i, j])

        threshold = self._find_optimal_similarity_threshold(
            similarity_values, similarity_matrix, n
        )
        self._auto_threshold = threshold

        print(
            f"[{self.__class__.__name__}] Auto-determined threshold: {threshold:.7f} to create {self.num_groups} groups"
        )

        adj_matrix = similarity_matrix >= threshold
        groups, index_groups = self._complete_linkage_grouping(
            adj_matrix, n, valid_indices
        )
        actual_groups = len(groups)

        print(
            f"[{self.__class__.__name__}] Created {actual_groups} groups (requested: {self.num_groups})"
        )

        if actual_groups > self.num_groups:
            groups, index_groups = self._merge_groups_to_target(
                groups, index_groups
            )

        return groups, index_groups

    def _find_optimal_similarity_threshold(
        self, similarity_values, similarity_matrix, n
    ):
        """Find similarity threshold using binary search."""
        sorted_similarities = sorted(
            [sim for sim in similarity_values if not np.isnan(sim)],
            reverse=True,
        )

        if not sorted_similarities:
            return 1.0

        low, high = 0, len(sorted_similarities) - 1
        best_threshold = sorted_similarities[0]

        while low <= high:
            mid = (low + high) // 2
            threshold = sorted_similarities[mid]
            adj_matrix = similarity_matrix >= threshold
            num_groups_found = self._count_groups(adj_matrix, n)

            if num_groups_found == self.num_groups:
                return threshold
            elif num_groups_found > self.num_groups:
                high = mid - 1
            else:
                best_threshold = threshold
                low = mid + 1

        return best_threshold

    def _count_groups(self, adj_matrix, n):
        """Count number of groups using complete linkage."""
        assigned = [False] * n
        num_groups = 0
        for i in range(n):
            if assigned[i]:
                continue
            current_group = [i]
            assigned[i] = True
            for j in range(i + 1, n):
                if assigned[j]:
                    continue
                can_join = all(adj_matrix[j, m] for m in current_group)
                if can_join:
                    current_group.append(j)
                    assigned[j] = True
            num_groups += 1
        return num_groups

    def _merge_groups_to_target(self, groups, index_groups):
        """Merge groups to reach target number."""
        while len(groups) > self.num_groups:
            min_idx = min(range(len(groups)), key=lambda i: len(groups[i]))
            largest_idx = max(
                range(len(groups)),
                key=lambda i: len(groups[i]) if i != min_idx else -1,
            )
            groups[largest_idx].extend(groups[min_idx])
            index_groups[largest_idx].extend(index_groups[min_idx])
            groups.pop(min_idx)
            index_groups.pop(min_idx)
        return groups, index_groups

    def _save_tanimoto_matrix(
        self,
        tanimoto_matrix: np.ndarray,
        filename: str,
        valid_indices: List[int],
        grouping_time: float = None,
        groups: List[List[Molecule]] = None,
        index_groups: List[List[int]] = None,
    ):
        """Save Tanimoto similarity matrix to Excel file."""
        n = len(self.molecules)

        if filename.endswith(".txt"):
            filename = filename[:-4] + ".xlsx"
        elif not filename.endswith(".xlsx"):
            filename = filename + ".xlsx"

        # Create full matrix with invalid molecules marked as NaN
        full_matrix = np.full((n, n), np.nan)
        for i, idx_i in enumerate(valid_indices):
            for j, idx_j in enumerate(valid_indices):
                full_matrix[idx_i, idx_j] = tanimoto_matrix[i, j]

        # Create DataFrame with labels
        if self.conformer_ids is not None and len(self.conformer_ids) == n:
            row_labels = self.conformer_ids
            col_labels = self.conformer_ids
        else:
            row_labels = [str(i + 1) for i in range(n)]
            col_labels = [str(j + 1) for j in range(n)]

        df = pd.DataFrame(full_matrix, index=row_labels, columns=col_labels)

        with pd.ExcelWriter(filename, engine="openpyxl") as writer:
            df.to_excel(
                writer,
                sheet_name="Tanimoto_Matrix",
                startrow=8,
                float_format="%.7f",
            )

            worksheet = writer.sheets["Tanimoto_Matrix"]

            row = 1
            worksheet[f"A{row}"] = (
                f"Full Tanimoto Similarity Matrix ({n}x{n}) - {self.__class__.__name__}"
            )
            row += 1
            worksheet[f"A{row}"] = f"Fingerprint Type: {self.fingerprint_type}"
            row += 1

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
                        f"Auto-determined Threshold: {self._auto_threshold:.7f}"
                    )
                    row += 1
            else:
                worksheet[f"A{row}"] = f"Threshold: {self.threshold}"
                row += 1

            worksheet[f"A{row}"] = f"Number of Processors: {self.num_procs}"
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

            # Add Groups sheet if groups are provided
            if groups is not None and index_groups is not None:
                self._write_groups_sheet(writer, groups, index_groups)

        logger.info(f"Tanimoto matrix saved to {filename}")

    def _write_groups_sheet(self, writer, groups, index_groups):
        """Write groups information to a separate sheet."""
        groups_data = []
        for i, indices in enumerate(index_groups):
            if self.conformer_ids is not None:
                member_labels = [self.conformer_ids[idx] for idx in indices]
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

    def __repr__(self):
        if self.num_groups is not None:
            return (
                f"{self.__class__.__name__}(num_groups={self.num_groups}, "
                f"num_procs={self.num_procs}, fingerprint_type={self.fingerprint_type})"
            )
        else:
            return (
                f"{self.__class__.__name__}(threshold={self.threshold}, "
                f"num_procs={self.num_procs}, fingerprint_type={self.fingerprint_type})"
            )
