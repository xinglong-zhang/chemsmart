"""
Fingerprint- and shape-similarity-based molecular grouping.

Binary/topological fingerprints use Tanimoto similarity, while USR/USRCAT
3D descriptors use RDKit's native USR similarity score.
"""

import logging
from multiprocessing.pool import ThreadPool
from typing import Iterable, List, Optional, Tuple

import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import rdDetermineBonds, rdMolDescriptors
from rdkit.Chem.rdFingerprintGenerator import GetRDKitFPGenerator

from chemsmart.io.molecules.structure import Molecule

from .base import MatrixGrouper

logger = logging.getLogger(__name__)


class TanimotoSimilarityGrouper(MatrixGrouper):
    """
    Groups molecules based on fingerprint or 3D-descriptor similarity using
    hierarchical complete-linkage clustering.

    This class supports different fingerprint types and uses hierarchical
    complete-linkage clustering to group structurally or shape-similar molecules.

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
        label: str = None,  # Label for output files
        ignore_hydrogens: bool = False,
        conformer_ids: List[str] = None,
        matrix_format: str = "xlsx",
        energy_type: str = "E",
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
            label (str): Label/name for output files. Defaults to None.
            ignore_hydrogens (bool): Whether to remove hydrogens before fingerprint
                calculation. Defaults to False.
            conformer_ids (list[str]): Custom IDs for each molecule (e.g., ['c1', 'c2']).
            matrix_format (str): Output format ('xlsx', 'csv', 'txt'). Defaults to 'xlsx'.
        """
        if threshold is None and num_groups is None:
            threshold = 0.9
        if threshold is not None and not (0.0 <= threshold <= 1.0):
            raise ValueError(
                "Tanimoto similarity threshold must be between 0.0 and 1.0."
            )

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

        self.ignore_hydrogens = ignore_hydrogens

        self.fingerprint_type = fingerprint_type.lower()

        # Convert molecules to list for indexing
        self.molecules = list(molecules)

        # 2D/topological fingerprints use RDKit-only bond perception, while
        # USR/USRCAT preserve the previous conformer preparation pathway.
        self.rdkit_molecules = []
        self.rdkit_molecule_indices = []
        tanimoto_skipped_indices = []
        for original_idx, mol in enumerate(self.molecules):
            rdkit_mol = self._molecule_to_rdkit(mol)
            if rdkit_mol is not None:
                self.rdkit_molecules.append(rdkit_mol)
                self.rdkit_molecule_indices.append(original_idx)
            else:
                tanimoto_skipped_indices.append(original_idx)

        self._tanimoto_skipped_indices = sorted(set(tanimoto_skipped_indices))
        self._rdkit_molecule_by_original_index = {
            original_idx: rdkit_mol
            for rdkit_mol, original_idx in zip(
                self.rdkit_molecules, self.rdkit_molecule_indices
            )
        }

    def calculate_tanimoto_pair(self, i: int, j: int) -> float:
        """Public API: calculate pairwise similarity between two molecules.

        Returns NaN when either molecule cannot be prepared/fingerprinted for
        the configured descriptor type.
        """
        i_idx, j_idx = self._validate_pair_indices(i, j)

        rdkit_mol_i = self._rdkit_molecule_by_original_index.get(i_idx)
        rdkit_mol_j = self._rdkit_molecule_by_original_index.get(j_idx)
        if rdkit_mol_i is None or rdkit_mol_j is None:
            logger.warning(
                f"[{self.__class__.__name__}] Cannot compute similarity for pair "
                f"({i_idx}, {j_idx}): one or both molecules failed RDKit preparation"
            )
            return float("nan")

        fp_i = self._get_fingerprint(rdkit_mol_i)
        fp_j = self._get_fingerprint(rdkit_mol_j)
        if fp_i is None or fp_j is None:
            logger.warning(
                f"[{self.__class__.__name__}] Cannot compute similarity for pair "
                f"({i_idx}, {j_idx}): fingerprint generation failed"
            )
            return float("nan")

        if self.fingerprint_type in {"usr", "usrcat"}:
            return float(rdMolDescriptors.GetUSRScore(fp_i, fp_j))

        return float(DataStructs.FingerprintSimilarity(fp_i, fp_j))

    def _molecule_to_rdkit(self, mol: Molecule) -> Optional[Chem.Mol]:
        """Build an RDKit molecule using preparation appropriate to the descriptor."""
        try:
            if self.fingerprint_type in {"usr", "usrcat"}:
                # Preserve the previous working preparation for 3D descriptors.
                # USR/USRCAT operate on the RDKit molecule/conformer produced by
                # CHEMSMART without re-perceiving bond orders from every XYZ conformer.
                rdkit_mol = mol.to_rdkit()
                if rdkit_mol is None:
                    return None

                if self.ignore_hydrogens:
                    try:
                        rdkit_mol = Chem.RemoveHs(rdkit_mol, sanitize=False)
                        Chem.SanitizeMol(
                            rdkit_mol,
                            sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL
                            ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE,
                        )
                    except Exception as e:
                        logger.warning(
                            f"Failed to remove hydrogens for {self.fingerprint_type}: {e}. "
                            "Using the original RDKit molecule."
                        )
                        rdkit_mol = mol.to_rdkit()

                return rdkit_mol

            xyz_lines = [str(mol.num_atoms), ""]
            xyz_lines.extend(
                f"{symbol} {x:.10f} {y:.10f} {z:.10f}"
                for symbol, (x, y, z) in zip(
                    mol.chemical_symbols, mol.positions
                )
            )
            xyz_block = "\n".join(xyz_lines) + "\n"

            rdkit_mol = Chem.MolFromXYZBlock(xyz_block)
            if rdkit_mol is None:
                return None

            charge = mol.charge
            if charge is None:
                charge = 0
            charge = int(charge)

            # 2D/topological fingerprints require a complete, sanitized
            # chemical graph perceived by RDKit from the input geometry.
            try:
                rdDetermineBonds.DetermineBonds(
                    rdkit_mol,
                    charge=charge,
                )
            except Exception as first_error:
                if not rdDetermineBonds.hueckelEnabled():
                    raise first_error

                rdkit_mol = Chem.MolFromXYZBlock(xyz_block)
                if rdkit_mol is None:
                    return None
                rdDetermineBonds.DetermineBonds(
                    rdkit_mol,
                    useHueckel=True,
                    charge=charge,
                )

            Chem.SanitizeMol(rdkit_mol)

            if self.ignore_hydrogens:
                rdkit_mol = Chem.RemoveHs(rdkit_mol, sanitize=False)

            return rdkit_mol
        except Exception as e:
            logger.warning(
                f"RDKit molecule preparation failed for {self.fingerprint_type}: {e}"
            )
            return None

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
                conf = rdkit_mol.GetConformer(0)
                return np.array(
                    rdMolDescriptors.GetUSR(rdkit_mol, confId=conf.GetId())
                )
            elif self.fingerprint_type == "usrcat":
                conf = rdkit_mol.GetConformer(0)
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
        Group molecules based on fingerprint or 3D-descriptor similarity using
        hierarchical complete-linkage clustering.

        Binary/topological fingerprints are compared with Tanimoto similarity,
        whereas USR/USRCAT descriptors are compared with RDKit's native USR
        similarity score. Similarities are converted internally to distances
        (1 - similarity) before hierarchical complete-linkage clustering.

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
        logger.info(
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
        valid_original_indices = [
            self.rdkit_molecule_indices[i] for i in valid_indices
        ]

        fingerprint_skipped_indices = [
            self.rdkit_molecule_indices[i]
            for i, fp in enumerate(fingerprints)
            if fp is None
        ]
        self._matrix_skipped_indices = sorted(
            set(self._tanimoto_skipped_indices)
            | set(fingerprint_skipped_indices)
        )
        if self.conformer_ids is not None:
            self._matrix_skipped_ids = [
                self.conformer_ids[idx] for idx in self._matrix_skipped_indices
            ]
        else:
            self._matrix_skipped_ids = [
                str(idx + 1) for idx in self._matrix_skipped_indices
            ]

        num_valid = len(valid_indices)

        if num_valid == 0:
            groups = []
            index_groups = []
            self._cached_groups = groups
            self._cached_group_indices = index_groups

            grouping_time = time.time() - grouping_start_time
            self.record(
                tanimoto_matrix=np.empty((0, 0), dtype=np.float32),
                valid_indices=[],
                grouping_time=grouping_time,
                groups=groups,
                index_groups=index_groups,
            )
            return groups, index_groups

        logger.info(
            f"[{self.__class__.__name__}] Computing pairwise similarities for {num_valid} valid molecules"
        )

        # Compute similarity matrix
        similarity_matrix = np.zeros((num_valid, num_valid), dtype=np.float32)

        pairs = [
            (i, j) for i in range(num_valid) for j in range(i + 1, num_valid)
        ]
        total_pairs = len(pairs)
        similarities = []
        next_progress = 10

        if self.fingerprint_type in {"usr", "usrcat"}:
            # USR/USRCAT are continuous 3D descriptors and should be compared
            # with RDKit's native USR similarity score, not Tanimoto.
            with ThreadPool(self.num_procs) as pool:
                for idx, sim in enumerate(
                    pool.imap(
                        lambda p: rdMolDescriptors.GetUSRScore(
                            valid_fps[p[0]], valid_fps[p[1]]
                        ),
                        pairs,
                    )
                ):
                    similarities.append(sim)
                    next_progress = self._log_progress(
                        idx + 1, total_pairs, next_progress
                    )
        else:
            # Binary/topological fingerprints use Tanimoto similarity.
            with ThreadPool(self.num_procs) as pool:
                for idx, sim in enumerate(
                    pool.imap(
                        lambda p: DataStructs.FingerprintSimilarity(
                            valid_fps[p[0]], valid_fps[p[1]]
                        ),
                        pairs,
                    )
                ):
                    similarities.append(sim)
                    next_progress = self._log_progress(
                        idx + 1, total_pairs, next_progress
                    )

        # Fill similarity matrix
        for (i, j), sim in zip(pairs, similarities):
            similarity_matrix[i, j] = similarity_matrix[j, i] = sim

        np.fill_diagonal(similarity_matrix, 1.0)

        # Choose grouping strategy based on parameters
        full_similarity_matrix = np.full(
            (n, n), np.nan, dtype=similarity_matrix.dtype
        )
        full_similarity_matrix[
            np.ix_(valid_original_indices, valid_original_indices)
        ] = similarity_matrix

        if self.num_groups is not None:
            groups, index_groups = self.group_by_num_groups(
                full_similarity_matrix
            )
        else:
            groups, index_groups = self.group_by_threshold(
                full_similarity_matrix
            )

        # Calculate total grouping time
        grouping_end_time = time.time()
        grouping_time = grouping_end_time - grouping_start_time

        # Cache results
        self._cached_groups = groups
        self._cached_group_indices = index_groups

        # Save full matrix through unified record entrypoint
        self.record(
            tanimoto_matrix=similarity_matrix,
            valid_indices=valid_original_indices,
            grouping_time=grouping_time,
            groups=groups,
            index_groups=index_groups,
        )

        return groups, index_groups

    def _prepare_similarity_output_matrix(
        self,
        matrix: np.ndarray,
    ) -> Tuple[np.ndarray, List[int], List[int], np.dtype]:
        """Validate full recorded similarity matrix and convert to distance."""
        similarity_matrix = np.asarray(matrix)
        if similarity_matrix.ndim != 2:
            raise ValueError(
                "Tanimoto similarity matrix must be 2-dimensional."
            )
        if similarity_matrix.shape != (
            len(self.molecules),
            len(self.molecules),
        ):
            raise ValueError(
                "Tanimoto similarity matrix dimensions must match the number "
                "of molecules."
            )
        if not np.issubdtype(similarity_matrix.dtype, np.number):
            raise ValueError("Tanimoto similarity matrix must be numeric.")

        valid_mask = np.isfinite(np.diag(similarity_matrix))
        valid_indices = np.flatnonzero(valid_mask).tolist()
        skipped_indices = np.flatnonzero(~valid_mask).tolist()

        for index in skipped_indices:
            if not (
                np.all(np.isnan(similarity_matrix[index, :]))
                and np.all(np.isnan(similarity_matrix[:, index]))
            ):
                raise ValueError(
                    "Skipped Tanimoto entries must be represented by complete "
                    "NaN rows and columns."
                )

        if not valid_indices:
            return (
                np.zeros((0, 0), dtype=float),
                [],
                skipped_indices,
                np.dtype(np.float64),
            )

        valid_similarity_matrix = similarity_matrix[
            np.ix_(valid_indices, valid_indices)
        ]
        if not np.all(np.isfinite(valid_similarity_matrix)):
            raise ValueError(
                "Valid Tanimoto similarity entries must all be finite."
            )
        if not np.allclose(valid_similarity_matrix, valid_similarity_matrix.T):
            raise ValueError("Tanimoto similarity matrix must be symmetric.")
        if not np.allclose(np.diag(valid_similarity_matrix), 1.0):
            raise ValueError(
                "Tanimoto similarity matrix diagonal must be one."
            )
        if np.any(valid_similarity_matrix < 0.0) or np.any(
            valid_similarity_matrix > 1.0
        ):
            raise ValueError(
                "Tanimoto similarity values must be between 0.0 and 1.0."
            )

        dtype = valid_similarity_matrix.dtype
        one = dtype.type(1.0)
        distance_matrix = one - valid_similarity_matrix
        np.fill_diagonal(distance_matrix, dtype.type(0.0))

        return distance_matrix, valid_indices, skipped_indices, np.dtype(dtype)

    def _normalize_output_matrix_for_threshold(
        self,
        matrix: np.ndarray,
        threshold: Optional[float],
    ) -> Tuple[np.ndarray, float, Optional[List[int]], List[int]]:
        """Normalize a recorded Tanimoto similarity matrix for clustering.

        The matrix must use the recorded full-size format: valid similarities
        are finite with a diagonal of one, while skipped molecules are
        represented by complete NaN rows and columns. ``threshold`` is a
        similarity threshold in the inclusive range [0, 1].
        """
        (
            distance_matrix,
            valid_indices,
            skipped_indices,
            distance_dtype,
        ) = self._prepare_similarity_output_matrix(matrix)

        effective_threshold = self._resolve_threshold(threshold)
        if effective_threshold > 1.0:
            raise ValueError(
                "Tanimoto similarity threshold must be between 0.0 and 1.0."
            )

        one = distance_dtype.type(1.0)
        distance_threshold = float(
            one - distance_dtype.type(effective_threshold)
        )

        return (
            distance_matrix,
            distance_threshold,
            valid_indices,
            skipped_indices,
        )

    def _normalize_output_matrix_for_num_groups(
        self,
        matrix: np.ndarray,
        num_groups: Optional[int],
    ) -> Tuple[np.ndarray, int, Optional[List[int]], List[int]]:
        """Normalize a recorded Tanimoto similarity matrix for -N grouping."""
        (
            distance_matrix,
            valid_indices,
            skipped_indices,
            distance_dtype,
        ) = self._prepare_similarity_output_matrix(matrix)
        self._num_groups_distance_dtype = distance_dtype

        return (
            distance_matrix,
            self._resolve_num_groups(num_groups),
            valid_indices,
            skipped_indices,
        )

    def _postprocess_auto_threshold_for_num_groups(self) -> None:
        """Convert auto threshold from distance back to similarity."""
        if self._auto_threshold is None:
            return
        dtype = getattr(self, "_num_groups_distance_dtype", np.float64)
        one = dtype.type(1.0)
        self._auto_threshold = float(one - dtype.type(self._auto_threshold))

    def _record_results(
        self,
        tanimoto_matrix: np.ndarray,
        valid_indices: List[int],
        grouping_time: float = None,
        groups: List[List[Molecule]] = None,
        index_groups: List[List[int]] = None,
    ):
        """Strategy-specific result writer for Tanimoto grouping."""
        n = len(self.molecules)

        # Create full matrix with invalid molecules marked as NaN
        full_matrix = np.full((n, n), np.nan)
        for i, idx_i in enumerate(valid_indices):
            for j, idx_j in enumerate(valid_indices):
                full_matrix[idx_i, idx_j] = tanimoto_matrix[i, j]

        # Use ResultsRecorder to save
        recorder = self._get_results_recorder()
        labels = recorder.get_labels(n)

        # Build header info
        header_info = [
            (
                "",
                f"Full Tanimoto Similarity Matrix ({n}x{n}) - {self.__class__.__name__}",
            ),
            ("Fingerprint Type", self.fingerprint_type),
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
                        f"{self._auto_threshold:.7f} (A value closer to 1 indicates greater similarity)",
                    )
                )
        else:
            header_info.append(
                (
                    "Threshold",
                    f"{self.threshold} (A value closer to 1 indicates greater similarity)",
                )
            )

        header_info.append(("Energy Type", self.energy_type))
        header_info.append(("Ignore Hydrogens", self.ignore_hydrogens))
        header_info.append(("Num Procs", self.num_procs))

        if grouping_time is not None:
            header_info.append(
                ("Grouping Time", f"{grouping_time:.2f} seconds")
            )

        self._append_input_usage_header(header_info)
        self._append_thermo_header(header_info)

        # Build sheets data using recorder's method
        sheets_data = {}
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
            matrix_data=("Tanimoto_Matrix", full_matrix, labels),
            suffix=suffix,
            startrow=len(header_info) + 2,
        )

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
