"""Molecular structure grouping using different grouping strategies."""

import logging
import multiprocessing
from abc import ABC, abstractmethod
from functools import partial
from typing import Dict, Iterable, List, Set, Tuple, Optional

import networkx as nx
import numpy as np
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core import Structure
from rdkit import Chem
from rdkit.Chem import DataStructs, rdMolHash
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components, reverse_cuthill_mckee

from chemsmart.io.molecules.structure import Molecule

logger = logging.getLogger(__name__)


class StructureGrouperConfig:
    """Configuration container for structure matching parameters."""

    def __init__(
        self, ltol: float = 0.1, stol: float = 0.18, angle_tol: float = 1
    ):
        self.ltol = ltol
        self.stol = stol
        self.angle_tol = angle_tol


class MoleculeGrouper(ABC):
    """Abstract base class for molecular structure grouping."""

    def __init__(self, molecules: Iterable[Molecule], num_procs: int = 1):
        self.molecules = list(molecules)
        self.num_procs = max(1, num_procs)
        self._validate_inputs()

    def _validate_inputs(self) -> None:
        """Validate input molecules."""
        if not isinstance(self.molecules, Iterable):
            raise TypeError("Molecules must be an iterable collection")
        if not all(isinstance(m, Molecule) for m in self.molecules):
            raise TypeError("All items must be Molecule instances")

    @abstractmethod
    def group(self) -> Tuple[List[List[Molecule]], List[List[int]]]:
        """Main grouping method to return grouped molecules and their indices."""
        pass

    def unique(self) -> List[Molecule]:
        """Get unique representative molecules."""
        groups, _ = self.group()
        return [group[0] for group in groups]


class RDKitFingerprintGrouper(MoleculeGrouper):
    """Groups molecules using RDKit fingerprint similarity."""

    def __init__(
        self,
        molecules: Iterable[Molecule],
        num_procs: int = 1,
        threshold: float = 0.8,
    ):
        super().__init__(molecules, num_procs)
        self.threshold = threshold

    def _get_fingerprint(
        self, mol: Molecule
    ) -> Optional[DataStructs.ExplicitBitVect]:
        """Generate RDKit fingerprint for a molecule."""
        try:
            rdkit_mol = mol.to_rdkit()
            return Chem.RDKFingerprint(rdkit_mol) if rdkit_mol else None
        except Exception as e:
            logger.warning(f"Fingerprint generation failed: {str(e)}")
            return None

    def group(self) -> Tuple[List[List[Molecule]], List[List[int]]]:
        """Group molecules based on fingerprint similarity."""
        with multiprocessing.Pool(self.num_procs) as pool:
            fps = [
                fp
                for fp in pool.map(self._get_fingerprint, self.molecules)
                if fp
            ]

        n = len(fps)
        similarity_matrix = np.zeros((n, n))
        range_indices = [(i, j) for i in range(n) for j in range(i + 1, n)]

        with multiprocessing.Pool(self.num_procs) as pool:
            results = pool.starmap(
                DataStructs.FingerprintSimilarity,
                [(fps[i], fps[j]) for i, j in range_indices],
            )

        for (i, j), sim in zip(range_indices, results):
            similarity_matrix[i][j] = similarity_matrix[j][i] = sim

        # Find connected components
        n_components, labels = connected_components(
            csr_matrix(similarity_matrix > self.threshold)
        )
        groups: List[List[Molecule]] = []
        indices: List[List[int]] = []
        for i in range(n_components):
            component = [self.molecules[j] for j in np.where(labels == i)[0]]
            groups.append(component)
            indices.append(list(np.where(labels == i)[0]))

        return groups, indices


class RDKitIsomorphismGrouper(MoleculeGrouper):
    """Group molecules using RDKit-based hashing and isomorphism checks."""

    def __init__(
        self,
        molecules: Iterable[Molecule],
        num_procs: int = 1,
        use_stereochemistry: bool = True,
        use_tautomers: bool = False,
    ):
        super().__init__(molecules, num_procs)
        self.use_stereochemistry = use_stereochemistry
        self.use_tautomers = use_tautomers

    def _get_mol_hash(self, mol: Molecule) -> Optional[str]:
        """Generate canonical hash for molecule."""
        try:
            rdkit_mol = mol.to_rdkit()
            if not rdkit_mol:
                return None

            hash_func = (
                rdMolHash.HashFunction.Tautomer
                if self.use_tautomers
                else (
                    rdMolHash.HashFunction.AnonymousGraph
                    if self.use_stereochemistry
                    else rdMolHash.HashFunction.MolFormula
                )
            )
            return rdMolHash.MolHash(rdkit_mol, hash_func)
        except Exception as e:
            logger.warning(f"Hash generation failed: {str(e)}")
            return None

    def group(self) -> Tuple[List[List[Molecule]], List[List[int]]]:
        """Group molecules using structural isomorphism."""
        with multiprocessing.Pool(self.num_procs) as pool:
            hashes = pool.map(self._get_mol_hash, self.molecules)

        groups: List[List[Molecule]] = []
        indices = list(range(len(self.molecules)))
        index_groups: List[List[int]] = []

        while indices:
            pivot_idx = indices[0]
            current_hash = hashes[pivot_idx]
            matches = [i for i in indices if hashes[i] == current_hash]

            # Verify isomorphism if needed
            if self.use_stereochemistry or self.use_tautomers:
                matches = [
                    i
                    for i in matches
                    if self._check_isomorphism(
                        self.molecules[pivot_idx], self.molecules[i]
                    )
                ]

            groups.append([self.molecules[i] for i in matches])
            index_groups.append(matches)
            indices = [i for i in indices if i not in matches]

        return groups, index_groups

    def _check_isomorphism(self, mol1: Molecule, mol2: Molecule) -> bool:
        """Check graph isomorphism considering stereochemistry."""
        try:
            return Chem.MolToInchiKey(mol1.to_rdkit()) == Chem.MolToInchiKey(
                mol2.to_rdkit()
            )
        except Exception as e:
            logger.warning(f"Isomorphism check failed: {str(e)}")
            return False


class RCMSimilarityGrouper(MoleculeGrouper):
    """Group molecules using Reverse Cuthill-McKee algorithm on similarity matrix."""

    def __init__(
        self,
        molecules: Iterable[Molecule],
        similarity_threshold: float = 0.9,
        num_procs: int = 1,
    ):
        super().__init__(molecules, num_procs)
        self.similarity_threshold = similarity_threshold

    def group(self) -> Tuple[List[List[Molecule]], List[List[int]]]:
        """Cluster molecules using RCM on similarity matrix."""
        with multiprocessing.Pool(self.num_procs) as pool:
            fps = pool.map(
                Chem.RDKFingerprint, (m.to_rdkit() for m in self.molecules)
            )

        n = len(self.molecules)
        similarity_matrix = np.zeros((n, n))
        indices = [(i, j) for i in range(n) for j in range(i, n)]

        with multiprocessing.Pool(self.num_procs) as pool:
            similarities = pool.starmap(
                DataStructs.FingerprintSimilarity,
                [(fps[i], fps[j]) for i, j in indices],
            )

        for (i, j), sim in zip(indices, similarities):
            similarity_matrix[i, j] = similarity_matrix[j, i] = sim

        # Apply RCM clustering
        sparse_matrix = csr_matrix(
            similarity_matrix > self.similarity_threshold
        )
        perm = reverse_cuthill_mckee(sparse_matrix, symmetric_mode=True)

        # Extract groups from permutation
        groups: List[List[Molecule]] = []
        index_groups: List[List[int]] = []
        current_group: List[int] = [perm[0]]

        for idx in perm[1:]:
            if (
                similarity_matrix[idx, current_group[0]]
                > self.similarity_threshold
            ):
                current_group.append(idx)
            else:
                groups.append([self.molecules[i] for i in current_group])
                index_groups.append(current_group)
                current_group = [idx]

        if current_group:
            groups.append([self.molecules[i] for i in current_group])
            index_groups.append(current_group)

        return groups, index_groups


class RMSDMoleculeGrouper(MoleculeGrouper):
    """Group molecules based on RMSD of atomic positions."""

    def __init__(
        self,
        molecules: Iterable[Molecule],
        rmsd_threshold: float = 0.5,
        num_procs: int = 1,
        align_molecules: bool = True,
    ):
        super().__init__(molecules, num_procs)
        self.rmsd_threshold = rmsd_threshold
        self.align_molecules = align_molecules

    def _kabsch_align(self, P: np.ndarray, Q: np.ndarray) -> np.ndarray:
        """Kabsch algorithm for optimal molecular alignment."""
        centroid_P = np.mean(P, axis=0)
        centroid_Q = np.mean(Q, axis=0)

        P_centered = P - centroid_P
        Q_centered = Q - centroid_Q

        H = P_centered.T @ Q_centered
        U, _, Vt = np.linalg.svd(H)
        d = np.sign(np.linalg.det(Vt.T @ U.T))
        R = Vt.T @ np.diag([1, 1, d]) @ U.T

        return (P_centered @ R) + centroid_Q

    def _calculate_rmsd(self, idx_pair: Tuple[int, int]) -> float:
        """Calculate RMSD between two molecules."""
        i, j = idx_pair
        mol1, mol2 = self.molecules[i], self.molecules[j]

        if mol1.num_atoms != mol2.num_atoms or sorted(
            mol1.chemical_symbols
        ) != sorted(mol2.chemical_symbols):
            return np.inf

        pos1 = mol1.positions
        pos2 = mol2.positions

        if self.align_molecules:
            pos1 = self._kabsch_align(pos1, pos2)

        return np.sqrt(np.mean(np.sum((pos1 - pos2) ** 2, axis=1)))

    def group(self) -> Tuple[List[List[Molecule]], List[List[int]]]:
        """Group molecules by geometric similarity."""
        n = len(self.molecules)
        indices = [(i, j) for i in range(n) for j in range(i + 1, n)]

        with multiprocessing.Pool(self.num_procs) as pool:
            rmsd_values = pool.map(self._calculate_rmsd, indices)

        # Build adjacency matrix
        adj_matrix = np.zeros((n, n), dtype=bool)
        for (i, j), rmsd in zip(indices, rmsd_values):
            if rmsd < self.rmsd_threshold:
                adj_matrix[i, j] = adj_matrix[j, i] = True

        # Find connected components
        _, labels = connected_components(csr_matrix(adj_matrix))
        groups: List[List[Molecule]] = []
        index_groups: List[List[int]] = []
        for label in np.unique(labels):
            mask = labels == label
            groups.append([self.molecules[i] for i in np.where(mask)[0]])
            index_groups.append(list(np.where(mask)[0]))

        return groups, index_groups


class PymatgenMoleculeGrouper(MoleculeGrouper):
    """Group molecules using pymatgen's StructureMatcher."""

    def __init__(
        self,
        molecules: Iterable[Molecule],
        config: StructureGrouperConfig = StructureGrouperConfig(),
        num_procs: int = 1,
    ):
        super().__init__(molecules, num_procs)
        self.config = config
        self.structures = [m.to_pymatgen() for m in self.molecules]

    def group(self) -> Tuple[List[List[Molecule]], List[List[int]]]:
        """Group using pymatgen's structure matching."""
        matcher = StructureMatcher(
            ltol=self.config.ltol,
            stol=self.config.stol,
            angle_tol=self.config.angle_tol,
            scale=False,
            primitive_cell=False,
            attempt_supercell=False,
        )

        groups = matcher.group_structures(self.structures)
        mol_groups = []
        index_groups = []
        structure_indices = {id(s): i for i, s in enumerate(self.structures)}

        for group in groups:
            indices = [structure_indices[id(s)] for s in group]
            mol_groups.append([self.molecules[i] for i in indices])
            index_groups.append(indices)

        return mol_groups, index_groups


class StructureGrouperFactory:
    """Factory for creating molecular grouping strategies."""

    _GROUPER_MAP = {
        "fingerprint": RDKitFingerprintGrouper,
        "isomorphism": RDKitIsomorphismGrouper,
        "rcm_similarity": RCMSimilarityGrouper,
        "rmsd": RMSDMoleculeGrouper,
        "pymatgen": PymatgenMoleculeGrouper,
    }

    @classmethod
    def create_grouper(
        cls,
        molecules: Iterable[Molecule],
        strategy: str = "fingerprint",
        strategy_params: Optional[Dict] = None,
        num_procs: int = 1,
    ) -> MoleculeGrouper:
        """
        Create a molecular grouper with specified strategy.

        Parameters:
            strategy: One of 'fingerprint', 'isomorphism', 'rcm_similarity', 'rmsd', 'pymatgen'
            strategy_params: Dictionary of parameters for the selected strategy
            num_procs: Number of processes for parallel processing
        """
        strategy_params = strategy_params or {}

        if strategy not in cls._GROUPER_MAP:
            valid = ", ".join(cls._GROUPER_MAP.keys())
            raise ValueError(
                f"Invalid strategy: {strategy}. Valid options: {valid}"
            )

        return cls._GROUPER_MAP[strategy](
            molecules=molecules, num_procs=num_procs, **strategy_params
        )
