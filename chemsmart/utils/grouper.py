"""Molecular structure grouping using different grouping strategies."""

import logging
import multiprocessing
from abc import ABC, abstractmethod
from multiprocessing.pool import ThreadPool
from typing import Iterable, List, Optional, Tuple

import networkx as nx
import numpy as np
from rdkit import Chem
from rdkit.Chem import DataStructs, rdMolHash
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components, reverse_cuthill_mckee

from chemsmart.io.molecules.structure import Molecule

logger = logging.getLogger(__name__)


class StructureGrouperConfig:
    """Configuration container for parameters for pymatgen StructureMatcher.
    Default values are for heterogenous systems, may need to check for molecules.
    """

    def __init__(self, ltol=0.1, stol=0.18, angle_tol=1):
        self.ltol = ltol
        self.stol = stol
        self.angle_tol = angle_tol


class MoleculeGrouper(ABC):
    """Abstract base class for molecular structure grouping.
    Specific type of base class that cannot be directly instantiated and
    designed to define a common interface that subclasses must implement"""

    def __init__(self, molecules: Iterable[Molecule], num_procs: int = 1):
        self.molecules = molecules
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
        self.threshold = threshold  # Similarity threshold

        # Convert Molecules to RDKit molecules
        self.molecules = [
            mol.to_rdkit() for mol in self.molecules if mol.to_rdkit()
        ]

    @staticmethod
    def _get_fingerprint(
        mol: Molecule,
    ) -> Optional[DataStructs.ExplicitBitVect]:
        """Generate RDKit fingerprint for a molecule."""
        try:
            rdkit_mol = mol.to_rdkit()
            return Chem.RDKFingerprint(rdkit_mol) if rdkit_mol else None
        except Exception as e:
            logger.warning(f"Fingerprint generation failed: {str(e)}")
            return None

    @staticmethod
    def _compute_similarity(idx1, idx2, fingerprints):
        """Compute similarity between two fingerprints.
        Uses TanimotoSimilarity metric from RDKit by default."""
        return (
            idx1,
            idx2,
            DataStructs.FingerprintSimilarity(
                fingerprints[idx1], fingerprints[idx2]
            ),
        )

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
                self._compute_similarity,
                [(i, j, fps) for i, j in range_indices],
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
    """Group molecules using RDKit-based hashing and isomorphism checks"""

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
        """Generate canonical hash for molecules"""
        try:
            rdkit_mol = mol.to_rdkit()
            if not rdkit_mol:
                return None
            # Choose hashing function based on requirements
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
    """Group molecules using Reverse Cuthill-McKee algorithm on similarity matrix.
    Utilize molecular fingerprinting to focus on connectivity and similarity,
    suitable when detailed structural matching is needed."""

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
        with ThreadPool(self.num_procs) as pool:
            # Use ThreadPool instead of Pool since the latter
            # cannot pickle non-python objects
            fps = pool.map(
                Chem.RDKFingerprint, (m.to_rdkit() for m in self.molecules)
            )

        n = len(self.molecules)
        similarity_matrix = np.zeros((n, n))
        indices = [(i, j) for i in range(n) for j in range(i, n)]

        with ThreadPool(self.num_procs) as pool:
            similarities = pool.starmap(
                Chem.DataStructs.TanimotoSimilarity,
                [(fps[idx[0]], fps[idx[1]]) for idx in indices],
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

    def _tanimoto_similarity(self, fp1, fp2) -> float:
        """Calculate Tanimoto similarity between fingerprints"""
        return Chem.DataStructs.TanimotoSimilarity(fp1, fp2)


class RMSDGrouper(MoleculeGrouper):
    """Group molecules based on RMSD (Root Mean Square Deviation) of atomic positions.
    Effective for precise 3D comparisons, ideal in contexts like crystallography or drug
    binding where exact spatial alignment is crucial."""

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

    def _check_similarity(self, mol1: Molecule, mol2: Molecule) -> bool:
        """Check if two molecules are geometrically similar"""
        if mol1.num_atoms != mol2.num_atoms:
            return False

        if sorted(mol1.chemical_symbols) != sorted(mol2.chemical_symbols):
            return False

        return self._calculate_rmsd(mol1, mol2) < self.rmsd_threshold

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
            logger.info("Aligning molecules using Kabsch algorithm.")
            pos1, pos2, _, _, _ = self._kabsch_align(pos1, pos2)

        return np.sqrt(np.mean(np.sum((pos1 - pos2) ** 2, axis=1)))

    def _kabsch_align(
        self, P: np.ndarray, Q: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Kabsch algorithm for molecular alignment"""
        # Center molecules
        assert P.shape == Q.shape, "Matrix dimensions must match"

        # Compute centroids
        centroid_P = np.mean(P, axis=0)
        centroid_Q = np.mean(Q, axis=0)

        # Optimal translation
        t = centroid_Q - centroid_P

        # Center the points
        p = P - centroid_P
        q = Q - centroid_Q

        # Compute the covariance matrix
        H = np.dot(p.T, q)

        # SVD
        U, S, Vt = np.linalg.svd(H)

        # Validate right-handed coordinate system
        if np.linalg.det(np.dot(Vt.T, U.T)) < 0.0:
            Vt[-1, :] *= -1.0

        # Optimal rotation
        R = np.dot(Vt.T, U.T)

        # rotate p
        p = np.dot(p, R.T)
        # RMSD
        rmsd = np.sqrt(np.sum(np.square(np.dot(p, R.T) - q)) / P.shape[0])

        return p, q, R, t, rmsd


class HybridMoleculeGrouper(MoleculeGrouper):
    """Hybrid grouping strategy combining multiple approaches.
    Combines geometric and substructure similarities, useful for comprehensive analysis
    requiring both bond arrangement and overall structure."""

    def __init__(
        self,
        molecules: Iterable[Molecule],
        num_procs: int = 1,
        strategies: List[Tuple[str, dict]] = None,
    ):
        super().__init__(molecules, num_procs)
        self.strategies = strategies or [
            ("formula", {}),
            ("connectivity", {"bond_cutoff": 1.5}),
            ("geometry", {"rmsd_threshold": 0.5}),
        ]

    def group(self) -> Tuple[List[List[Molecule]], List[List[int]]]:
        current_groups = [[m] for m in self.molecules]

        for strategy, params in self.strategies:
            current_groups = self._apply_strategy(
                current_groups, strategy, params
            )

        return current_groups, [list(range(len(g))) for g in current_groups]

    def _apply_strategy(
        self, groups: List[List["Molecule"]], strategy: str, params: dict
    ) -> List[List["Molecule"]]:
        new_groups = []
        for group in groups:
            if len(group) == 1:
                new_groups.append(group)
                continue

            if strategy == "formula":
                subgrouper = FormulaGrouper(group, **params)
            elif strategy == "connectivity":
                subgrouper = ConnectivityGrouper(group, **params)
            elif strategy == "rmsd":
                subgrouper = RMSDGrouper(group, **params)
            else:
                raise ValueError(f"Unknown strategy: {strategy}")

            subgroups, _ = subgrouper.group()
            new_groups.extend(subgroups)

        return new_groups


class FormulaGrouper(MoleculeGrouper):
    """Group by chemical formula.
    Ideal for grouping molecules based solely on their chemical formula, making it suitable
    when elemental composition is the primary concern."""

    def group(self):
        formula_groups = {}
        for idx, mol in enumerate(self.molecules):
            formula = mol.get_chemical_formula()
            if formula not in formula_groups:
                formula_groups[formula] = []
            formula_groups[formula].append((mol, idx))

        mol_groups = []
        idx_groups = []
        for formula, group in formula_groups.items():
            mols, indices = zip(*group)
            mol_groups.append(list(mols))
            idx_groups.append(list(indices))

        return mol_groups, idx_groups


class ConnectivityGrouper(MoleculeGrouper):
    """Group by molecular connectivity.
    Useful for identifying molecules with similar bond arrangements and structures, useful
    in scenarios where molecular connectivity is critical."""

    def __init__(self, molecules: List[Molecule], bond_cutoff: float = 1.5):
        super().__init__(molecules)
        self.bond_cutoff = bond_cutoff

    def group(self):
        groups = []
        remaining = list(enumerate(self.molecules))

        while remaining:
            pivot_idx, pivot_mol = remaining.pop(0)
            current_group = [pivot_mol]
            current_indices = [pivot_idx]

            pivot_graph = pivot_mol.to_graph(
                bond_cutoff_buffer=self.bond_cutoff
            )

            to_remove = []
            for i, (idx, mol) in enumerate(remaining):
                mol_graph = mol.to_graph(bond_cutoff_buffer=self.bond_cutoff)
                if self._are_isomorphic(pivot_graph, mol_graph):
                    current_group.append(mol)
                    current_indices.append(idx)
                    to_remove.append(i)

            for i in reversed(to_remove):
                remaining.pop(i)

            groups.append((current_group, current_indices))

        mol_groups = [g[0] for g in groups]
        idx_groups = [g[1] for g in groups]
        return mol_groups, idx_groups

    def _are_isomorphic(self, g1: nx.Graph, g2: nx.Graph) -> bool:
        """Check if two molecular graphs are isomorphic"""
        return nx.is_isomorphic(
            g1,
            g2,
            node_match=lambda a, b: a["element"] == b["element"],
            edge_match=lambda a, b: a["bond_order"] == b["bond_order"],
        )


class StructureGrouperFactory:
    @staticmethod
    def create_grouper(structures, strategy="rdkit", num_procs=1):
        groupers = {
            "fingerprint": RDKitFingerprintGrouper,
            "isomorphism": RDKitIsomorphismGrouper,
            "rcm": RCMSimilarityGrouper,
            "rmsd": RMSDGrouper,
            "hybrid": HybridMoleculeGrouper,
            "formula": FormulaGrouper,
            "connectivity": ConnectivityGrouper,
        }
        if strategy in groupers:
            logger.info(f"Using {strategy} grouping strategy.")
            return groupers[strategy](structures, num_procs)
        raise ValueError(f"Unknown grouping strategy: {strategy}")
