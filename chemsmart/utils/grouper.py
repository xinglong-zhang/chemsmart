"""Molecular structure grouping using different grouping strategies."""

import logging
import multiprocessing
from abc import ABC, abstractmethod
from typing import Iterable, List, Tuple

import networkx as nx
import numpy as np
from pymatgen.analysis.structure_matcher import StructureMatcher
from rdkit import Chem
from rdkit.Chem import DataStructs, rdMolHash
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import reverse_cuthill_mckee

from chemsmart.io.molecules.structure import Molecule

logger = logging.getLogger(__name__)


class StructureGrouperConfig:
    def __init__(self, ltol=0.1, stol=0.18, angle_tol=1):
        self.ltol = ltol
        self.stol = stol
        self.angle_tol = angle_tol


class MoleculeGrouper(ABC):
    """Abstract base class for molecular structure grouping.
    Specific type of base class that cannot be directly instantiated and
    designed to define a common interface that subclasses must implement"""

    def __init__(self, molecules, num_procs=1):
        self.molecules = molecules
        self.num_procs = max(1, num_procs)

        self._validate_inputs()

    def _validate_inputs(self):
        if not isinstance(self.molecules, Iterable):
            raise TypeError("Molecules must be an iterable collection")

    @abstractmethod
    def group(self):
        pass

    def unique(self) -> List:
        groups, _ = self.group()
        return [group[0] for group in groups]


class RDKitFingerprintGrouper(MoleculeGrouper):
    """Groups molecules using RDKit fingerprint similarity."""

    def __init__(
        self,
        molecules: List["Molecule"],
        num_procs: int = 1,
        threshold: float = 0.8,
    ):
        super().__init__(molecules, num_procs)
        self.threshold = threshold  # Similarity threshold

        # Convert Molecules to RDKit molecules
        self.molecules = [
            mol.to_rdkit() for mol in self.molecules if mol.to_rdkit()
        ]

    def _compute_similarity(self, idx1, idx2, fingerprints):
        return (
            idx1,
            idx2,
            DataStructs.FingerprintSimilarity(
                fingerprints[idx1], fingerprints[idx2]
            ),
        )

    def group(self):
        fingerprints = [Chem.RDKFingerprint(mol) for mol in self.molecules]
        indices = [
            (i, j)
            for i in range(len(fingerprints))
            for j in range(i + 1, len(fingerprints))
        ]

        with multiprocessing.Pool(self.num_procs) as pool:
            similarities = pool.starmap(
                self._compute_similarity,
                [(i, j, fingerprints) for i, j in indices],
            )

        G = nx.Graph()
        G.add_nodes_from(range(len(self.molecules)))
        for i, j, sim in similarities:
            if sim > self.threshold:
                G.add_edge(i, j)

        return list(nx.connected_components(G))

    def unique(self) -> List:
        # overrides parent method
        return self.group()


class RDKitIsomorphismGrouper(MoleculeGrouper):
    """Group molecules using RDKit-based hashing and isomorphism checks"""

    def __init__(
        self,
        molecules: List["Molecule"],
        num_procs: int = 1,
        use_stereochemistry: bool = True,
        use_tautomers: bool = False,
    ):
        super().__init__(molecules, num_procs)
        self.use_stereochemistry = use_stereochemistry
        self.use_tautomers = use_tautomers

    def _get_mol_hash(self, mol: Molecule) -> str:
        """Generate canonical hash for molecules"""
        rdkit_mol = mol.to_rdkit()

        # Choose hashing function based on requirements
        if self.use_tautomers:
            hash_func = rdMolHash.HashFunction.Tautomer
        elif self.use_stereochemistry:
            hash_func = rdMolHash.HashFunction.AnonymousGraph
        else:
            hash_func = rdMolHash.HashFunction.MolFormula

        return rdMolHash.MolHash(rdkit_mol, hash_func)

    def group(self) -> Tuple[List[List["Molecule"]], List[List[int]]]:
        groups = []
        indices = list(range(len(self.molecules)))
        hashes = [self._get_mol_hash(m) for m in self.molecules]

        while indices:
            pivot_idx = indices[0]
            current_hash = hashes[pivot_idx]

            # Find matches
            matches = [i for i in indices if hashes[i] == current_hash]

            # Verify isomorphism for non-formula hashes
            if self.use_stereochemistry or self.use_tautomers:
                matches = [
                    i
                    for i in matches
                    if self._check_isomorphism(
                        self.molecules[pivot_idx], self.molecules[i]
                    )
                ]

            groups.append([self.molecules[i] for i in matches])
            indices = [i for i in indices if i not in matches]

        return groups, [list(range(len(g))) for g in groups]

    def _check_isomorphism(self, mol1: Molecule, mol2: Molecule) -> bool:
        """Check graph isomorphism considering stereochemistry"""
        return Chem.MolToInchiKey(mol1.to_rdkit()) == Chem.MolToInchiKey(
            mol2.to_rdkit()
        )


class RCMAdjacencyGrouper(MoleculeGrouper):
    """Group molecules using Reverse Cuthill-McKee algorithm on adjacency matrix.
    Utilize molecular fingerprinting to focus on connectivity and similarity,
    suitable when detailed structural matching is needed."""

    def __init__(self, molecules, num_procs=1):
        super().__init__(molecules, num_procs)
        self.molecules = [
            mol.to_rdkit() for mol in self.molecules if mol.to_rdkit()
        ]

    def _get_adjacency_matrix(self, mol):
        num_atoms = mol.GetNumAtoms()
        A = np.zeros((num_atoms, num_atoms))
        for bond in mol.GetBonds():
            i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            A[i, j] = A[j, i] = 1
        return A

    def group(self):
        with multiprocessing.Pool(self.num_procs) as pool:
            adj_matrices = pool.map(self._get_adjacency_matrix, self.molecules)

        return [
            reverse_cuthill_mckee(csr_matrix(A), symmetric_mode=True)
            for A in adj_matrices
        ]

    def unique(self) -> List:
        # overrides parent method
        return self.group()


class RCMSimilarityGrouper(MoleculeGrouper):
    """Group molecules using Reverse Cuthill-McKee algorithm on similarity matrix.
    Utilize molecular fingerprinting to focus on connectivity and similarity,
    suitable when detailed structural matching is needed."""

    def __init__(
        self, molecules: List["Molecule"], similarity_threshold: float = 0.9
    ):
        super().__init__(molecules)
        self.similarity_threshold = similarity_threshold

    def group(self) -> Tuple[List[List["Molecule"]], List[List[int]]]:
        sim_matrix = self._create_similarity_matrix()
        return self._cluster_molecules(sim_matrix)

    def _create_similarity_matrix(self) -> np.ndarray:
        """Create similarity matrix using molecular fingerprints"""
        n = len(self.molecules)
        fps = [self._get_fingerprint(m) for m in self.molecules]

        matrix = np.zeros((n, n))
        for i in range(n):
            for j in range(i, n):
                similarity = self._tanimoto_similarity(fps[i], fps[j])
                matrix[i, j] = similarity
                matrix[j, i] = similarity

        return matrix > self.similarity_threshold

    def _get_fingerprint(self, mol: Molecule):
        """Get RDKit fingerprint for similarity comparison"""
        return Chem.RDKFingerprint(mol.to_rdkit())

    def _tanimoto_similarity(self, fp1, fp2) -> float:
        """Calculate Tanimoto similarity between fingerprints"""
        return Chem.DataStructs.TanimotoSimilarity(fp1, fp2)

    def _cluster_molecules(self, sim_matrix: np.ndarray):
        """Cluster using RCM algorithm"""
        sparse_matrix = csr_matrix(sim_matrix)
        perm = reverse_cuthill_mckee(sparse_matrix, symmetric_mode=True)

        # Reorder matrix and extract blocks
        reordered = sparse_matrix[perm][:, perm]
        n = reordered.shape[0]

        groups = []
        current_group = [perm[0]]
        for i in range(1, n):
            if reordered[i, current_group[0]]:
                current_group.append(perm[i])
            else:
                groups.append(sorted(current_group))
                current_group = [perm[i]]

        if current_group:
            groups.append(sorted(current_group))

        # Convert to molecule groups
        mol_groups = [[self.molecules[i] for i in group] for group in groups]
        return mol_groups, groups


class RMSDMoleculeGrouper(MoleculeGrouper):
    """Group molecules based on RMSD (Root Mean Square Deviation) of atomic positions.
    Effective for precise 3D comparisons, ideal in contexts like crystallography or drug
    binding where exact spatial alignment is crucial."""

    def __init__(
        self,
        molecules: List["Molecule"],
        num_procs: int = 1,
        rmsd_threshold: float = 0.5,
        align_molecules: bool = True,
    ):
        super().__init__(molecules, num_procs)
        self.rmsd_threshold = rmsd_threshold
        self.align_molecules = align_molecules

    def group(self) -> Tuple[List[List["Molecule"]], List[List[int]]]:
        groups = []
        remaining = list(enumerate(self.molecules))

        while remaining:
            pivot_idx, pivot_mol = remaining.pop(0)
            current_group = [pivot_mol]
            current_indices = [pivot_idx]

            to_remove = []
            for i, (idx, mol) in enumerate(remaining):
                #     with multiprocessing.Pool(self.num_procs) as pool:
                #         distances = pool.starmap(
                #             self._calculate_rmsd,
                #             [(self.molecules[i], self.molecules[j]) for i, j in indices],
                #         )
                #         print(distances)
                if self._check_similarity(pivot_mol, mol):
                    current_group.append(mol)
                    current_indices.append(idx)
                    to_remove.append(i)

            # Remove processed molecules
            for i in reversed(to_remove):
                remaining.pop(i)

            groups.append((current_group, current_indices))

        # Unpack results
        mol_groups = [g[0] for g in groups]
        idx_groups = [g[1] for g in groups]
        return mol_groups, idx_groups

    def _check_similarity(self, mol1: Molecule, mol2: Molecule) -> bool:
        """Check if two molecules are geometrically similar"""
        if mol1.num_atoms != mol2.num_atoms:
            return False

        if sorted(mol1.chemical_symbols) != sorted(mol2.chemical_symbols):
            return False

        return self._calculate_rmsd(mol1, mol2) < self.rmsd_threshold

    def _calculate_rmsd(self, mol1: Molecule, mol2: Molecule) -> float:
        """Calculate RMSD between two molecules"""
        pos1 = mol1.positions
        pos2 = mol2.positions

        if self.align_molecules:
            logger.info("Aligning molecules using Kabsch algorithm.")
            pos1, pos2, _, _, _ = self._kabsch_align(pos1, pos2)
            # print(f"pos1 {pos1}")
            # print(f"pos2 {pos2}")
            # print(f"pos1 - pos2 {(pos1 - pos2)**2}")
            # print(np.sum((pos1 - pos2) ** 2, axis=1))
            # print(np.sum(np.sum((pos1 - pos2) ** 2, axis=1))/6)
            # print(np.mean(np.sum((pos1 - pos2) ** 2, axis=1)))

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


    # def __init__(self, molecules, rmsd_threshold=0.5, num_procs=1):
    #     super().__init__(molecules, num_procs)
    #     self.rmsd_threshold = rmsd_threshold
    #
    # def _calculate_rmsd(self, mol1, mol2):
    #     """Calculate RMSD between two molecules.
    #     Apply Kabsch algorithm to align the two molecules and calculate RMSD.
    #     """
    #     # Center molecules
    #     pos1, pos2 = mol1.positions - mol1.positions.mean(
    #         axis=0
    #     ), mol2.positions - mol2.positions.mean(axis=0)
    #     # Compute covariance matrix
    #     H = pos1.T @ pos2
    #     # SVD
    #     U, _, Vt = svd(H)
    #     # Rotation matrix
    #     R = Vt.T @ U.T
    #     if np.linalg.det(R) < 0:
    #         Vt[-1, :] *= -1
    #         R = Vt.T @ U.T
    #     # Apply rotation
    #     pos1_aligned = pos1 @ R
    #     return np.sqrt(np.mean(np.sum((pos1_aligned - pos2) ** 2, axis=1)))
    #
    # def group(self):
    #     indices = [
    #         (i, j)
    #         for i in range(len(self.molecules))
    #         for j in range(i + 1, len(self.molecules))
    #     ]
    #
    #     with multiprocessing.Pool(self.num_procs) as pool:
    #         distances = pool.starmap(
    #             self._calculate_rmsd,
    #             [(self.molecules[i], self.molecules[j]) for i, j in indices],
    #         )
    #         print(distances)
    #
    #     G = nx.Graph()
    #     G.add_nodes_from(range(len(self.molecules)))
    #     for (i, j), rmsd in zip(indices, distances):
    #         if rmsd < self.rmsd_threshold:
    #             G.add_edge(i, j)
    #
    #     return list(nx.connected_components(G))
    #
    # def unique(self) -> List:
    #     # overrides parent method
    #     return self.group()


class PymatgenMoleculeGrouper(MoleculeGrouper):
    """Group molecules using pymatgen's StructureMatcher.
    Useful for analyzing crystal structures based on unit cell parameters and
    composition, relevant in material science.
    Cells are turned off herein."""

    def __init__(self, molecules, num_procs=1):
        super().__init__(molecules, num_procs)
        self.config = StructureGrouperConfig()
        self.molecules = [
            mol.to_pymatgen() for mol in self.molecules if mol.to_pymatgen()
        ]

    def _match_structure(self, idx1, idx2, matcher):
        return (
            idx1,
            idx2,
            matcher.fit(self.molecules[idx1], self.molecules[idx2]),
        )

    def group(self):
        matcher = StructureMatcher(
            ltol=self.config.ltol,
            stol=self.config.stol,
            angle_tol=self.config.angle_tol,
            scale=False,
            primitive_cell=False,
            attempt_supercell=False,
        )
        indices = [
            (i, j)
            for i in range(len(self.molecules))
            for j in range(i + 1, len(self.molecules))
        ]

        with multiprocessing.Pool(self.num_procs) as pool:
            matches = pool.starmap(
                self._match_structure, [(i, j, matcher) for i, j in indices]
            )

        grouped = []
        for i, j, is_match in matches:
            if is_match:
                for group in grouped:
                    if i in group or j in group:
                        group.update([i, j])
                        break
                else:
                    grouped.append({i, j})

        return [list(group) for group in grouped]

    def unique(self) -> List:
        # overrides parent method
        return self.group()


class HybridMoleculeGrouper(MoleculeGrouper):
    """Hybrid grouping strategy combining multiple approaches.
    Combines geometric and substructure similarities, useful for comprehensive analysis
    requiring both bond arrangement and overall structure."""

    def __init__(
        self,
        molecules: List["Molecule"],
        strategies: List[Tuple[str, dict]] = None,
    ):
        super().__init__(molecules)
        self.strategies = strategies or [
            ("formula", {}),
            ("connectivity", {"bond_cutoff": 1.5}),
            ("geometry", {"rmsd_threshold": 0.5}),
        ]

    def group(self) -> Tuple[List[List["Molecule"]], List[List[int]]]:
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
            elif strategy == "geometry":
                subgrouper = GeometryGrouper(group, **params)
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

    def __init__(self, molecules: List["Molecule"], bond_cutoff: float = 1.5):
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


class GeometryGrouper(RMSDMoleculeGrouper):
    """Specialized geometry-based grouper (inherits from RMSDMoleculeGrouper).
    Best for grouping based on 3D structural similarity, essential in applications like drug
    design where precise molecular shape matters."""

    pass


class StructureGrouperFactory:
    @staticmethod
    def create_grouper(structures, strategy="rdkit", num_procs=1):
        groupers = {
            "fingerprint": RDKitFingerprintGrouper,
            "isomorphism": RDKitIsomorphismGrouper,
            "rcm_adjacency": RCMAdjacencyGrouper,
            "rcm_similarity": RCMSimilarityGrouper,
            "rmsd": RMSDMoleculeGrouper,
            "pymatgen": PymatgenMoleculeGrouper,
            "hybrid": HybridMoleculeGrouper,
            "formula": FormulaGrouper,
            "connectivity": ConnectivityGrouper,
            "geometry": GeometryGrouper,
        }
        if strategy in groupers:
            logger.info(f"Using {strategy} grouping strategy.")
            return groupers[strategy](structures, num_procs)
        raise ValueError(f"Unknown grouping strategy: {strategy}")
