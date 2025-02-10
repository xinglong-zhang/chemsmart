"""Molecular structure grouping using different grouping strategies."""

import logging
import multiprocessing
import pickle
from abc import ABC, abstractmethod
from multiprocessing import RawArray, shared_memory
from multiprocessing.pool import ThreadPool
from typing import Iterable, List, Optional, Tuple

import networkx as nx
import numpy as np
from joblib import Parallel, delayed  # More efficient parallelization
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator, rdMolHash
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

from chemsmart.io.molecules.structure import Molecule
from chemsmart.utils.utils import kabsch_align

logger = logging.getLogger(__name__)


def to_graph_wrapper(
    mol: Molecule, bond_cutoff_buffer: float = 0.0, adjust_H: bool = True
) -> nx.Graph:
    """Global Helper function to call to_graph()."""
    return mol.to_graph(
        bond_cutoff_buffer=bond_cutoff_buffer, adjust_H=adjust_H
    )


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


class RMSDGrouper(MoleculeGrouper):
    """Group molecules based on RMSD (Root Mean Square Deviation) of atomic positions.
    Effective for precise 3D comparisons, ideal in contexts like crystallography or drug
    binding where exact spatial alignment is crucial.
    """

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
        # Cache sorted chemical symbols as sets for faster comparison
        self._chemical_symbol_sets = [
            set(mol.chemical_symbols) for mol in molecules
        ]

    def group(self) -> Tuple[List[List[Molecule]], List[List[int]]]:
        """Group molecules by geometric similarity."""
        n = len(self.molecules)
        indices = [(i, j) for i in range(n) for j in range(i + 1, n)]

        # Use map instead of imap_unordered for better parallelism
        with multiprocessing.Pool(self.num_procs) as pool:
            rmsd_values = pool.map(self._calculate_rmsd, indices)

        # Build adjacency matrix
        adj_matrix = np.zeros((n, n), dtype=bool)
        for (i, j), rmsd in zip(indices, rmsd_values):
            if rmsd < self.rmsd_threshold:
                adj_matrix[i, j] = adj_matrix[j, i] = True

        # Find connected components
        _, labels = connected_components(csr_matrix(adj_matrix))

        # Use np.unique(labels) approach for better memory efficiency
        unique_labels = np.unique(labels)
        groups = [
            [self.molecules[i] for i in np.where(labels == label)[0]]
            for label in unique_labels
        ]
        index_groups = [
            list(np.where(labels == label)[0]) for label in unique_labels
        ]

        return groups, index_groups

    def _calculate_rmsd(self, idx_pair: Tuple[int, int]) -> float:
        """Calculate RMSD between two molecules."""
        i, j = idx_pair
        mol1, mol2 = self.molecules[i], self.molecules[j]

        if (
            mol1.num_atoms != mol2.num_atoms
            or self._chemical_symbol_sets[i] != self._chemical_symbol_sets[j]
        ):
            return np.inf

        pos1, pos2 = mol1.positions, mol2.positions

        if self.align_molecules:
            logger.debug("Aligning molecules using Kabsch algorithm.")
            pos1, pos2, _, _, _ = kabsch_align(pos1, pos2)

        return np.sqrt(np.mean(np.sum((pos1 - pos2) ** 2, axis=1)))


class RMSDGrouperSharedMemory(MoleculeGrouper):
    """Group molecules based on RMSD using shared memory with minimal locking."""

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
        """Group molecules using shared memory with optimized parallelism."""
        n = len(self.molecules)
        num_atoms = self.molecules[0].positions.shape[0]

        # 🧠 **1️⃣ Create Shared Memory (RawArray - Faster, Less Locking)**
        shared_pos = RawArray("d", n * num_atoms * 3)  # 'd' -> float64

        # Convert RawArray into numpy view
        pos_np = np.frombuffer(shared_pos, dtype=np.float64).reshape(
            n, num_atoms, 3
        )

        # Copy molecular positions into shared memory (only once!)
        for i, mol in enumerate(self.molecules):
            pos_np[i] = mol.positions

        # 🏃‍♂️ **2️⃣ Run Parallel RMSD Calculation Using Explicit Shared Memory**
        indices = [(i, j) for i in range(n) for j in range(i + 1, n)]
        with multiprocessing.Pool(
            self.num_procs,
            initializer=self._init_worker,
            initargs=(shared_pos, (n, num_atoms, 3)),
        ) as pool:
            rmsd_values = pool.map(self._calculate_rmsd, indices)

        # 🏗️ **3️⃣ Construct Adjacency Matrix for Clustering**
        adj_matrix = np.zeros((n, n), dtype=bool)
        for (i, j), rmsd in zip(indices, rmsd_values):
            if rmsd < self.rmsd_threshold:
                adj_matrix[i, j] = adj_matrix[j, i] = True

        # 🔗 **4️⃣ Find Connected Components (Groups)**
        _, labels = connected_components(csr_matrix(adj_matrix))
        groups, index_groups = [], []
        for label in np.unique(labels):
            mask = labels == label
            groups.append([self.molecules[i] for i in np.where(mask)[0]])
            index_groups.append(list(np.where(mask)[0]))

        return groups, index_groups

    @staticmethod
    def _init_worker(shared_pos, pos_shape):
        """Worker process initializer to attach shared memory."""
        global shared_positions
        shared_positions = np.frombuffer(shared_pos, dtype=np.float64).reshape(
            pos_shape
        )

    def _calculate_rmsd(self, idx_pair: Tuple[int, int]) -> float:
        """Calculate RMSD efficiently using local copies of shared memory."""
        i, j = idx_pair

        # ✅ **Read from Shared Memory ONCE (No repeated locking)**
        pos1 = np.array(shared_positions[i])  # Copying reduces lock contention
        pos2 = np.array(shared_positions[j])

        if pos1.shape != pos2.shape:
            return np.inf

        if self.align_molecules:
            pos1, pos2, _, _, _ = kabsch_align(pos1, pos2)

        return np.sqrt(np.mean(np.sum((pos1 - pos2) ** 2, axis=1)))


class TanimotoSimilarityGrouper(MoleculeGrouper):
    """Group molecules using molecular fingerprinting and Tanimoto similarity.
    Tanimoto similarity is a measure of how similar two molecular fingerprints are,
    ranging from 0 (completely different) to 1 (identical).
    Default = 0.9 ensures molecules have a strong structural resemblance while
    allowing minor variations.

    Threshold	Effect	Use Case
    0.95 - 1.0	Very strict: Only almost identical molecules are grouped.
                            Ideal for highly similar molecules (e.g., stereoisomers).
    0.80 - 0.95	Moderately strict: Groups structurally similar molecules.
                            Useful for clustering molecules with minor functional group differences.
    0.50 - 0.80	More relaxed: Groups molecules with broad structural similarities.
                            Good for structural analogs or scaffold-based grouping.
    < 0.50	Very lenient: Even molecules with weak similarity are grouped.
                            Not recommended unless looking for very broad chemical families.
    """

    def __init__(
        self,
        molecules: Iterable[Molecule],
        similarity_threshold: float = 0.9,
        num_procs: int = 1,
    ):
        super().__init__(molecules, num_procs)
        self.similarity_threshold = similarity_threshold

    def _compute_fingerprint(
        self, mol: Molecule
    ) -> Optional[Tuple[int, DataStructs.ExplicitBitVect]]:
        """Compute RDKit fingerprint for a molecule, returning (index, fingerprint)."""
        rdkit_mol = mol.to_rdkit()
        if rdkit_mol:
            return Chem.RDKFingerprint(rdkit_mol)
        return None

    def group(self) -> Tuple[List[List[Molecule]], List[List[int]]]:
        """Cluster molecules using similarity-based connected components clustering."""
        # Compute fingerprints in parallel using ThreadPool (RDKit objects are not pickleable)
        with ThreadPool(self.num_procs) as pool:
            fingerprints = pool.map(self._compute_fingerprint, self.molecules)

        # Filter out molecules that failed to generate fingerprints
        valid_indices = [
            i for i, fp in enumerate(fingerprints) if fp is not None
        ]
        valid_fps = [fingerprints[i] for i in valid_indices]
        num_valid = len(valid_indices)

        if num_valid == 0:
            return [], []  # No valid molecules

        # Compute pairwise similarity matrix in parallel
        similarity_matrix = np.zeros((num_valid, num_valid), dtype=np.float32)
        pairs = [
            (i, j) for i in range(num_valid) for j in range(i + 1, num_valid)
        ]

        with ThreadPool(self.num_procs) as pool:
            similarities = pool.starmap(
                DataStructs.TanimotoSimilarity,
                [(valid_fps[i], valid_fps[j]) for i, j in pairs],
            )

        # Fill the similarity matrix
        for (i, j), sim in zip(pairs, similarities):
            similarity_matrix[i, j] = similarity_matrix[j, i] = sim

        # Apply thresholding to create an adjacency matrix
        adj_matrix = csr_matrix(similarity_matrix >= self.similarity_threshold)

        # Use connected components to find molecule groups
        _, labels = connected_components(adj_matrix)

        # Group molecules by their component labels
        unique_labels = np.unique(labels)
        mol_groups = [
            [
                self.molecules[valid_indices[i]]
                for i in np.where(labels == label)[0]
            ]
            for label in unique_labels
        ]
        idx_groups = [
            list(np.array(valid_indices)[np.where(labels == label)[0]])
            for label in unique_labels
        ]

        return mol_groups, idx_groups


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

        # Convert to RDKit molecules and filter out invalid ones
        self.molecules = [mol for mol in molecules if mol.to_rdkit()]
        self.rdkit_molecules = [mol.to_rdkit() for mol in self.molecules]

    @staticmethod
    def _get_fingerprint(
        rdkit_mol: Chem.Mol,
    ) -> Optional[DataStructs.ExplicitBitVect]:
        """Generate RDKit fingerprint for a molecule."""
        try:
            return rdFingerprintGenerator.GetRDKitFPGenerator().GetFingerprint(
                rdkit_mol
            )
        except Exception as e:
            logger.warning(f"Fingerprint generation failed: {str(e)}")
            return None

    @staticmethod
    def _compute_similarity(args):
        """Compute similarity between two fingerprints using Tanimoto similarity."""
        idx1, idx2, fps = args
        return (
            idx1,
            idx2,
            DataStructs.FingerprintSimilarity(fps[idx1], fps[idx2]),
        )

    def group(self) -> Tuple[List[List[Molecule]], List[List[int]]]:
        """Group molecules based on fingerprint similarity."""
        with multiprocessing.Pool(self.num_procs) as pool:
            fps = pool.map(self._get_fingerprint, self.rdkit_molecules)

        # Remove None fingerprints and adjust molecule lists accordingly
        valid_indices = [i for i, fp in enumerate(fps) if fp is not None]
        self.molecules = [self.molecules[i] for i in valid_indices]
        fps = [fps[i] for i in valid_indices]

        n = len(fps)
        if n == 0:
            return [], []

        # Compute similarity matrix
        range_indices = [
            (i, j, fps) for i in range(n) for j in range(i + 1, n)
        ]
        with multiprocessing.Pool(self.num_procs) as pool:
            results = pool.map(self._compute_similarity, range_indices)

        # Construct similarity matrix
        similarity_matrix = np.zeros((n, n))
        for i, j, sim in results:
            similarity_matrix[i, j] = similarity_matrix[j, i] = sim

        # Identify connected components
        graph = csr_matrix(similarity_matrix >= self.threshold)
        n_components, labels = connected_components(graph)

        # Group molecules based on component labels
        groups = [[] for _ in range(n_components)]
        indices = [[] for _ in range(n_components)]
        for idx, label in enumerate(labels):
            groups[label].append(self.molecules[idx])
            indices[label].append(valid_indices[idx])

        return groups, indices


class RDKitIsomorphismGrouper(MoleculeGrouper):
    """Group molecules using RDKit-based hashing and isomorphism checks.
    Seems very expensive for large molecules and unoptimized, yet."""

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
    """Group molecules based on molecular connectivity (graph isomorphism).
    Efficient for recognizing similar bond arrangements in large datasets.
    """

    def __init__(
        self,
        molecules: Iterable[Molecule],
        num_procs: int = 1,
        bond_cutoff_buffer: float = 0.0,
        adjust_H: bool = True,
    ):
        super().__init__(molecules, num_procs)
        self.bond_cutoff_buffer = bond_cutoff_buffer
        self.adjust_H = adjust_H

    def _are_isomorphic(self, g1: nx.Graph, g2: nx.Graph) -> bool:
        """Check if two molecular graphs are isomorphic."""
        return nx.is_isomorphic(
            g1,
            g2,
            node_match=lambda a, b: a["element"] == b["element"],
            edge_match=lambda a, b: a["bond_order"] == b["bond_order"],
        )

    def _check_isomorphism(
        self, idx_pair: Tuple[int, int]
    ) -> Tuple[int, int, bool]:
        """Multiprocessing-compatible function to check graph isomorphism."""
        i, j = idx_pair
        return i, j, self._are_isomorphic(self.graphs[i], self.graphs[j])

    def group(self) -> Tuple[List[List[Molecule]], List[List[int]]]:
        """Group molecules by connectivity using parallel isomorphism checks."""
        n = len(self.molecules)

        # Convert molecules to graphs in parallel
        with multiprocessing.Pool(self.num_procs) as pool:
            self.graphs = pool.starmap(
                to_graph_wrapper,
                [
                    (mol, self.bond_cutoff_buffer, self.adjust_H)
                    for mol in self.molecules
                ],
            )

        # Compute pairwise isomorphism in parallel
        indices = [(i, j) for i in range(n) for j in range(i + 1, n)]
        with multiprocessing.Pool(self.num_procs) as pool:
            isomorphic_pairs = pool.map(self._check_isomorphism, indices)

        # Build adjacency matrix
        adj_matrix = np.zeros((n, n), dtype=bool)
        for i, j, is_iso in isomorphic_pairs:
            if is_iso:
                adj_matrix[i, j] = adj_matrix[j, i] = True

        # Find connected components (groups of isomorphic molecules)
        _, labels = connected_components(csr_matrix(adj_matrix))

        unique_labels = np.unique(labels)
        groups = [
            [self.molecules[i] for i in np.where(labels == label)[0]]
            for label in unique_labels
        ]
        index_groups = [
            list(np.where(labels == label)[0]) for label in unique_labels
        ]

        return groups, index_groups


class ConnectivityGrouperSharedMemory(MoleculeGrouper):
    """Group molecules based on molecular connectivity."""

    def __init__(
        self,
        molecules: Iterable[Molecule],
        num_procs: int = 1,
        bond_cutoff_buffer: float = 0.0,
        adjust_H: bool = True,
    ):
        super().__init__(molecules, num_procs)
        self.bond_cutoff_buffer = bond_cutoff_buffer
        self.adjust_H = adjust_H

    def _are_isomorphic(self, g1: nx.Graph, g2: nx.Graph) -> bool:
        """Check if two molecular graphs are isomorphic."""
        return nx.is_isomorphic(
            g1,
            g2,
            node_match=lambda a, b: a["element"] == b["element"],
            edge_match=lambda a, b: a["bond_order"] == b["bond_order"],
        )

    def _check_isomorphism(
        self, pivot_graph_bytes, mol_graph_bytes, idx: int
    ) -> Tuple[int, bool]:
        """Helper function for multiprocessing: Checks isomorphism for a molecule using shared memory."""
        pivot_graph = pickle.loads(pivot_graph_bytes)
        mol_graph = pickle.loads(mol_graph_bytes)
        return idx, self._are_isomorphic(pivot_graph, mol_graph)

    def _convert_to_graphs(self):
        """Convert molecules to graphs in parallel and store in shared memory."""
        with multiprocessing.Pool(self.num_procs) as pool:
            graphs = pool.starmap(
                to_graph_wrapper,
                [
                    (mol, self.bond_cutoff_buffer, self.adjust_H)
                    for mol in self.molecules
                ],
            )

        # Serialize graphs using pickle and store in shared memory
        graph_bytes = [pickle.dumps(graph) for graph in graphs]

        # Store in shared NumPy array
        shared_array = np.array(graph_bytes, dtype=object)
        shm = shared_memory.SharedMemory(create=True, size=shared_array.nbytes)
        np.ndarray(
            shared_array.shape, dtype=shared_array.dtype, buffer=shm.buf
        )[:] = shared_array[:]

        return shm, shared_array.shape, shared_array.dtype

    def group(self) -> Tuple[List[List[Molecule]], List[List[int]]]:
        """Group molecules based on molecular connectivity using multiprocessing and shared memory."""
        groups = []
        shm, shape, dtype = self._convert_to_graphs()
        graphs = np.ndarray(shape, dtype=dtype, buffer=shm.buf)

        remaining = list(enumerate(zip(self.molecules, graphs)))

        while remaining:
            pivot_idx, (pivot_mol, pivot_graph_bytes) = remaining.pop(0)

            # Parallel isomorphism check
            results = Parallel(n_jobs=self.num_procs, backend="loky")(
                delayed(self._check_isomorphism)(
                    pivot_graph_bytes, g_bytes, idx
                )
                for idx, (_, g_bytes) in remaining
            )

            # Collect isomorphic molecules
            to_remove = {idx for idx, is_iso in results if is_iso}
            current_group = [pivot_mol] + [
                mol for idx, (mol, _) in remaining if idx in to_remove
            ]
            current_indices = [pivot_idx] + [
                idx for idx, _ in remaining if idx in to_remove
            ]

            remaining = [
                (idx, (mol, g_bytes))
                for idx, (mol, g_bytes) in remaining
                if idx not in to_remove
            ]
            groups.append((current_group, current_indices))

        shm.close()
        shm.unlink()  # Free shared memory

        mol_groups = [g[0] for g in groups]
        idx_groups = [g[1] for g in groups]
        return mol_groups, idx_groups


class StructureGrouperFactory:
    @staticmethod
    def create(structures, strategy="rdkit", num_procs=1):
        groupers = {
            "rmsd": RMSDGrouper,
            "tanimoto": TanimotoSimilarityGrouper,
            "fingerprint": RDKitFingerprintGrouper,
            "isomorphism": RDKitIsomorphismGrouper,
            "formula": FormulaGrouper,
            "connectivity": ConnectivityGrouper,
        }
        if strategy in groupers:
            logger.info(f"Using {strategy} grouping strategy.")
            return groupers[strategy](structures, num_procs)
        raise ValueError(f"Unknown grouping strategy: {strategy}")
