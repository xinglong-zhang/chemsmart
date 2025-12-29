"""
Molecular structure grouping using multiple strategies.

Algorithms for grouping molecules by geometric (RMSD), topological
(connectivity/graph isomorphism), and chemical (fingerprint similarity)
criteria. Many implementations support parallel execution and shared-memory
optimizations for large datasets.

Available strategies:
- RMSD: Root Mean Square Deviation of atomic positions
- Tanimoto: Fingerprint-based chemical similarity
- Connectivity: Graph isomorphism-based grouping (NetworkX)
- Formula: Chemical formula-based grouping
- Isomorphism: RDKit-based hashing/isomorphism

Key classes include:
- MoleculeGrouper: Abstract base class for all groupers
- RMSDGrouper: Abstract parent class for RMSD-based groupers
- BasicRMSDGrouper: Standard RMSD calculation using Euclidean distance
- HungarianRMSDGrouper: RMSD grouping with optimal atom assignment
- RMSDGrouperSharedMemory: RMSD grouping with shared-memory optimization
- TanimotoSimilarityGrouper: Chemical fingerprint similarity
- RDKitIsomorphismGrouper: RDKit hashing and isomorphism grouping
- ConnectivityGrouper: Molecular connectivity (graph isomorphism)
- ConnectivityGrouperSharedMemory: Connectivity grouping with shared memory
- FormulaGrouper: Chemical formula-based grouping
- StructureGrouperFactory: Factory for creating grouper instances
"""

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
from networkx.algorithms import isomorphism
from rdkit import Chem, DataStructs
from rdkit.Chem import rdMolHash
from rdkit.Chem.rdFingerprintGenerator import GetRDKitFPGenerator
from scipy.optimize import linear_sum_assignment  # for Hungarian algorithm
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial.distance import cdist  # for Hungarian algorithm cost matrix

from chemsmart.io.molecules.structure import Molecule
from chemsmart.utils.utils import kabsch_align

logger = logging.getLogger(__name__)


def to_graph_wrapper(
    mol: Molecule, bond_cutoff_buffer: float = 0.0, adjust_H: bool = True
) -> nx.Graph:
    """
    Global helper function to call Molecule.to_graph() for multiprocessing.

    Provides a picklable wrapper for the Molecule.to_graph() method that
    can be used with multiprocessing pools.

    Args:
        mol (Molecule): Molecule instance to convert to graph.
        bond_cutoff_buffer (float): Buffer for bond cutoff distance.
        adjust_H (bool): Whether to adjust hydrogen bond detection.

    Returns:
        networkx.Graph: Molecular graph representation.
    """
    return mol.to_graph(
        bond_cutoff_buffer=bond_cutoff_buffer, adjust_H=adjust_H
    )


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
    """

    def __init__(self, molecules: Iterable[Molecule], num_procs: int = 1):
        """
        Initialize the molecular grouper.

        Args:
            molecules (Iterable[Molecule]): Collection of molecules to group.
            num_procs (int): Number of processes for parallel computation.
                Defaults to 1.
        """
        self.molecules = molecules
        self.num_procs = int(max(1, num_procs))

        self._validate_inputs()

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

    def unique(self) -> List[Molecule]:
        """
        Get unique representative molecules from each group.

        Returns the first molecule from each group as a representative
        of that structural family.

        Returns:
            List[Molecule]: List of unique representative molecules.
        """
        groups, _ = self.group()
        return [group[0] for group in groups]


class RMSDGrouper(MoleculeGrouper):
    """
    Abstract base class for RMSD-based molecular grouping.

    Groups molecules based on geometric similarity of atomic positions using
    various RMSD calculation methods. This base class provides common
    functionality for all RMSD-based groupers while allowing subclasses
    to implement specific RMSD calculation algorithms.

    Follows the same design pattern as JobRunner - provides the core
    grouping logic while allowing subclasses to implement specific
    RMSD calculation methods via _calculate_rmsd_core().

    Attributes:
        molecules (Iterable[Molecule]): Inherited; collection of molecules to
            group.
        num_procs (int): Inherited; number of worker processes/threads.
        threshold (float): RMSD threshold for grouping molecules.
        align_molecules (bool): Whether to align molecules before RMSD calculation.
        ignore_hydrogens (bool): Whether to exclude hydrogen atoms from RMSD.
    """

    def __init__(
        self,
        molecules: Iterable[Molecule],
        threshold=None,  # RMSD threshold for grouping
        num_procs: int = 1,
        align_molecules: bool = True,
        ignore_hydrogens: bool = False,  # Option to ignore H atoms for grouping
    ):
        """
        Initialize RMSD-based molecular grouper.

        Args:
            molecules (Iterable[Molecule]): Collection of molecules to group.
            threshold (float): RMSD threshold for grouping. Defaults to 0.5.
            num_procs (int): Number of processes for parallel computation.
            align_molecules (bool): Whether to align molecules using Kabsch
                algorithm before RMSD calculation. Defaults to True.
            ignore_hydrogens (bool): Whether to exclude hydrogen atoms from
                RMSD calculation. Defaults to False.
        """
        super().__init__(molecules, num_procs)
        if threshold is None:
            threshold = 0.5
        self.threshold = threshold  # RMSD threshold for grouping
        self.align_molecules = align_molecules
        self.ignore_hydrogens = ignore_hydrogens
        # Cache sorted chemical symbols as sets for faster comparison
        self._chemical_symbol_sets = [
            set(mol.chemical_symbols) for mol in molecules
        ]

    def _get_heavy_atoms(self, mol: Molecule) -> Tuple[np.ndarray, List[str]]:
        """
        Extract heavy atoms (non-hydrogen) if ignore_hydrogens is enabled.

        Args:
            mol (Molecule): Molecule to process.

        Returns:
            Tuple[np.ndarray, List[str]]: Tuple containing positions array
                and chemical symbols list (filtered or full based on settings).
        """
        if self.ignore_hydrogens:
            non_h_indices = [
                i for i, sym in enumerate(mol.chemical_symbols) if sym != "H"
            ]
            return mol.positions[non_h_indices], [
                mol.chemical_symbols[i] for i in non_h_indices
            ]
        return (
            mol.positions,
            mol.chemical_symbols,
        )  # Use all atoms if flag is False

    @abstractmethod
    def _calculate_rmsd(self, idx_pair: Tuple[int, int]) -> float:
        """Calculate RMSD between two molecules. Must be implemented by subclasses."""
        pass

    def group(self) -> Tuple[List[List[Molecule]], List[List[int]]]:
        """
        Group molecules by geometric similarity using RMSD clustering.

        Computes pairwise RMSD values between all molecules and groups
        those within the specified threshold using connected components
        clustering.

        Returns:
            Tuple[List[List[Molecule]], List[List[int]]]: Tuple containing:
                - List of molecule groups (each group is a list of molecules)
                - List of index groups (corresponding indices for each group)
        """
        n = len(self.molecules)
        indices = [(i, j) for i in range(n) for j in range(i + 1, n)]

        # Use map instead of imap_unordered for better parallelism
        with multiprocessing.Pool(self.num_procs) as pool:
            rmsd_values = pool.map(self._calculate_rmsd, indices)

        # Build adjacency matrix
        adj_matrix = np.zeros((n, n), dtype=bool)
        for (i, j), rmsd in zip(indices, rmsd_values):
            if rmsd < self.threshold:
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

    def __repr__(self):
        return (
            f"{self.__class__.__name__}(threshold={self.threshold}, "
            f"num_procs={self.num_procs}, align_molecules={self.align_molecules}, "
            f"ignore_hydrogens={self.ignore_hydrogens})"
        )


class BasicRMSDGrouper(RMSDGrouper):
    """
    Basic RMSD grouper using standard Euclidean distance calculation.

    Implements the most straightforward RMSD calculation using the standard
    formula: sqrt(mean(sum((pos1 - pos2)^2))). This is the classic RMSD
    implementation that compares atomic positions directly.
    """

    def _calculate_rmsd(self, idx_pair: Tuple[int, int]) -> float:
        """Calculate RMSD between two molecules without atom reordering."""
        i, j = idx_pair
        mol1, mol2 = self.molecules[i], self.molecules[j]

        if self.ignore_hydrogens:
            pos1, symbols1 = self._get_heavy_atoms(mol1)
            pos2, symbols2 = self._get_heavy_atoms(mol2)
        else:
            pos1, symbols1 = mol1.positions, list(mol1.chemical_symbols)
            pos2, symbols2 = mol2.positions, list(mol2.chemical_symbols)

        # Quick check for compatibility
        if len(symbols1) != len(symbols2) or sorted(symbols1) != sorted(
            symbols2
        ):
            return np.inf

        if self.align_molecules:
            logger.debug("Aligning molecules using Kabsch algorithm.")
            pos1, pos2, _, _, _ = kabsch_align(pos1, pos2)
        rmsd = np.sqrt(np.mean(np.sum((pos1 - pos2) ** 2, axis=1)))
        return rmsd


class HungarianRMSDGrouper(RMSDGrouper):
    """
    Hungarian RMSD grouper for optimal atom assignment.

    Uses the Hungarian algorithm (Kuhn-Munkres algorithm) to find the optimal
    assignment of atoms between two molecules that minimizes RMSD. This approach
    handles cases where atom ordering might differ between molecules of the same
    chemical structure.

    The Hungarian algorithm ensures that atoms of the same element type are
    optimally paired to minimize the sum of squared distances, resulting in
    more accurate RMSD calculations for molecules with permuted atom arrangements.

    Supports additional parameters for spyrmsd compatibility:
    - center: Center coordinates at origin before calculation
    - minimize: Use optimal superposition (Kabsch alignment)
    """

    def __init__(
        self,
        molecules,
        threshold=None,
        num_procs: int = 1,
        align_molecules: bool = True,
        ignore_hydrogens: bool = False,
        center: bool = False,
    ):
        """
        Initialize Hungarian RMSD grouper.

        Args:
            molecules: Collection of molecules to group.
            threshold (float): RMSD threshold for grouping.
            num_procs (int): Number of processes for parallel computation.
            align_molecules (bool): Whether to align molecules (legacy parameter).
            ignore_hydrogens (bool): Whether to exclude hydrogen atoms.
            center (bool): Center coordinates at origin before calculation.
            minimize (bool): Use optimal superposition (Kabsch alignment).
        """
        super().__init__(
            molecules, threshold, num_procs, align_molecules, ignore_hydrogens
        )
        self.center = center

    def _calculate_rmsd(self, idx_pair):
        i, j = idx_pair
        mol1, mol2 = self.molecules[i], self.molecules[j]

        if self.ignore_hydrogens:
            pos1, symbols1 = self._get_heavy_atoms(mol1)
            pos2, symbols2 = self._get_heavy_atoms(mol2)
        else:
            pos1, symbols1 = mol1.positions, list(mol1.chemical_symbols)
            pos2, symbols2 = mol2.positions, list(mol2.chemical_symbols)

        if len(symbols1) != len(symbols2) or sorted(symbols1) != sorted(
            symbols2
        ):
            return np.inf

        # Center coordinates if requested
        if self.center:
            pos1 = pos1 - np.mean(pos1, axis=0)
            pos2 = pos2 - np.mean(pos2, axis=0)

        # Use Hungarian algorithm for optimal atom matching
        elements = sorted(set(symbols1))
        matched_idx1 = []
        matched_idx2 = []
        for elem in elements:
            idxs1 = [k for k, s in enumerate(symbols1) if s == elem]
            idxs2 = [k for k, s in enumerate(symbols2) if s == elem]

            if len(idxs1) == 1 and len(idxs2) == 1:
                # For single atoms, direct match
                matched_idx1.extend(idxs1)
                matched_idx2.extend(idxs2)
            else:
                # For multiple atoms of same type, use Hungarian algorithm
                pos1_elem = pos1[idxs1]
                pos2_elem = pos2[idxs2]
                # Use scipy's cdist for efficient distance matrix calculation
                dist_matrix = cdist(pos1_elem, pos2_elem, metric="sqeuclidean")
                row_ind, col_ind = linear_sum_assignment(dist_matrix)
                matched_idx1.extend([idxs1[r] for r in row_ind])
                matched_idx2.extend([idxs2[c] for c in col_ind])

        pos1_matched = pos1[matched_idx1]
        pos2_matched = pos2[matched_idx2]

        # Apply alignment based on minimize parameter or legacy align_molecules
        if self.align_molecules:
            logger.debug("Aligning molecules using Kabsch algorithm.")
            pos1_matched, pos2_matched, _, _, _ = kabsch_align(
                pos1_matched, pos2_matched
            )

        rmsd = np.sqrt(
            np.mean(np.sum((pos1_matched - pos2_matched) ** 2, axis=1))
        )
        return rmsd


class SpyRMSDGrouper(RMSDGrouper):
    """
    SpyRMSD grouper with graph isomorphism symmetry correction.

    Uses graph isomorphism for advanced RMSD calculations with symmetry correction.
    This approach handles molecular symmetries and atom permutations more
    comprehensively than basic Hungarian assignment.

    Supports optimal superposition and can handle cases where molecular graphs
    need to be compared for finding the best atom mapping.
    """

    def __init__(
        self,
        molecules,
        threshold=None,
        num_procs: int = 1,
        align_molecules: bool = True,
        ignore_hydrogens: bool = False,
        center: bool = False,  # Center molecules
        cache: bool = True,  # Cache graph isomorphisms
    ):
        """
        Initialize SpyRMSD-based molecular grouper.

        Args:
            molecules: Collection of molecules to group.
            threshold (float): RMSD threshold for grouping. Defaults to 0.5.
            num_procs (int): Number of processes for parallel computation.
            align_molecules (bool): Whether to align molecules (legacy parameter).
            ignore_hydrogens (bool): Whether to exclude hydrogen atoms.
            center (bool): Center molecules at origin before calculation.
            cache (bool): Cache graph isomorphisms for efficiency.
        """
        super().__init__(
            molecules, threshold, num_procs, align_molecules, ignore_hydrogens
        )
        self.center = center
        self.cache = cache

    @staticmethod
    def _center_coordinates(pos: np.ndarray) -> np.ndarray:
        """
        Calculate the center of geometry (centroid) of atomic coordinates.

        Args:
            pos (np.ndarray): Array of atomic coordinates (N x 3).

        Returns:
            np.ndarray: Center of geometry as a 3D vector.
        """
        return np.mean(pos, axis=0)

    @staticmethod
    def _symbol_to_atomicnum(symbol: str) -> int:
        """
        Convert element symbol to atomic number.

        Args:
            symbol (str): Element symbol (e.g., 'H', 'C', 'O').

        Returns:
            int: Atomic number, or 0 if unknown.
        """
        symbol_map = {
            "H": 1,
            "He": 2,
            "Li": 3,
            "Be": 4,
            "B": 5,
            "C": 6,
            "N": 7,
            "O": 8,
            "F": 9,
            "Ne": 10,
            "Na": 11,
            "Mg": 12,
            "Al": 13,
            "Si": 14,
            "P": 15,
            "S": 16,
            "Cl": 17,
            "Ar": 18,
            "K": 19,
            "Ca": 20,
            "Br": 35,
            "I": 53,
        }
        return symbol_map.get(symbol, 0)

    def _symmrmsd(
        self,
        pos1: np.ndarray,
        pos2: np.ndarray,
        symbols1: list,
        symbols2: list,
        adj_matrix1: np.ndarray,
        adj_matrix2: np.ndarray,
    ) -> float:
        """
        Calculate symmetry-corrected RMSD using graph isomorphism.

        Accounts for molecular symmetry by finding all possible atom mappings
        through graph isomorphism and selecting the one that gives minimum RMSD.

        Args:
            pos1 (np.ndarray): First set of coordinates (N x 3).
            pos2 (np.ndarray): Second set of coordinates (N x 3).
            symbols1 (list): Chemical symbols for first structure.
            symbols2 (list): Chemical symbols for second structure.
            adj_matrix1 (np.ndarray): Adjacency matrix for first structure.
            adj_matrix2 (np.ndarray): Adjacency matrix for second structure.

        Returns:
            float: Symmetry-corrected RMSD value.

        Notes:
            This is a simplified implementation. For complex molecules with
            high symmetry, this may be computationally expensive. Falls back
            to Hungarian RMSD if graphs are not isomorphic.
        """

        # Verify same number of atoms
        if len(pos1) != len(pos2):
            return np.inf

        # Verify same atomic composition
        if sorted(symbols1) != sorted(symbols2):
            return np.inf

        # Center coordinates if requested
        if self.center:
            pos1 = pos1 - self._center_coordinates(pos1)
            pos2 = pos2 - self._center_coordinates(pos2)

        # Create NetworkX graphs from adjacency matrices
        G1 = nx.from_numpy_array(adj_matrix1)
        G2 = nx.from_numpy_array(adj_matrix2)

        # Add atomic numbers as node attributes
        for i, symbol in enumerate(symbols1):
            G1.nodes[i]["element"] = self._symbol_to_atomicnum(symbol)
        for i, symbol in enumerate(symbols2):
            G2.nodes[i]["element"] = self._symbol_to_atomicnum(symbol)

        # Find all graph isomorphisms (atom mappings)
        node_match = isomorphism.categorical_node_match("element", 0)
        GM = isomorphism.GraphMatcher(G1, G2, node_match=node_match)

        # Try all possible isomorphisms and find minimum RMSD
        min_rmsd = np.inf

        if GM.is_isomorphic():
            for mapping in GM.isomorphisms_iter():
                # Reorder pos2 according to this mapping
                reordered_pos2 = np.array(
                    [pos2[mapping[i]] for i in range(len(pos2))]
                )

                # Calculate RMSD for this mapping
                if self.align_molecules:
                    _, _, _, _, rmsd = kabsch_align(pos1, reordered_pos2)
                else:
                    diff = pos1 - reordered_pos2
                    rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))

                if rmsd < min_rmsd:
                    min_rmsd = rmsd

            return min_rmsd
        else:
            # If not isomorphic, return infinity (structures are fundamentally different)
            return np.inf

    def _calculate_rmsd(self, idx_pair):
        """Calculate RMSD using symmetry correction with graph isomorphism."""

        i, j = idx_pair
        mol1, mol2 = self.molecules[i], self.molecules[j]

        if self.ignore_hydrogens:
            pos1, symbols1 = self._get_heavy_atoms(mol1)
            pos2, symbols2 = self._get_heavy_atoms(mol2)
        else:
            pos1, symbols1 = mol1.positions, list(mol1.chemical_symbols)
            pos2, symbols2 = mol2.positions, list(mol2.chemical_symbols)

        # Quick compatibility check
        if len(symbols1) != len(symbols2) or sorted(symbols1) != sorted(
            symbols2
        ):
            return np.inf

        # Use symmetry-corrected RMSD if molecules have connectivity
        # Try to get molecular graphs for symmetry correction
        adj_matrix1 = (
            mol1.adjacency_matrix
            if hasattr(mol1, "adjacency_matrix")
            else None
        )
        adj_matrix2 = (
            mol2.adjacency_matrix
            if hasattr(mol2, "adjacency_matrix")
            else None
        )

        if adj_matrix1 is not None and adj_matrix2 is not None:
            rmsd = self._symmrmsd(
                pos1=pos1,
                pos2=pos2,
                symbols1=symbols1,
                symbols2=symbols2,
                adj_matrix1=adj_matrix1,
                adj_matrix2=adj_matrix2,
            )
        else:
            return np.inf
        return rmsd


class RMSDGrouperSharedMemory(MoleculeGrouper):
    """
    Group molecules based on RMSD using shared memory optimization.

    Optimized version of RMSDGrouper that uses shared memory to reduce
    data copying overhead in multiprocessing scenarios. Provides faster
    computation for large datasets by minimizing memory allocation and
    inter-process communication costs.

    Attributes:
        molecules (Iterable[Molecule]): Inherited; collection of molecules to
            group.
        num_procs (int): Inherited; number of worker processes.
        threshold (float): RMSD threshold for grouping molecules.
        align_molecules (bool): Whether to align molecules before RMSD
            calculation.
    """

    def __init__(
        self,
        molecules: Iterable[Molecule],
        threshold: float = 0.5,  # RMSD threshold for grouping
        num_procs: int = 1,
        align_molecules: bool = True,
    ):
        """
        Initialize RMSD grouper with shared memory optimization.

        Args:
            molecules (Iterable[Molecule]): Collection of molecules to group.
            threshold (float): RMSD threshold for grouping. Defaults to 0.5.
            num_procs (int): Number of processes for parallel computation.
            align_molecules (bool): Whether to align molecules using Kabsch
                algorithm before RMSD calculation. Defaults to True.
        """
        super().__init__(molecules, num_procs)
        self.threshold = threshold
        self.align_molecules = align_molecules

    def group(self) -> Tuple[List[List[Molecule]], List[List[int]]]:
        """
        Group molecules using shared memory with optimized parallelism.

        Uses RawArray shared memory to minimize data copying between processes.
        Molecular positions are stored once in shared memory and accessed
        by worker processes for RMSD calculations.

        Returns:
            Tuple[List[List[Molecule]], List[List[int]]]: Tuple containing:
                - List of molecule groups (each group is a list of molecules)
                - List of index groups (corresponding indices for each group)
        """
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
            if rmsd < self.threshold:
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
        """
        Initialize worker process with shared memory access.

        Sets up global shared memory access for worker processes,
        allowing them to read molecular positions without data copying.

        Args:
            shared_pos: RawArray containing shared position data.
            pos_shape (tuple): Shape tuple for reshaping the shared array.
        """
        global shared_positions
        shared_positions = np.frombuffer(shared_pos, dtype=np.float64).reshape(
            pos_shape
        )

    def _calculate_rmsd(self, idx_pair: Tuple[int, int]) -> float:
        """
        Calculate RMSD efficiently using shared memory.

        Computes RMSD between two molecules by reading their positions
        from shared memory and creating local copies to reduce lock
        contention during computation.

        Args:
            idx_pair (Tuple[int, int]): Pair of molecule indices to compare.

        Returns:
            float: RMSD value between the two molecules, or np.inf if
                   shapes don't match.
        """
        i, j = idx_pair

        # Read from Shared Memory ONCE (No repeated locking)
        pos1 = np.array(shared_positions[i])  # Copying reduces lock contention
        pos2 = np.array(shared_positions[j])

        if pos1.shape != pos2.shape:
            return np.inf

        if self.align_molecules:
            pos1, pos2, _, _, _ = kabsch_align(pos1, pos2)

        return np.sqrt(np.mean(np.sum((pos1 - pos2) ** 2, axis=1)))


class TanimotoSimilarityGrouper(MoleculeGrouper):
    """
    Groups molecules based on fingerprint similarity using Tanimoto coefficient.

    This class supports different fingerprint types and uses connected components
    clustering to group structurally similar molecules.

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
        threshold=None,  # Tanimoto similarity threshold
        num_procs: int = 1,
        use_rdkit_fp: bool = True,  # Allows switching between RDKit FP and RDKFingerprint
    ):
        """
        Initialize Tanimoto similarity-based molecular grouper.

        Args:
            molecules (Iterable[Molecule]): Collection of molecules to group.
            threshold (float): Tanimoto similarity threshold. Defaults to 0.9.
            num_procs (int): Number of processes for parallel computation.
            use_rdkit_fp (bool): Whether to use RDKit fingerprints (True) or
                RDKFingerprint method (False). Defaults to True.
        """
        super().__init__(molecules, num_procs)
        if threshold is None:
            threshold = 0.9
        self.threshold = threshold
        self.use_rdkit_fp = use_rdkit_fp  # Choose fingerprinting method

        # Convert valid molecules to RDKit format
        self.rdkit_molecules = [
            mol.to_rdkit() for mol in molecules if mol.to_rdkit()
        ]
        self.valid_molecules = [mol for mol in molecules if mol.to_rdkit()]

    def _get_fingerprint(
        self, rdkit_mol: Chem.Mol
    ) -> Optional[DataStructs.ExplicitBitVect]:
        """
        Generate an RDKit fingerprint for a molecule.

        Creates molecular fingerprints using either RDKit fingerprint
        generator or RDKFingerprint method based on configuration.

        Args:
            rdkit_mol (Chem.Mol): RDKit molecule object.

        Returns:
            Optional[DataStructs.ExplicitBitVect]: Molecular fingerprint
                or None if generation fails.
        """
        try:
            if self.use_rdkit_fp:
                return GetRDKitFPGenerator().GetFingerprint(rdkit_mol)
            else:
                return Chem.RDKFingerprint(rdkit_mol)  # Alternative method
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

        # Compute similarity matrix
        similarity_matrix = np.zeros((num_valid, num_valid), dtype=np.float32)
        pairs = [
            (i, j) for i in range(num_valid) for j in range(i + 1, num_valid)
        ]

        with ThreadPool(self.num_procs) as pool:
            similarities = pool.starmap(
                DataStructs.FingerprintSimilarity,
                [(valid_fps[i], valid_fps[j]) for i, j in pairs],
            )

        # Fill similarity matrix
        for (i, j), sim in zip(pairs, similarities):
            similarity_matrix[i, j] = similarity_matrix[j, i] = sim

        # Apply threshold and create adjacency matrix
        adj_matrix = csr_matrix(similarity_matrix >= self.threshold)

        # Use connected components clustering
        _, labels = connected_components(adj_matrix)

        # Build molecule groups
        unique_labels = np.unique(labels)
        mol_groups = [
            [
                self.valid_molecules[valid_indices[i]]
                for i in np.where(labels == label)[0]
            ]
            for label in unique_labels
        ]
        idx_groups = [
            list(np.array(valid_indices)[np.where(labels == label)[0]])
            for label in unique_labels
        ]

        return mol_groups, idx_groups


class RDKitIsomorphismGrouper(MoleculeGrouper):
    """
    Group molecules using RDKit hashing with optional isomorphism checks.

    First clusters molecules by RDKit molecular hash (choice depends on
    options), then optionally verifies equivalence using InChIKey equality
    as a lightweight isomorphism proxy. This can be computationally
    expensive for large sets.

    Hashing choices (see `_get_mol_hash`):
    - If `use_tautomers` is True: `rdMolHash.HashFunction.Tautomer`.
    - Else if `use_stereochemistry` is True: `rdMolHash.HashFunction.AnonymousGraph`.
    - Else: `rdMolHash.HashFunction.MolFormula`.

    Attributes:
        molecules (Iterable[Molecule]): Inherited; collection of molecules to
            group.
        num_procs (int): Inherited; number of worker processes.
        use_stereochemistry (bool): Whether to consider stereochemistry.
        use_tautomers (bool): Whether to consider tautomeric forms.
    """

    def __init__(
        self,
        molecules: Iterable[Molecule],
        num_procs: int = 1,
        use_stereochemistry: bool = True,
        use_tautomers: bool = False,
    ):
        """
        Initialize RDKit isomorphism-based molecular grouper.

        Args:
            molecules (Iterable[Molecule]): Collection of molecules to group.
            num_procs (int): Number of processes for parallel computation.
            use_stereochemistry (bool): Whether to consider stereochemical
                differences in grouping. Defaults to True.
            use_tautomers (bool): Whether to consider tautomeric forms
                as equivalent. Defaults to False.
        """
        super().__init__(molecules, num_procs)
        self.use_stereochemistry = use_stereochemistry
        self.use_tautomers = use_tautomers

    def _get_mol_hash(self, mol: Molecule) -> Optional[str]:
        """
        Generate canonical hash for molecular structure identification.

        Creates a canonical hash string for the molecule using RDKit's
        hashing functions. Hash type depends on stereochemistry and
        tautomer configuration settings.

        Args:
            mol (Molecule): Molecule to generate hash for.

        Returns:
            Optional[str]: Canonical hash string, or None if generation fails.
        """
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
        """
        Group molecules using structural isomorphism detection.

        Uses RDKit molecular hashing for initial grouping, followed by
        detailed isomorphism checks when stereochemistry or tautomer
        considerations are enabled.

        Returns:
            Tuple[List[List[Molecule]], List[List[int]]]: Tuple containing:
                - List of molecule groups (each group is a list of molecules)
                - List of index groups (corresponding indices for each group)
        """
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
        """
        Check equivalence via InChIKey (isomorphism proxy).

        Compares the RDKit InChIKey strings of both molecules. This is a
        lightweight proxy for isomorphism when stereochemistry or tautomer
        checks are enabled; it does not perform an explicit graph
        isomorphism test.

        Args:
            mol1 (Molecule): First molecule to compare.
            mol2 (Molecule): Second molecule to compare.

        Returns:
            bool: True if InChIKeys match; False otherwise (including on failure).
        """
        try:
            return Chem.MolToInchiKey(mol1.to_rdkit()) == Chem.MolToInchiKey(
                mol2.to_rdkit()
            )
        except Exception as e:
            logger.warning(f"Isomorphism check failed: {str(e)}")
            return False


class FormulaGrouper(MoleculeGrouper):
    """
    Group molecules by chemical formula.

    Groups molecules based solely on their chemical formula composition,
    making it suitable when elemental composition is the primary concern.
    Ideal for initial filtering and broad chemical classification.
    """

    def group(self):
        """
        Group molecules by chemical formula composition.

        Creates groups based on identical chemical formulas, regardless
        of structural or stereochemical differences. Each group contains
        molecules with the same elemental composition.

        Returns:
            Tuple[List[List[Molecule]], List[List[int]]]: Tuple containing:
                - List of molecule groups (each group is a list of molecules)
                - List of index groups (corresponding indices for each group)
        """
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
    """
    Group molecules based on molecular connectivity (graph isomorphism).

    Groups molecules by analyzing their bond connectivity patterns using
    graph isomorphism. Efficient for recognizing similar bond arrangements
    in large datasets regardless of 3D spatial configuration.

    Attributes:
        molecules (Iterable[Molecule]): Inherited; collection of molecules to
            group.
        num_procs (int): Inherited; number of worker processes.
        threshold (float): Buffer for bond cutoff distance.
        adjust_H (bool): Whether to adjust hydrogen bond detection.
    """

    def __init__(
        self,
        molecules: Iterable[Molecule],
        num_procs: int = 1,
        threshold=None,  # Buffer for bond cutoff
        adjust_H: bool = True,
    ):
        """
        Initialize connectivity-based molecular grouper.

        Args:
            molecules (Iterable[Molecule]): Collection of molecules to group.
            num_procs (int): Number of processes for parallel computation.
            threshold (float): Buffer for bond cutoff distance. Defaults to 0.0.
            adjust_H (bool): Whether to adjust hydrogen bond detection.
                Defaults to True.
        """
        super().__init__(molecules, num_procs)
        if threshold is None:
            threshold = 0.0
        self.threshold = threshold  # Buffer for bond cutoff
        self.adjust_H = adjust_H

    def _are_isomorphic(self, g1: nx.Graph, g2: nx.Graph) -> bool:
        """
        Check if two molecular graphs are isomorphic (NetworkX).

        Uses `networkx.is_isomorphic` with attribute-aware matching:
        - Nodes must have equal `element` values.
        - Edges must have equal `bond_order` values.

        Args:
            g1 (nx.Graph): First molecular graph.
            g2 (nx.Graph): Second molecular graph.

        Returns:
            bool: True if graphs are isomorphic, False otherwise.
        """
        return nx.is_isomorphic(
            g1,
            g2,
            node_match=lambda a, b: a["element"] == b["element"],
            edge_match=lambda a, b: a["bond_order"] == b["bond_order"],
        )

    def _check_isomorphism(
        self, idx_pair: Tuple[int, int]
    ) -> Tuple[int, int, bool]:
        """
        Check graph isomorphism between two molecules for multiprocessing.

        Multiprocessing-compatible function that checks whether two
        molecular graphs are isomorphic based on their connectivity
        patterns and atomic properties.

        Args:
            idx_pair (Tuple[int, int]): Pair of molecule indices to compare.

        Returns:
            Tuple[int, int, bool]: Original indices and isomorphism result.
        """
        i, j = idx_pair
        return i, j, self._are_isomorphic(self.graphs[i], self.graphs[j])

    def group(self) -> Tuple[List[List[Molecule]], List[List[int]]]:
        """
        Group molecules by connectivity using parallel isomorphism checks.

        Converts molecules to graph representations and performs pairwise
        isomorphism checks in parallel. Uses connected components clustering
        to identify groups of structurally equivalent molecules.

        Returns:
            Tuple[List[List[Molecule]], List[List[int]]]: Tuple containing:
                - List of molecule groups (each group is a list of molecules)
                - List of index groups (corresponding indices for each group)
        """
        n = len(self.molecules)

        # Convert molecules to graphs in parallel
        with multiprocessing.Pool(self.num_procs) as pool:
            self.graphs = pool.starmap(
                to_graph_wrapper,
                [
                    (mol, self.threshold, self.adjust_H)
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
    """
    Group molecules based on molecular connectivity using shared memory.

    Optimized version of ConnectivityGrouper that uses shared memory
    for storing molecular graph data to reduce memory overhead in
    multiprocessing scenarios. Particularly useful for large molecular
    datasets where memory efficiency is critical.

    Attributes:
        molecules (Iterable[Molecule]): Inherited; collection of molecules to
            group.
        num_procs (int): Inherited; number of worker processes.
        threshold (float): Buffer for bond cutoff distance.
        adjust_H (bool): Whether to adjust hydrogen bond detection.
    """

    def __init__(
        self,
        molecules: Iterable[Molecule],
        num_procs: int = 1,
        threshold: float = 0.0,  # Buffer for bond cutoff
        adjust_H: bool = True,
    ):
        """
        Initialize connectivity grouper with shared memory optimization.

        Args:
            molecules (Iterable[Molecule]): Collection of molecules to group.
            num_procs (int): Number of processes for parallel computation.
            threshold (float): Buffer for bond cutoff distance. Defaults to 0.0.
            adjust_H (bool): Whether to adjust hydrogen bond detection.
                Defaults to True.
        """
        super().__init__(molecules, num_procs)
        self.threshold = threshold
        self.adjust_H = adjust_H

    def _are_isomorphic(self, g1: nx.Graph, g2: nx.Graph) -> bool:
        """
        Check if two molecular graphs are isomorphic (NetworkX).

        Uses `networkx.is_isomorphic` with attribute-aware matching:
        - Nodes must have equal `element` values.
        - Edges must have equal `bond_order` values.

        Args:
            g1 (nx.Graph): First molecular graph.
            g2 (nx.Graph): Second molecular graph.

        Returns:
            bool: True if graphs are isomorphic, False otherwise.
        """
        return nx.is_isomorphic(
            g1,
            g2,
            node_match=lambda a, b: a["element"] == b["element"],
            edge_match=lambda a, b: a["bond_order"] == b["bond_order"],
        )

    def _check_isomorphism(
        self, pivot_graph_bytes, mol_graph_bytes, idx: int
    ) -> Tuple[int, bool]:
        """
        Check equivalence for multiprocessing using serialized graphs.

        Deserializes NetworkX graphs from bytes and tests connectivity
        isomorphism against a pivot graph.

        Args:
            pivot_graph_bytes (bytes): Pickled pivot molecular graph.
            mol_graph_bytes (bytes): Pickled molecular graph to compare.
            idx (int): Index of the molecule being compared.

        Returns:
            Tuple[int, bool]: (idx, is_isomorphic) result.
        """
        pivot_graph = pickle.loads(pivot_graph_bytes)
        mol_graph = pickle.loads(mol_graph_bytes)
        return idx, self._are_isomorphic(pivot_graph, mol_graph)

    def _convert_to_graphs(self):
        """
        Convert molecules to graphs and store in shared memory.

        Converts all molecules to NetworkX graph representations in parallel,
        then serializes them using pickle and stores in shared memory for
        efficient access by worker processes.

        Returns:
            Tuple[shared_memory.SharedMemory, tuple, np.dtype]: Shared memory
            handle, array shape, and dtype for reconstructing the object array.
        """
        with multiprocessing.Pool(self.num_procs) as pool:
            graphs = pool.starmap(
                to_graph_wrapper,
                [
                    (mol, self.threshold, self.adjust_H)
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
        """
        Group molecules using connectivity with shared memory optimization.

        Groups molecules based on molecular connectivity using multiprocessing
        and shared memory for graph storage. Uses iterative comparison with
        a pivot molecule approach to identify structurally equivalent groups.

        Returns:
            Tuple[List[List[Molecule]], List[List[int]]]: Tuple containing:
                - List of molecule groups (each group is a list of molecules)
                - List of index groups (corresponding indices for each group)
        """
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
    """
    Factory for creating molecular grouper instances.

    Provides a unified entry point to construct groupers by name. Supported
    strategies (case-insensitive):
    - "rmsd": RMSDGrouper
    - "tanimoto" or "fingerprint": TanimotoSimilarityGrouper
    - "isomorphism" or "rdkit": RDKitIsomorphismGrouper
    - "formula": FormulaGrouper
    - "connectivity": ConnectivityGrouper

    Additional keyword arguments are forwarded to the specific grouper
    constructors (e.g., thresholds or flags).
    """

    @staticmethod
    def create(
        structures,
        strategy="rmsd",
        num_procs=1,
        threshold=None,
        ignore_hydrogens=None,
        **kwargs,
    ):
        """
        Create a molecular grouper instance by strategy name.

        Args:
            structures: Iterable of `Molecule` to group.
            strategy (str): One of "rmsd", "tanimoto"/"fingerprint",
                "isomorphism"/"rdkit", "formula", or "connectivity".
                Defaults to "rdkit" (alias of RDKitIsomorphismGrouper).
            num_procs (int): Number of workers for parallel computation.
            **kwargs: Extra options forwarded to the grouper constructor
                (e.g., `threshold`, `align_molecules`, `adjust_H`).

        Returns:
            MoleculeGrouper: An instance of the selected grouper subclass.

        Raises:
            ValueError: If `strategy` is not a supported name.
        """
        groupers = {
            "rmsd": BasicRMSDGrouper,
            "hrmsd": HungarianRMSDGrouper,
            "spyrmsd": SpyRMSDGrouper,
            "tanimoto": TanimotoSimilarityGrouper,
            "isomorphism": RDKitIsomorphismGrouper,
            "formula": FormulaGrouper,
            "connectivity": ConnectivityGrouper,
        }

        threshold_supported = {
            "rmsd",
            "spyrmsd",
            "hrmsd",
            "tanimoto",
            "connectivity",
        }
        if strategy in groupers:
            logger.info(f"Using {strategy} grouping strategy.")
            if strategy in threshold_supported:
                if strategy in ["rmsd", "hrmsd", "spyrmsd"]:
                    return groupers[strategy](
                        structures,
                        threshold=threshold,
                        num_procs=num_procs,
                        ignore_hydrogens=ignore_hydrogens,
                        **kwargs,
                    )
                return groupers[strategy](
                    structures,
                    threshold=threshold,
                    num_procs=num_procs,
                    **kwargs,
                )
            else:
                return groupers[strategy](
                    structures,
                    num_procs=num_procs,
                    **kwargs,
                )
        raise ValueError(f"Unknown grouping strategy: {strategy}")
