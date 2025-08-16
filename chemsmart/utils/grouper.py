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
from rdkit.Chem import rdMolHash
from rdkit.Chem.rdFingerprintGenerator import GetRDKitFPGenerator
from scipy.optimize import linear_sum_assignment  # Hungarian algorithm
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
        self.num_procs = int(max(1, num_procs))

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
    """Base class for RMSD-based molecular grouping.

    This abstract class provides common functionality for RMSD-based molecular grouping,
    including threshold management, heavy atom filtering, and parallel processing setup.
    """

    def __init__(
        self,
        molecules: Iterable[Molecule],
        threshold: float = 1.5,  # RMSD threshold for grouping
        num_procs: int = 1,
        align_molecules: bool = True,
        ignore_hydrogens: bool = False,  # option to ignore H atoms for grouping
    ):
        super().__init__(molecules, num_procs)
        self.threshold = threshold
        self.align_molecules = align_molecules
        self.ignore_hydrogens = ignore_hydrogens

    def _get_heavy_atoms(self, mol):
        """Extract heavy atom positions and symbols."""
        heavy_indices = [i for i, s in enumerate(mol.symbols) if s != "H"]
        positions = mol.positions[heavy_indices]
        symbols = [mol.symbols[i] for i in heavy_indices]
        return positions, symbols

    @abstractmethod
    def _calculate_rmsd(self, idx_pair: Tuple[int, int]) -> float:
        """Calculate RMSD between two molecules. Must be implemented by subclasses."""
        pass

    def group(self) -> Tuple[List[List[Molecule]], List[List[int]]]:
        """Group molecules by geometric similarity using connected components."""
        n = len(self.molecules)
        indices = [(i, j) for i in range(n) for j in range(i + 1, n)]

        # Use map for better parallelism
        with multiprocessing.Pool(self.num_procs) as pool:
            rmsd_values = pool.map(self._calculate_rmsd, indices)

        # Build adjacency matrix
        adj_matrix = np.zeros((n, n), dtype=bool)
        for (i, j), rmsd in zip(indices, rmsd_values):
            if rmsd < self.threshold:
                adj_matrix[i, j] = adj_matrix[j, i] = True

        # Find connected components
        _, labels = connected_components(csr_matrix(adj_matrix))

        # Group molecules and indices
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


class RMSDGrouperSymmetric(RMSDGrouper):
    """Groups molecules based on RMSD (Root Mean Square Deviation) with symmetry consideration.

    This version handles molecular symmetries by finding optimal atom correspondence
    using the Hungarian algorithm before alignment.
    References are as follows: http://dx.doi.org/10.1021/ci400534h,
    https://doi.org/10.1186/s13321-020-00455-2.
    It's designed for scenarios where
    atom ordering might differ between conformers, such as CREST-generated structures
    or molecules with equivalent atoms in different arrangements.
    Use "rsmd" or "rsmd_symmetric" to utilize this default RMSDGrouper.

    Attributes:
        threshold (float): RMSD threshold for grouping molecules (Å).
        num_procs (int): Number of processes for parallel computation.
        align_molecules (bool): Whether to perform Kabsch alignment.
        ignore_hydrogens (bool): Whether to exclude hydrogen atoms from RMSD calculation.

    Performance Notes:
        - Hungarian algorithm is computationally expensive for large molecules
        - Parallel processing recommended for multiple conformers
        - Memory usage scales with O(n²) for pairwise comparisons

    Best suited for:
        - CREST conformer analysis
        - Molecules with symmetrical arrangements
        - Cases where atom ordering is inconsistent
        - High-accuracy structural comparisons

    Guideline for threshold values:
        Molecule type / context	                    Typical RMSD threshold for grouping
        Small rigid molecules                        0.5–1.0 Å
        (≤20 heavy atoms, little flexibility)	     – distinguishes small conformational changes clearly

        Flexible drug-like molecules                 1.0–2.0 Å
        (rotatable bonds, ~20–50 heavy atoms)	     – allows grouping of conformers differing mainly in torsions

        Large biomolecules or peptides	             2.0–3.0 Å for same fold;
                                                     ≥3 Å often indicates different conformations/folds

        Coarse clustering for diverse sets	         2.0–4.0 Å
                                                     – groups into broad families rather than fine differences
    """

    def __init__(
        self,
        molecules: Iterable[Molecule],
        threshold: float = 1.5,  # RMSD threshold for grouping
        num_procs: int = 1,
        align_molecules: bool = True,
        ignore_hydrogens: bool = False,  # option to ignore H atoms for grouping
        use_fallback: bool = True,
    ):
        super().__init__(
            molecules, threshold, num_procs, align_molecules, ignore_hydrogens
        )
        self.use_fallback = use_fallback
        # Cache sorted chemical symbols as sets for faster comparison
        self._chemical_symbol_sets = [
            set(mol.chemical_symbols) for mol in molecules
        ]

    def _calculate_rmsd(self, idx_pair):
        i, j = idx_pair
        mol1, mol2 = self.molecules[i], self.molecules[j]

        if self.ignore_hydrogens:
            pos1, symbols1 = self._get_heavy_atoms(mol1)
            pos2, symbols2 = self._get_heavy_atoms(mol2)
        else:
            pos1, symbols1 = mol1.positions, list(mol1.symbols)
            pos2, symbols2 = mol2.positions, list(mol2.symbols)

        if len(symbols1) != len(symbols2) or sorted(symbols1) != sorted(
            symbols2
        ):
            return np.inf

        if (
            self.use_fallback
            and not self.ignore_hydrogens
            and self.align_molecules
        ):
            _, _, _, _, rmsd_direct = kabsch_align(pos1, pos2)
            # If RMSD is very small, the molecular ordering is likely correct
            if (
                rmsd_direct < 1e-6
            ):  # Nearly zero, molecules are identical conformations
                return rmsd_direct

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
                dist_matrix = np.zeros((len(idxs1), len(idxs2)))
                for a, idx1 in enumerate(idxs1):
                    for b, idx2 in enumerate(idxs2):
                        dist_matrix[a, b] = np.sum(
                            (pos1[idx1] - pos2[idx2]) ** 2
                        )
                row_ind, col_ind = linear_sum_assignment(dist_matrix)
                matched_idx1.extend([idxs1[r] for r in row_ind])
                matched_idx2.extend([idxs2[c] for c in col_ind])

        pos1_matched = pos1[matched_idx1]
        pos2_matched = pos2[matched_idx2]

        # Apply alignment if requested
        if self.align_molecules:
            pos1_matched, pos2_matched, _, _, _ = kabsch_align(
                pos1_matched, pos2_matched
            )

        rmsd = np.sqrt(
            np.mean(np.sum((pos1_matched - pos2_matched) ** 2, axis=1))
        )
        return rmsd


class RMSDGrouperSimple(RMSDGrouper):
    """Simple RMSD grouper without symmetry consideration for fast molecular comparison.

    This lightweight version assumes atoms are in the same order and directly applies
    Kabsch alignment without atom reordering. It's optimized for scenarios where
    molecular conformations maintain consistent atom ordering, such as:
    - Rotated/translated versions of the same molecule
    - MD trajectory frames with preserved atom indices
    - Molecules generated from the same starting structure
    Use "rsmd_simple" to utilize this RMSDGrouper.

    Best suited for:
        - Molecular dynamics trajectory clustering
        - Conformer screening with preserved atom order
        - High-throughput virtual screening
        - Cases where speed is prioritized over symmetry handling

    Warning:
        Although it is much faster than RMSDGrouperSymmetric，it will produce incorrect results
        for symmetric molecules or when conformers have mismatched atom ordering.
    """

    def _calculate_rmsd(self, idx_pair: Tuple[int, int]) -> float:
        """Calculate RMSD between two molecules without atom reordering."""
        i, j = idx_pair
        mol1, mol2 = self.molecules[i], self.molecules[j]

        if self.ignore_hydrogens:
            pos1, symbols1 = self._get_heavy_atoms(mol1)
            pos2, symbols2 = self._get_heavy_atoms(mol2)
        else:
            pos1, symbols1 = mol1.positions, list(mol1.symbols)
            pos2, symbols2 = mol2.positions, list(mol2.symbols)

        # Quick check for compatibility
        if len(symbols1) != len(symbols2) or sorted(symbols1) != sorted(
            symbols2
        ):
            return np.inf

        # Direct alignment without atom reordering - assumes same atom ordering
        if self.align_molecules:
            _, _, _, _, rmsd = kabsch_align(pos1, pos2)
            return rmsd
        else:
            # No alignment, direct RMSD calculation
            return np.sqrt(np.mean(np.sum((pos1 - pos2) ** 2, axis=1)))


class RMSDGrouperSharedMemory(MoleculeGrouper):
    """Group molecules based on RMSD using shared memory with minimal locking.
    Uses symmetric RMSD grouper.

    This implementation uses multiprocessing.RawArray for shared memory to minimize
    locking overhead and improve parallel performance for large molecule sets.
    """

    def __init__(
        self,
        molecules: Iterable[Molecule],
        threshold: float = 1.5,
        num_procs: int = 1,
        align_molecules: bool = True,
        ignore_hydrogens: bool = False,
    ):
        super().__init__(molecules, num_procs)
        self.threshold = threshold
        self.align_molecules = align_molecules
        self.ignore_hydrogens = ignore_hydrogens
        # Store symbols list for each molecule
        self.symbols_list = [list(mol.symbols) for mol in self.molecules]

    def group(self) -> Tuple[List[List[Molecule]], List[List[int]]]:
        """Group molecules using shared memory with optimized parallelism."""
        n = len(self.molecules)
        if n == 0:
            return [], []

        # Check if all molecules have the same number of atoms
        num_atoms = self.molecules[0].positions.shape[0]
        if not all(
            mol.positions.shape[0] == num_atoms for mol in self.molecules
        ):
            raise ValueError(
                "All molecules must have the same number of atoms for shared memory approach"
            )

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
            initargs=(shared_pos, (n, num_atoms, 3), self.symbols_list),
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
    def _init_worker(shared_pos, pos_shape, symbols_list):
        """Worker process initializer to attach shared memory."""
        global shared_positions, shared_symbols_list
        n, num_atoms, _ = pos_shape
        shared_positions = np.frombuffer(shared_pos, dtype=np.float64).reshape(
            n, num_atoms, 3
        )
        shared_symbols_list = symbols_list

    def _calculate_rmsd(self, idx_pair: Tuple[int, int]) -> float:
        """Calculate RMSD efficiently using local copies of shared memory."""
        i, j = idx_pair

        pos1 = np.array(shared_positions[i])
        pos2 = np.array(shared_positions[j])

        symbols1 = shared_symbols_list[i]
        symbols2 = shared_symbols_list[j]

        # Apply hydrogen filtering if requested
        if self.ignore_hydrogens:
            heavy_indices1 = [k for k, s in enumerate(symbols1) if s != "H"]
            heavy_indices2 = [k for k, s in enumerate(symbols2) if s != "H"]
            pos1 = pos1[heavy_indices1]
            pos2 = pos2[heavy_indices2]
            symbols1 = [symbols1[k] for k in heavy_indices1]
            symbols2 = [symbols2[k] for k in heavy_indices2]

        # Quick compatibility check
        if len(symbols1) != len(symbols2) or sorted(symbols1) != sorted(
            symbols2
        ):
            return np.inf

        matched_idx1 = []
        matched_idx2 = []

        for elem in sorted(set(symbols1)):
            idxs1 = [k for k, s in enumerate(symbols1) if s == elem]
            idxs2 = [k for k, s in enumerate(symbols2) if s == elem]

            if len(idxs1) == 1 and len(idxs2) == 1:
                matched_idx1.extend(idxs1)
                matched_idx2.extend(idxs2)
            else:
                dist_matrix = np.sum(
                    (pos1[idxs1, None, :] - pos2[None, idxs2, :]) ** 2, axis=2
                )
                row_ind, col_ind = linear_sum_assignment(dist_matrix)
                matched_idx1.extend([idxs1[r] for r in row_ind])
                matched_idx2.extend([idxs2[c] for c in col_ind])

        pos1_matched = pos1[matched_idx1]
        pos2_matched = pos2[matched_idx2]

        if self.align_molecules:
            pos1_matched, pos2_matched, _, _, _ = kabsch_align(
                pos1_matched, pos2_matched
            )

        rmsd = np.sqrt(
            np.mean(np.sum((pos1_matched - pos2_matched) ** 2, axis=1))
        )
        return rmsd

    def __repr__(self):
        return (
            f"{self.__class__.__name__}(threshold={self.threshold}, "
            f"num_procs={self.num_procs}, align_molecules={self.align_molecules}, "
            f"ignore_hydrogens={self.ignore_hydrogens})"
        )


class TanimotoSimilarityGrouper(MoleculeGrouper):
    """Groups molecules based on fingerprint similarity using Tanimoto coefficient.

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
        threshold: float = 0.9,  # Tanimoto similarity threshold
        num_procs: int = 1,
        use_rdkit_fp: bool = True,  # Allows switching between RDKit FP and RDKFingerprint
    ):
        super().__init__(molecules, num_procs)
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
        """Generate an RDKit fingerprint for a molecule."""
        try:
            if self.use_rdkit_fp:
                return GetRDKitFPGenerator().GetFingerprint(rdkit_mol)
            else:
                return Chem.RDKFingerprint(rdkit_mol)  # Alternative method
        except Exception as e:
            logger.warning(f"Fingerprint generation failed: {str(e)}")
            return None

    def group(self) -> Tuple[List[List[Molecule]], List[List[int]]]:
        """Groups molecules based on Tanimoto similarity using connected components clustering."""
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
        threshold: float = 0.0,  # Buffer for bond cutoff
        adjust_H: bool = True,
    ):
        super().__init__(molecules, num_procs)
        self.threshold = threshold  # Buffer for bond cutoff
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
    """Group molecules based on molecular connectivity."""

    def __init__(
        self,
        molecules: Iterable[Molecule],
        num_procs: int = 1,
        threshold: float = 0.0,  # Buffer for bond cutoff
        adjust_H: bool = True,
    ):
        super().__init__(molecules, num_procs)
        self.threshold = threshold
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
    def create(
        structures, strategy="rmsd", num_procs=1, threshold=5.0, **kwargs
    ):
        groupers = {
            "rmsd": RMSDGrouperSymmetric,  # Use symmetric version as default
            "rmsd_simple": RMSDGrouperSimple,  # Simple version without symmetry
            "rmsd_symmetric": RMSDGrouperSymmetric,  # Explicit symmetric version
            "tanimoto": TanimotoSimilarityGrouper,
            "isomorphism": RDKitIsomorphismGrouper,
            "formula": FormulaGrouper,
            "connectivity": ConnectivityGrouper,
        }
        if strategy in groupers:
            logger.info(f"Using {strategy} grouping strategy.")
            if strategy.startswith("rmsd"):
                return groupers[strategy](
                    structures,
                    threshold=threshold,
                    num_procs=num_procs,
                    **kwargs,
                )
            else:
                return groupers[strategy](
                    structures, num_procs=num_procs, **kwargs
                )
        raise ValueError(f"Unknown grouping strategy: {strategy}")
