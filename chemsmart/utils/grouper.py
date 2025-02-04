"""Molecular structure grouping using structural similarity matching."""

from typing import List

import networkx as nx
import numpy as np

from chemsmart.io.molecules.structure import Molecule


import logging
import multiprocessing
from functools import partial
from typing import Iterable, List, Tuple

import numpy as np
from ase import Atoms
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.io.ase import AseAtomsAdaptor
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import reverse_cuthill_mckee

logger = logging.getLogger(__name__)


class StructureGrouperConfig:
    """Configuration container for structure matching parameters."""

    def __init__(self, ltol=0.1, stol=0.18, angle_tol=1):
        self.ltol = ltol
        self.stol = stol
        self.angle_tol = angle_tol


class BaseGrouper:
    """Base class for molecular structure grouping."""

    def __init__(self, objects, comparator, num_procs=1):
        self.objects = objects
        self.comparator = comparator
        self.num_procs = max(1, num_procs)
        self._validate_inputs()

    def _validate_inputs(self):
        if not callable(self.comparator):
            raise ValueError("Comparator must be a callable function")
        if not isinstance(self.objects, Iterable):
            raise TypeError("Objects must be an iterable collection")

    def group(self) -> Tuple[List[List[Atoms]], List[List[int]]]:
        """Main grouping interface."""
        raise NotImplementedError

    def unique(self) -> List[Atoms]:
        """Get unique representative structures."""
        groups, _ = self.group()
        return [group[0] for group in groups]


class MatrixGrouper(BaseGrouper):
    """Group structures using similarity matrix approach."""

    def group(self):
        similarity_matrix = self._create_similarity_matrix()
        return self._cluster_similar_structures(similarity_matrix)

    def _create_similarity_matrix(self) -> np.ndarray:
        """Create boolean similarity matrix using parallel processing."""
        n = len(self.objects)
        matrix = np.eye(n, dtype=bool)

        with multiprocessing.Pool(self.num_procs) as pool:
            for i in range(n):
                row = pool.map(
                    partial(self.comparator, self.objects[i]),
                    self.objects[i:]
                )
                matrix[i, i:] = row
                matrix[i:, i] = row

        return matrix | matrix.T  # Ensure symmetry

    def _cluster_similar_structures(self, matrix: np.ndarray):
        """Cluster structures using RCM algorithm."""
        sparse_matrix = csr_matrix(matrix)
        permutation = reverse_cuthill_mckee(sparse_matrix, symmetric_mode=True)
        reordered = sparse_matrix[permutation[:, None], permutation]
        return self._extract_groups(reordered.toarray(), permutation)

    def _extract_groups(self, matrix: np.ndarray, permutation: np.ndarray):
        """Extract similarity groups from ordered matrix."""
        groups = []
        current_group = [permutation[0]]

        for idx in permutation[1:]:
            if matrix[idx, current_group[0]]:
                current_group.append(idx)
            else:
                groups.append(sorted(current_group))
                current_group = [idx]

        if current_group:
            groups.append(sorted(current_group))

        structure_groups = [[self.objects[i] for i in group] for group in groups]
        return structure_groups, groups


class SequentialGrouper(BaseGrouper):
    """Group structures using sequential greedy approach."""

    def group(self):
        indices = list(range(len(self.objects)))
        groups = []

        while indices:
            pivot_idx = indices[0]
            matches = self._find_matches(pivot_idx, indices)
            groups.append(matches)
            indices = [i for i in indices if i not in matches]

        structure_groups = [[self.objects[i] for i in group] for group in groups]
        return structure_groups, groups

    def _find_matches(self, pivot_idx: int, candidates: List[int]) -> List[int]:
        """Find all matches for a pivot structure using parallel processing."""
        with multiprocessing.Pool(self.num_procs) as pool:
            results = pool.map(
                partial(self.comparator, self.objects[pivot_idx]),
                (self.objects[i] for i in candidates)
            )
        return [candidates[i] for i, match in enumerate(results) if match]


class StructureComparator:
    """Handles structure comparison logic."""

    def __init__(self, config: StructureGrouperConfig):
        self.matcher = StructureMatcher(
            ltol=config.ltol,
            stol=config.stol,
            angle_tol=config.angle_tol,
            scale=False,
            primitive_cell = False,  # Disable primitive cell reduction
            attempt_supercell = False,  # Disable supercell matching
        )

    def __call__(self, struct1, struct2) -> bool:
        """Main comparison method."""
        return self._structures_match(struct1, struct2)

    def _structures_match(self, struct1, struct2) -> bool:
        """Core structure matching logic."""
        try:
            return self.matcher.fit(struct1, struct2)
        except Exception as e:
            logger.warning(f"Structure matching failed: {str(e)}")
            return False


class StructureGrouperFactory:
    """Factory for creating appropriate grouper instances."""

    @staticmethod
    def create_grouper(
            structures: Iterable[Atoms],
            config: StructureGrouperConfig,
            strategy: str = "matrix",
            num_procs: int = 1
    ) -> BaseGrouper:
        adaptor = AseAtomsAdaptor()
        pmg_structures = [adaptor.get_structure(a) for a in structures]
        comparator = StructureComparator(config)

        if strategy == "matrix":
            return MatrixGrouper(pmg_structures, comparator, num_procs)
        if strategy == "sequential":
            return SequentialGrouper(pmg_structures, comparator, num_procs)
        raise ValueError(f"Unknown grouping strategy: {strategy}")


class MolecularGrouper:
    def __init__(self, rmsd_threshold=0.2, bond_cutoff_buffer=0.3):
        self.rmsd_threshold = rmsd_threshold
        self.bond_cutoff_buffer = bond_cutoff_buffer

    def _calculate_rmsd(self, mol1: Molecule, mol2: Molecule) -> float:
        """Calculate RMSD between two Molecule objects."""
        pos1 = mol1.positions
        pos2 = mol2.positions

        # Center molecules
        pos1_centered = pos1 - np.mean(pos1, axis=0)
        pos2_centered = pos2 - np.mean(pos2, axis=0)

        # Calculate RMSD
        return np.sqrt(
            np.mean(np.sum((pos1_centered - pos2_centered) ** 2, axis=1))
        )

    def _are_similar(self, mol1: Molecule, mol2: Molecule) -> bool:
        """Check if two molecules are similar using composition, connectivity, and geometry."""
        # Check basic composition
        if sorted(mol1.chemical_symbols) != sorted(mol2.chemical_symbols):
            return False

        # Check graph isomorphism
        graph1 = mol1.to_graph(bond_cutoff_buffer=self.bond_cutoff_buffer)
        graph2 = mol2.to_graph(bond_cutoff_buffer=self.bond_cutoff_buffer)
        if not nx.is_isomorphic(
            graph1,
            graph2,
            node_match=lambda a, b: a["element"] == b["element"],
        ):
            return False

        # Check geometric similarity
        return self._calculate_rmsd(mol1, mol2) < self.rmsd_threshold

    def group_molecules(
        self, molecules: List[Molecule]
    ) -> List[List[Molecule]]:
        """Group similar molecules using a greedy algorithm."""
        groups = []
        remaining = set(range(len(molecules)))

        while remaining:
            pivot_idx = remaining.pop()
            group = [molecules[pivot_idx]]
            group_indices = [pivot_idx]

            to_remove = set()
            for idx in remaining:
                if self._are_similar(molecules[pivot_idx], molecules[idx]):
                    group.append(molecules[idx])
                    group_indices.append(idx)
                    to_remove.add(idx)

            remaining -= to_remove
            groups.append(group)

        return groups

    def get_unique_representatives(
        self, molecules: List[Molecule]
    ) -> List[Molecule]:
        """Get unique representative molecules from each group."""
        groups = self.group_molecules(molecules)
        return [group[0] for group in groups]
