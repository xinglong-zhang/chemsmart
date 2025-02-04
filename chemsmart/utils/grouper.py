from typing import List

import networkx as nx
import numpy as np

from chemsmart.io.molecules.structure import Molecule


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
