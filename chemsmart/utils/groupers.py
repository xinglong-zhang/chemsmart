from typing import List, Tuple

import networkx as nx
import numpy as np
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core import Structure
from rdkit import Chem
from rdkit.Chem import rdMolHash
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import reverse_cuthill_mckee

from chemsmart.io.molecules.structure import Molecule


class MoleculeGrouper:
    """Base class for molecular grouping strategies."""

    def __init__(self, molecules: List["Molecule"], **kwargs):
        self.molecules = molecules
        self._validate_molecules()

    def _validate_molecules(self):
        if not all(isinstance(m, Molecule) for m in self.molecules):
            raise TypeError("All inputs must be Molecule instances")

    def group(self) -> Tuple[List[List["Molecule"]], List[List[int]]]:
        """Main grouping method to be implemented by subclasses"""
        raise NotImplementedError

    def unique(self) -> List["Molecule"]:
        """Get unique representatives from groups"""
        groups, _ = self.group()
        return [group[0] for group in groups]


class RDKitMoleculeGrouper(MoleculeGrouper):
    """Group molecules using RDKit-based hashing and isomorphism checks"""

    def __init__(
        self,
        molecules: List["Molecule"],
        use_stereochemistry: bool = True,
        use_tautomers: bool = False,
    ):
        super().__init__(molecules)
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


class RCMMoleculeGrouper(MoleculeGrouper):
    """Group molecules using Reverse Cuthill-McKee algorithm on similarity matrix"""

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
    """Group molecules based on geometric similarity using RMSD"""

    def __init__(
        self,
        molecules: List["Molecule"],
        rmsd_threshold: float = 0.5,
        align_molecules: bool = True,
    ):
        super().__init__(molecules)
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
            pos1, pos2 = self._kabsch_align(pos1, pos2)

        return np.sqrt(np.mean(np.sum((pos1 - pos2) ** 2, axis=1)))

    def _kabsch_align(
        self, P: np.ndarray, Q: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Kabsch algorithm for molecular alignment"""
        # Center molecules
        P_centered = P - np.mean(P, axis=0)
        Q_centered = Q - np.mean(Q, axis=0)

        # Compute covariance matrix
        H = P_centered.T @ Q_centered

        # SVD
        U, S, Vt = np.linalg.svd(H)
        d = np.sign(np.linalg.det(Vt.T @ U.T))

        # Rotation matrix
        R = Vt.T @ np.diag([1, 1, d]) @ U.T

        # Apply rotation
        aligned_P = P_centered @ R
        return aligned_P, Q_centered


class PMGMoleculeGrouper(MoleculeGrouper):
    """Group molecules using pymatgen's StructureMatcher"""

    def __init__(
        self,
        molecules: List["Molecule"],
        ltol: float = 0.2,
        stol: float = 0.3,
        angle_tol: float = 5,
        primitive_cell: bool = False,
    ):
        super().__init__(molecules)
        self.matcher = StructureMatcher(
            ltol=ltol,
            stol=stol,
            angle_tol=angle_tol,
            primitive_cell=primitive_cell,
            scale=False,
            attempt_supercell=False,
        )

    def group(self) -> Tuple[List[List["Molecule"]], List[List[int]]]:
        structures = [self._molecule_to_pmg_struct(m) for m in self.molecules]
        groups = self.matcher.group_structures(structures)

        # Convert back to molecule groups
        mol_groups = []
        idx_groups = []
        for group in groups:
            indices = [structures.index(s) for s in group]
            mol_groups.append([self.molecules[i] for i in indices])
            idx_groups.append(indices)

        return mol_groups, idx_groups

    def _molecule_to_pmg_struct(self, mol: Molecule) -> Structure:
        """Convert Molecule to pymatgen Structure"""
        return Structure(
            lattice=np.eye(3) * 20,  # Large box for molecular systems
            species=mol.chemical_symbols,
            coords=mol.positions,
            coords_are_cartesian=True,
        )


class HybridMoleculeGrouper(MoleculeGrouper):
    """Hybrid grouping strategy combining multiple approaches"""

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
    """Group by chemical formula"""

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
    """Group by molecular connectivity"""

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
    """Specialized geometry-based grouper (inherits from RMSDMoleculeGrouper)"""

    pass
