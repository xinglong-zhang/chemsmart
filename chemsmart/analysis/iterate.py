import logging

import networkx as nx
import numpy as np
from ase.data import atomic_numbers as ase_atomic_numbers
from ase.data import covalent_radii
from rdkit import Chem
from scipy.optimize import minimize

from chemsmart.io.molecules.structure import Molecule
from chemsmart.utils.periodictable import PeriodicTable

logger = logging.getLogger(__name__)

# Default buffer for bond cutoff calculations in Angstroms
DEFAULT_BUFFER = 0.3

# PeriodicTable instance for symbol/atomic number conversion
pt = PeriodicTable()


class SkeletonPreprocessor:
    """Preprocessor to prepare skeleton molecule by removing substituent at link position.
    
    This class handles two scenarios:
    1. User provides skeleton_indices: Keep only atoms at specified indices
    2. User doesn't provide skeleton_indices: Auto-detect and remove substituent at link position
    """
    
    def __init__(
            self,
            molecule: Molecule,
            link_index: int,
            skeleton_indices: list[int] | None = None,
        ):
        """
        Initialize SkeletonPreprocessor.
        
        Parameters
        ----------
        molecule : Molecule
            Input molecule object
        link_index : int
            Index of the link atom (1-based, will be converted to 0-based internally)
        skeleton_indices : list[int] | None
            Indices of atoms belonging to skeleton (1-based).
            If None, auto-detection will be used.
        """
        self.molecule = molecule
        # Convert 1-based to 0-based index
        self.link_index = link_index - 1
        # Convert skeleton_indices to 0-based if provided
        self.skeleton_indices = None
        if skeleton_indices is not None:
            self.skeleton_indices = [i - 1 for i in skeleton_indices]
    
    def run(self) -> Molecule:
        """
        Execute preprocessing to get clean skeleton molecule.
        
        Returns
        -------
        Molecule
            Skeleton molecule with substituent removed.
        """
        # First check if link_index atom has available bonding position
        if self._has_available_bonding_position():
            logger.info(
                f"Link atom at index {self.link_index + 1} (1-based) has available bonding position. "
                "No substituent removal needed."
            )
            return self.molecule
        
        if self.skeleton_indices is not None:
            # Mode 1: User provided skeleton indices
            # Check if link_index is in skeleton_indices
            if self.link_index not in self.skeleton_indices:
                logger.warning(
                    f"Link atom at index {self.link_index + 1} (1-based) is not in skeleton_indices. "
                    "Falling back to auto-detection mode."
                )
                return self._run_auto_detect()
            
            # Find substituent branches that don't contain skeleton atoms
            substituent_branches = self._find_non_skeleton_branches()
            
            if len(substituent_branches) == 0:
                # No substituent to remove at this link_index, return molecule as-is
                return self.molecule
            elif len(substituent_branches) == 1:
                # Exactly one substituent branch - remove only this branch
                # Keep all other atoms (including substituents on other link atoms)
                substituent_indices = substituent_branches[0]
                keep_indices = self._get_complement_indices(substituent_indices)
                return self._extract_by_indices(keep_indices)
            else:
                # Multiple substituent branches - ambiguous, fall back to auto-detection
                logger.warning(
                    f"Found {len(substituent_branches)} branches without skeleton atoms. "
                    "Falling back to auto-detection mode."
                )
                return self._run_auto_detect()
        else:
            # Mode 2: Auto-detect substituent
            return self._run_auto_detect()
    
    def _has_available_bonding_position(self) -> bool:
        """
        Check if the link atom has an available bonding position.
        
        This is determined by comparing the actual number of neighbors
        with the expected maximum bonding capacity of the element.
        
        Returns
        -------
        bool
            True if link atom has available position for new bond.
        """
        # Build molecular graph
        graph = self.molecule.to_graph()
        
        # Get current number of neighbors (bonds)
        current_bonds = len(list(graph.neighbors(self.link_index)))
        
        # Get element symbol
        element = self.molecule.chemical_symbols[self.link_index]
        
        # Get expected maximum bonding capacity
        max_bonds = self._get_max_bonding_capacity(element)
        
        return current_bonds < max_bonds
    
    @staticmethod
    def _get_max_bonding_capacity(element: str) -> int:
        """
        Get the maximum bonding capacity for an element.
        
        Uses RDKit's periodic table to get the default valence.
        
        Parameters
        ----------
        element : str
            Element symbol
        
        Returns
        -------
        int
            Maximum number of bonds the element can form
        """
        
        periodic_table = Chem.GetPeriodicTable()
        atomic_num = periodic_table.GetAtomicNumber(element)
        
        # GetDefaultValence returns a tuple of possible valences, take the max
        default_valence = periodic_table.GetDefaultValence(atomic_num)
        
        if isinstance(default_valence, tuple):
            return max(default_valence)
        return default_valence
    
    def _run_auto_detect(self) -> Molecule:
        """
        Run auto-detection mode to find and remove substituent.
        
        Returns
        -------
        Molecule
            Skeleton molecule with substituent removed.
        """
        substituent_indices = self.detect_substituent()
        skeleton_indices = self._get_complement_indices(substituent_indices)
        return self._extract_by_indices(skeleton_indices)
    
    def _find_non_skeleton_branches(self) -> list[list[int]]:
        """
        Find all branches from link_index that don't contain any skeleton atoms.
        
        Uses DFS traversal from link_index to explore each neighbor branch.
        A branch is considered a "non-skeleton branch" (substituent) if none of
        its atoms are in skeleton_indices.
        
        Returns
        -------
        list[list[int]]
            List of branches (each branch is a list of atom indices, 0-based)
            that don't contain skeleton atoms.
        """
        # Build molecular graph
        graph = self.molecule.to_graph()
        
        # Get neighbors of link atom
        neighbors = list(graph.neighbors(self.link_index))
        
        if len(neighbors) == 0:
            logger.warning(f"Link atom at index {self.link_index} has no neighbors")
            return []
        
        skeleton_set = set(self.skeleton_indices)
        non_skeleton_branches = []
        
        for neighbor in neighbors:
            # DFS to collect all atoms in this branch (excluding link_index)
            branch_atoms = self._dfs_collect_branch(graph, neighbor, self.link_index)
            
            # Check if any atom in this branch is in skeleton_indices
            has_skeleton_atom = any(atom in skeleton_set for atom in branch_atoms)
            
            if not has_skeleton_atom:
                # This branch has no skeleton atoms - it's a substituent
                non_skeleton_branches.append(branch_atoms)
        
        return non_skeleton_branches
    
    def _dfs_collect_branch(
        self, 
        graph: nx.Graph, 
        start: int, 
        excluded: int
    ) -> list[int]:
        """
        Collect all atoms in a branch using DFS, excluding the starting point's parent.
        
        Parameters
        ----------
        graph : nx.Graph
            Molecular connectivity graph
        start : int
            Starting atom index (0-based)
        excluded : int
            Atom index to exclude (the link atom, 0-based)
        
        Returns
        -------
        list[int]
            All atom indices in this branch (0-based)
        """
        visited = set()
        stack = [start]
        branch_atoms = []
        
        while stack:
            node = stack.pop()
            if node in visited or node == excluded:
                continue
            visited.add(node)
            branch_atoms.append(node)
            
            # Add unvisited neighbors to stack
            for neighbor in graph.neighbors(node):
                if neighbor not in visited and neighbor != excluded:
                    stack.append(neighbor)
        
        return branch_atoms
    
    def detect_substituent(self) -> list[int]:
        """
        Auto-detect substituent atoms connected to link atom.
        
        The substituent is identified as the smallest connected component
        when the bond between link_atom and its neighbor is broken.
        
        Returns
        -------
        list[int]
            Indices of atoms belonging to the substituent (0-based).
        """
        # Build molecular graph
        graph = self.molecule.to_graph()
        
        # Get neighbors of link atom
        neighbors = list(graph.neighbors(self.link_index))
        
        if len(neighbors) == 0:
            logger.warning(f"Link atom at index {self.link_index} has no neighbors")
            return []
        
        # For each neighbor, calculate the size of the component if bond is broken
        # The substituent is the smallest component
        min_component_size = float('inf')
        substituent_indices = []
        
        for neighbor in neighbors:
            # Temporarily remove the edge
            graph_copy = graph.copy()
            graph_copy.remove_edge(self.link_index, neighbor)
            
            # Find connected components
            components = list(nx.connected_components(graph_copy))
            
            # Find which component contains the neighbor (not the link atom)
            for component in components:
                if neighbor in component and self.link_index not in component:
                    # This component is a potential substituent
                    if len(component) < min_component_size:
                        min_component_size = len(component)
                        substituent_indices = list(component)
                    break
        
        return substituent_indices
    
    def _get_complement_indices(self, exclude_indices: list[int]) -> list[int]:
        """
        Get indices of all atoms except those in exclude_indices.
        
        Parameters
        ----------
        exclude_indices : list[int]
            Indices to exclude (0-based)
        
        Returns
        -------
        list[int]
            Complement indices (0-based)
        """
        exclude_set = set(exclude_indices)
        return [i for i in range(len(self.molecule)) if i not in exclude_set]
    
    def _extract_by_indices(self, indices: list[int]) -> Molecule:
        """
        Extract a subset of atoms from molecule by indices.
        
        Parameters
        ----------
        indices : list[int]
            Indices of atoms to keep (0-based)
        
        Returns
        -------
        Molecule
            New molecule containing only specified atoms
        """
        if not indices:
            raise ValueError("Cannot create molecule with no atoms")
        
        # Sort indices to maintain order
        sorted_indices = sorted(indices)
        
        # Extract symbols
        symbols = [self.molecule.chemical_symbols[i] for i in sorted_indices]
        
        # Extract positions
        positions = self.molecule.positions[sorted_indices]
        
        # Extract frozen_atoms if present
        frozen_atoms = None
        if self.molecule.frozen_atoms is not None:
            frozen_atoms = [self.molecule.frozen_atoms[i] for i in sorted_indices]
        
        return Molecule(
            symbols=symbols,
            positions=positions,
            charge=self.molecule.charge,
            multiplicity=self.molecule.multiplicity,
            frozen_atoms=frozen_atoms,
        )
    
    def get_new_link_index(self) -> int:
        """
        Get the new link index in the processed skeleton molecule (1-based).
        
        After extracting skeleton atoms, the original link_index may change.
        This method returns the new index.
        
        Returns
        -------
        int
            New link index (1-based)
        """
        if self.skeleton_indices is not None:
            indices = sorted(self.skeleton_indices)
        else:
            substituent_indices = self.detect_substituent()
            indices = sorted(self._get_complement_indices(substituent_indices))
        
        # Find position of link_index in sorted indices
        try:
            new_index_0based = indices.index(self.link_index)
            return new_index_0based + 1  # Convert to 1-based
        except ValueError:
            raise ValueError(
                f"Link atom at index {self.link_index + 1} (1-based) is not in skeleton indices"
            )


class SubstituentPreprocessor:
    """Preprocessor to prepare substituent molecule by removing atom/group at link position if needed.
    
    This class checks if the link atom has an available bonding position.
    If not, it auto-detects and removes the smallest substituent group at the link position.
    
    Unlike SkeletonPreprocessor, SubstituentPreprocessor does not have a "skeleton_indices" concept.
    It only operates in auto-detect mode when preprocessing is needed.
    """
    
    def __init__(
            self,
            molecule: Molecule,
            link_index: int,
        ):
        """
        Initialize SubstituentPreprocessor.
        
        Parameters
        ----------
        molecule : Molecule
            Input substituent molecule object
        link_index : int
            Index of the link atom (1-based, will be converted to 0-based internally)
        """
        self.molecule = molecule
        # Convert 1-based to 0-based index
        self.link_index = link_index - 1
    
    def run(self) -> Molecule:
        """
        Execute preprocessing to get clean substituent molecule.
        
        If the link atom has an available bonding position, the molecule is
        returned as-is. Otherwise, the smallest connected group at the link
        position is removed.
        
        Returns
        -------
        Molecule
            Substituent molecule ready for attachment.
        """
        # First check if link_index atom has available bonding position
        if self._has_available_bonding_position():
            logger.info(
                f"Substituent link atom at index {self.link_index + 1} (1-based) has available bonding position. "
                "No removal needed."
            )
            return self.molecule
        
        # Auto-detect and remove the smallest group at link position
        logger.info(
            f"Substituent link atom at index {self.link_index + 1} (1-based) has no available bonding position. "
            "Auto-detecting and removing group."
        )
        return self._run_auto_detect()
    
    def _has_available_bonding_position(self) -> bool:
        """
        Check if the link atom has an available bonding position.
        
        This is determined by comparing the actual number of neighbors
        with the expected maximum bonding capacity of the element.
        
        Returns
        -------
        bool
            True if link atom has available position for new bond.
        """
        # Build molecular graph
        graph = self.molecule.to_graph()
        
        # Get current number of neighbors (bonds)
        current_bonds = len(list(graph.neighbors(self.link_index)))
        
        # Get element symbol
        element = self.molecule.chemical_symbols[self.link_index]
        
        # Get expected maximum bonding capacity
        max_bonds = self._get_max_bonding_capacity(element)
        
        return current_bonds < max_bonds
    
    @staticmethod
    def _get_max_bonding_capacity(element: str) -> int:
        """
        Get the maximum bonding capacity for an element.
        
        Uses RDKit's periodic table to get the default valence.
        
        Parameters
        ----------
        element : str
            Element symbol
        
        Returns
        -------
        int
            Maximum number of bonds the element can form
        """
        periodic_table = Chem.GetPeriodicTable()
        atomic_num = periodic_table.GetAtomicNumber(element)
        
        # GetDefaultValence returns a tuple of possible valences, take the max
        default_valence = periodic_table.GetDefaultValence(atomic_num)
        
        if isinstance(default_valence, tuple):
            return max(default_valence)
        return default_valence
    
    def _run_auto_detect(self) -> Molecule:
        """
        Run auto-detection mode to find and remove the smallest group at link position.
        
        Returns
        -------
        Molecule
            Substituent molecule with group removed.
        """
        removed_indices = self._detect_group_to_remove()
        keep_indices = self._get_complement_indices(removed_indices)
        return self._extract_by_indices(keep_indices)
    
    def _detect_group_to_remove(self) -> list[int]:
        """
        Auto-detect the smallest group connected to link atom to remove.
        
        The group to remove is identified as the smallest connected component
        when the bond between link_atom and its neighbor is broken.
        
        Returns
        -------
        list[int]
            Indices of atoms to remove (0-based).
        """
        # Build molecular graph
        graph = self.molecule.to_graph()
        
        # Get neighbors of link atom
        neighbors = list(graph.neighbors(self.link_index))
        
        if len(neighbors) == 0:
            logger.warning(f"Substituent link atom at index {self.link_index} has no neighbors")
            return []
        
        # For each neighbor, calculate the size of the component if bond is broken
        # The group to remove is the smallest component (not containing link atom)
        min_component_size = float('inf')
        remove_indices = []
        
        for neighbor in neighbors:
            # Temporarily remove the edge
            graph_copy = graph.copy()
            graph_copy.remove_edge(self.link_index, neighbor)
            
            # Find connected components
            components = list(nx.connected_components(graph_copy))
            
            # Find which component contains the neighbor (not the link atom)
            for component in components:
                if neighbor in component and self.link_index not in component:
                    # This component is a potential group to remove
                    if len(component) < min_component_size:
                        min_component_size = len(component)
                        remove_indices = list(component)
                    break
        
        return remove_indices
    
    def _get_complement_indices(self, exclude_indices: list[int]) -> list[int]:
        """
        Get indices of all atoms except those in exclude_indices.
        
        Parameters
        ----------
        exclude_indices : list[int]
            Indices to exclude (0-based)
        
        Returns
        -------
        list[int]
            Complement indices (0-based)
        """
        exclude_set = set(exclude_indices)
        return [i for i in range(len(self.molecule)) if i not in exclude_set]
    
    def _extract_by_indices(self, indices: list[int]) -> Molecule:
        """
        Extract a subset of atoms from molecule by indices.
        
        Parameters
        ----------
        indices : list[int]
            Indices of atoms to keep (0-based)
        
        Returns
        -------
        Molecule
            New molecule containing only specified atoms
        """
        if not indices:
            raise ValueError("Cannot create substituent molecule with no atoms")
        
        # Sort indices to maintain order
        sorted_indices = sorted(indices)
        
        # Extract symbols
        symbols = [self.molecule.chemical_symbols[i] for i in sorted_indices]
        
        # Extract positions
        positions = self.molecule.positions[sorted_indices]
        
        # Extract frozen_atoms if present
        frozen_atoms = None
        if self.molecule.frozen_atoms is not None:
            frozen_atoms = [self.molecule.frozen_atoms[i] for i in sorted_indices]
        
        return Molecule(
            symbols=symbols,
            positions=positions,
            charge=self.molecule.charge,
            multiplicity=self.molecule.multiplicity,
            frozen_atoms=frozen_atoms,
        )
    
    def get_new_link_index(self) -> int:
        """
        Get the new link index in the processed substituent molecule (1-based).
        
        After removing atoms, the original link_index may change.
        This method returns the new index.
        
        Returns
        -------
        int
            New link index (1-based)
        """
        if self._has_available_bonding_position():
            # No removal happened, index unchanged
            return self.link_index + 1
        
        removed_indices = self._detect_group_to_remove()
        indices = sorted(self._get_complement_indices(removed_indices))
        
        # Find position of link_index in sorted indices
        try:
            new_index_0based = indices.index(self.link_index)
            return new_index_0based + 1  # Convert to 1-based
        except ValueError:
            raise ValueError(
                f"Link atom at index {self.link_index + 1} (1-based) was removed during preprocessing"
            )


class IterateAnalyzer:
    """Analyzer for a pair of skeleton and substituent in Iterate job, find the optimal position to place the substituent on skeleton."""
    
    def __init__(
            self,
            skeleton: Molecule,
            substituent: Molecule,
            skeleton_link_index: int,
            substituent_link_index: int,
            buffer: float = DEFAULT_BUFFER,
            algorithm: str = 'lagrange_multipliers',
        ):
        """
        Initialize IterateAnalyzer.
        
        Parameters
        ----------
        skeleton : Molecule
            Skeleton molecule object
        substituent : Molecule
            Substituent molecule object
        skeleton_link_index : int
            Index of the link atom in skeleton (1-based, will be converted to 0-based internally)
        substituent_link_index : int
            Index of the link atom in substituent (1-based, will be converted to 0-based internally)
        buffer : float
            Buffer for min distance constraint (default: 0.3 Å)
        algorithm : str
            Optimization algorithm to use (default: 'lagrange_multipliers')
            Supported: 'lagrange_multipliers'
        """
        self.skeleton = skeleton
        self.substituent = substituent
        # Convert 1-based to 0-based index
        self.skeleton_link_index = skeleton_link_index - 1
        self.substituent_link_index = substituent_link_index - 1
        self.buffer = buffer
        self.algorithm = algorithm
    
    def run(self) -> Molecule:
        """
        Execute the iterate analysis to find optimal substituent position.
        
        Returns
        -------
        Molecule
            Combined molecule with skeleton and optimally positioned substituent.
        """
        # Convert Molecule to np.ndarray [atomic_number, x, y, z]
        skeleton_arr = self._molecule_to_array(self.skeleton)
        sub_arr = self._molecule_to_array(self.substituent)
        
        # Find optimal position for substituent
        sub_optimal_arr = self._find_optimal_position(
            skeleton_arr,
            sub_arr,
            self.skeleton_link_index,
            self.substituent_link_index,
            self.buffer,
            self.algorithm
        )
        
        # Update substituent positions with optimized positions
        self.substituent.positions = sub_optimal_arr[:, 1:4]
        
        # Combine skeleton and optimized substituent
        combined_mol = self._combine_molecules(self.skeleton, self.substituent)
        
        return combined_mol
    
    @staticmethod
    def _update_molecule_positions(mol: Molecule, new_positions: np.ndarray) -> Molecule:
        """
        Create a new Molecule with updated positions, preserving other attributes.
        
        Parameters
        ----------
        mol : Molecule
            Original molecule object
        new_positions : np.ndarray
            New positions array, shape (n, 3)
        
        Returns
        -------
        Molecule
            New molecule with updated positions
        """
        return Molecule(
            symbols=mol.symbols,
            positions=new_positions,
            charge=mol.charge,
            multiplicity=mol.multiplicity,
            frozen_atoms=mol.frozen_atoms,
            energy=mol.energy,
        )
    
    @staticmethod
    def _combine_molecules(mol1: Molecule, mol2: Molecule) -> Molecule:
        """
        Combine two Molecule objects into one.
        
        Parameters
        ----------
        mol1 : Molecule
            First molecule (skeleton)
        mol2 : Molecule
            Second molecule (substituent)
        
        Returns
        -------
        Molecule
            Combined molecule
        
        Notes
        -----
        - symbols and positions are concatenated
        - charge and multiplicity are inherited from mol1 (skeleton)
        - frozen_atoms are merged (mol2 indices are offset by mol1 atom count)
        """
        # Combine symbols
        combined_symbols = list(mol1.chemical_symbols) + list(mol2.chemical_symbols)
        
        # Combine positions
        combined_positions = np.vstack((mol1.positions, mol2.positions))
        
        # Merge frozen_atoms (offset mol2 indices)
        combined_frozen = None
        if mol1.frozen_atoms is not None or mol2.frozen_atoms is not None:
            n1 = len(mol1)
            frozen1 = mol1.frozen_atoms if mol1.frozen_atoms is not None else [0] * n1
            frozen2 = mol2.frozen_atoms if mol2.frozen_atoms is not None else [0] * len(mol2)
            combined_frozen = list(frozen1) + list(frozen2)
        
        return Molecule(
            symbols=combined_symbols,
            positions=combined_positions,
            charge=mol1.charge,
            multiplicity=mol1.multiplicity,
            frozen_atoms=combined_frozen,
        )
    
    @staticmethod
    def _molecule_to_array(mol: Molecule) -> np.ndarray:
        """
        Convert Molecule object to numpy array.
        
        Parameters
        ----------
        mol : Molecule
            Molecule object
        
        Returns
        -------
        np.ndarray
            Array of shape (n, 4), each row is [atomic_number, x, y, z]
        """
        n_atoms = len(mol)
        arr = np.zeros((n_atoms, 4), dtype=np.float64)
        
        for i, symbol in enumerate(mol.chemical_symbols):
            arr[i, 0] = pt.to_atomic_number(symbol)
        
        arr[:, 1:4] = mol.positions
        
        return arr
        
    @staticmethod
    def _calc_relative_coords(mol: np.ndarray, base_index: int) -> np.ndarray:
        """
        Calculate relative coordinates of all atoms with respect to a base atom.
        
        Parameters
        ----------
        mol : np.ndarray
            Molecule array, shape (n, 4), each row is [atomic_number, x, y, z]
        base_index : int
            Index of the base atom (0-based)
        
        Returns
        -------
        np.ndarray
            Array of shape (n, 4), same format as input, but with coordinates
            relative to the base atom. The base atom will have coordinates (0, 0, 0).
        """
        n = mol.shape[0]
        if base_index < 0 or base_index >= n:
            raise IndexError(f"base_index {base_index} out of range [0, {n-1}]")
        
        result = mol.copy()
        base_coords = mol[base_index, 1:4]
        result[:, 1:4] = mol[:, 1:4] - base_coords
        
        return result

    @staticmethod
    def _find_optimal_position(
        skeleton_coord: np.ndarray,
        sub_coord: np.ndarray,
        skeleton_link_index: int,
        sub_link_index: int,
        buffer: float = DEFAULT_BUFFER,
        algorithm: str = 'lagrange_multipliers'
    ) -> np.ndarray:
        """
        Find the optimal position for a substituent molecule (sub) to be attached to a 
        skeleton molecule, such that the total distance sum is maximized while avoiding
        atomic overlaps.
        
        The optimization is performed by finding the optimal position for the link atom
        in sub (sub_link), and then placing all other atoms in sub according to their
        relative positions to sub_link.
        
        Parameters
        ----------
        skeleton_coord : np.ndarray
            Skeleton molecule, shape (n, 4), each row is [atomic_number, x, y, z]
        sub_coord : np.ndarray
            Substituent molecule, shape (m, 4), each row is [atomic_number, x, y, z]
        skeleton_link_index : int
            Index of the link atom in skeleton (0-based)
        sub_link_index : int
            Index of the link atom in sub (0-based)
        buffer : float
            Buffer for min distance constraint (default: 0.3 Å)
        algorithm : str
            Optimization algorithm to use (default: 'lagrange_multipliers')
            Supported: 'lagrange_multipliers'
        
        Returns
        -------
        np.ndarray
            Optimal position of the substituent molecule, shape (m, 4),
            each row is [atomic_number, x, y, z]
        
        Constraints
        -----------
        1. All atoms in sub maintain their relative positions to sub_link
        2. Distance(sub_link, skeleton_link) = bond_dist (equality constraint)
        3. Distance(sub_i, skeleton_j) >= min_dist for all i, j (inequality constraints)
        
        Objective
        ---------
        Maximize: sum of all pairwise distances between sub atoms and skeleton atoms
        """
        n_skeleton = skeleton_coord.shape[0]
        n_sub = sub_coord.shape[0]
        
        if skeleton_link_index < 0 or skeleton_link_index >= n_skeleton:
            raise IndexError(f"skeleton_link_index {skeleton_link_index} out of range [0, {n_skeleton-1}]")
        if sub_link_index < 0 or sub_link_index >= n_sub:
            raise IndexError(f"sub_link_index {sub_link_index} out of range [0, {n_sub-1}]")
        
        # Calculate relative coordinates of sub atoms with respect to sub_link
        sub_relative = IterateAnalyzer._calc_relative_coords(sub_coord, sub_link_index)
        relative_offsets = sub_relative[:, 1:4]  # shape (m, 3)
        
        # Get element indices
        sub_elements = sub_coord[:, 0].astype(int)
        skeleton_elements = skeleton_coord[:, 0].astype(int)
        
        # Get skeleton coordinates
        skeleton_coords = skeleton_coord[:, 1:4]  # shape (n, 3)
        skeleton_link_coords = skeleton_coords[skeleton_link_index]
        
        # Calculate bond distance for equality constraint (sub_link to skeleton_link)
        sub_link_element = sub_elements[sub_link_index]
        skeleton_link_element = skeleton_elements[skeleton_link_index]
        bond_dist = covalent_radii[sub_link_element] + covalent_radii[skeleton_link_element]
        
        # Calculate minimum distance matrix for inequality constraints
        sub_radii = np.array([covalent_radii[z] for z in sub_elements])  # (m,)
        skeleton_radii = np.array([covalent_radii[z] for z in skeleton_elements])  # (n,)
        min_dist_matrix = sub_radii[:, np.newaxis] + skeleton_radii[np.newaxis, :] + buffer  # (m, n)
        
        # Select algorithm
        if algorithm == 'lagrange_multipliers':
            optimal_sub = IterateAnalyzer._optimize_lagrange(
                sub_coord, skeleton_coords, skeleton_link_coords,
                relative_offsets, bond_dist, min_dist_matrix, sub_link_index
            )
        else:
            raise ValueError(f"Unknown algorithm: {algorithm}. Supported: 'lagrange_multipliers'")
        
        return optimal_sub

    @staticmethod
    def _optimize_lagrange(
        sub_coord: np.ndarray,
        skeleton_coords: np.ndarray,
        skeleton_link_coords: np.ndarray,
        relative_offsets: np.ndarray,
        bond_dist: float,
        min_dist_matrix: np.ndarray,
        sub_link_index: int
    ) -> np.ndarray:
        """
        Lagrange multiplier optimization using SLSQP.
        
        Parameters
        ----------
        sub_coord : np.ndarray
            Substituent molecule, shape (m, 4)
        skeleton_coords : np.ndarray
            Skeleton coordinates, shape (n, 3)
        skeleton_link_coords : np.ndarray
            Skeleton link atom coordinates, shape (3,)
        relative_offsets : np.ndarray
            Relative offsets of sub atoms from sub_link, shape (m, 3)
        bond_dist : float
            Bond distance for equality constraint
        min_dist_matrix : np.ndarray
            Minimum distance matrix for inequality constraints, shape (m, n)
        sub_link_index : int
            Index of the link atom in sub (0-based)
        
        Returns
        -------
        np.ndarray
            Optimal position of the substituent molecule, shape (m, 4)
        """
        def objective(x):
            """Negative sum of all pairwise distances (to maximize via minimization)"""
            sub_positions = x + relative_offsets  # (m, 3)
            diff = sub_positions[:, np.newaxis, :] - skeleton_coords[np.newaxis, :, :]  # (m, n, 3)
            distances = np.linalg.norm(diff, axis=2)  # (m, n)
            return -np.sum(distances)
        
        def objective_gradient(x):
            """Gradient of negative objective function"""
            sub_positions = x + relative_offsets  # (m, 3)
            diff = sub_positions[:, np.newaxis, :] - skeleton_coords[np.newaxis, :, :]  # (m, n, 3)
            distances = np.linalg.norm(diff, axis=2, keepdims=True)  # (m, n, 1)
            distances = np.maximum(distances, 1e-10)
            normalized_diff = diff / distances  # (m, n, 3)
            grad = -np.sum(normalized_diff, axis=(0, 1))  # (3,)
            return grad
        
        def eq_constraint(x):
            """Equality constraint: ||x - skeleton_link||^2 - bond_dist^2 = 0"""
            diff = x - skeleton_link_coords
            return np.dot(diff, diff) - bond_dist**2
        
        def eq_constraint_gradient(x):
            """Gradient of equality constraint"""
            return 2 * (x - skeleton_link_coords)
        
        def ineq_constraint(x):
            """Inequality constraints: ||sub_i - skeleton_j||^2 - min_dist[i,j]^2 >= 0"""
            sub_positions = x + relative_offsets  # (m, 3)
            diff = sub_positions[:, np.newaxis, :] - skeleton_coords[np.newaxis, :, :]  # (m, n, 3)
            dist_sq = np.sum(diff**2, axis=2)  # (m, n)
            constraints = dist_sq - min_dist_matrix**2  # (m, n)
            return constraints.flatten()  # (m * n,)
        
        def ineq_constraint_gradient(x):
            """Gradient of inequality constraints"""
            sub_positions = x + relative_offsets  # (m, 3)
            diff = sub_positions[:, np.newaxis, :] - skeleton_coords[np.newaxis, :, :]  # (m, n, 3)
            grad = 2 * diff.reshape(-1, 3)  # (m * n, 3)
            return grad
        
        # Initial guess: project onto constraint sphere
        x_init = sub_coord[sub_link_index, 1:4].copy()
        init_diff = x_init - skeleton_link_coords
        init_dist = np.linalg.norm(init_diff)
        
        if init_dist < 1e-6:
            # Default direction: opposite to the center of mass of skeleton
            skeleton_com = np.mean(skeleton_coords, axis=0)
            direction = skeleton_link_coords - skeleton_com
            dir_norm = np.linalg.norm(direction)
            if dir_norm < 1e-6:
                direction = np.array([1.0, 0.0, 0.0])
            else:
                direction = direction / dir_norm
            x_init = skeleton_link_coords + direction * bond_dist
        else:
            x_init = skeleton_link_coords + (init_diff / init_dist) * bond_dist
        
        # Define constraints for scipy
        constraints = [
            {'type': 'eq', 'fun': eq_constraint, 'jac': eq_constraint_gradient},
            {'type': 'ineq', 'fun': ineq_constraint, 'jac': ineq_constraint_gradient}
        ]
        
        # Solve using SLSQP
        result = minimize(
            objective,
            x_init,
            method='SLSQP',
            jac=objective_gradient,
            constraints=constraints,
            options={'ftol': 1e-9, 'maxiter': 1000}
        )
        
        # Extract optimal positions
        optimal_sub_link = result.x  # (3,)
        optimal_sub_positions = optimal_sub_link + relative_offsets  # (m, 3)
        
        # Build result array
        optimal_sub = np.zeros_like(sub_coord)
        optimal_sub[:, 0] = sub_coord[:, 0]  # Keep atomic numbers
        optimal_sub[:, 1:4] = optimal_sub_positions
        
        return optimal_sub
    

if __name__ == "__main__":
    test_path = '/Users/wanglewen/Desktop/Project/test/ase'
    skeleton_path = f"{test_path}/skeleton.xyz"
    sub1_path = f"{test_path}/sub1.xyz"
    skeleton_index = 5  # 1-based index
    sub_index = 1       # 1-based index
    
    # Load molecules from xyz files
    skeleton = Molecule.from_filepath(skeleton_path)
    sub1 = Molecule.from_filepath(sub1_path)
    
    print(f"Skeleton: {skeleton.num_atoms} atoms")
    print(f"Substituent: {sub1.num_atoms} atoms")
    
    # Create analyzer and run
    analyzer = IterateAnalyzer(
        skeleton=skeleton,
        substituent=sub1,
        skeleton_link_index=skeleton_index,
        substituent_link_index=sub_index,
    )
    
    combined = analyzer.run()
    
    print(f"Combined: {combined.num_atoms} atoms")
    print(f"Combined formula: {combined.chemical_formula}")
    
    # Write output
    output_path = f"{test_path}/skeleton_with_sub1_test.xyz"
    with open(output_path, 'w') as f:
        combined.write_coordinates(f)
    print(f"Output written to: {output_path}")