import logging

import networkx as nx
from rdkit import Chem

from chemsmart.io.molecules.structure import Molecule

logger = logging.getLogger(__name__)


class BasePreprocessor:
    """Base preprocessor class with shared logic
    for Skeleton and Substituent preprocessors."""

    def __init__(self, molecule: Molecule, link_index: int):
        self.molecule = molecule
        # Convert 1-based to 0-based index
        self.link_index = link_index - 1
        # Store the final indices used to create the processed molecule
        self._final_indices = None

    def _has_available_bonding_position(self) -> bool:
        """
        Check if the link atom has an available bonding position.
        """
        # Build molecular graph
        graph = self.molecule.to_graph()

        # Calculate current total bond order
        current_bond_order_sum = 0.0
        for neighbor in graph.neighbors(self.link_index):
            edge_data = graph.get_edge_data(self.link_index, neighbor)
            # Default to 1.0 if bond_order is missing
            # (shouldn't happen with to_graph)
            bond_order = edge_data.get("bond_order", 1.0)
            current_bond_order_sum += bond_order

        # Get element symbol
        element = self.molecule.chemical_symbols[self.link_index]

        # Get expected maximum bonding capacity
        max_bonds = self._get_max_bonding_capacity(element)

        # Use a small epsilon for float comparison
        return current_bond_order_sum < (max_bonds - 0.1)

    @staticmethod
    def _get_max_bonding_capacity(element: str) -> int:
        """
        Get the maximum bonding capacity for an element.
        """
        periodic_table = Chem.GetPeriodicTable()
        atomic_num = periodic_table.GetAtomicNumber(element)

        # GetDefaultValence returns a tuple of possible valences, take the max
        default_valence = periodic_table.GetDefaultValence(atomic_num)

        if isinstance(default_valence, tuple):
            return max(default_valence)
        return default_valence

    def detect_substituent(self) -> list[int]:
        """
        Auto-detect substituent/group connected to link atom.
        Identifies the smallest connected component when the bond is broken.
        """
        # Build molecular graph
        graph = self.molecule.to_graph()

        # Get neighbors of link atom
        neighbors = list(graph.neighbors(self.link_index))

        if len(neighbors) == 0:
            logger.warning(
                f"Link atom at index {self.link_index} has no neighbors"
            )
            return []

        # For each neighbor, calculate the size
        # of the component if bond is broken.
        # The substituent/group to remove is the smallest component.
        min_component_size = float("inf")
        target_indices = []

        for neighbor in neighbors:
            # Temporarily remove the edge
            graph_copy = graph.copy()
            graph_copy.remove_edge(self.link_index, neighbor)

            # Find connected components
            components = list(nx.connected_components(graph_copy))

            # Find which component contains the neighbor (not the link atom)
            for component in components:
                if neighbor in component and self.link_index not in component:
                    if len(component) < min_component_size:
                        min_component_size = len(component)
                        target_indices = list(component)
                    break

        return target_indices

    def _get_complement_indices(self, exclude_indices: list[int]) -> list[int]:
        """
        Get indices of all atoms except those in exclude_indices.
        """
        exclude_set = set(exclude_indices)
        return [i for i in range(len(self.molecule)) if i not in exclude_set]

    def _extract_by_indices(self, indices: list[int]) -> Molecule:
        """
        Extract a subset of atoms from molecule by indices.
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
            frozen_atoms = [
                self.molecule.frozen_atoms[i] for i in sorted_indices
            ]

        return Molecule(
            symbols=symbols,
            positions=positions,
            charge=self.molecule.charge,
            multiplicity=self.molecule.multiplicity,
            frozen_atoms=frozen_atoms,
        )

    def _run_auto_detect(self) -> Molecule:
        """
        Run auto-detection mode to find and remove
        the smallest group at link position.

        Returns
        -------
        Molecule
            Processed molecule with group removed.
        """
        removed_indices = self.detect_substituent()
        keep_indices = self._get_complement_indices(removed_indices)
        self._final_indices = sorted(keep_indices)
        return self._extract_by_indices(keep_indices)

    def get_new_link_index(self) -> int:
        """
        Get the new link index in the processed molecule (1-based).

        After removing atoms, the original link_index may change.
        This method returns the new index.

        Returns
        -------
        int
            New link index (1-based)
        """
        if self._final_indices is None:
            indices = self._get_fallback_indices()
        else:
            indices = self._final_indices

        # Find position of link_index in sorted indices
        try:
            new_index_0based = indices.index(self.link_index)
            return new_index_0based + 1  # Convert to 1-based
        except ValueError:
            raise ValueError(
                f"Link atom at index {self.link_index + 1} (1-based) was removed during preprocessing"
            )

    def _get_fallback_indices(self) -> list[int]:
        """
        Get indices to use for new link index
        calculation if run() wasn't called.
        Default implementation: Auto-detect.
        """
        substituent_indices = self.detect_substituent()
        return sorted(self._get_complement_indices(substituent_indices))


class SkeletonPreprocessor(BasePreprocessor):
    """Preprocessor to prepare skeleton molecule
    by removing substituent at link position.

    This class handles two scenarios:
    1. User provides skeleton_indices: Keep only atoms at specified indices
    2. User doesn't provide skeleton_indices:
    Auto-detect and remove substituent at link position
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
            Index of the link atom (1-based, will
            be converted to 0-based internally)
        skeleton_indices : list[int] | None
            Indices of atoms belonging to skeleton (1-based).
            If None, auto-detection will be used.
        """
        super().__init__(molecule, link_index)
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
            logger.debug(
                f"Link atom at index {self.link_index + 1} (1-based) has available bonding position. "
                "No substituent removal needed."
            )
            # No removal, so final indices are all atoms
            self._final_indices = list(range(len(self.molecule)))
            return self.molecule

        if self.skeleton_indices is not None:
            # Mode 1: User provided skeleton indices
            # Check if link_index is in skeleton_indices
            if self.link_index not in self.skeleton_indices:
                raise ValueError(
                    f"Link atom at index {self.link_index + 1} (1-based) must be included in skeleton_indices."
                )

            # Find substituent branches that don't contain skeleton atoms
            substituent_branches = self._find_non_skeleton_branches()

            if len(substituent_branches) == 0:
                # No substituent to remove at this
                # link_index, return molecule as-is
                self._final_indices = list(range(len(self.molecule)))
                return self.molecule
            elif len(substituent_branches) == 1:
                # Exactly one substituent branch - remove only this branch
                # Keep all other atoms (including
                # substituents on other link atoms)
                substituent_indices = substituent_branches[0]
                keep_indices = self._get_complement_indices(
                    substituent_indices
                )
                self._final_indices = sorted(keep_indices)
                return self._extract_by_indices(keep_indices)
            else:
                # Multiple substituent branches -
                # ambiguous, fall back to auto-detection
                logger.warning(
                    f"Found {len(substituent_branches)} branches without skeleton atoms. "
                    "Falling back to auto-detection mode."
                )
                return self._run_auto_detect()
        else:
            # Mode 2: Auto-detect substituent
            return self._run_auto_detect()

    def _find_non_skeleton_branches(self) -> list[list[int]]:
        """
        Find all branches from link_index
        that don't contain any skeleton atoms.

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
            logger.warning(
                f"Link atom at index {self.link_index} has no neighbors"
            )
            return []

        skeleton_set = set(self.skeleton_indices)
        non_skeleton_branches = []

        for neighbor in neighbors:
            # DFS to collect all atoms in this branch (excluding link_index)
            branch_atoms = self._dfs_collect_branch(
                graph, neighbor, self.link_index
            )

            # Check if any atom in this branch is in skeleton_indices
            has_skeleton_atom = any(
                atom in skeleton_set for atom in branch_atoms
            )

            if not has_skeleton_atom:
                # This branch has no skeleton atoms - it's a substituent
                non_skeleton_branches.append(branch_atoms)

        return non_skeleton_branches

    def _dfs_collect_branch(
        self, graph: nx.Graph, start: int, excluded: int
    ) -> list[int]:
        """
        Collect all atoms in a branch using DFS,
        excluding the starting point's parent.

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

    def _get_fallback_indices(self) -> list[int]:
        if self.skeleton_indices is not None:
            return sorted(self.skeleton_indices)
        return super()._get_fallback_indices()


class SubstituentPreprocessor(BasePreprocessor):
    """Preprocessor to prepare substituent molecule
    by removing atom/group at link position if needed.

    This class checks if the link atom has an available bonding position.
    If not, it auto-detects and removes the
    smallest substituent group at the link position.

    Unlike SkeletonPreprocessor, SubstituentPreprocessor
    does not have a "skeleton_indices" concept.
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
            Index of the link atom (1-based, will
            be converted to 0-based internally)
        """
        super().__init__(molecule, link_index)

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
            logger.debug(
                f"Substituent link atom at index {self.link_index + 1} (1-based) has available bonding position. "
                "No removal needed."
            )
            self._final_indices = list(range(len(self.molecule)))
            return self.molecule

        # Auto-detect and remove the smallest group at link position
        logger.debug(
            f"Substituent link atom at index {self.link_index + 1} (1-based) has no available bonding position. "
            "Auto-detecting and removing group."
        )
        return self._run_auto_detect()

    def _get_fallback_indices(self) -> list[int]:
        if self._has_available_bonding_position():
            return list(range(len(self.molecule)))
        return super()._get_fallback_indices()
