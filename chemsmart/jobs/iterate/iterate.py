import logging

import networkx as nx
from rdkit import Chem

from chemsmart.io.molecules.structure import Molecule

logger = logging.getLogger(__name__)


class BasePreprocessor:
    """Shared single- and multi-site molecular preprocessing framework.

    A single attachment site is represented by a one-element ``link_indices``
    list. For any number of sites, removal decisions are made against the
    original molecule and applied together so atom indices are remapped only
    once. Subclasses customize only how a removable branch is selected.

    Parameters
    ----------
    molecule : Molecule
        Original molecule to preprocess.
    link_indices : list[int]
        One or more 1-based attachment atom indices.
    """

    def __init__(self, molecule: Molecule, link_indices: list[int]):
        self.molecule = molecule
        self.link_indices = [index - 1 for index in link_indices]
        self._validate_link_indices()

    def _validate_link_indices(self) -> None:
        """Validate the original 0-based attachment indices."""
        if not self.link_indices:
            raise ValueError("At least one link index is required.")
        if len(set(self.link_indices)) != len(self.link_indices):
            raise ValueError("Link indices must be unique.")

        atom_count = len(self.molecule)
        invalid = [
            index + 1
            for index in self.link_indices
            if index < 0 or index >= atom_count
        ]
        if invalid:
            raise ValueError(
                f"Link indices {invalid} are out of bounds for a molecule "
                f"with {atom_count} atoms."
            )

    def run(self) -> tuple[Molecule, dict[int, int]]:
        """Remove one branch per saturated attachment site in one pass.

        Returns
        -------
        tuple[Molecule, dict[int, int]]
            Processed molecule and a mapping from original 1-based atom
            indices to processed 1-based atom indices.
        """
        graph = self.molecule.to_graph()
        link_set = set(self.link_indices)
        remove_indices: set[int] = set()

        for link_index in self.link_indices:
            if self._has_available_bonding_position(graph, link_index):
                continue

            branch = self._select_removal_branch(graph, link_index)
            other_links = (link_set - {link_index}).intersection(branch)
            if other_links:
                formatted = sorted(index + 1 for index in other_links)
                raise ValueError(
                    f"Removing a branch at link index {link_index + 1} "
                    f"would remove other attachment link indices {formatted}."
                )

            overlap = remove_indices.intersection(branch)
            if overlap:
                formatted = sorted(index + 1 for index in overlap)
                raise ValueError(
                    "Removal branches for different attachment sites "
                    f"overlap at atom indices {formatted}."
                )
            remove_indices.update(branch)

        keep_indices = sorted(set(range(len(self.molecule))) - remove_indices)
        if not keep_indices:
            raise ValueError(
                "All atoms would be removed during preprocessing."
            )

        processed = self._extract_molecule(keep_indices)
        index_map = {
            original + 1: new + 1 for new, original in enumerate(keep_indices)
        }
        return processed, index_map

    def _select_removal_branch(
        self, graph: nx.Graph, link_index: int
    ) -> list[int]:
        """Select the smallest detachable branch at a saturated link atom."""
        branch = self._detect_smallest_branch(graph, link_index)
        if not branch:
            raise ValueError(
                f"Link atom at index {link_index + 1} is saturated, but no "
                "removable branch was found."
            )
        return branch

    def _has_available_bonding_position(
        self, graph: nx.Graph, link_index: int
    ) -> bool:
        """Return whether a link atom can accept another single bond."""
        bond_order_sum = sum(
            graph.get_edge_data(link_index, neighbor).get("bond_order", 1.0)
            for neighbor in graph.neighbors(link_index)
        )
        element = self.molecule.chemical_symbols[link_index]
        return bond_order_sum < (self._get_max_bonding_capacity(element) - 0.1)

    def _detect_smallest_branch(
        self, graph: nx.Graph, link_index: int
    ) -> list[int]:
        """Find the smallest branch disconnected by cutting one link bond."""
        neighbors = sorted(graph.neighbors(link_index))
        if not neighbors:
            logger.warning(
                "Link atom at index %s (1-based) has no neighbors.",
                link_index + 1,
            )
            return []

        branches: list[list[int]] = []
        for neighbor in neighbors:
            graph_copy = graph.copy()
            graph_copy.remove_edge(link_index, neighbor)
            component = nx.node_connected_component(graph_copy, neighbor)
            if link_index not in component:
                branches.append(sorted(component))

        if not branches:
            return []
        return min(
            branches,
            key=lambda branch: self._branch_sort_key(
                graph, link_index, branch
            ),
        )

    def _branch_sort_key(self, graph, link_index, branch_atom_indices):
        """Rank by size, actual bonded start's atomic number, then index."""
        branch_atom_count = len(branch_atom_indices)
        branch_start_index = next(
            neighbor_index
            for neighbor_index in graph.neighbors(link_index)
            if neighbor_index in branch_atom_indices
        )
        branch_start_atomic_number = Chem.GetPeriodicTable().GetAtomicNumber(
            self.molecule.chemical_symbols[branch_start_index]
        )
        return (
            branch_atom_count,
            branch_start_atomic_number,
            branch_start_index,
        )

    def _extract_molecule(self, indices: list[int]) -> Molecule:
        """Return a molecule containing selected 0-based atom indices."""
        if not indices:
            raise ValueError("Cannot create molecule with no atoms.")

        indices = sorted(indices)
        frozen_atoms = None
        if self.molecule.frozen_atoms is not None:
            frozen_atoms = [self.molecule.frozen_atoms[i] for i in indices]

        return Molecule(
            symbols=[self.molecule.chemical_symbols[i] for i in indices],
            positions=self.molecule.positions[indices],
            charge=self.molecule.charge,
            multiplicity=self.molecule.multiplicity,
            frozen_atoms=frozen_atoms,
        )

    @staticmethod
    def _get_max_bonding_capacity(element: str) -> int:
        """Return RDKit's default maximum valence for an element."""
        periodic_table = Chem.GetPeriodicTable()
        atomic_number = periodic_table.GetAtomicNumber(element)
        default_valence = periodic_table.GetDefaultValence(atomic_number)
        if isinstance(default_valence, tuple):
            return max(default_valence)
        return default_valence


class SkeletonPreprocessor(BasePreprocessor):
    """Prepare a skeleton for one or more simultaneous attachments.

    When ``skeleton_indices`` are supplied, a removable branch may not contain
    a skeleton-core atom. If several non-skeleton branches are available at a
    saturated site, exactly one is selected deterministically.

    Parameters
    ----------
    molecule : Molecule
        Original skeleton molecule.
    link_indices : list[int]
        One or more 1-based attachment atom indices.
    skeleton_indices : list[int] or None
        1-based indices defining the skeleton core.
    """

    def __init__(
        self,
        molecule: Molecule,
        link_indices: list[int],
        skeleton_indices: list[int] | None = None,
    ):
        super().__init__(molecule, link_indices)
        self.skeleton_indices = (
            None
            if skeleton_indices is None
            else [index - 1 for index in skeleton_indices]
        )
        self._validate_skeleton_indices()

    def _validate_skeleton_indices(self) -> None:
        """Validate skeleton-core indices and their attachment membership."""
        if self.skeleton_indices is None:
            return

        atom_count = len(self.molecule)
        invalid = [
            index + 1
            for index in self.skeleton_indices
            if index < 0 or index >= atom_count
        ]
        if invalid:
            raise ValueError(
                f"Skeleton indices {invalid} are out of bounds for a molecule "
                f"with {atom_count} atoms."
            )

        skeleton_set = set(self.skeleton_indices)
        missing_links = [
            index + 1
            for index in self.link_indices
            if index not in skeleton_set
        ]
        if missing_links:
            raise ValueError(
                f"Skeleton link indices {missing_links} must be included in "
                "skeleton_indices."
            )

    def _select_removal_branch(
        self, graph: nx.Graph, link_index: int
    ) -> list[int]:
        """Select one branch that contains no skeleton-core atom."""
        if self.skeleton_indices is None:
            return super()._select_removal_branch(graph, link_index)

        candidates = self._find_non_skeleton_branches(graph, link_index)
        if not candidates:
            raise ValueError(
                f"Skeleton link atom at index {link_index + 1} is saturated, "
                "but no removable non-skeleton branch was found."
            )

        selected = min(
            candidates,
            key=lambda branch: self._branch_sort_key(
                graph, link_index, branch
            ),
        )
        if len(candidates) > 1:
            logger.debug(
                "Found %s removable branches at skeleton link index %s; "
                "selected atom indices %s.",
                len(candidates),
                link_index + 1,
                [index + 1 for index in selected],
            )
        return selected

    def _find_non_skeleton_branches(
        self, graph: nx.Graph, link_index: int
    ) -> list[list[int]]:
        """Return branches containing no skeleton-core atom."""
        skeleton_set = set(self.skeleton_indices or [])
        branches: list[list[int]] = []
        for neighbor in sorted(graph.neighbors(link_index)):
            graph_copy = graph.copy()
            graph_copy.remove_edge(link_index, neighbor)
            branch = sorted(nx.node_connected_component(graph_copy, neighbor))
            if link_index not in branch and not skeleton_set.intersection(
                branch
            ):
                branches.append(branch)
        return branches

    @staticmethod
    def _collect_branch(
        graph: nx.Graph, start: int, excluded: int
    ) -> list[int]:
        """Collect a branch without traversing through its attachment atom."""
        visited: set[int] = set()
        stack = [start]
        while stack:
            atom_index = stack.pop()
            if atom_index in visited or atom_index == excluded:
                continue
            visited.add(atom_index)
            stack.extend(
                neighbor
                for neighbor in graph.neighbors(atom_index)
                if neighbor not in visited and neighbor != excluded
            )
        return sorted(visited)


class SubstituentPreprocessor(BasePreprocessor):
    """Prepare one substituent using the shared one-site preprocessing path."""

    def __init__(
        self,
        molecule: Molecule,
        link_index: int,
        remove_branch_start: int | None = None,
        skip_cleanup: bool = False,
        label: str = "unnamed",
    ):
        self.label = label
        self.remove_branch_start = remove_branch_start
        self.skip_cleanup = skip_cleanup
        super().__init__(molecule, [link_index])
        context = f"Substituent '{label}' (link atom {link_index})"
        if skip_cleanup and remove_branch_start is not None:
            raise ValueError(
                f"{context}: skip_cleanup conflicts with remove_branch_start "
                f"{remove_branch_start}."
            )
        if remove_branch_start is not None and (
            type(remove_branch_start) is not int
            or not 1 <= remove_branch_start <= len(molecule)
            or remove_branch_start == link_index
        ):
            raise ValueError(
                f"{context}: invalid remove_branch_start {remove_branch_start}; "
                "must be an original 1-based atom index different from the link atom."
            )

    def run(self) -> tuple[Molecule, dict[int, int]]:
        if self.skip_cleanup:
            return self.molecule.copy(), {
                i: i for i in range(1, len(self.molecule) + 1)
            }
        try:
            return super().run()
        except ValueError as exc:
            raise ValueError(f"Substituent '{self.label}': {exc}") from exc

    def _has_available_bonding_position(self, graph, link_index):
        if self.remove_branch_start is not None:
            return False
        return super()._has_available_bonding_position(graph, link_index)

    def _select_removal_branch(self, graph, link_index):
        if self.remove_branch_start is None:
            return super()._select_removal_branch(graph, link_index)
        start = self.remove_branch_start - 1
        context = (
            f"link atom {link_index + 1}, remove_branch_start {start + 1}"
        )
        if not graph.has_edge(link_index, start):
            raise ValueError(f"{context}: atoms are not bonded.")
        graph_copy = graph.copy()
        graph_copy.remove_edge(link_index, start)
        branch = nx.node_connected_component(graph_copy, start)
        if link_index in branch:
            raise ValueError(
                f"{context}: cutting this bond cannot separate the branch "
                "because an alternative path connects it to the link atom."
            )
        return sorted(branch)
