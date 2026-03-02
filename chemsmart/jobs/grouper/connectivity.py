"""
Connectivity-based molecular grouping.

Groups molecules by molecular connectivity (graph isomorphism).
"""

import logging
from typing import Iterable, List, Tuple

import networkx as nx
import pandas as pd
from joblib import Parallel, delayed

from chemsmart.io.molecules.structure import Molecule
from chemsmart.utils.utils import to_graph_wrapper

from .base import MoleculeGrouper

logger = logging.getLogger(__name__)


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
        adjust_H (bool): Whether to adjust hydrogen bond detection.
        ignore_hydrogens (bool): Whether to exclude hydrogen atoms from comparison.
    """

    def __init__(
        self,
        molecules: Iterable[Molecule],
        num_procs: int = 1,
        adjust_H: bool = True,
        ignore_hydrogens: bool = False,
        label: str = None,
        conformer_ids: List[str] = None,
        output_format: str = "xlsx",
        **kwargs,
    ):
        """
        Initialize connectivity-based molecular grouper.

        Args:
            molecules (Iterable[Molecule]): Collection of molecules to group.
            num_procs (int): Number of processes for parallel computation.
            adjust_H (bool): Whether to adjust hydrogen bond detection.
                Defaults to True.
            ignore_hydrogens (bool): Whether to exclude hydrogen atoms from
                graph comparison. Defaults to False.
            label (str): Label/name for output files. Defaults to None.
            conformer_ids (list[str]): Custom IDs for each molecule (e.g., ['c1', 'c2']).
            output_format (str): Output format ('xlsx', 'csv', 'txt'). Defaults to 'xlsx'.
        """
        super().__init__(
            molecules,
            num_procs,
            label=label,
            conformer_ids=conformer_ids,
            output_format=output_format,
        )
        self.adjust_H = adjust_H
        self.ignore_hydrogens = ignore_hydrogens

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
                - List of molecule groups (same connectivity)
                - List of index groups (corresponding indices for each group)
        """
        import time

        grouping_start_time = time.time()

        molecules_list = list(self.molecules)
        n = len(molecules_list)

        logger.info(
            f"[{self.__class__.__name__}] Converting {n} molecules to graphs..."
        )

        # Parallel graph conversion (use fixed bond_cutoff_buffer=0.0)
        self.graphs = Parallel(n_jobs=self.num_procs)(
            delayed(to_graph_wrapper)(mol, 0.0, self.adjust_H)
            for mol in molecules_list
        )

        # Remove hydrogens from graphs if requested
        if self.ignore_hydrogens:
            for g in self.graphs:
                h_nodes = [
                    node
                    for node, data in g.nodes(data=True)
                    if data.get("element") == "H"
                ]
                g.remove_nodes_from(h_nodes)

        logger.info(
            f"[{self.__class__.__name__}] Performing pairwise isomorphism checks..."
        )

        # Generate all pairs
        indices = [(i, j) for i in range(n) for j in range(i + 1, n)]
        total_pairs = len(indices)

        # Use Union-Find for efficient grouping
        parent = list(range(n))

        def find(x):
            if parent[x] != x:
                parent[x] = find(parent[x])
            return parent[x]

        def union(x, y):
            px, py = find(x), find(y)
            if px != py:
                parent[px] = py

        # Check isomorphism for all pairs
        for idx, (i, j) in enumerate(indices):
            if self._are_isomorphic(self.graphs[i], self.graphs[j]):
                union(i, j)
            if (idx + 1) % 100 == 0 or idx + 1 == total_pairs:
                logger.debug(f"Checked {idx + 1}/{total_pairs} pairs")

        # Build groups from Union-Find
        group_map = {}
        for i in range(n):
            root = find(i)
            if root not in group_map:
                group_map[root] = []
            group_map[root].append(i)

        # Convert to output format
        groups = []
        index_groups = []
        for indices_list in group_map.values():
            indices_list.sort()  # Sort by index (energy order)
            groups.append([molecules_list[i] for i in indices_list])
            index_groups.append(indices_list)

        # Sort groups by first index
        sorted_pairs = sorted(zip(index_groups, groups), key=lambda x: x[0][0])
        index_groups = [p[0] for p in sorted_pairs]
        groups = [p[1] for p in sorted_pairs]

        grouping_time = time.time() - grouping_start_time

        # Save results to Excel
        self._save_connectivity_results(groups, index_groups, grouping_time)

        # Cache results
        self._cached_groups = groups
        self._cached_group_indices = index_groups

        logger.info(
            f"[{self.__class__.__name__}] Found {len(groups)} connectivity groups"
        )

        return groups, index_groups

    def _save_connectivity_results(
        self,
        groups: List[List[Molecule]],
        index_groups: List[List[int]],
        grouping_time: float = None,
    ):
        """Save connectivity grouping results to file using ResultsRecorder."""
        n = sum(len(g) for g in groups)

        # Build header info
        header_info = [
            ("", f"Connectivity Grouping Results - {self.__class__.__name__}"),
            ("Total Molecules", n),
            ("adjust H", self.adjust_H),
            ("Ignore Hydrogens", self.ignore_hydrogens),
            ("Num Procs", self.num_procs),
        ]

        if grouping_time is not None:
            header_info.append(
                ("Grouping Time", f"{grouping_time:.2f} seconds")
            )

        # Use ResultsRecorder to save
        recorder = self._get_results_recorder()

        # Build summary data for main sheet
        summary_data = []
        for i, group in enumerate(groups):
            summary_data.append(
                {
                    "Group": i + 1,
                    "Members": len(group),
                }
            )
        summary_df = pd.DataFrame(summary_data)

        # Build sheets data using recorder's method
        sheets_data = {
            "Connectivity_Groups": summary_df,
            "Groups": recorder.build_groups_dataframe(index_groups, n),
        }

        # For connectivity, we write summary_df with header info
        # Need custom handling since it's not a matrix
        recorder.record_results(
            grouper_name=self.__class__.__name__,
            header_info=header_info,
            sheets_data=sheets_data,
            matrix_data=None,
            suffix=None,
            startrow=7,
        )

    def __repr__(self):
        return (
            f"{self.__class__.__name__}("
            f"num_procs={self.num_procs}, ignore_hydrogens={self.ignore_hydrogens})"
        )
