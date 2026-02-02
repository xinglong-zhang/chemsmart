"""
Grouper job runners for molecular structure clustering.

This module implements the GrouperJobRunner that executes grouping jobs.
All grouper algorithm implementations have been split into separate files:
- rmsd.py: RMSD-based groupers
- tanimoto.py: Tanimoto similarity grouper
- tfd.py: Torsion fingerprint grouper
- formula.py: Formula grouper
- connectivity.py: Connectivity grouper
- isomorphism.py: RDKit isomorphism grouper
"""

import logging
import os
from abc import ABC, abstractmethod
from typing import Iterable, List, Tuple

import networkx as nx

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.runner import JobRunner

logger = logging.getLogger(__name__)


class GrouperJobRunner(JobRunner):
    """Job runner for molecular grouping/clustering jobs."""

    JOBTYPES = ["grouper"]
    PROGRAM = "grouper"
    FAKE = False
    SCRATCH = False

    def __init__(
        self, server, scratch=None, fake=False, scratch_dir=None, **kwargs
    ):
        if scratch is None:
            scratch = self.SCRATCH
        super().__init__(
            server=server,
            scratch=scratch,
            scratch_dir=scratch_dir,
            fake=fake,
            **kwargs,
        )

    def _prerun(self, job):
        if not os.path.exists(job.output_dir):
            os.makedirs(job.output_dir)

    def _get_command(self, job):
        return None

    def _get_executable(self):
        return None

    def _create_process(self, job, command, env):
        try:
            grouper = self._create_grouper(job)
            groups, group_indices = grouper.group()
            job._grouper = grouper
            job._groups = groups
            job._group_indices = group_indices
            self._write_outputs(job, groups, group_indices)
            return 0
        except Exception as e:
            logger.error(f"Error during grouping: {str(e)}")
            with open(job.errfile, "w") as err:
                err.write(f"Error during grouping: {str(e)}\n")
            return 1

    def _get_grouper_classes(self):
        """Lazy load grouper classes to avoid circular imports."""
        from .connectivity import ConnectivityGrouper
        from .energy import EnergyGrouper
        from .formula import FormulaGrouper
        from .isomorphism import RDKitIsomorphismGrouper
        from .rmsd import (
            BasicRMSDGrouper,
            HungarianRMSDGrouper,
            IRMSDGrouper,
            PymolRMSDGrouper,
            SpyRMSDGrouper,
        )
        from .tanimoto import TanimotoSimilarityGrouper
        from .tfd import TorsionFingerprintGrouper

        return {
            "rmsd": BasicRMSDGrouper,
            "hrmsd": HungarianRMSDGrouper,
            "spyrmsd": SpyRMSDGrouper,
            "irmsd": IRMSDGrouper,
            "pymolrmsd": PymolRMSDGrouper,
            "tanimoto": TanimotoSimilarityGrouper,
            "torsion": TorsionFingerprintGrouper,
            "isomorphism": RDKitIsomorphismGrouper,
            "formula": FormulaGrouper,
            "connectivity": ConnectivityGrouper,
            "energy": EnergyGrouper,
        }

    def _create_grouper(self, job):
        """Create appropriate grouper instance based on job strategy."""
        strategy = job.grouping_strategy
        grouper_classes = self._get_grouper_classes()

        if strategy not in grouper_classes:
            raise ValueError(f"Unknown grouping strategy: {strategy}")

        grouper_cls = grouper_classes[strategy]

        # Common kwargs
        kwargs = {
            "molecules": job.molecules,
            "num_procs": job.num_procs,
            "label": job.label,
            "conformer_ids": job.conformer_ids,  # Pass custom conformer IDs
        }

        # Strategy-specific kwargs
        if strategy in ["rmsd", "hrmsd", "spyrmsd", "irmsd", "pymolrmsd"]:
            kwargs["threshold"] = job.threshold
            kwargs["num_groups"] = job.num_groups
            kwargs["ignore_hydrogens"] = job.ignore_hydrogens
        elif strategy in ["tanimoto", "torsion", "energy"]:
            kwargs["threshold"] = job.threshold
            kwargs["num_groups"] = job.num_groups
        elif strategy == "connectivity":
            kwargs["threshold"] = job.threshold
        # isomorphism and formula don't need threshold

        # Add any additional kwargs from job
        kwargs.update(job.grouper_kwargs)

        return grouper_cls(**kwargs)

    def _write_outputs(self, job, groups, group_indices):
        """
        Write grouping results following the same logic as utils/grouper.py unique() method.

        - Creates per-group xyz files with molecules sorted by energy
        - Appends Groups summary to existing RMSD matrix excel file
        """
        unique_molecules = []
        conformer_ids = (
            job.conformer_ids
        )  # Get custom conformer IDs if provided

        # Determine file prefix
        file_prefix = f"{job.label}_group" if job.label else "group"

        for i, (group, indices) in enumerate(zip(groups, group_indices)):
            # Create tuples of (molecule, original_index) for tracking
            mol_index_pairs = list(zip(group, indices))

            # Filter molecules that have energy information and sort by energy
            molecules_with_energy = [
                (mol, idx)
                for mol, idx in mol_index_pairs
                if mol.energy is not None
            ]
            molecules_without_energy = [
                (mol, idx)
                for mol, idx in mol_index_pairs
                if mol.energy is None
            ]

            # Sort molecules with energy by energy (ascending - lowest first)
            if molecules_with_energy:
                sorted_pairs = sorted(
                    molecules_with_energy, key=lambda pair: pair[0].energy
                )
                # Add molecules without energy at the end
                sorted_pairs.extend(molecules_without_energy)
            else:
                # If no molecules have energy, use original group order
                sorted_pairs = mol_index_pairs

            # Write group XYZ file with all molecules sorted by energy
            group_filename = os.path.join(
                job.output_dir, f"{file_prefix}_{i+1}.xyz"
            )
            with open(group_filename, "w") as f:
                for j, (mol, original_idx) in enumerate(sorted_pairs):
                    # Write the molecule coordinates
                    f.write(f"{mol.num_atoms}\n")

                    # Determine original index label (use conformer_id if available)
                    if conformer_ids is not None:
                        original_label = conformer_ids[original_idx]
                    else:
                        original_label = str(original_idx + 1)

                    # Create comment line with energy info and original molecule index
                    if mol.energy is not None:
                        comment = f"Group {i+1} Molecule {j+1} Original_Index: {original_label} Energy(Hartree): {mol.energy:.8f}"
                    else:
                        comment = f"Group {i+1} Molecule {j+1} Original_Index: {original_label} Energy: N/A"

                    f.write(f"{comment}\n")

                    # Write coordinates
                    for symbol, position in zip(
                        mol.chemical_symbols, mol.positions
                    ):
                        f.write(
                            f"{symbol:2s} {position[0]:15.10f} {position[1]:15.10f} {position[2]:15.10f}\n"
                        )

            logger.info(
                f"Written group {i+1} with {len(sorted_pairs)} molecules to {group_filename}"
            )

            # Add the lowest energy molecule (first in sorted pairs) as representative
            unique_molecules.append(sorted_pairs[0][0])

        logger.info(
            f"Generated {len(groups)} group XYZ files in {job.output_dir}"
        )

    def _run(self, process, **kwargs):
        pass

    def _postrun(self, job, **kwargs):
        pass


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
        label (str): Label/name for this grouping task (used in output filenames).
        conformer_ids (list[str]): Optional custom IDs for each molecule (e.g., ['c1', 'c2']).
    """

    def __init__(
        self,
        molecules: Iterable[Molecule],
        num_procs: int = 1,
        label: str = None,
        conformer_ids: List[str] = None,
    ):
        """
        Initialize the molecular grouper.

        Args:
            molecules (Iterable[Molecule]): Collection of molecules to group.
            num_procs (int): Number of processes for parallel computation.
                Defaults to 1.
            label (str): Label/name for this grouping task. Used in output folder
                and file names. Defaults to None.
            conformer_ids (list[str]): Optional custom IDs for each molecule (e.g., ['c1', 'c2']).
                If provided, these are used as labels in matrix output instead of numeric indices.
        """
        self.molecules = molecules
        self.num_procs = int(max(1, num_procs))
        self.label = label
        self.conformer_ids = conformer_ids

        # Cache for avoiding repeated grouping calculations
        self._cached_groups = None
        self._cached_group_indices = None

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

    def unique(
        self, output_dir: str = ".", prefix: str = "group"
    ) -> List[Molecule]:
        """
        Get unique representative molecules from each group.

        Returns the lowest energy molecule from each group as a representative
        of that structural family. Also generates XYZ files for each group,
        sorted by energy, in a dedicated subfolder.

        Args:
            output_dir (str): Base directory for output. Default is current directory.
            prefix (str): Prefix for output XYZ files. Default is "group".

        Returns:
            List[Molecule]: List of unique representative molecules (lowest energy from each group).
        """
        import os

        # Use cached results if available, otherwise compute and cache
        if (
            self._cached_groups is not None
            and self._cached_group_indices is not None
        ):
            print(f"[{self.__class__.__name__}] Using cached grouping results")
            groups, group_indices = (
                self._cached_groups,
                self._cached_group_indices,
            )
        else:
            print(
                f"[{self.__class__.__name__}] Computing groups for unique method"
            )
            groups, group_indices = self.group()
            # Cache the results
            self._cached_groups = groups
            self._cached_group_indices = group_indices

        unique_molecules = []

        # Create dedicated subfolder for XYZ files (include label if provided)
        if self.label:
            result_folder = f"{self.label}_group_result"
        else:
            result_folder = "group_result"
        full_output_path = os.path.join(output_dir, result_folder)
        os.makedirs(full_output_path, exist_ok=True)

        logger.info(f"Creating XYZ files in folder: {full_output_path}")

        # Determine the file prefix (include label if provided)
        if self.label:
            file_prefix = f"{self.label}_{prefix}"
        else:
            file_prefix = prefix

        for i, (group, indices) in enumerate(zip(groups, group_indices)):
            # Create tuples of (molecule, original_index) for tracking
            mol_index_pairs = list(zip(group, indices))

            # Filter molecules that have energy information and sort by energy
            molecules_with_energy = [
                (mol, idx)
                for mol, idx in mol_index_pairs
                if mol.energy is not None
            ]
            molecules_without_energy = [
                (mol, idx)
                for mol, idx in mol_index_pairs
                if mol.energy is None
            ]

            # Sort molecules with energy by energy (ascending - lowest first)
            if molecules_with_energy:
                sorted_pairs = sorted(
                    molecules_with_energy, key=lambda pair: pair[0].energy
                )
                # Add molecules without energy at the end
                sorted_pairs.extend(molecules_without_energy)
            else:
                # If no molecules have energy, use original group order
                sorted_pairs = mol_index_pairs

            # Write group XYZ file with all molecules sorted by energy
            group_filename = os.path.join(
                full_output_path, f"{file_prefix}_{i+1}.xyz"
            )
            with open(group_filename, "w") as f:
                for j, (mol, original_idx) in enumerate(sorted_pairs):
                    # Write the molecule coordinates
                    f.write(f"{mol.num_atoms}\n")

                    # Create comment line with energy info and original molecule index
                    if mol.energy is not None:
                        comment = f"Group {i+1} Member {j+1} Original_Index: {original_idx+1} Energy(Hartree): {mol.energy:.8f}"
                    else:
                        comment = f"Group {i+1} Member {j+1} Original_Index: {original_idx+1} Energy: N/A"

                    f.write(f"{comment}\n")

                    # Write coordinates
                    for symbol, position in zip(
                        mol.chemical_symbols, mol.positions
                    ):
                        f.write(
                            f"{symbol:2s} {position[0]:15.10f} {position[1]:15.10f} {position[2]:15.10f}\n"
                        )

            logger.info(
                f"Written group {i+1} with {len(sorted_pairs)} molecules to {group_filename}"
            )

            # Add the lowest energy molecule (first in sorted pairs) as representative
            unique_molecules.append(sorted_pairs[0][0])

        logger.info(
            f"Generated {len(groups)} group XYZ files in {full_output_path}"
        )

        return unique_molecules
