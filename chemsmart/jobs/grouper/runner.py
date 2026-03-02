"""
Grouper job runners for molecular structure clustering.

This module implements the GrouperJobRunner that executes grouping jobs.
All grouper algorithm implementations have been split into separate files:
- base.py: MoleculeGrouper base class and ResultsRecorder
- rmsd.py: RMSD-based groupers
- tanimoto.py: Tanimoto similarity grouper
- tfd.py: Torsion fingerprint grouper
- formula.py: Formula grouper
- connectivity.py: Connectivity grouper
- isomorphism.py: RDKit isomorphism grouper
"""

import logging
import os

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
            "output_format": job.output_format,  # Pass output format
        }

        # Strategy-specific kwargs
        if strategy in ["rmsd", "hrmsd", "spyrmsd", "irmsd", "pymolrmsd"]:
            kwargs["threshold"] = job.threshold
            kwargs["num_groups"] = job.num_groups
            kwargs["ignore_hydrogens"] = job.ignore_hydrogens
        elif strategy == "torsion":
            kwargs["threshold"] = job.threshold
            kwargs["num_groups"] = job.num_groups
            kwargs["ignore_hydrogens"] = job.ignore_hydrogens
        elif strategy in ["tanimoto", "isomorphism"]:
            kwargs["threshold"] = (
                job.threshold if strategy == "tanimoto" else None
            )
            kwargs["num_groups"] = (
                job.num_groups if strategy == "tanimoto" else None
            )
            kwargs["ignore_hydrogens"] = job.ignore_hydrogens
        elif strategy == "energy":
            kwargs["threshold"] = job.threshold
            kwargs["num_groups"] = job.num_groups
        elif strategy == "connectivity":
            kwargs["ignore_hydrogens"] = job.ignore_hydrogens
        # formula doesn't need threshold or ignore_hydrogens

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
