"""
Grouper base job for molecular structure clustering.

Defines the core `GrouperJob` used by higher-level grouping strategy jobs.
It provides common file paths, output utilities, and factory methods
to create grouping jobs from molecular structures.

This follows the same pattern as mol/job.py (PyMOLJob).
"""

import logging
import os
from typing import List, Optional

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.job import Job

logger = logging.getLogger(__name__)


class GrouperJob(Job):
    """
    Base class for molecular structure grouping jobs.

    Provides core functionality for creating and managing grouping
    jobs including structure loading, grouping execution, and output
    management for molecular clustering tasks.

    Attributes:
        PROGRAM (str): Program identifier ('grouper').
        molecules (list[Molecule]): Molecules to be grouped.
        grouping_strategy (str): Grouping strategy to use.
        threshold (float): Similarity threshold for grouping.
        num_groups (int): Target number of groups (alternative to threshold).
        ignore_hydrogens (bool): Whether to ignore hydrogen atoms.
        num_procs (int): Number of parallel processes.
        label (str): Job identifier used for file naming and outputs.
    """

    TYPE = "grouper"
    PROGRAM = "grouper"

    def __init__(
        self,
        molecules: List[Molecule],
        grouping_strategy: str,
        threshold: Optional[float] = None,
        num_groups: Optional[int] = None,
        label: Optional[str] = None,
        jobrunner=None,
        ignore_hydrogens: bool = False,
        num_procs: int = 1,
        skip_completed: bool = True,
        conformer_ids: Optional[List[str]] = None,
        output_format: str = "xlsx",
        **kwargs,
    ):
        """
        Initialize a Grouper job.

        Args:
            molecules (list[Molecule]): List of Molecule objects to group.
            grouping_strategy (str): Grouping strategy (rmsd, hrmsd, irmsd, tfd, tanimoto, etc.)
            threshold (float, optional): Similarity threshold for grouping.
            num_groups (int, optional): Target number of groups (alternative to threshold).
            label (str, optional): Label for the job.
            jobrunner (JobRunner, optional): Job execution handler.
            ignore_hydrogens (bool): Whether to ignore hydrogens.
            num_procs (int): Number of processors for parallel calculation.
            skip_completed (bool): If True, skip completed jobs.
            conformer_ids (list[str], optional): Custom IDs for each molecule (e.g., ['c1', 'c2', 'c3']).
                If provided, these are used instead of numeric indices for matrix labels and output.
            output_format (str): Output file format ('xlsx', 'csv', 'txt'). Default is 'xlsx'.
            **kwargs: Additional strategy-specific arguments.
        """
        if not isinstance(molecules, list) or len(molecules) < 2:
            raise ValueError(
                "Molecules must be a list of at least 2 Molecule objects."
            )

        # Initialize parent with first molecule
        super().__init__(
            molecule=molecules[0],
            label=label or "grouper",
            jobrunner=jobrunner,
            skip_completed=skip_completed,
        )

        self.molecules = molecules
        self.grouping_strategy = grouping_strategy
        self.threshold = threshold
        self.num_groups = num_groups
        self.ignore_hydrogens = ignore_hydrogens
        self.num_procs = num_procs
        self.conformer_ids = conformer_ids
        self.output_format = output_format
        self.grouper_kwargs = kwargs

        # Results storage (populated after grouping)
        self._grouper = None
        self._groups = None
        self._group_indices = None

    @property
    def job_basename(self):
        """Get the base name for job-related files."""
        return self._get_job_basename()

    def _get_job_basename(self):
        """Internal method to derive the job base name."""
        return self.label

    @property
    def num_molecules(self) -> int:
        """Get total number of molecules."""
        return len(self.molecules)

    @property
    def output_dir(self) -> str:
        """Get output directory path."""
        if self.label:
            return os.path.join(self.folder, f"{self.label}_group_result")
        return os.path.join(self.folder, "group_result")

    @property
    def outputfile(self) -> str:
        """Get the path to the main output file (unique structures)."""
        return os.path.join(self.output_dir, f"{self.label}_unique.xyz")

    @property
    def logfile(self) -> str:
        """Get the path to the log file."""
        return os.path.join(self.folder, f"log.{self.job_basename}")

    @property
    def errfile(self) -> str:
        """Get the path to the error file."""
        return os.path.join(self.folder, f"{self.job_basename}.err")

    def _job_is_complete(self) -> bool:
        """Check if the grouper job has completed successfully.

        Checks if output directory exists.
        """
        return os.path.exists(self.output_dir)

    def is_complete(self) -> bool:
        """Check if grouping job is complete."""
        return self._job_is_complete()

    def _run(self):
        """Execute the grouper job using the configured job runner."""
        self.jobrunner.run(self)

    @classmethod
    def from_filename(
        cls,
        filename: str,
        grouping_strategy: str,
        index: str = ":",
        label: Optional[str] = None,
        jobrunner=None,
        threshold: Optional[float] = None,
        num_groups: Optional[int] = None,
        ignore_hydrogens: bool = False,
        num_procs: int = 1,
        **kwargs,
    ):
        """
        Create a Grouper job from a molecular structure file.

        Args:
            filename: Path to the molecular structure file.
            grouping_strategy: Grouping strategy to use.
            index: Molecule index selection string (default: ":" for all).
            label: Job identifier (default: derived from filename).
            jobrunner: Job execution runner (default: auto-created).
            threshold: Similarity threshold for grouping.
            num_groups: Target number of groups.
            ignore_hydrogens: Whether to ignore hydrogens.
            num_procs: Number of processors.
            **kwargs: Additional arguments for job configuration.

        Returns:
            GrouperJob: Configured grouper job instance.
        """
        from chemsmart.jobs.runner import JobRunner
        from chemsmart.utils.utils import string2index_1based

        logger.info(f"Reading molecules from file: {filename}")
        molecules = Molecule.from_filepath(
            filepath=filename, index=":", return_list=True
        )
        logger.info(f"Num of molecules read: {len(molecules)}")

        if label is None:
            label = os.path.basename(filename).split(".")[0]

        # Apply index selection
        molecules = molecules[string2index_1based(index)]
        logger.info(f"Num of molecules to use: {len(molecules)}")

        if len(molecules) < 2:
            raise ValueError(
                f"Need at least 2 molecules for grouping, got {len(molecules)}"
            )

        # Create jobrunner if not provided
        if jobrunner is None:
            jobrunner = JobRunner.from_job(
                cls(
                    molecules=molecules,
                    grouping_strategy=grouping_strategy,
                    label=label,
                    jobrunner=None,
                    threshold=threshold,
                    num_groups=num_groups,
                    ignore_hydrogens=ignore_hydrogens,
                    num_procs=num_procs,
                    **kwargs,
                ),
                server=kwargs.get("server"),
                scratch=kwargs.get("scratch"),
                fake=kwargs.get("fake", False),
            )

        return cls(
            molecules=molecules,
            grouping_strategy=grouping_strategy,
            label=label,
            jobrunner=jobrunner,
            threshold=threshold,
            num_groups=num_groups,
            ignore_hydrogens=ignore_hydrogens,
            num_procs=num_procs,
            **kwargs,
        )
