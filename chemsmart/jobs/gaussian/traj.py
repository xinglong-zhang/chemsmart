"""
Gaussian trajectory calculation job implementation.

This module provides the GaussianTrajJob class for processing
molecular trajectory data with Gaussian calculations. It handles
structure grouping, energy sorting, and selective optimization
of unique conformations from trajectory datasets.
"""

import logging
from functools import cached_property

import numpy as np

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.gaussian.job import GaussianGeneralJob, GaussianJob
from chemsmart.utils.grouper import StructureGrouperFactory

logger = logging.getLogger(__name__)


class GaussianTrajJob(GaussianJob):
    """
    Gaussian job class for trajectory structure processing.

    Processes molecular trajectory data by grouping similar structures,
    sorting by energy, and running Gaussian calculations on unique
    conformations. Supports various grouping strategies and selective
    processing of trajectory subsets.

    The workflow includes:
    1. Extract structures from trajectory end portion
    2. Group similar structures to identify unique conformations (default: no grouping)
    3. Sort structures by energy for prioritization
    4. Run Gaussian calculations on selected unique structures

    Note:
        By default, no structure grouping is performed (grouping_strategy=None),
        meaning all trajectory structures will be treated as unique and potentially
        processed. To enable structure grouping and identify truly unique conformations,
        specify a grouping strategy (e.g., 'rmsd', 'tanimoto') via the CLI -g option.

    Attributes:
        TYPE (str): Job type identifier ('g16traj').
        molecules (list[Molecule]): Processed trajectory structures from the
            last portion of the input set, determined by
            `proportion_to_opt`.
        grouper: Structure grouping utility for identifying unique
            conformations.
        num_structures_to_run (int): Number of structures to process.
        grouping_strategy (str): Method for structure grouping.
        num_procs (int): Number of processes for execution
            (currently not used).
        proportion_to_opt (float): Fraction of the trajectory tail used to
            build the `molecules` subset.
        skip_completed (bool): If True, completed jobs are not rerun.
        proportion_structures_to_use (float): Initializer input that sets
            `proportion_to_opt` and controls the trajectory fraction used.
    """

    TYPE = "g16traj"

    def __init__(
        self,
        molecules,
        settings,
        label,
        jobrunner,
        num_structures_to_run=None,
        grouping_strategy=None,
        num_procs=1,
        proportion_structures_to_use=0.1,
        skip_completed=True,
        **kwargs,
    ):
        """
        Initialize a Gaussian trajectory processing calculation.

        Sets up trajectory processing with structure grouping and
        selective optimization. Validates input and configures
        processing parameters for efficient trajectory analysis.

        Args:
            molecules (list): List of Molecule objects from trajectory.
            settings (GaussianJobSettings): Calculation configuration.
            label (str): Base label for structure jobs.
            jobrunner: Job execution handler.
            num_structures_to_run (int, optional): Number of structures
                to process. If None, process all unique structures.
            grouping_strategy (str): Structure grouping method.
            num_procs (int): Number of processes for parallel execution
                (currently not implemented).
            proportion_structures_to_use (float): Fraction of trajectory
                end to use for structure extraction (default: 0.1).
            skip_completed (bool): Skip already completed jobs.
            **kwargs: Additional keyword arguments for parent class
                and structure grouping.

        Raises:
            ValueError: If molecules is not a list, is empty, or contains
                non-Molecule objects.
        """
        if not isinstance(molecules, list):
            raise ValueError("Molecules must be a list of Molecule objects.")
        if len(molecules) == 0:
            raise ValueError("Molecules list cannot be empty.")
        if not all(isinstance(mol, Molecule) for mol in molecules):
            raise ValueError("All molecules must be instances of Molecule.")

        super().__init__(
            molecule=molecules[0],  # Use the first molecule as a placeholder
            # has to be instantiated due to parent class check
            settings=settings,
            label=label,
            jobrunner=jobrunner,
            skip_completed=skip_completed,
            **kwargs,
        )
        self.num_structures_to_run = num_structures_to_run
        self.grouping_strategy = grouping_strategy
        self.num_procs = num_procs

        # proportion of the traj from last portion to obtain structures
        self.proportion_to_opt = proportion_structures_to_use
        last_num_structures = int(
            round(len(molecules) * proportion_structures_to_use, 1)
        )
        self.molecules = molecules[-last_num_structures:]
        if grouping_strategy is not None:
            self.grouper = StructureGrouperFactory.create(
                self.molecules, strategy=self.grouping_strategy, **kwargs
            )
            self.grouper.group()
        else:
            self.grouper = None

    @cached_property
    def num_structures(self):
        """
        Get the number of processed structures from trajectory.

        Returns:
            int: Number of structures extracted from trajectory
                end portion.
        """
        return len(self.molecules)

    @cached_property
    def all_energies(self):
        """
        Get energies of all processed structures.

        Returns:
            list: Energy values for all trajectory structures.
        """
        return [structure.energy for structure in self.molecules]

    @cached_property
    def energies_indices(self):
        """
        Get indices of structures sorted by ascending energy.

        Provides index mapping for energy-based structure sorting,
        useful for identifying lowest-energy conformations.

        Returns:
            numpy.ndarray: Indices of structures in ascending
                energy order.
        """
        return np.argsort(self.all_energies)

    @cached_property
    def sorted_energies(self):
        """
        Get energies sorted in ascending order.

        Returns:
            list: Structure energies sorted from lowest to highest.
        """
        return [self.all_energies[i] for i in self.energies_indices]

    @cached_property
    def sorted_molecules(self):
        """
        Get molecules sorted by ascending energy.

        Returns:
            list: Molecule objects sorted from lowest to highest energy.
        """
        return [self.molecules[i] for i in self.energies_indices]

    @property
    def unique_structures(self):
        """
        Get unique structures after grouping similar conformations.

        Uses the configured grouping strategy to identify structurally
        distinct conformations from the trajectory data.

        Returns:
            list: Unique Molecule objects after grouping.
        """
        if self.grouper is None:
            return self.molecules
        return self.grouper.unique()

    @property
    def unique_structures_energies(self):
        """
        Get energies of unique structures.

        Returns:
            list: Energy values for each unique structure after grouping.
        """
        return [structure.energy for structure in self.unique_structures]

    @property
    def num_unique_structures(self):
        """
        Get the number of unique structures after grouping.

        Returns:
            int: Count of structurally distinct conformations.
        """
        if self.grouper is None:
            return len(self.molecules)
        return len(self.unique_structures)

    def _prepare_all_jobs(self):
        """
        Create Gaussian jobs for all unique structures.

        Generates GaussianGeneralJob objects for each unique structure
        identified by the grouping algorithm. Logs structure information
        for debugging and tracking purposes.

        Returns:
            list: GaussianGeneralJob objects for unique structures.
        """
        jobs = []
        logger.debug(
            f"Number of structures used for optimization: "
            f"{self.num_unique_structures}\n"
        )
        logger.debug(f"Unique structures: {self.unique_structures}")
        logger.debug(
            f"Unique structures energies: {self.unique_structures_energies}"
        )
        for i in range(self.num_unique_structures):
            label = f"{self.label}_c{i + 1}"  # 1-indexed for structures
            jobs += [
                GaussianGeneralJob(
                    molecule=self.unique_structures[i],
                    settings=self.settings,
                    label=label,
                    jobrunner=self.jobrunner,
                )
            ]
        return jobs

    @property
    def all_structures_run_jobs(self):
        """
        Get all prepared structure calculation jobs.

        Returns:
            list: All GaussianGeneralJob objects for trajectory
                structure processing.
        """
        return self._prepare_all_jobs()

    @property
    def last_run_job_index(self):
        """
        Get the index of the last completed job.

        Tracks progress through the structure job list for resuming
        interrupted calculations.

        Returns:
            int: Index of last finished job, or total number if
                all jobs are complete.
        """
        return self._check_last_finished_job_index()

    @property
    def incomplete_structure_run_jobs(self):
        """
        Get incomplete structure calculation jobs.

        Filters the job list to return only those that have not
        completed successfully, useful for selective resubmission.

        Returns:
            list: Incomplete GaussianGeneralJob objects.
        """
        return [
            job
            for job in self.all_structures_run_jobs
            if not job.is_complete()
        ]

    def _check_last_finished_job_index(self):
        """
        Find the index of the last completed job in sequence.

        Iterates through structure jobs to identify progress and
        determine where to resume if needed.

        Returns:
            int: Index of last finished job, or total number of unique
                structures if all jobs are complete.
        """
        for i, job in enumerate(self.all_structures_run_jobs):
            if not job.is_complete():
                return i
        # If all complete
        return self.num_unique_structures

    def _run_all_jobs(self):
        """
        Execute structure calculation jobs based on configuration.

        Runs either all available jobs or a specified subset based
        on the num_structures_to_run setting. Handles both complete
        processing and selective structure optimization.
        """
        if self.num_structures_to_run is None:
            # run all jobs if num_structures_to_run is not specified
            jobs_to_run = self.all_structures_run_jobs
        else:
            jobs_to_run = self.incomplete_structure_run_jobs[
                : self.num_structures_to_run
            ]
        for job in jobs_to_run:
            job.run()

    def _run(self):
        """
        Execute the trajectory structure processing workflow.

        Main execution method that initiates all configured structure
        calculations. Called internally by the job runner framework.
        """
        self._run_all_jobs()

    def is_complete(self):
        """
        Check if all trajectory structure jobs are complete.

        Returns:
            bool: True if all required structure calculations have
                finished successfully, False otherwise.
        """
        return self._run_all_structure_set_jobs_are_complete()

    def _run_all_structure_set_jobs_are_complete(self):
        """
        Verify completion status of all required structure jobs.

        Checks completion based on whether all jobs or a subset
        (num_structures_to_run) should be processed.

        Returns:
            bool: True if all required jobs are complete,
                False otherwise.
        """
        if self.num_structures_to_run is None:
            return all(
                job.is_complete() for job in self.all_structures_run_jobs
            )
        return all(
            job.is_complete()
            for job in self.all_structures_run_jobs[
                : self.num_structures_to_run
            ]
        )
