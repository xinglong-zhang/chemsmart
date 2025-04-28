import logging
from functools import cached_property

import numpy as np

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.gaussian.job import GaussianGeneralJob, GaussianJob
from chemsmart.utils.grouper import StructureGrouperFactory

logger = logging.getLogger(__name__)


class GaussianTrajJob(GaussianJob):
    """Gaussian job for running a set of structures.
    Args:
        molecules (list): List of Molecule objects to process.
        settings: Gaussian job settings.
        label (str): Base label for the job.
        jobrunner: JobRunner instance for executing the jobs.
        num_structures_to_run (int, optional): Number of structures to run. If None, run all.
        grouping_strategy (str): Strategy for grouping structures (default: "rmsd").
        num_procs (int): Number of processes for parallel execution (not implemented).
        proportion_structures_to_use (float): Proportion of structures to use from the end of the list.
        skip_completed (bool): Skip completed jobs if True.
        **kwargs: Additional keyword arguments for GaussianJob.
    """

    TYPE = "g16traj"

    def __init__(
        self,
        molecules,
        settings,
        label,
        jobrunner,
        num_structures_to_run=None,
        grouping_strategy="rmsd",
        num_procs=1,
        proportion_structures_to_use=0.1,
        skip_completed=True,
        **kwargs,
    ):
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
        self.grouper = StructureGrouperFactory.create(
            self.molecules, strategy=self.grouping_strategy, **kwargs
        )

    @cached_property
    def num_structures(self):
        return len(self.molecules)

    @cached_property
    def all_energies(self):
        return [structure.energy for structure in self.molecules]

    @cached_property
    def energies_indices(self):
        """List of energy indices in ascending energies order."""
        return np.argsort(self.all_energies)

    @cached_property
    def sorted_energies(self):
        """List of sorted energies in ascending order."""
        return [self.all_energies[i] for i in self.energies_indices]

    @cached_property
    def sorted_molecules(self):
        """List of molecules sorted by ascending energies."""
        return [self.molecules[i] for i in self.energies_indices]

    @property
    def unique_structures(self):
        return self.grouper.unique()

    @property
    def unique_structures_energies(self):
        return [structure.energy for structure in self.unique_structures]

    @property
    def num_unique_structures(self):
        return len(self.unique_structures)

    def _prepare_all_jobs(self):
        jobs = []
        logger.debug(
            f"Number of structures used for optimization: {self.num_unique_structures}\n"
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
        return self._prepare_all_jobs()

    @property
    def last_run_job_index(self):
        return self._check_last_finished_job_index()

    @property
    def incomplete_structure_run_jobs(self):
        return [
            job
            for job in self.all_structures_run_jobs
            if not job.is_complete()
        ]

    def _check_last_finished_job_index(self):
        for i, job in enumerate(self.all_structures_run_jobs):
            if not job.is_complete():
                return i
        # If all complete
        return self.num_unique_structures

    def _run_all_jobs(self):
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
        self._run_all_jobs()

    def is_complete(self):
        return self._run_all_structure_set_jobs_are_complete()

    def _run_all_structure_set_jobs_are_complete(self):
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
