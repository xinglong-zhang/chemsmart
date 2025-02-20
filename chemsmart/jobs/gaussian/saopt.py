import logging
from functools import cached_property

import numpy as np

from chemsmart.jobs.gaussian.job import GaussianGeneralJob, GaussianJob
from chemsmart.utils.grouper import StructureGrouperFactory

logger = logging.getLogger(__name__)


class GaussianSAOptJob(GaussianJob):
    TYPE = "g16saopt"

    def __init__(
        self,
        molecules,
        settings,
        label,
        num_structures_to_run=None,
        grouping_strategy="rmsd",
        num_procs=1,
        proportion_structures_to_use=0.1,
        skip_completed=True,
        **kwargs,
    ):
        if not isinstance(molecules, list) and len(molecules) == 0:
            raise ValueError("Molecules must be a list of Molecule objects.")

        super().__init__(
            molecule=molecules[0],
            settings=settings,
            label=label,
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
                )
            ]
        return jobs

    @property
    def all_structures_opt_jobs(self):
        return self._prepare_all_jobs()

    @property
    def last_run_job_index(self):
        return self._check_last_finished_job_index()

    @property
    def incomplete_structure_opt_jobs(self):
        return [
            job
            for job in self.all_structures_opt_jobs
            if not job.is_complete()
        ]

    def _check_last_finished_job_index(self):
        for i, job in enumerate(self.all_structures_opt_jobs):
            if not job.is_complete():
                return i
        # If all complete
        return self.num_unique_structures

    def _run_all_jobs(self, jobrunner):
        if self.num_structures_to_run is None:
            # run all jobs if num_structures_to_run is not specified
            jobs_to_run = self.all_structures_opt_jobs
        else:
            jobs_to_run = self.incomplete_structure_opt_jobs[
                : self.num_structures_to_run
            ]
        for job in jobs_to_run:
            job.run(jobrunner=jobrunner)

    def _run(self, jobrunner):
        self._run_all_jobs(jobrunner=jobrunner)

    def is_complete(self):
        return self._run_all_sa_opt_jobs_are_complete()

    def _run_all_sa_opt_jobs_are_complete(self):
        if self.num_structures_to_run is None:
            return all(
                job.is_complete() for job in self.all_structures_opt_jobs
            )
        return all(
            job.is_complete()
            for job in self.all_structures_opt_jobs[
                : self.num_structures_to_run
            ]
        )
