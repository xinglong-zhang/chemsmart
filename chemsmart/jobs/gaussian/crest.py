import logging

from chemsmart.jobs.gaussian.job import GaussianGeneralJob, GaussianJob

logger = logging.getLogger(__name__)


class GaussianCrestJob(GaussianJob):
    """Object representing a Gaussian CREST job that takes in a list of molecules object."""

    TYPE = "g16crest"

    def __init__(
        self,
        molecules,
        settings=None,
        label=None,
        jobrunner=None,
        num_confs_to_run=None,
        grouping_strategy=None,
        skip_completed=True,
        **kwargs,
    ):
        if not isinstance(molecules, list) and len(molecules) == 0:
            raise ValueError("Molecules must be a list of Molecule objects.")

        super().__init__(
            molecule=molecules[0],
            settings=settings,
            label=label,
            jobrunner=jobrunner,
            skip_completed=skip_completed,
            **kwargs,
        )

        if num_confs_to_run is None:
            num_confs_to_run = len(molecules)

        self.num_confs_to_opt = num_confs_to_run

        # if grouping strategy is provided, set the grouper
        # and carry out the grouping before running the group of molecules
        if grouping_strategy is not None:
            logger.info(f"Using grouping strategy: {grouping_strategy}")
            from chemsmart.utils.grouper import StructureGrouperFactory

            logger.info(f"Total structures to group: {len(molecules)}")
            grouper = StructureGrouperFactory.create(
                molecules, strategy=grouping_strategy, **kwargs
            )
            grouper.group()
            unique_molecules = grouper.unique()
            self.grouper = grouper
            logger.debug(f"Grouping strategy: {grouper.__repr__()}")
            logger.info(f"Number of unique groups: {len(unique_molecules)}")
            logger.info(f"Unique molecules: {unique_molecules}")
            self.all_conformers = unique_molecules

        else:
            # if no grouping strategy is provided, use all molecules as conformers
            self.grouper = None
            self.all_conformers = molecules

    @property
    def num_conformers(self):
        return len(self.all_conformers)

    @property
    def last_run_job_index(self):
        return self._check_last_finished_job_index()

    @property
    def all_conformers_opt_jobs(self):
        return self._prepare_all_jobs()

    @property
    def incomplete_conformers_opt_jobs(self):
        return [
            job
            for job in self.all_conformers_opt_jobs
            if not job.is_complete()
        ]

    def _check_last_finished_job_index(self):
        for i, job in enumerate(self.all_conformers_opt_jobs):
            if not job.is_complete():
                return i

        # If all complete
        return self.num_conformers

    def _prepare_all_jobs(self):
        jobs = []
        for i in range(self.num_conformers):
            label = f"{self.label}_c{i + 1}"  # 1-indexed for conformers
            jobs += [
                GaussianGeneralJob(
                    molecule=self.all_conformers[i],
                    settings=self.settings,
                    label=label,
                    jobrunner=self.jobrunner,
                    skip_completed=self.skip_completed,
                )
            ]
        return jobs

    def _run_all_jobs(self):
        for job in self.all_conformers_opt_jobs[: self.num_confs_to_opt]:
            job.run()

    def _run(self):
        self._run_all_jobs()

    def is_complete(self):
        return self._run_all_crest_opt_jobs_are_complete()

    def _run_all_crest_opt_jobs_are_complete(self):
        return all(
            job.is_complete()
            for job in self.all_conformers_opt_jobs[: self.num_confs_to_opt]
        )
