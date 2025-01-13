from chemsmart.jobs.gaussian.job import GaussianGeneralJob, GaussianJob


class GaussianCrestOptJob(GaussianJob):
    TYPE = "g16crestopt"

    def __init__(
        self, folder, atoms, settings, num_confs_to_opt=None, **kwargs
    ):
        super().__init__(
            folder=folder, atoms=atoms, settings=settings, **kwargs
        )

        if num_confs_to_opt is None:
            num_confs_to_opt = len(atoms)

        self.all_conformers = atoms
        self.num_confs_to_opt = num_confs_to_opt

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
                    folder=self.folder,
                    atoms=self.all_conformers[i],
                    settings=self.settings,
                    label=label,
                )
            ]
        return jobs

    def _run_all_jobs(self, jobrunner, queue_manager=None):
        for job in self.all_conformers_opt_jobs[: self.num_confs_to_opt]:
            job.run(jobrunner=jobrunner, queue_manager=queue_manager)

    def _run(self, jobrunner, queue_manager=None, **kwargs):
        self._run_all_jobs(jobrunner=jobrunner, queue_manager=queue_manager)

    def is_complete(self):
        return self._run_all_crest_opt_jobs_are_complete()

    def _run_all_crest_opt_jobs_are_complete(self):
        return all(
            job.is_complete()
            for job in self.all_conformers_opt_jobs[: self.num_confs_to_opt]
        )
