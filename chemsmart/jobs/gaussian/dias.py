from chemsmart.utils.utils import get_list_from_string_range
from chemsmart.jobs.gaussian.job import GaussianGeneralJob, GaussianJob


class GaussianDIASJob(GaussianJob):
    TYPE = "g16dias"

    def __init__(
        self,
        molecules,
        settings,
        label,
        fragment_indices,
        every_n_points,
        **kwargs,
    ):
        super().__init__(
            molecule=molecules,
            settings=settings,
            label=label,
            **kwargs,
        )
        self.all_molecules = molecules  # alone IRC coordinate
        self.fragment_indices = fragment_indices
        self.every_n_points = every_n_points

    @property
    def num_images(self):
        return len(self.all_molecules)

    def _fragment_structure(self, molecule):
        fragment1_indices = get_list_from_string_range(self.fragment_indices)
        fragment2_indices = [
            i + 1
            for i in range(len(molecule))
            if i + 1 not in fragment1_indices
        ]

        fragment1_atoms = molecule[fragment1_indices]
        fragment2_atoms = molecule[fragment2_indices]
        return fragment1_atoms, fragment2_atoms

    @property
    def fragment1_atoms(self):
        return [
            self._fragment_structure(molecule=molecule)[0]
            for molecule in self.all_molecules
        ]

    @property
    def fragment2_atoms(self):
        return [
            self._fragment_structure(molecule=molecule)[1]
            for molecule in self.all_molecules
        ]

    def _sample_images(self, images):
        """Samples provided images every self.every_n_points."""
        filtered_images = images[0 :: self.every_n_points]
        if (self.num_images - 1) / self.every_n_points != 0:
            filtered_images.append(images[-1])
        return filtered_images

    @property
    def fragment1_jobs(self):
        images = self._sample_images(self.fragment1_atoms)
        jobs = []
        for i, molecule in enumerate(images):
            label = f"{self.label}_p{i}_f1"
            jobs += [
                GaussianGeneralJob(
                    molecule=molecule,
                    settings=self.settings,
                    label=label,
                    skip_completed=self.skip_completed,
                )
            ]
        return jobs

    @property
    def fragment2_jobs(self):
        images = self._sample_images(self.fragment2_atoms)
        jobs = []
        for i, molecule in enumerate(images):
            label = f"{self.label}_p{i}_f2"
            jobs += [
                GaussianGeneralJob(
                    molecule=molecule,
                    settings=self.settings,
                    label=label,
                    skip_completed=self.skip_completed,
                )
            ]
        return jobs

    @property
    def all_molecules_jobs(self):
        images = self._sample_images(self.all_molecules)
        jobs = []
        for i, molecule in enumerate(images):
            label = f"{self.label}_p{i}"
            jobs += [
                GaussianGeneralJob(
                    molecule=molecule,
                    settings=self.settings,
                    label=label,
                    skip_completed=self.skip_completed,
                )
            ]
        return jobs

    def _run_all_molecules_jobs(self, jobrunner):
        for job in self.all_molecules_jobs:
            job.run(jobrunner=jobrunner)

    def _run_fragment1_jobs(self, jobrunner):
        for job in self.fragment1_jobs:
            job.run(jobrunner=jobrunner)

    def _run_fragment2_jobs(self, jobrunner):
        for job in self.fragment2_jobs:
            job.run(jobrunner=jobrunner)

    def _run(self, jobrunner, **kwargs):
        self._run_all_molecules_jobs(jobrunner=jobrunner)
        self._run_fragment1_jobs(jobrunner=jobrunner)
        self._run_fragment2_jobs(jobrunner=jobrunner)

    def is_complete(self):
        return (
            self._run_all_molecules_jobs_are_complete()
            and self._run_fragment1_jobs_are_complete()
            and self._run_fragment2_jobs_are_complete()
        )

    def _run_all_molecules_jobs_are_complete(self):
        return all(job.is_complete() for job in self.all_molecules_jobs)

    def _run_fragment1_jobs_are_complete(self):
        return all(job.is_complete() for job in self.fragment1_jobs)

    def _run_fragment2_jobs_are_complete(self):
        return all(job.is_complete() for job in self.fragment2_jobs)
