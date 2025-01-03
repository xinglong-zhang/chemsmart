from pyatoms.io.gaussian.utils import get_list_from_string_range
from pyatoms.jobs.gaussian.job import GaussianGeneralJob, GaussianJob


class GaussianDIASJob(GaussianJob):
    TYPE = 'g16dias'

    def __init__(self, folder, atoms, settings, fragment_indices, every_n_points, **kwargs):
        super().__init__(folder=folder, atoms=atoms, settings=settings, **kwargs)
        self.all_atoms = atoms
        self.fragment_indices = fragment_indices
        self.every_n_points = every_n_points

    @property
    def num_images(self):
        return len(self.all_atoms)

    def _fragment_structure(self, atoms):
        fragment1_indices = get_list_from_string_range(self.fragment_indices)
        fragment2_indices = [i for i in range(len(atoms)) if i not in fragment1_indices]

        fragment1_atoms = atoms[fragment1_indices]
        fragment2_atoms = atoms[fragment2_indices]
        return fragment1_atoms, fragment2_atoms

    @property
    def fragment1_atoms(self):
        return [self._fragment_structure(atoms=atoms)[0] for atoms in self.all_atoms]

    @property
    def fragment2_atoms(self):
        return [self._fragment_structure(atoms=atoms)[1] for atoms in self.all_atoms]

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
        for i, atoms in enumerate(images):
            label = f'{self.label}_p{i}_f1'
            jobs += [GaussianGeneralJob(folder=self.folder, atoms=atoms, settings=self.settings, label=label)]
        return jobs

    @property
    def fragment2_jobs(self):
        images = self._sample_images(self.fragment2_atoms)
        jobs = []
        for i, atoms in enumerate(images):
            label = f'{self.label}_p{i}_f2'
            jobs += [GaussianGeneralJob(folder=self.folder, atoms=atoms, settings=self.settings, label=label)]
        return jobs

    @property
    def all_atoms_jobs(self):
        images = self._sample_images(self.all_atoms)
        jobs = []
        for i, atoms in enumerate(images):
            label = f'{self.label}_p{i}'
            jobs += [GaussianGeneralJob(folder=self.folder, atoms=atoms, settings=self.settings, label=label)]
        return jobs

    def _run_all_atoms_jobs(self, jobrunner, queue_manager=None):
        for job in self.all_atoms_jobs:
            job.run(jobrunner=jobrunner, queue_manager=queue_manager)

    def _run_fragment1_jobs(self, jobrunner, queue_manager=None):
        for job in self.fragment1_jobs:
            job.run(jobrunner=jobrunner, queue_manager=queue_manager)

    def _run_fragment2_jobs(self, jobrunner, queue_manager=None):
        for job in self.fragment2_jobs:
            job.run(jobrunner=jobrunner, queue_manager=queue_manager)

    def _run(self, jobrunner, queue_manager=None, **kwargs):
        self._run_all_atoms_jobs(jobrunner=jobrunner, queue_manager=queue_manager)
        self._run_fragment1_jobs(jobrunner=jobrunner, queue_manager=queue_manager)
        self._run_fragment2_jobs(jobrunner=jobrunner, queue_manager=queue_manager)

    def is_complete(self):
        return (
            self._run_all_atoms_jobs_are_complete()
            and self._run_fragment1_jobs_are_complete()
            and self._run_fragment2_jobs_are_complete()
        )

    def _run_all_atoms_jobs_are_complete(self):
        return all(job.is_complete() for job in self.all_atoms_jobs)

    def _run_fragment1_jobs_are_complete(self):
        return all(job.is_complete() for job in self.fragment1_jobs)

    def _run_fragment2_jobs_are_complete(self):
        return all(job.is_complete() for job in self.fragment2_jobs)
