from pyatoms.jobs.gaussian.job import GaussianJob


class GaussianNCIJob(GaussianJob):
    TYPE = 'g16nci'

    def __init__(self, folder, atoms, settings, **kwargs):
        super().__init__(folder=folder, atoms=atoms, settings=settings, **kwargs)
