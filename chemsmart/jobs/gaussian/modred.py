from pyatoms.jobs.gaussian.job import GaussianJob


class GaussianModredundantJob(GaussianJob):
    TYPE = 'g16modred'

    def __init__(self, folder, atoms, settings, **kwargs):
        super().__init__(folder=folder, atoms=atoms, settings=settings, **kwargs)
