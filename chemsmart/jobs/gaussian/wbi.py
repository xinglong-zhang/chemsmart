from pyatoms.jobs.gaussian.job import GaussianJob


class GaussianWBIJob(GaussianJob):
    TYPE = 'g16wbi'

    def __init__(self, folder, atoms, settings, **kwargs):
        super().__init__(folder=folder, atoms=atoms, settings=settings, **kwargs)
