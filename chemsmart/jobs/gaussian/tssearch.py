from pyatoms.jobs.gaussian.job import GaussianJob


class GaussianTSJob(GaussianJob):
    TYPE = "g16ts"

    def __init__(self, folder, atoms, settings, **kwargs):
        super().__init__(
            folder=folder, atoms=atoms, settings=settings, **kwargs
        )
