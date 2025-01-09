from pyatoms.jobs.gaussian.job import GaussianJob


class GaussianCustomJob(GaussianJob):
    TYPE = "g16job"

    def __init__(self, folder, atoms, settings, **kwargs):
        super().__init__(
            folder=folder, atoms=atoms, settings=settings, **kwargs
        )
