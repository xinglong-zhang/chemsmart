from pyatoms.jobs.gaussian.job import GaussianJob


class GaussianSinglePointJob(GaussianJob):
    TYPE = "g16sp"

    def __init__(self, folder, atoms, settings, **kwargs):
        super().__init__(
            folder=folder, atoms=atoms, settings=settings, **kwargs
        )
