from pyatoms.jobs.gaussian.job import GaussianJob


class GaussianPESScanJob(GaussianJob):
    TYPE = "g16scan"

    def __init__(self, folder, atoms, settings, **kwargs):
        super().__init__(
            folder=folder, atoms=atoms, settings=settings, **kwargs
        )
