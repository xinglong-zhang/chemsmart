from chemsmart.jobs.gaussian.job import GaussianJob


class GaussianQMMMJob(GaussianJob):
    TYPE = "g16qmmm"

    def __init__(self, molecule, settings, label, **kwargs):
        super().__init__(
            molecule=molecule, settings=settings, label=label, **kwargs
        )
