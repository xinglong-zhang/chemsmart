from chemsmart.jobs.gaussian.job import GaussianJob


class GaussianModredJob(GaussianJob):
    TYPE = "g16modred"

    def __init__(self, molecule, settings, label, **kwargs):
        super().__init__(
            molecule=molecule, settings=settings, label=label, **kwargs
        )
