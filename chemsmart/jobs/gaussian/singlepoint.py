from chemsmart.jobs.gaussian.job import GaussianJob


class GaussianSinglePointJob(GaussianJob):
    TYPE = "g16sp"

    def __init__(self, molecule, settings, label, **kwargs):
        super().__init__(
            molecule=molecule, settings=settings, label=label, **kwargs
        )