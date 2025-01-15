from chemsmart.jobs.gaussian.job import GaussianJob


class GaussianTSJob(GaussianJob):
    TYPE = "g16ts"

    def __init__(self, molecule, settings, label, **kwargs):
        super().__init__(
            molecule=molecule, settings=settings, label=label, **kwargs
        )