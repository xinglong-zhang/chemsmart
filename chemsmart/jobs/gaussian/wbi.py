from chemsmart.jobs.gaussian.job import GaussianJob


class GaussianWBIJob(GaussianJob):
    TYPE = "g16wbi"

    def __init__(self, molecule, settings, label, **kwargs):
        super().__init__(
            molecule=molecule, settings=settings, label=label, **kwargs
        )
