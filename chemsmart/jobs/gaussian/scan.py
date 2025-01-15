from chemsmart.jobs.gaussian.job import GaussianJob


class GaussianScanJob(GaussianJob):
    TYPE = "g16scan"

    def __init__(self, molecule, settings, label, **kwargs):
        super().__init__(
            molecule=molecule, settings=settings, label=label, **kwargs
        )
