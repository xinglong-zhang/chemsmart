from chemsmart.jobs.gaussian.job import GaussianJob


class GaussianOptJob(GaussianJob):
    TYPE = "g16opt"

    def __init__(self, molecule, settings, label, **kwargs):
        super().__init__(
            molecule=molecule, settings=settings, label=label, **kwargs
        )
