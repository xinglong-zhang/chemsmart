from chemsmart.jobs.gaussian.job import GaussianJob


class GaussianNCIJob(GaussianJob):
    TYPE = "g16nci"

    def __init__(self, molecule, settings, label, **kwargs):
        super().__init__(
            molecule=molecule, settings=settings, label=label, **kwargs
        )

        self.settings.freq = False  # turn off freq calc for NCI
