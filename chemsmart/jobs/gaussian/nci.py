from chemsmart.jobs.gaussian.job import GaussianJob
from chemsmart.jobs.runner import JobRunner


class GaussianNCIJob(GaussianJob):
    TYPE = "g16nci"

    def __init__(self, molecule, settings, label, jobrunner=None, **kwargs):
        super().__init__(
            molecule=molecule,
            settings=settings,
            label=label,
            jobrunner=JobRunner,
            **kwargs,
        )

        self.settings.freq = False  # turn off freq calc for NCI
