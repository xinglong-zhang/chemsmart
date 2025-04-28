from chemsmart.jobs.gaussian.job import GaussianJob


class GaussianModredJob(GaussianJob):
    TYPE = "g16modred"

    def __init__(self, molecule, settings, label, jobrunner=None, **kwargs):
        super().__init__(
            molecule=molecule,
            settings=settings,
            label=label,
            jobrunner=jobrunner,
            **kwargs,
        )
