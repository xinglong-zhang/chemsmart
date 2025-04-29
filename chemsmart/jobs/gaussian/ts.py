from chemsmart.jobs.gaussian.job import GaussianJob


class GaussianTSJob(GaussianJob):
    TYPE = "g16ts"

    def __init__(self, molecule, settings, label, jobrunner=None, **kwargs):
        super().__init__(
            molecule=molecule,
            settings=settings,
            label=label,
            jobrunner=jobrunner,
            **kwargs,
        )
