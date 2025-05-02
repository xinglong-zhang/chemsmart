from chemsmart.jobs.gaussian.job import GaussianJob


class GaussianOptJob(GaussianJob):
    TYPE = "g16opt"

    def __init__(self, molecule, settings, label, jobrunner=None, **kwargs):
        super().__init__(
            molecule=molecule,
            settings=settings,
            label=label,
            jobrunner=jobrunner,
            **kwargs,
        )
