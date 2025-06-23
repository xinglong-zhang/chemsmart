from chemsmart.jobs.gaussian.job import GaussianJob


class GaussianCustomJob(GaussianJob):
    TYPE = "g16job"

    def __init__(
        self, molecule, settings=None, label=None, jobrunner=None, **kwargs
    ):
        super().__init__(
            molecule=molecule,
            settings=settings,
            label=label,
            jobrunner=jobrunner,
            **kwargs,
        )
