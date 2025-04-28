from chemsmart.jobs.gaussian.job import GaussianJob


class GaussianWBIJob(GaussianJob):
    TYPE = "g16wbi"

    def __init__(self, molecule, settings, label, jobrunner, **kwargs):
        super().__init__(
            molecule=molecule,
            settings=settings,
            label=label,
            jobrunner=jobrunner,
            **kwargs,
        )

        self.settings.freq = False  # turn off freq calc for WBI
