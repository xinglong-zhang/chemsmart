from chemsmart.jobs.xtb.job import XTBJob


class XTBOptJob(XTBJob):
    TYPE = "xtbopt"

    def __init__(self, molecule, settings, label, jobrunner=None, **kwargs):
        super().__init__(
            molecule=molecule,
            settings=settings,
            label=label,
            jobrunner=jobrunner,
            **kwargs,
        )
