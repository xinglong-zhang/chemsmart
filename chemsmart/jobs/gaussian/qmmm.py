import os

from chemsmart.jobs.gaussian.job import GaussianJob


class GaussianQMMMJob(GaussianJob):
    TYPE = "g16qmmm"

    def __init__(self, molecule, settings, label, **kwargs):
        super().__init__(
            molecule=molecule, settings=settings, label=label, **kwargs
        )

    def _output(self):
        """Return a Gaussian16QMMMOutput parser for this job's logfile."""
        if not os.path.exists(self.outputfile):
            return None

        try:
            from chemsmart.io.gaussian.output import Gaussian16QMMMOutput

            return Gaussian16QMMMOutput(filename=self.outputfile)
        except AttributeError:
            return None
