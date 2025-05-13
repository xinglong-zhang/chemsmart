import logging

from chemsmart.jobs.gaussian.job import GaussianJob

logger = logging.getLogger(__name__)


class GaussianRESPJob(GaussianJob):
    TYPE = "g16resp"

    def __init__(self, molecule, settings, label, jobrunner=None, **kwargs):
        super().__init__(
            molecule=molecule,
            settings=settings,
            label=label,
            jobrunner=jobrunner,
            **kwargs,
        )
