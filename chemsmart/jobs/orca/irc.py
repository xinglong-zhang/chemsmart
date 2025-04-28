import logging

from chemsmart.jobs.orca.job import ORCAJob

logger = logging.getLogger(__name__)


class ORCAIRCJob(ORCAJob):
    TYPE = "orcairc"

    def __init__(self, molecule, settings, label, jobrunner=None, **kwargs):
        super().__init__(
            molecule=molecule,
            settings=settings,
            label=label,
            jobrunner=jobrunner,
            **kwargs,
        )
