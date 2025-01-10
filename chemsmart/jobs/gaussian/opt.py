import logging
from chemsmart.jobs.gaussian.job import GaussianJob

logger = logging.getLogger(__name__)


class GaussianGeomOptJob(GaussianJob):
    TYPE = "g16opt"

    def __init__(self, molecule, settings, label, **kwargs):
        super().__init__(
            molecule=molecule, settings=settings, label=label, **kwargs
        )
