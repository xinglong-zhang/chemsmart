"""PySCF single point job."""

import logging

from chemsmart.jobs.pyscf.job import PySCFJob

logger = logging.getLogger(__name__)


class PySCFSinglePointJob(PySCFJob):
    """Single point energy at a fixed geometry.

    Attributes:
        TYPE (str): Job type identifier ('pyscf_sp').
    """

    TYPE = "pyscf_sp"

    def __init__(self, molecule, settings, label, jobrunner=None, **kwargs):
        super().__init__(
            molecule=molecule,
            settings=settings,
            label=label,
            jobrunner=jobrunner,
            **kwargs,
        )

    # Stages come from the resolved settings (``PySCFJob.stages``): the SCF
    # alone, or ``scf, corr`` for a correlated method.  ``freq=True`` is
    # rejected during settings validation rather than silently ignored; use
    # an explicit ``hess`` node when a fixed-geometry Hessian is intended.
