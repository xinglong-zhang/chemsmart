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

    @property
    def stages(self):
        """Return the ordered stage list the generated script executes.

        A single point converges the SCF and stops. ``freq`` is ignored here
        even if the project YAML sets it, because a Hessian at a non-optimised
        geometry is not a meaningful frequency calculation -- use ``hess``
        explicitly if that is really wanted.
        """
        return ["scf"]
