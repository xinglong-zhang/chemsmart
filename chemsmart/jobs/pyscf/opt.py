"""PySCF geometry optimisation job."""

import logging

from chemsmart.jobs.pyscf.job import PySCFJob

logger = logging.getLogger(__name__)


class PySCFOptJob(PySCFJob):
    """Geometry optimisation as one explicit workflow node.

    Attributes:
        TYPE (str): Job type identifier ('pyscf_opt').
    """

    TYPE = "pyscf_opt"

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

        A Hessian is deliberately not appended in-process.  The host binds the
        validated optimized-geometry artifact to an explicit ``hess`` node,
        keeping stage identity and evidence consistent across programs.
        """
        return ["scf", "opt"]
