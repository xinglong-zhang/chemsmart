"""PySCF transition-state search job."""

import logging

from chemsmart.jobs.pyscf.job import PySCFJob

logger = logging.getLogger(__name__)


class PySCFTSJob(PySCFJob):
    """First-order saddle-point search, optionally followed by a Hessian.

    Attributes:
        TYPE (str): Job type identifier ('pyscf_ts').
    """

    TYPE = "pyscf_ts"

    @property
    def stages(self):
        """Return the ordered stage list the generated script executes.

        geomeTRIC's P-RFO search starts from the analytic Hessian of the
        guess. With ``freq: True`` the Hessian at the result counts its
        imaginary modes, and a structure that is not a first-order saddle
        point fails the search instead of being reported as one.
        """
        stages = ["scf", "ts"]
        if self.settings.freq:
            stages.append("hess")
        return stages
