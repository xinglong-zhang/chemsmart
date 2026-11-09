"""PySCF intrinsic reaction coordinate job."""

import logging

from chemsmart.jobs.pyscf.job import PySCFJob

logger = logging.getLogger(__name__)


class PySCFIRCJob(PySCFJob):
    """Follow the reaction path down from a transition state.

    Attributes:
        TYPE (str): Job type identifier ('pyscf_irc').
    """

    TYPE = "pyscf_irc"
    JOBTYPE = "irc"

    @property
    def stages(self):
        """Return the ordered stage list the generated script executes.

        The ``irc`` stage computes the analytic Hessian of the starting
        structure itself, to confirm it is a first-order saddle point and to
        seed geomeTRIC, so ``freq`` adds no stage here.
        """
        return ["scf", "irc"]
