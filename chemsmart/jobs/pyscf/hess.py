"""PySCF Hessian / harmonic frequency job."""

import logging

from chemsmart.jobs.pyscf.job import PySCFJob

logger = logging.getLogger(__name__)


class PySCFHessJob(PySCFJob):
    """Analytic Hessian and harmonic frequencies at a fixed geometry.

    Attributes:
        TYPE (str): Job type identifier ('pyscf_hess').
    """

    TYPE = "pyscf_hess"

    @property
    def stages(self):
        """Return the ordered stage list the generated script executes.

        The geometry is taken as given and is *not* re-optimised. Frequencies
        are only physically meaningful at a stationary point of the same
        method and basis, so the caller is responsible for supplying an
        optimised structure.
        """
        return ["scf", "hess"]
