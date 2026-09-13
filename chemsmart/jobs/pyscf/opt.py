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

    # Stages come from the resolved settings (``PySCFJob.stages``): a
    # ground-state optimisation runs ``scf, opt``; an optimisation on an
    # excited root re-evaluates the spectrum at the reached geometry
    # (``scf, opt, td``); a correlated optimisation computes the method's
    # components there (``scf, opt, corr``).  A Hessian is deliberately not
    # appended in-process: the host binds the validated optimized-geometry
    # artifact to an explicit ``hess`` node, keeping stage identity and
    # evidence consistent across programs.
