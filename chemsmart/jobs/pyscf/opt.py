"""PySCF geometry optimisation job."""

import logging

from chemsmart.jobs.pyscf.job import PySCFJob

logger = logging.getLogger(__name__)


class PySCFOptJob(PySCFJob):
    """Geometry optimisation, optionally followed by a Hessian.

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

        When the project YAML sets ``freq: True`` the Hessian runs in the
        *same* process, against the already-converged mean-field object.
        Splitting it into a second ChemSmart job would re-converge the SCF
        from scratch and re-pay density-fitting factorisation and grid
        construction -- the one place PySCF's in-process nature is a genuine
        advantage over the file-based engines.
        """
        stages = ["scf", "opt"]
        if self.settings.freq:
            stages.append("hess")
        return stages
