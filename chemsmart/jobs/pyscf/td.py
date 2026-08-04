"""Preview-only PySCF TDA/TDDFT job."""

from chemsmart.jobs.pyscf.job import PySCFJob


class PySCFTDJob(PySCFJob):
    """Closed-shell gas-phase singlet vertical-excitation preview."""

    TYPE = "pyscf_td"

    @property
    def stages(self):
        """Return the intended ground-state and response stages."""

        return ["scf", "td"]
