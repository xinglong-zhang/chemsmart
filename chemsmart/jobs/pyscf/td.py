"""PySCF TDA/TDDFT vertical-excitation job."""

from chemsmart.jobs.pyscf.job import PySCFJob


class PySCFTDJob(PySCFJob):
    """Vertical excitations of the supplied geometry.

    A Kohn-Sham reference (closed shell: singlet or triplet manifold; open
    shell: the one unrestricted manifold) followed by a TDA or TDDFT
    response stage.  The geometry is fixed; roots are ascending indices
    within the manifold at this geometry.  Stages: ``scf, td``.
    """

    TYPE = "pyscf_td"
