"""
PyMOL NBO perturbation visualization job implementation.
"""

from chemsmart.jobs.mol.job import PyMOLJob


class PyMOLNBOJob(PyMOLJob):
    """
    PyMOL job for visualizing NBO second-order perturbation interactions.
    """

    TYPE = "pymol_nbo"

    def __init__(
        self,
        molecule,
        label,
        analysis_filename,
        threshold=10.0,
        max_interactions=20,
        nbo_basename=None,
        **kwargs,
    ):
        super().__init__(molecule=molecule, label=label, **kwargs)
        self.analysis_filename = analysis_filename
        self.threshold = threshold
        self.max_interactions = max_interactions

        if self.max_interactions <= 0:
            raise ValueError("--max-interactions must be greater than 0.")

        if nbo_basename is None:
            nbo_basename = f"{self.label}_nbo"
        self.nbo_basename = nbo_basename

    def _get_job_basename(self):
        return self.nbo_basename
