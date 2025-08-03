import os.path

from chemsmart.jobs.mol.job import PyMOLJob

class PyMOLAlignJob(PyMOLJob):
    TYPE = "pymol_align"

    def __init__(
        self,
        molecule,
        label,
        jobrunner = None,
        **kwargs,
    ):
        super().__init__(
            molecule=molecule,
            label=label,
            jobrunner=jobrunner,
            **kwargs,
        )

