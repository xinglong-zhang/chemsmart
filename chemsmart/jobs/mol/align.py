
from chemsmart.jobs.mol.job import PyMOLJob


class PyMOLAlignJob(PyMOLJob):
    TYPE = "pymol_align"

    def __init__(
        self,
        molecule,
        label,
        jobrunner=None,
        **kwargs,
    ):
        self.xyz_absolute_paths = []
        self.mol_names = []
        super().__init__(
            molecule=molecule,
            label=label,
            jobrunner=jobrunner,
            **kwargs,
        )
