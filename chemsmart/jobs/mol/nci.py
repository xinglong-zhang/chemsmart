from chemsmart.jobs.mol.job import PyMOLJob


class PyMOLNCIJob(PyMOLJob):
    TYPE = "pymol_nci"

    def __init__(
        self,
        molecule,
        label,
        isosurface=0.5,
        range=1.0,
        binary=False,
        **kwargs,
    ):
        super().__init__(
            molecule=molecule,
            label=label,
            **kwargs,
        )
        self.isosurface = isosurface
        self.range = range
        self.binary = binary
