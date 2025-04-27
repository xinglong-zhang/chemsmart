from chemsmart.jobs.mol.job import PyMOLJob


class PyMOLMOJob(PyMOLJob):
    TYPE = "pymol_mo"

    def __init__(
        self,
        molecule,
        label,
        number=None,
        homo=None,
        lumo=None,
        **kwargs,
    ):
        super().__init__(
            molecule=molecule,
            label=label,
            **kwargs,
        )
        self.number = number
        self.homo = homo
        self.lumo = lumo
