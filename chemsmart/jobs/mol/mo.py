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
        mo_basename=None,
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

        if self.number:
            mo_basename = f"{self.label}_MO{self.number}"
        if self.homo:
            mo_basename = f"{self.label}_HOMO"
        if self.lumo:
            mo_basename = f"{self.label}_LUMO"

        assert mo_basename, (
            "Molecular orbitals should be specified!\n"
            "Please specify MO number, or HOMO or LUMO to plot."
        )

        self.mo_basename = mo_basename
