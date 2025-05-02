import os

from chemsmart.jobs.mol.job import PyMOLJob


class PyMOLNCIJob(PyMOLJob):
    TYPE = "pymol_nci"

    def __init__(
        self,
        molecule,
        label,
        isosurface=0.5,
        color_range=1.0,
        binary=False,
        intermediate=False,
        nci_basename=None,
        **kwargs,
    ):
        super().__init__(
            molecule=molecule,
            label=label,
            **kwargs,
        )
        self.isosurface = isosurface
        self.color_range = color_range
        self.binary = binary
        self.intermediate = intermediate

        if nci_basename is None:
            nci_basename = self.label

        if self.binary:
            nci_basename += "_binary"
        if self.intermediate:
            nci_basename += "_intermediate"

        self.nci_basename = nci_basename

    def _job_is_complete(self):
        """PyMOL MO job is complete if the corresponding .pse file exists."""
        return os.path.exists(f"{self.nci_basename}.pse")
