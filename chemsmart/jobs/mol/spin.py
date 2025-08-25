import os.path

from chemsmart.jobs.mol.job import PyMOLJob


class PyMOLSpinJob(PyMOLJob):
    TYPE = "pymol_spin"

    def __init__(
        self,
        molecule,
        label,
        spin_basename=None,
        npts=80,
        **kwargs,
    ):
        """
        PyMOL spin density job.
        Args:
            npts: Number of points per side in the cube.
            A value of 0 selects default of 80^3 points distributed evenly over a rectangular grid.
            Positive values of npts specify the number of points per “side”;
            e.g., 100 specified a grid of 1,000,000 (100^3) points.
        """
        super().__init__(
            molecule=molecule,
            label=label,
            **kwargs,
        )

        if spin_basename is None:
            spin_basename = f"{self.label}_spin"

        self.spin_basename = spin_basename
        self.npts = npts

    def _job_is_complete(self):
        """PyMOL MO job is complete if the corresponding .pse file exists."""
        return os.path.exists(f"{self.spin_basename}.pse")
