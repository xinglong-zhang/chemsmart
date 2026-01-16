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

    def _get_job_basename(self):
        """
        Internal method to derive the job base name.
        Job specific implementation that overrides parent class method.

        Returns:
            str: Base name derived from the label.
        """
        return self.label
