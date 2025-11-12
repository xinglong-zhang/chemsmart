from chemsmart.jobs.mol.job import PyMOLJob


class PyMOLAlignJob(PyMOLJob):
    TYPE = "pymol_align"

    def __init__(
        self,
        molecule,
        label,
        align_basename=None,
        use_raw_label=False,
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

        if align_basename is None:
            if use_raw_label:
                # Use the label directly without modification
                align_basename = self.label
            else:
                # Use the original logic with molecule count suffix
                n_molecules = (
                    len(molecule) if isinstance(molecule, list) else 1
                )
                if n_molecules > 2:
                    align_basename = (
                        f"{self.label}_and_{n_molecules-1}_molecules_align"
                    )
                else:
                    align_basename = f"{self.label}_and_1_molecule_align"

        self.align_basename = align_basename

    def _get_job_basename(self):
        """
        Internal method to derive the job base name.
        Job specific implementation that overrides parent class method.

        Returns:
            str: Base name derived from the align basename.
        """
        return self.align_basename
