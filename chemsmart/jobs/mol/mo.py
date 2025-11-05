"""
PyMOL molecular orbital (MO) visualization job implementation.

Defines the PyMOLMOJob used to visualize molecular orbitals (HOMO, LUMO,
or a specific orbital index) via isosurfaces and coloring. Works in
conjunction with the PyMOLMOJobRunner, which generates cube files and
drives the PyMOL rendering workflow.
"""

from chemsmart.jobs.mol.job import PyMOLJob


class PyMOLMOJob(PyMOLJob):
    """
    PyMOL job for molecular orbital visualization.

    Specialized job class for visualizing molecular orbitals using
    PyMOL, supporting HOMO, LUMO, or specific orbital number
    visualization with isosurface generation and coloring.

    Attributes:
        TYPE (str): Job type identifier ('pymol_mo').
        molecule: Molecule object to visualize.
        label (str): Job identifier used for file naming and outputs.
        number (int | None): Specific MO number to visualize.
        homo (bool | None): Whether to visualize the HOMO.
        lumo (bool | None): Whether to visualize the LUMO.
        mo_basename (str): Basename for MO-related artifacts (cube, pml, pse).
        jobrunner (JobRunner): Execution backend for running the job.
        skip_completed (bool): If True, completed jobs are not rerun.
    """

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
        """
        Initialize a PyMOL molecular orbital visualization job.

        Sets up the job with molecular structure and orbital specification.
        At least one of `number`, `homo`, or `lumo` must be provided. If
        multiple are provided, the precedence is number < HOMO < LUMO,
        i.e., LUMO overrides HOMO which overrides number.

        Args:
            molecule: Molecule object to visualize.
            label: Job identifier string.
            number (int, optional): Specific molecular orbital number.
            homo (bool, optional): Visualize HOMO if True.
            lumo (bool, optional): Visualize LUMO if True.
            mo_basename (str, optional): Base name for output files
                (auto-generated if not provided).
            **kwargs: Additional arguments passed to parent PyMOLJob.

        Raises:
            AssertionError: If no orbital specification is provided.
        """
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

    def _get_job_basename(self):
        """
        Internal method to derive the job base name.
        Job specific implementation that overrides parent class method.

        Returns:
            str: Base name derived from the job label.
        """
        return self.mo_basename
