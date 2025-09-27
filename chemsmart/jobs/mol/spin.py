"""
PyMOL Spin Density Visualization Jobs Module.

This module provides specialized PyMOL jobs for visualizing molecular spin
density distributions. Handles cube file generation and PyMOL session
creation for spin density analysis.
"""

import os.path

from chemsmart.jobs.mol.job import PyMOLJob


class PyMOLSpinJob(PyMOLJob):
    """
    PyMOL job for spin density visualization.

    Specialized PyMOL job class for creating spin density visualizations
    from quantum chemistry calculations. Generates cube files and PyMOL
    sessions for analyzing electron spin distribution in molecules.

    Attributes:
        TYPE (str): Job type identifier ('pymol_spin').
        molecule: Molecular structure object to visualize.
        label (str): Job identifier used for file naming and outputs.
        spin_basename (str): Basename used for spin density outputs
            (e.g., cube and session files).
        npts (int): Grid points per side used for cube generation.
        jobrunner (JobRunner): Execution backend for running the job.
        skip_completed (bool): If True, completed jobs are not rerun.
    """

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
        Initialize PyMOL spin density visualization job.

        Creates a PyMOL job for generating spin density visualizations
        from quantum chemistry calculations. Configures cube file
        generation parameters and output file naming.

        Args:
            molecule: Molecular structure object for visualization.
            label: Unique identifier label for the job.
            spin_basename: Base name for spin density output files.
                Defaults to "{label}_spin" if not provided.
            npts: Number of grid points per side for cube generation
                (default: 80). Positive values specify points per side
                (e.g., 100 => 100^3 total grid points).
            **kwargs: Additional keyword arguments passed to parent class.
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
        """
        Job is complete if the spin density .pse session exists.
        """
        return os.path.exists(f"{self.spin_basename}.pse")

    @property
    def outputfile(self):
        """
        Get the path to the spin density PyMOL session file output.

        Returns:
            str: Absolute path to the spin density PSE session file.
        """
        outputfile = self.spin_basename + ".pse"
        return os.path.join(self.folder, outputfile)
