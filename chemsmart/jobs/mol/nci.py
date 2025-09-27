"""
PyMOL non-covalent interaction (NCI) visualization job implementation.

Defines the PyMOLNCIJob used to visualize NCI isosurfaces from density
and reduced gradient data. Works together with PyMOLNCIJobRunner, which
loads cube files and executes the PyMOL commands to render the analysis.
"""

import os

from chemsmart.jobs.mol.job import PyMOLJob


class PyMOLNCIJob(PyMOLJob):
    """
    PyMOL job for Non-Covalent Interactions (NCI) analysis visualization.

    Specialized job class for visualizing non-covalent interactions
    using PyMOL with density and gradient cube files, supporting
    different analysis modes and customizable visualization parameters.

    Attributes:
        TYPE (str): Job type identifier ('pymol_nci').
        molecule: Molecule object to analyze.
        label (str): Job identifier used for file naming and outputs.
        isosurface_value (float): Isosurface value for NCI visualization.
        color_range (float): Range used for coloring the interaction map.
        binary (bool): Whether to use binary NCI analysis mode.
        intermediate (bool): Whether to use intermediate NCI mode.
        nci_basename (str): Basename for NCI artifacts (pml/pse naming).
        jobrunner (JobRunner): Execution backend for running the job.
        skip_completed (bool): If True, completed jobs are not rerun.
    """

    TYPE = "pymol_nci"

    def __init__(
        self,
        molecule,
        label,
        isosurface_value,
        color_range,
        binary=False,
        intermediate=False,
        nci_basename=None,
        **kwargs,
    ):
        """
        Initialize a PyMOL NCI analysis visualization job.

        Sets up the job for non-covalent interaction analysis with
        customizable isosurface levels, color ranges, and analysis
        modes for comprehensive interaction visualization.

        Args:
            molecule: Molecule object to analyze.
            label (str): Job identifier string.
            isosurface_value (float): Isosurface level for NCI visualization (default 0.5).
            color_range (float): Color range for interaction mapping (default 1.0).
            binary (bool): Use binary NCI analysis mode (default False).
            intermediate (bool): Use intermediate NCI analysis mode (default False).
            nci_basename (str, optional): Base name for output files (auto-generated).
            **kwargs: Additional arguments passed to parent PyMOLJob.
        """
        super().__init__(
            molecule=molecule,
            label=label,
            **kwargs,
        )
        # set defaults
        if isosurface_value is None:
            isosurface_value = 0.5
        if color_range is None:
            color_range = 1.0
        self.isosurface_value = isosurface_value
        self.color_range = color_range
        self.binary = binary
        self.intermediate = intermediate

        if nci_basename is None:
            nci_basename = f"{self.label}_nci"

        if self.binary:
            nci_basename += "_binary"
        if self.intermediate:
            nci_basename += "_intermediate"

        self.nci_basename = nci_basename

    def _job_is_complete(self):
        """
        Check if the PyMOL NCI analysis job has completed.

        Determines job completion by checking for the existence of
        the NCI analysis PyMOL session file.

        Returns:
            bool: True if the NCI PSE file exists, False otherwise.
        """
        return os.path.exists(f"{self.nci_basename}.pse")

    @property
    def outputfile(self):
        """
        Get the path to the NCI PyMOL session file output.

        Returns:
            str: Absolute path to the NCI PSE session file.
        """
        outputfile = self.nci_basename + ".pse"
        return os.path.join(self.folder, outputfile)
