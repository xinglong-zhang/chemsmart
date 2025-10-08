"""
PyMOL molecular visualization job implementations.

This module defines basic PyMOL visualization jobs used to generate
static molecular images with configurable styling, labeling, and
rendering options. Jobs in this module rely on the shared PyMOLJob
infrastructure for file management and execution via job runners.
"""

from chemsmart.jobs.mol.job import PyMOLJob


class PyMOLVisualizationJob(PyMOLJob):
    """
    PyMOL job for basic molecular visualization.

    Provides standard molecular visualization capabilities using PyMOL
    for creating static molecular images with customizable styling,
    labeling, and rendering options for publication-quality figures.

    Attributes:
        TYPE (str): Job type identifier ('pymol_visualization').
        molecule: Molecule object to visualize.
        label (str): Job identifier used for file naming and outputs.
        jobrunner (JobRunner): Execution backend for running the job.
        skip_completed (bool): If True, completed jobs are not rerun.
    """

    TYPE = "pymol_visualization"

    def __init__(
        self,
        molecule,
        label,
        jobrunner=None,
        **kwargs,
    ):
        """
        Initialize a PyMOL visualization job.

        Sets up a basic molecular visualization job with standard
        PyMOL rendering capabilities for creating static molecular
        images with professional styling.

        Args:
            molecule: Molecule object to visualize.
            label: Job identifier string.
            jobrunner: Job execution runner (default: None).
            **kwargs: Additional arguments passed to parent PyMOLJob.
        """
        super().__init__(
            molecule=molecule,
            label=label,
            jobrunner=jobrunner,
            **kwargs,
        )


class PyMOLHybridVisualizationJob(PyMOLVisualizationJob):
    """Pymol job for hybrid visualization
    Default settings for hybrid visualization of molecules,
    with selected groups highlighted in different colors,
    background in faded color and surface style
    - Group 1: cbap color scheme
    - Group 2: cbac color scheme
    - Group 3: cbay color scheme
    - Group 4: cbag color scheme"""

    TYPE = "pymol_hybrid_visualization"

    def __init__(
        self,
        group1,
        group2=None,
        group3=None,
        group4=None,
        group1_color=None,
        group2_color=None,
        group3_color=None,
        group4_color=None,
        stick_radius=None,
        surface_transcparency=None,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.group1 = group1
        self.group2 = group2
        self.group3 = group3
        self.group4 = group4
        self.group1_color = group1_color
        self.group2_color = group2_color
        self.group3_color = group3_color
        self.group4_color = group4_color
        self.stick_radius = stick_radius
        self.surface_transcparency = surface_transcparency
