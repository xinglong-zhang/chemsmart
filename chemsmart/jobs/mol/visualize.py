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
        color1=None,
        color2=None,
        color3=None,
        color4=None,
        stick_radius=None,
        surface_color=None,
        surface_transparency=None,
        new_color_carbon=None,
        new_color_nitrogen=None,
        new_color_oxygen=None,
        new_color_phosphorus=None,
        new_color_sulfur=None,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.group1 = group1
        self.group2 = group2
        self.group3 = group3
        self.group4 = group4
        self.color1 = color1
        self.color2 = color2
        self.color3 = color3
        self.color4 = color4
        self.stick_radius = stick_radius
        self.surface_color = surface_color
        self.surface_transparency = surface_transparency
        self.new_color_carbon = new_color_carbon
        self.new_color_nitrogen = new_color_nitrogen
        self.new_color_oxygen = new_color_oxygen
        self.new_color_phosphorus = new_color_phosphorus
        self.new_color_sulfur = new_color_sulfur
