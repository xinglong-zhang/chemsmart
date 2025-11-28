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
    """
    PyMOL job for hybrid molecular visualization.

    Extends :class:`PyMOLVisualizationJob` to provide advanced hybrid
    visualization capabilities. This mode selectively highlights user-defined
    atom groups with distinct color schemes while rendering the remainder of
    the molecule in a subdued (faded) background style. Additional options
    enable customization of surface appearance and atom-type recoloring,
    allowing the generation of publication-quality hybrid representations.

    Features:
    - Accepts an arbitrary number of highlight groups.
    - Assign independent colors to each group (optional).
    - Render background atoms using a faded color palette.
    - Optionally customize surface color and transparency.
    - Override default atomic colors (C, N, O, S, P) with user-specified RGB values.

    Command-line integration (matching CLI behavior) should provide:
    - `groups`: a list of group specifications (repeatable `--group`).
    - `colors`: optional list of colors, one per group (repeatable `--color`).
    """

    TYPE = "pymol_hybrid_visualization"

    def __init__(
        self,
        groups,
        colors=None,
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
        """
        Initialize a hybrid visualization job that supports an arbitrary number
        of highlight groups.

        Args:
            groups (Iterable): Sequence of group specifications (parsed from CLI).
            colors (Iterable, optional): Sequence of colors corresponding to groups.
            stick_radius (float, optional): Stick radius for rendering.
            surface_color (str or tuple, optional): Surface color override.
            surface_transparency (float, optional): Surface transparency override.
            new_color_* (tuple, optional): RGB triplets for element recoloring.
            **kwargs: Additional arguments passed to parent PyMOLJob.
        """
        super().__init__(**kwargs)

        # Normalize to lists
        self.groups = list(groups) if groups is not None else []
        self.colors = list(colors) if colors is not None else []

        # Provide compatibility with PyMOLHybridVisualizationJobRunner
        # by populating dynamic attributes group1, group2, ... and color1, color2, ...
        for idx, group in enumerate(self.groups, start=1):
            setattr(self, f"group{idx}", group)
        for idx, color in enumerate(self.colors, start=1):
            setattr(self, f"color{idx}", color)

        # Expose number of groups for convenience
        self.group_count = len(self.groups)

        self.stick_radius = stick_radius
        self.surface_color = surface_color
        self.surface_transparency = surface_transparency
        self.new_color_carbon = new_color_carbon
        self.new_color_nitrogen = new_color_nitrogen
        self.new_color_oxygen = new_color_oxygen
        self.new_color_phosphorus = new_color_phosphorus
        self.new_color_sulfur = new_color_sulfur
