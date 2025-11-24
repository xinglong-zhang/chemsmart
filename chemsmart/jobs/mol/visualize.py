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

    **Hybrid Visualization Features**
    - Highlight multiple atom groups using user-defined selections.
    - Assign independent colors to each group.
    - Render background atoms using a faded color palette.
    - Optionally customize surface color and transparency.
    - Override default atomic colors (C, N, O, S, P) with user-specified RGB values.

    **Default group color schemes**
    - Group 1: ``cbap``
    - Group 2: ``cbac``
    - Group 3: ``cbay``
    - Group 4: ``cbag``

    **Command-Line Options**
    - ``--hybrid``: Activate hybrid visualization mode.
    - ``-g, --group``: Atom indices to include in a highlight group.
      Accepts ranges (e.g., ``1-5``) or comma-separated lists
      (e.g., ``6,7,8``). Repeatable for multiple groups.
    - ``-C, --color``: Color for each group. Must correspond to the number
      of ``-g`` occurrences.
    - ``--surface-color``: Custom surface color for the molecule.
    - ``--surface-transparency``: Custom transparency for the molecular surface.
    - ``--new-color-carbon``: Override carbon atom color using an RGB triplet.
    - ``--new-color-nitrogen``: Override nitrogen atom color using an RGB triplet.
    - ``--new-color-oxygen``: Override oxygen atom color using an RGB triplet.
    - ``--new-color-sulfur``: Override sulfur atom color using an RGB triplet.
    - ``--new-color-phosphorus``: Override phosphorus atom color using an RGB triplet.

    Attributes:
        TYPE (str): Job type identifier (``'pymol_hybrid_visualization'``).
        molecule: Molecule object to visualize.
        label (str): Job identifier used for file naming and outputs.
        jobrunner (JobRunner): Execution backend for running the job.
        skip_completed (bool): If True, completed jobs are not rerun.

    """

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
