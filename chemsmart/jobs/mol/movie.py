"""
PyMOL molecular movie job implementation.

Defines the PyMOLMovieJob used to generate rotating molecular animations
by producing frame sequences suitable for video encoding. Relies on the
shared PyMOLJob base class for file management and execution.
"""

from chemsmart.jobs.mol.job import PyMOLJob


class PyMOLMovieJob(PyMOLJob):
    """
    PyMOL job for creating molecular animation movies.

    Specialized job class for generating rotating molecular animations
    using PyMOL, creating frame sequences that can be converted to
    video files for dynamic molecular visualization presentations.

    Attributes:
        TYPE (str): Job type identifier ('pymol_movie').
        molecule: Molecule object to animate.
        label (str): Job identifier used for file naming (suffix '_movie' added).
        jobrunner (JobRunner): Execution backend for running the job.
        overwrite (bool): Overwrite an existing MP4 when post-processing.
        skip_completed (bool): If True, completed jobs are not rerun.
    """

    TYPE = "pymol_movie"

    def __init__(
        self,
        molecule,
        label,
        jobrunner=None,
        overwrite=False,
        **kwargs,
    ):
        """
        Initialize a PyMOL movie generation job.

        Sets up the job for creating animated molecular visualizations
        with rotation effects and optional ray tracing for high-quality
        output videos.

        Args:
            molecule: Molecule object to animate.
            label: Job identifier string.
            jobrunner: Job execution runner (default: None).
            overwrite: Whether to overwrite existing movie files (default: False).
            **kwargs: Additional arguments passed to parent PyMOLJob.

        Note:
            The internal `label` used for outputs is suffixed with `_movie` to
            distinguish movie artifacts from other visualizations.
        """
        super().__init__(
            molecule=molecule,
            label=label,
            jobrunner=jobrunner,
            **kwargs,
        )
        self.overwrite = overwrite

        self.label = f"{label}_movie"
