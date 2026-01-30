"""
xTB single point job implementation.

This module provides the XTBSinglePointJob class for performing
single point energy calculations using xTB.
"""

from chemsmart.jobs.xtb.job import XTBJob


class XTBSinglePointJob(XTBJob):
    """
    xTB single point energy calculation job.

    This class handles single point energy calculations at fixed geometries
    without geometry optimization.

    Attributes:
        TYPE (str): Job type identifier ('xtbsp').
        molecule (Molecule): Molecular structure for calculation.
        settings (XTBJobSettings): Calculation configuration options.
        label (str): Job identifier used for file naming.
        jobrunner (JobRunner): Execution backend that runs the job.
        skip_completed (bool): If True, completed jobs are not rerun.
    """

    TYPE = "xtbsp"

    def __init__(self, molecule, settings, label, jobrunner=None, **kwargs):
        """
        Initialize XTBSinglePointJob.

        Args:
            molecule: Molecule object for the calculation
            settings: XTBJobSettings instance
            label: Job label for identification
            jobrunner: Job runner instance
            **kwargs: Additional keyword arguments
        """
        super().__init__(
            molecule=molecule,
            settings=settings,
            label=label,
            jobrunner=jobrunner,
            **kwargs,
        )
