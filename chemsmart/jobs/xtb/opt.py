"""
xTB geometry optimization job implementation.

This module provides the XTBOptJob class for performing
molecular geometry optimization calculations using xTB.
"""

from chemsmart.jobs.xtb.job import XTBJob


class XTBOptJob(XTBJob):
    """
    xTB geometry optimization job.

    This class handles molecular geometry optimizations to find energy
    minima on the potential energy surface.

    Attributes:
        TYPE (str): Job type identifier ('xtbopt').
        molecule (Molecule): Molecular structure to optimize.
        settings (XTBJobSettings): Calculation configuration options.
        label (str): Job identifier used for file naming.
        jobrunner (JobRunner): Execution backend that runs the job.
        skip_completed (bool): If True, completed jobs are not rerun.
    """

    TYPE = "xtbopt"

    def __init__(self, molecule, settings, label, jobrunner=None, **kwargs):
        """
        Initialize XTBOptJob.

        Args:
            molecule: Molecule object for the optimization
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
