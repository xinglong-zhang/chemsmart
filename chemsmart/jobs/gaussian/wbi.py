"""
Gaussian Wiberg Bond Index (WBI) calculation job implementation.

This module provides the GaussianWBIJob class for performing
Wiberg Bond Index calculations using Gaussian. WBI calculations
analyze bond orders and electronic structure properties of
molecular systems.
"""

from chemsmart.jobs.gaussian.job import GaussianJob


class GaussianWBIJob(GaussianJob):
    """
    Gaussian job class for Wiberg Bond Index (WBI) calculations.

    Performs Wiberg Bond Index calculations to analyze bond orders
    and electronic structure properties. These calculations provide
    quantitative measures of bond strength and character in molecular
    systems through natural bond orbital (NBO) analysis.

    WBI calculations are typically used for:
    - Bond order analysis in organic and inorganic systems
    - Electronic structure characterization
    - Validation of bonding models
    - Studying bond order changes in reactions

    Attributes:
        TYPE (str): Job type identifier ('g16wbi').
        molecule (Molecule): Molecular structure used for the calculation.
        settings (GaussianJobSettings): Calculation configuration; frequency
            analysis (`freq`) is disabled for WBI runs.
        label (str): Job identifier label used for file naming.
        jobrunner (JobRunner): Execution backend that runs the job.
        skip_completed (bool): If True, completed jobs are not rerun.
    """

    TYPE = "g16wbi"

    def __init__(self, molecule, settings, label, jobrunner, **kwargs):
        """
        Initialize a Gaussian Wiberg Bond Index calculation.

        Sets up WBI calculation with appropriate settings. Frequency
        calculations are automatically disabled as they are not
        needed for WBI analysis.

        Args:
            molecule (Molecule): Molecular structure for WBI calculation.
            settings (GaussianJobSettings): Calculation configuration.
            label (str): Job identifier label.
            jobrunner: Job execution handler.
            **kwargs: Additional keyword arguments for parent class.
        """
        super().__init__(
            molecule=molecule,
            settings=settings,
            label=label,
            jobrunner=jobrunner,
            **kwargs,
        )

        self.settings.freq = False  # turn off freq calc for WBI
