"""
Gaussian Time-Dependent Density Functional Theory (TDDFT) job implementation.

This module provides the GaussianTDDFTJob class for performing
TDDFT calculations using Gaussian. TDDFT calculations compute
electronic excitation energies, oscillator strengths, and excited
state properties for studying UV-Vis spectra and photochemistry.
"""

from typing import Type

from chemsmart.jobs.gaussian.job import GaussianJob
from chemsmart.jobs.gaussian.settings import GaussianTDDFTJobSettings


class GaussianTDDFTJob(GaussianJob):
    """
    Gaussian job class for Time-Dependent DFT calculations.

    Performs TDDFT calculations to compute electronic excitation
    energies and properties. Used for predicting UV-Vis absorption
    spectra, emission spectra, and studying excited state processes.

    TDDFT provides access to vertical excitation energies, oscillator
    strengths, and transition properties for photochemical studies.

    Attributes:
        TYPE (str): Job type identifier ('g16td').
        molecule (Molecule): Molecular structure for TDDFT calculation.
        settings (GaussianTDDFTJobSettings): TDDFT configuration options.
        label (str): Job identifier used for file naming.
        jobrunner (JobRunner): Execution backend that runs the job.
        skip_completed (bool): If True, completed jobs are not rerun.
    """

    TYPE = "g16td"

    def __init__(self, molecule, settings, label, jobrunner=None, **kwargs):
        """
        Initialize a Gaussian TDDFT calculation job.

        Sets up a TDDFT calculation with the specified molecular
        structure and excited state calculation settings.

        Args:
            molecule (Molecule): Molecular structure for TDDFT calculation.
            settings (GaussianTDDFTJobSettings): TDDFT configuration.
            label (str): Job identifier for file naming.
            jobrunner (JobRunner, optional): Job execution handler.
            **kwargs: Additional keyword arguments for parent class.
        """
        super().__init__(
            molecule=molecule,
            settings=settings,
            label=label,
            jobrunner=jobrunner,
            **kwargs,
        )

    @classmethod
    def settings_class(cls) -> Type[GaussianTDDFTJobSettings]:
        """
        Get the settings class used by this TDDFT job type.

        Returns the appropriate settings class for configuring
        TDDFT-specific parameters and calculation options.

        Returns:
            Type[GaussianTDDFTJobSettings]: Settings class for TDDFT jobs.
        """
        return GaussianTDDFTJobSettings
