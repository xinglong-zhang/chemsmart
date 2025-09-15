"""
Gaussian link job implementation for multi-step calculations.

This module provides the GaussianLinkJob class for performing
multi-step Gaussian calculations using link directives. Link jobs
allow for sequential calculations with different settings or
methods in a single input file, useful for complex workflows.
"""

from typing import Type

from chemsmart.jobs.gaussian.job import GaussianJob
from chemsmart.jobs.gaussian.settings import GaussianLinkJobSettings


class GaussianLinkJob(GaussianJob):
    """
    Gaussian job class for multi-step link calculations.
    
    Performs multi-step calculations using Gaussian link directives
    to chain multiple calculation steps together. Useful for complex
    workflows like optimization followed by frequency analysis,
    or method validation studies.
    
    Link jobs allow sharing of data between calculation steps
    and can optimize computational efficiency for related calculations.
    
    Attributes:
        TYPE (str): Job type identifier ('g16link').
        molecule (Molecule): Molecular structure used across link steps.
        settings (GaussianLinkJobSettings): Configuration for multi-step
            link calculations and their sequence.
        label (str): Job identifier used for file naming.
        jobrunner (JobRunner): Execution backend that runs the job.
        skip_completed (bool): If True, completed jobs are not rerun.
    """
    TYPE = "g16link"

    def __init__(self, molecule, settings, label, jobrunner=None, **kwargs):
        """
        Initialize a Gaussian link calculation job.
        
        Sets up a multi-step link calculation with the specified
        molecular structure and sequential calculation settings.
        
        Args:
            molecule (Molecule): Molecular structure for calculations.
            settings (GaussianLinkJobSettings): Link job configuration.
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
    def settings_class(cls) -> Type[GaussianLinkJobSettings]:
        """
        Get the settings class used by this link job type.
        
        Returns the appropriate settings class for configuring
        multi-step link job parameters and calculation sequences.
        
        Returns:
            Type[GaussianLinkJobSettings]: Settings class for link jobs.
        """
        return GaussianLinkJobSettings
