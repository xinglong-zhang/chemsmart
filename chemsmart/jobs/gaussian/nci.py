"""
Gaussian Non-Covalent Interaction (NCI) analysis job implementation.

This module provides the GaussianNCIJob class for performing
NCI analysis calculations using Gaussian. NCI analysis identifies
and visualizes non-covalent interactions in molecular systems
including hydrogen bonds, van der Waals interactions, and steric clashes.
"""

from chemsmart.jobs.gaussian.job import GaussianJob
from chemsmart.jobs.runner import JobRunner


class GaussianNCIJob(GaussianJob):
    """
    Gaussian job class for Non-Covalent Interaction analysis.
    
    Performs NCI analysis to identify and characterize non-covalent
    interactions in molecular systems. The analysis provides insights
    into hydrogen bonding, van der Waals interactions, and steric
    effects through reduced density gradient calculations.
    
    NCI analysis generates data for visualization of interaction
    regions and strength assessment of non-covalent contacts.
    
    Attributes:
        TYPE (str): Job type identifier ('g16nci').
        molecule (Molecule): Molecular structure for NCI analysis.
        settings (GaussianJobSettings): Calculation configuration options.
        label (str): Job identifier used for file naming.
        jobrunner (JobRunner): Execution backend that runs the job.
        skip_completed (bool): If True, completed jobs are not rerun.
    """
    TYPE = "g16nci"

    def __init__(self, molecule, settings, label, jobrunner=None, **kwargs):
        """
        Initialize a Gaussian NCI analysis job.
        
        Sets up an NCI analysis calculation with the specified
        molecular structure. Automatically disables frequency
        calculations as they are not needed for NCI analysis.
        
        Args:
            molecule (Molecule): Molecular structure for NCI analysis.
            settings (GaussianJobSettings): Calculation configuration.
            label (str): Job identifier for file naming.
            jobrunner (JobRunner, optional): Job execution handler.
            **kwargs: Additional keyword arguments for parent class
                (e.g., `skip_completed`).
        """
        super().__init__(
            molecule=molecule,
            settings=settings,
            label=label,
            jobrunner=JobRunner,
            **kwargs,
        )

        self.settings.freq = False  # turn off freq calc for NCI
