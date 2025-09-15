"""
Gaussian Non-Covalent Interaction (NCI) job implementation.

This module provides the GaussianNCIJob class for running Gaussian
job that creates the .wfn file required for NCIPLOT program to 
generate the associated dens.cube and grad.cube files. 
NCI analysis identifies and visualizes non-covalent interactions 
in molecular systems including hydrogen bonds, van der Waals 
interactions, and steric clashes.
"""

from chemsmart.jobs.gaussian.job import GaussianJob
from chemsmart.jobs.runner import JobRunner


class GaussianNCIJob(GaussianJob):
    """
    Gaussian job class for Non-Covalent Interaction (NCI) cube files
    generation.
    
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
        
        Sets up an NCI .wfn file calculation with the specified
        molecular structure. 
        
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
