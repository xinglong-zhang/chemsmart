"""
Gaussian geometry optimization job implementation.

This module provides the GaussianOptJob class for performing
molecular geometry optimization calculations using Gaussian.
Handles optimization of molecular structures to find energy
minima or transition states.
"""

from chemsmart.jobs.gaussian.job import GaussianJob


class GaussianOptJob(GaussianJob):
    """
    Gaussian job class for geometry optimization calculations.
    
    Performs molecular geometry optimization to find stationary
    points on the potential energy surface. Can be used for
    finding energy minima, transition states, or other critical
    points depending on the optimization settings.
    
    Attributes:
        TYPE (str): Job type identifier ('g16opt').
        molecule (Molecule): Molecular structure to optimize.
        settings (GaussianJobSettings): Calculation configuration options.
        label (str): Job identifier used for file naming.
        jobrunner (JobRunner): Execution backend that runs the job.
        skip_completed (bool): If True, completed jobs are not rerun.
    """
    TYPE = "g16opt"

    def __init__(self, molecule, settings, label, jobrunner=None, **kwargs):
        """
        Initialize a Gaussian optimization job.
        
        Sets up a geometry optimization calculation with the
        specified molecular structure and calculation settings.
        
        Args:
            molecule (Molecule): Molecular structure to optimize.
            settings (GaussianJobSettings): Calculation configuration.
            label (str): Job identifier for file naming.
            jobrunner (JobRunner, optional): Job execution handler.
            **kwargs: Additional keyword arguments for parent class.
            
        Raises:
            ValueError: If `settings` is not a GaussianJobSettings instance
                or `molecule` is not a Molecule (validated by base class).
        """
        super().__init__(
            molecule=molecule,
            settings=settings,
            label=label,
            jobrunner=jobrunner,
            **kwargs,
        )
