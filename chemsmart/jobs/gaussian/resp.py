"""
Gaussian RESP (Restrained Electrostatic Potential) charge calculation job.

This module provides the GaussianRESPJob class for performing
RESP charge calculations using Gaussian. RESP charges are derived
from electrostatic potential fitting and are commonly used in
molecular dynamics simulations and force field development.
"""

import logging

from chemsmart.jobs.gaussian.job import GaussianJob

logger = logging.getLogger(__name__)


class GaussianRESPJob(GaussianJob):
    """
    Gaussian job class for RESP charge calculations.

    Performs RESP (Restrained Electrostatic Potential) charge
    calculations to derive atomic partial charges fitted to
    the molecular electrostatic potential. These charges are
    particularly useful for molecular dynamics simulations.

    RESP charges provide a good balance between accuracy and
    transferability for force field applications.

    Attributes:
        TYPE (str): Job type identifier ('g16resp').
        molecule (Molecule): Molecular structure for RESP calculation.
        settings (GaussianJobSettings): Calculation configuration options.
        label (str): Job identifier used for file naming.
        jobrunner (JobRunner): Execution backend that runs the job.
        skip_completed (bool): If True, completed jobs are not rerun.
    """

    TYPE = "g16resp"

    def __init__(self, molecule, settings, label, jobrunner=None, **kwargs):
        """
        Initialize a Gaussian RESP charge calculation job.

        Sets up a RESP charge calculation with the specified
        molecular structure and calculation settings.

        Args:
            molecule (Molecule): Molecular structure for RESP calculation.
            settings (GaussianJobSettings): Calculation configuration.
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
