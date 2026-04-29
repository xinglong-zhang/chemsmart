"""
GROMACS job definitions.
"""

import logging
import os

from chemsmart.jobs.job import Job

logger = logging.getLogger(__name__)


class GromacsJob(Job):
    """
    Base class for GROMACS jobs.
    """

    PROGRAM = "gromacs"
    TYPE = "gmx"

    def __init__(self, molecule=None, settings=None, label=None, jobrunner=None, **kwargs):
        super().__init__(
            molecule=molecule,
            label=label,
            jobrunner=jobrunner,
            **kwargs,
        )
        self.molecule = molecule
        self.settings = settings

        if label is None:
            label = "gromacs_job"
        self.label = label

    @property
    def inputfile(self):
        """
        First version: use .gro as the main structure input placeholder.
        Later this may be generalized.
        """
        return os.path.join(self.folder, f"{self.label}.gro")

    @property
    def outputfile(self):
        return os.path.join(self.folder, f"{self.label}.out")

    @property
    def errfile(self):
        return os.path.join(self.folder, f"{self.label}.err")

    def _run(self, **kwargs):
        logger.info(f"Running GromacsJob {self} with jobrunner {self.jobrunner}")
        self.jobrunner.run(self, **kwargs)


class GromacsEMJob(GromacsJob):
    """
    Energy minimization job for GROMACS.
    """

    TYPE = "gmxem"

    def __init__(self, molecule=None, settings=None, label=None, jobrunner=None, **kwargs):
        super().__init__(
            molecule=molecule,
            settings=settings,
            label=label,
            jobrunner=jobrunner,
            **kwargs,
        )
