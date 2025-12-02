import logging
import os
from typing import Optional, Type

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.job import Job
from chemsmart.jobs.iterate.settings import IterateJobSettings

logger = logging.getLogger(__name__)


class IterateJob(Job):
    """
    Job for generating molecular structures by attaching substituents to skeletons.
    
    This job holds configuration from settings (skeleton_list, substituent_list)
    and delegates all execution logic to IterateJobRunner.
    """

    PROGRAM = "Iterate"
    TYPE = "iterate"

    def __init__(
            self, 
            settings=None,
            jobrunner=None,
            nprocs=1,
            timeout=120,
            outputfile="iterate_out",
            **kwargs,
            ):
        """
        Initialize an IterateJob instance.
        
        Parameters
        ----------
        settings : IterateJobSettings, optional
            Job settings containing skeleton_list, substituent_list, algorithm
        jobrunner : IterateJobRunner, optional
            Job runner instance for executing combinations
        nprocs : int, optional
            Number of processes for parallel execution. Default 1.
        timeout : int, optional
            Timeout in seconds for each worker process. Default 120 (2 minutes).
        outputfile : str, optional
            Output filename (without .xyz extension). Default is 'iterate_out'.
        **kwargs
            Additional keyword arguments
        """
        logger.debug("Initializing IterateJob")
        # IterateJob doesn't use a single molecule, pass None and a default label
        super().__init__(
            molecule=None,
            label="iterate_job",
            jobrunner=jobrunner,
            skip_completed=False,  # IterateJob doesn't support skip_completed
            **kwargs
        )
        
        self.settings = settings if settings is not None else IterateJobSettings()
        self.nprocs = nprocs
        self.timeout = timeout
        self._outputfile = outputfile
    
    @classmethod
    def settings_class(cls) -> Type[IterateJobSettings]:
        """Return the settings class for iterate jobs."""
        return IterateJobSettings
    
    @property
    def outputfile(self) -> str:
        """
        Get the path to the output xyz file.
        
        Returns
        -------
        str
            Path to the output file with .xyz extension
        """
        # Add .xyz extension if not already present
        if self._outputfile.endswith('.xyz'):
            return self._outputfile
        return f"{self._outputfile}.xyz"
    
    def run(self) -> str:
        """
        Run the iterate job by delegating to IterateJobRunner.
        
        The runner handles all execution logic including:
        - Loading molecules from settings
        - Generating combinations
        - Running combinations with multiprocessing
        - Writing results to output xyz file
        
        Returns
        -------
        str
            Path to the output xyz file
        """
        logger.info(f"Running IterateJob with {self.nprocs} process(es)")
        
        # Use jobrunner if available, otherwise create a temporary one
        if self.jobrunner is not None:
            runner = self.jobrunner
        else:
            from chemsmart.jobs.iterate.runner import IterateJobRunner
            runner = IterateJobRunner()
        
        # Delegate all execution to runner
        runner.run(self)
        
        return self.outputfile
