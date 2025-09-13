"""
Gaussian Intrinsic Reaction Coordinate (IRC) job implementation.

This module provides the GaussianIRCJob class for performing
IRC calculations using Gaussian. IRC calculations trace reaction
paths from transition states to reactants and products by following
the steepest descent path on the potential energy surface.

The implementation manages both forward and reverse IRC calculations
automatically, providing complete reaction path characterization.
"""

import logging
from copy import deepcopy
from typing import Type

from chemsmart.jobs.gaussian.job import GaussianGeneralJob, GaussianJob
from chemsmart.jobs.gaussian.settings import GaussianIRCJobSettings

logger = logging.getLogger(__name__)


class GaussianIRCJob(GaussianJob):
    """
    Gaussian job class for Intrinsic Reaction Coordinate calculations.
    
    Performs IRC calculations to trace reaction paths from transition
    states to reactants and products. Automatically manages both forward
    and reverse IRC calculations to provide complete reaction path
    characterization.
    
    IRC calculations follow the steepest descent path on the potential
    energy surface in mass-weighted coordinates, providing the minimum
    energy pathway connecting transition state to stable products.
    
    Attributes:
        TYPE (str): Job type identifier ('g16irc').
        molecule (Molecule): Transition state structure used as the IRC start.
        settings (GaussianIRCJobSettings): IRC-specific configuration options.
        label (str): Job identifier used for file naming.
        jobrunner (JobRunner): Execution backend that runs the job.
        skip_completed (bool): If True, completed jobs are not rerun.
    """
    TYPE = "g16irc"

    def __init__(
        self, molecule, settings=None, label=None, jobrunner=None, **kwargs
    ):
        """
        Initialize a Gaussian IRC calculation job.
        
        Sets up an IRC calculation with the specified transition state
        structure and calculation settings. Automatically disables
        frequency calculations for the IRC steps.
        
        Args:
            molecule (Molecule): Transition state structure for IRC.
            settings (GaussianIRCJobSettings, optional): IRC configuration.
            label (str, optional): Job identifier for file naming.
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
        self.settings = settings
        self.settings.freq = False  # turn off freq calc for IRC

    @classmethod
    def settings_class(cls) -> Type[GaussianIRCJobSettings]:
        """
        Get the settings class used by this IRC job type.
        
        Returns the appropriate settings class for configuring
        IRC-specific parameters and calculation options.
        
        Returns:
            Type[GaussianIRCJobSettings]: Settings class for IRC jobs.
        """
        return GaussianIRCJobSettings

    def _ircf_job(self):
        """
        Create forward IRC job configuration.
        
        Sets up the forward IRC calculation that follows the reaction
        coordinate from transition state toward products.
        
        Returns:
            GaussianGeneralJob: Configured forward IRC job.
        """
        # create IRCf job:
        ircf_settings = deepcopy(self.settings)
        label = self.label
        if self.label is not None:
            label = self.label + "f"
            if self.settings.flat_irc:
                label += "_flat"

        ircf_settings.job_type = "ircf"
        return GaussianGeneralJob(
            molecule=self.molecule,
            settings=ircf_settings,
            label=label,
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )

    def _ircr_job(self):
        """
        Create reverse IRC job configuration.
        
        Sets up the reverse IRC calculation that follows the reaction
        coordinate from transition state toward reactants.
        
        Returns:
            GaussianGeneralJob: Configured reverse IRC job.
        """
        # create IRCr job:
        ircr_settings = deepcopy(self.settings)
        label = self.label
        if self.label is not None:
            label = self.label + "r"
            if self.settings.flat_irc:
                label += "_flat"
        ircr_settings.job_type = "ircr"
        return GaussianGeneralJob(
            molecule=self.molecule,
            settings=ircr_settings,
            label=label,
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )

    def _run_forward(self):
        """
        Execute the forward IRC calculation.
        
        Runs the forward IRC job that traces the reaction path
        from transition state toward products.
        """
        logger.debug(
            f"Running forward IRC job: {self._ircf_job().settings.job_type}"
        )
        self._ircf_job().run()

    def _run_reverse(self):
        """
        Execute the reverse IRC calculation.
        
        Runs the reverse IRC job that traces the reaction path
        from transition state toward reactants.
        """
        logger.debug(
            f"Running reverse IRC job: {self._ircr_job().settings.job_type}"
        )
        self._ircr_job().run()

    def _run(self, **kwargs):
        """
        Execute both forward and reverse IRC calculations.
        
        Orchestrates the complete IRC calculation by running
        both forward and reverse IRC jobs sequentially.
        
        Args:
            **kwargs: Additional keyword arguments (currently unused).
        """
        self._run_forward()
        self._run_reverse()

    def _job_is_complete(self):
        """
        Check if the complete IRC calculation is finished.
        
        Determines completion status by verifying that both forward
        and reverse IRC calculations have completed successfully.
        
        Returns:
            bool: True if both IRC directions are complete, False otherwise.
        """
        return (
            self._run_forward_is_complete() and self._run_reverse_is_complete()
        )

    def _run_forward_is_complete(self):
        """
        Check if the forward IRC calculation is complete.
        
        Returns:
            bool: True if forward IRC is complete, False otherwise.
        """
        return self._ircf_job().is_complete()

    def _run_reverse_is_complete(self):
        """
        Check if the reverse IRC calculation is complete.
        
        Returns:
            bool: True if reverse IRC is complete, False otherwise.
        """
        return self._ircr_job().is_complete()

    def backup_files(self, backup_chk=False):
        """
        Create backup copies of IRC input and output files.
        
        Backs up all files from both forward and reverse IRC calculations
        to preserve important calculation data.
        
        Args:
            backup_chk (bool): Whether to backup checkpoint files.
        """
        self.backup_file(self._ircf_job().inputfile)
        self.backup_file(self._ircr_job().inputfile)
        self.backup_file(self._ircf_job().outputfile)
        self.backup_file(self._ircr_job().outputfile)
        if backup_chk:
            self.backup_file(self._ircf_job().chkfile)
            self.backup_file(self._ircr_job().chkfile)
