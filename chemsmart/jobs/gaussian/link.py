"""
Gaussian link job implementation for multi-step calculations.

This module provides the GaussianLinkJob class for performing
multi-step Gaussian calculations using link directives. Link jobs
allow for sequential calculations with different settings or
methods in a single input file, useful for complex workflows.
"""

import logging
from copy import deepcopy
from typing import Type

from chemsmart.jobs.gaussian.job import GaussianJob
from chemsmart.jobs.gaussian.settings import GaussianLinkJobSettings

logger = logging.getLogger(__name__)


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

        For IRC jobs, automatically handles the creation of separate
        forward and reverse calculations.

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

    def _is_irc_job(self):
        """
        Check if this is an IRC link job.

        Returns:
            bool: True if job type is IRC, False otherwise.
        """
        return self.settings.job_type == "irc"

    def _ircf_link_job(self):
        """
        Create forward IRC link job configuration.

        Sets up the forward IRC link calculation that follows the reaction
        coordinate from transition state toward products.

        Returns:
            GaussianLinkJob: Configured forward IRC link job.
        """
        if not self._is_irc_job():
            raise ValueError("This method is only for IRC link jobs")

        # Create IRCf link job
        ircf_settings = deepcopy(self.settings)
        label = self.label
        if self.label is not None:
            # Append ircf_link to the original label
            label = f"{self.label}_f"
            if hasattr(self.settings, "flat_irc") and self.settings.flat_irc:
                label += "_flat"

        # Set job type for forward IRC
        ircf_settings.job_type = "ircf"

        return GaussianLinkJob(
            molecule=self.molecule,
            settings=ircf_settings,
            label=label,
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )

    def _ircr_link_job(self):
        """
        Create reverse IRC link job configuration.

        Sets up the reverse IRC link calculation that follows the reaction
        coordinate from transition state toward reactants.

        Returns:
            GaussianLinkJob: Configured reverse IRC link job.
        """
        if not self._is_irc_job():
            raise ValueError("This method is only for IRC link jobs")

        # Create IRCr link job
        ircr_settings = deepcopy(self.settings)
        label = self.label
        if self.label is not None:
            # Append ircr_link to the original label
            label = f"{self.label}_r"
            if hasattr(self.settings, "flat_irc") and self.settings.flat_irc:
                label += "_flat"

        # Set job type for reverse IRC
        ircr_settings.job_type = "ircr"

        return GaussianLinkJob(
            molecule=self.molecule,
            settings=ircr_settings,
            label=label,
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )

    def _run_forward_irc(self):
        """
        Execute the forward IRC link calculation.
        """
        logger.debug(
            f"Running forward IRC link job: {self._ircf_link_job().settings.job_type}"
        )
        self._ircf_link_job().run()

    def _run_reverse_irc(self):
        """
        Execute the reverse IRC link calculation.
        """
        logger.debug(
            f"Running reverse IRC link job: {self._ircr_link_job().settings.job_type}"
        )
        self._ircr_link_job().run()

    def _run(self, **kwargs):
        """
        Execute the link calculation.

        For IRC jobs, runs both forward and reverse calculations.
        For other jobs, runs the standard link calculation.

        Args:
            **kwargs: Additional keyword arguments.
        """
        if self._is_irc_job():
            # For IRC jobs, run both directions
            self._run_forward_irc()
            self._run_reverse_irc()
        else:
            # For regular link jobs, use parent implementation
            super()._run(**kwargs)

    def _job_is_complete(self):
        """
        Check if the link calculation is complete.

        For IRC jobs, both forward and reverse must be complete.
        For other jobs, uses standard completion check.

        Returns:
            bool: True if job is complete, False otherwise.
        """
        if self._is_irc_job():
            return (
                self._ircf_link_job().is_complete()
                and self._ircr_link_job().is_complete()
            )
        else:
            return super()._job_is_complete()

    def backup_files(self, backup_chk=False):
        """
        Create backup copies of link input and output files.

        For IRC jobs, backs up files from both forward and reverse calculations.

        Args:
            backup_chk (bool): Whether to backup checkpoint files.
        """
        if self._is_irc_job():
            # Backup files for both IRC directions
            self.backup_file(self._ircf_link_job().inputfile)
            self.backup_file(self._ircr_link_job().inputfile)
            self.backup_file(self._ircf_link_job().outputfile)
            self.backup_file(self._ircr_link_job().outputfile)
            if backup_chk:
                self.backup_file(self._ircf_link_job().chkfile)
                self.backup_file(self._ircr_link_job().chkfile)
        else:
            # Use parent backup method for regular link jobs
            super().backup_files(backup_chk)
