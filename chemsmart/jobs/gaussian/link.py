"""
Gaussian link job implementation for multi-step calculations.

This module provides the GaussianLinkJob class for performing
multi-step Gaussian calculations using link directives. Link jobs
allow for sequential calculations with different settings or
methods in a single input file, useful for complex workflows.
"""

import logging
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
        """Check if this is an IRC link job."""
        return self.settings.jobtype == "irc"

    def _create_irc_subjob(self, direction):
        """
        Create IRC subjob for specified direction.

        Args:
            direction (str): Either 'f' for forward or 'r' for reverse

        Returns:
            GaussianLinkJob: Configured IRC subjob
        """
        if not self._is_irc_job():
            raise ValueError("This method is only for IRC link jobs")

        # Create IRC subjob
        settings_dict = self.settings.__dict__.copy()
        settings_dict["jobtype"] = f"irc{direction}"
        irc_settings = GaussianLinkJobSettings(**settings_dict)

        label = self.label
        if self.label is not None:
            if "irc_link" in label:
                new_irc_part = f"irc{direction}"
                if (
                    hasattr(self.settings, "flat_irc")
                    and self.settings.flat_irc
                ):
                    new_irc_part += "_flat"
                label = label.replace("irc_link", f"{new_irc_part}_link")
            else:
                label = f"{self.label}_{direction}"
                if (
                    hasattr(self.settings, "flat_irc")
                    and self.settings.flat_irc
                ):
                    label += "_flat"

        return GaussianLinkJob(
            molecule=self.molecule,
            settings=irc_settings,
            label=label,
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )

    def _ircf_link_job(self):
        """Create forward IRC link job."""
        return self._create_irc_subjob("f")

    def _ircr_link_job(self):
        """Create reverse IRC link job."""
        return self._create_irc_subjob("r")

    def _get_irc_jobs(self):
        """
        Get required IRC jobs based on direction setting.

        Returns:
            list: List of IRC jobs to run/check
        """
        if not self._is_irc_job():
            return []

        if (
            hasattr(self.settings, "direction")
            and self.settings.direction == "forward"
        ):
            return [self._ircf_link_job()]
        elif (
            hasattr(self.settings, "direction")
            and self.settings.direction == "reverse"
        ):
            return [self._ircr_link_job()]
        else:
            return [self._ircf_link_job(), self._ircr_link_job()]

    def _run(self, **kwargs):
        """Execute the link calculation."""
        if self._is_irc_job():
            irc_jobs = self._get_irc_jobs()
            direction_names = {1: "forward", 2: "both forward and reverse"}
            logger.info(
                f"Running {direction_names.get(len(irc_jobs), 'reverse')} IRC link calculation"
            )

            # Check if jobs should be run in serial based on jobrunner flag
            if self.jobrunner and self.jobrunner.run_in_serial:
                logger.info("Running IRC jobs in serial mode (one after another)")
                for job in irc_jobs:
                    logger.debug(f"Running IRC job: {job.settings.jobtype}")
                    job.run()
            else:
                logger.info("Running IRC jobs using default behavior")
                for job in irc_jobs:
                    logger.debug(f"Running IRC job: {job.settings.jobtype}")
                    job.run()
        else:
            super()._run(**kwargs)

    def _job_is_complete(self):
        """Check if the link calculation is complete."""
        if self._is_irc_job():
            irc_jobs = self._get_irc_jobs()
            return all(job._job_is_complete() for job in irc_jobs)
        else:
            return super()._job_is_complete()

    def backup_files(self, backup_chk=False):
        """Create backup copies of link files."""
        if self._is_irc_job():
            irc_jobs = self._get_irc_jobs()
            for job in irc_jobs:
                self.backup_file(job.inputfile)
                self.backup_file(job.outputfile)
                if backup_chk:
                    self.backup_file(job.chkfile)
        else:
            super()._backup_files(backup_chk)
