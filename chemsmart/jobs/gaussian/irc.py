import logging
from copy import deepcopy
from typing import Type

from chemsmart.jobs.gaussian.job import GaussianGeneralJob, GaussianJob
from chemsmart.jobs.gaussian.settings import GaussianIRCJobSettings

logger = logging.getLogger(__name__)


class GaussianIRCJob(GaussianJob):
    TYPE = "g16irc"

    def __init__(
        self, molecule, settings=None, label=None, jobrunner=None, **kwargs
    ):
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
        return GaussianIRCJobSettings

    def _ircf_job(self):
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
        logger.debug(
            f"Running forward IRC job: {self._ircf_job().settings.job_type}"
        )
        self._ircf_job().run()

    def _run_reverse(self):
        logger.debug(
            f"Running reverse IRC job: {self._ircr_job().settings.job_type}"
        )
        self._ircr_job().run()

    def _run(self, **kwargs):
        self._run_forward()
        self._run_reverse()

    def _job_is_complete(self):
        """private method for checking that the IRC job is complete.
        Job is complete when both ircf and ircr jobs are complete#"""
        return (
            self._run_forward_is_complete() and self._run_reverse_is_complete()
        )

    def _run_forward_is_complete(self):
        return self._ircf_job().is_complete()

    def _run_reverse_is_complete(self):
        return self._ircr_job().is_complete()

    def backup_files(self, backup_chk=False):
        self.backup_file(self._ircf_job().inputfile)
        self.backup_file(self._ircr_job().inputfile)
        self.backup_file(self._ircf_job().outputfile)
        self.backup_file(self._ircr_job().outputfile)
        if backup_chk:
            self.backup_file(self._ircf_job().chkfile)
            self.backup_file(self._ircr_job().chkfile)
