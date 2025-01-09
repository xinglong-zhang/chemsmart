from copy import deepcopy

from pyatoms.jobs.gaussian.job import GaussianGeneralJob, GaussianJob
from pyatoms.jobs.gaussian.settings import GaussianIRCJobSettings


class GaussianIRCJob(GaussianJob):
    TYPE = "g16irc"
    _SETTINGS_CLS = GaussianIRCJobSettings

    def __init__(self, folder, atoms, settings, **kwargs):
        super().__init__(
            folder=folder, atoms=atoms, settings=settings, **kwargs
        )
        self.settings = settings

    def _ircf_job(self):
        # create IRCf job:
        ircf_settings = deepcopy(self.settings)
        if self.label is not None:
            label = self.label + "f"
            if self.settings.flat_irc:
                label += "_flat"

        ircf_settings.job_type = "ircf"
        return GaussianGeneralJob(
            folder=self.folder,
            atoms=self.atoms,
            settings=ircf_settings,
            label=label,
        )

    def _ircr_job(self):
        # create IRCr job:
        ircr_settings = deepcopy(self.settings)
        if self.label is not None:
            label = self.label + "r"
            if self.settings.flat_irc:
                label += "_flat"
        ircr_settings.job_type = "ircr"
        return GaussianGeneralJob(
            folder=self.folder,
            atoms=self.atoms,
            settings=ircr_settings,
            label=label,
        )

    def _run_forward(self, jobrunner, queue_manager=None):
        self._ircf_job().run(jobrunner=jobrunner, queue_manager=queue_manager)

    def _run_reverse(self, jobrunner, queue_manager=None):
        self._ircr_job().run(jobrunner=jobrunner, queue_manager=queue_manager)

    def _run(self, jobrunner, queue_manager=None, **kwargs):
        self._run_forward(jobrunner=jobrunner, queue_manager=queue_manager)
        self._run_reverse(jobrunner=jobrunner, queue_manager=queue_manager)

    def _job_is_complete(self):
        # private method for checking that the IRC job is complete
        # job is complete when both ircf and ircr jobs are complete
        return (
            self._run_forward_is_complete() and self._run_reverse_is_complete()
        )

    def _run_forward_is_complete(self):
        return self._ircf_job().is_complete()

    def _run_reverse_is_complete(self):
        return self._ircr_job().is_complete()

    def backup_files(self, backup_chk=False):
        folder = self._create_backup_folder_name()
        self.backup_file(self._ircf_job().inputfile, folder=folder)
        self.backup_file(self._ircr_job().inputfile, folder=folder)
        self.backup_file(self._ircf_job().outputfile, folder=folder)
        self.backup_file(self._ircr_job().outputfile, folder=folder)
        if backup_chk:
            self.backup_file(self._ircf_job().chkfile, folder=folder)
            self.backup_file(self._ircr_job().chkfile, folder=folder)
