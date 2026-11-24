"""Gaussian fragment-guess jobs, for parallel or antiparallel fragment spins."""

from chemsmart.jobs.gaussian.link import GaussianLinkJob
from chemsmart.jobs.gaussian.settings import GaussianFragJobSettings


class GaussianFragJob(GaussianLinkJob):
    """One Gaussian execution containing fragment guess and continuation steps."""

    def __init__(
        self,
        molecule,
        settings,
        label=None,
        jobrunner=None,
        *,
        submission_index=None,
        **kwargs,
    ):
        if not isinstance(settings, GaussianFragJobSettings):
            raise ValueError("Fragment jobs require GaussianFragJobSettings.")
        super().__init__(molecule, settings, label, jobrunner, **kwargs)
        self.submission_index = submission_index
        self.settings.validate_fragments(self.molecule)

    @classmethod
    def settings_class(cls):
        return GaussianFragJobSettings

    def _create_irc_subjob(self, direction):
        """Retain fragment settings when the shared link lifecycle splits IRC."""
        if not self._is_irc_job():
            raise ValueError("This method is only for IRC fragment jobs")
        settings = self.settings.copy()
        settings.jobtype = f"irc{direction}"
        label = self.label
        if settings.direction is None:
            suffix = "_frag" if label.endswith("_irc_frag") else ""
            label = label.removesuffix(suffix)
            if not label.endswith("_irc"):
                label += "_irc"
            label += direction
            if settings.flat_irc:
                label += "_flat"
            label += suffix
        job = GaussianFragJob(
            molecule=self.molecule,
            settings=settings,
            label=label,
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )
        job.folder = self.folder
        return job
