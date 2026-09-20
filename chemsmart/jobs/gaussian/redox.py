"""Gaussian redox exchange job: Opt (Ox, Red) → Ref Opt → SP → Ref SP."""

from chemsmart.jobs.gaussian.job import GaussianJob
from chemsmart.jobs.gaussian.opt import GaussianOptJob
from chemsmart.jobs.gaussian.settings import GaussianJobSettings
from chemsmart.jobs.gaussian.singlepoint import GaussianSinglePointJob
from chemsmart.jobs.redox import RedoxChainMixin, RedoxJobSettingsMixin


class GaussianRedoxJobSettings(RedoxJobSettingsMixin, GaussianJobSettings):
    """Gaussian settings for exchange redox calculations."""

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        if not self.title:
            self.title = "Gaussian redox calculation job"


class GaussianRedoxJob(RedoxChainMixin, GaussianJob):
    """Gaussian job for the dual-level exchange redox cycle.

    Phases: Opt (Ox, Red) → Ref Opt → SP → Ref SP. The oxidized target
    comes from ``molecule``; the reduced target uses the same geometry
    or ``settings.red_file`` with charge ``ox − n``. The oxidized reference
    comes from ``ref_ox_file`` or the registry; the reduced reference uses
    the same geometry (or ``ref_red_file``) with charge ``ox − n``.
    """

    TYPE = "g16redox"
    _opt_job_class = GaussianOptJob
    _sp_job_class = GaussianSinglePointJob

    def __init__(self, molecule, settings=None, **kwargs):
        if not isinstance(settings, GaussianRedoxJobSettings):
            raise ValueError(
                "Settings must be instance of GaussianRedoxJobSettings "
                f"for {self.__class__.__name__}, but is {settings} instead!"
            )
        super().__init__(molecule=molecule, settings=settings, **kwargs)
