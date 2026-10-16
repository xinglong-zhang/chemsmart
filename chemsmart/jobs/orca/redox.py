"""ORCA redox exchange job: Opt (Ox, Red) → Ref Opt → SP → Ref SP."""

from chemsmart.jobs.orca.job import ORCAJob
from chemsmart.jobs.orca.opt import ORCAOptJob
from chemsmart.jobs.orca.settings import ORCAJobSettings
from chemsmart.jobs.orca.singlepoint import ORCASinglePointJob
from chemsmart.jobs.redox import RedoxChainMixin, RedoxJobSettingsMixin


class ORCARedoxJobSettings(RedoxJobSettingsMixin, ORCAJobSettings):
    """ORCA settings for exchange redox calculations."""

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        if not self.title:
            self.title = "ORCA redox calculation job"


class ORCARedoxJob(RedoxChainMixin, ORCAJob):
    """ORCA job for the dual-level exchange redox cycle.

    Phases: Opt (Ox, Red) → Ref Opt → SP → Ref SP. The oxidized target
    comes from ``molecule``; the reduced target uses the same geometry
    or ``settings.red_file`` with charge ``ox − n``. The oxidized reference
    comes from ``ref_ox_file`` or the registry; the reduced reference uses
    the same geometry (or ``ref_red_file``) with charge ``ox − n``.
    """

    TYPE = "orcaredox"
    _opt_job_class = ORCAOptJob
    _sp_job_class = ORCASinglePointJob

    def __init__(self, molecule, settings=None, **kwargs):
        if not isinstance(settings, ORCARedoxJobSettings):
            raise ValueError(
                "Settings must be instance of ORCARedoxJobSettings "
                f"for {self.__class__.__name__}, but is {settings} instead!"
            )
        super().__init__(molecule=molecule, settings=settings, **kwargs)
