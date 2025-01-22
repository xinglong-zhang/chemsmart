from typing import Type

from chemsmart.jobs.gaussian.job import GaussianJob
from chemsmart.jobs.gaussian.settings import GaussianLinkJobSettings


class GaussianLinkJob(GaussianJob):
    TYPE = "g16link"

    def __init__(self, molecule, settings, label, **kwargs):
        super().__init__(
            molecule=molecule, settings=settings, label=label, **kwargs
        )

    @classmethod
    def settings_class(cls) -> Type[GaussianLinkJobSettings]:
        return GaussianLinkJobSettings
