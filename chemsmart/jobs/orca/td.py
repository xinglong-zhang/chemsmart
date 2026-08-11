"""ORCA fixed-geometry CIS/TD-DFT job."""

from typing import Type

from chemsmart.jobs.orca.job import ORCAJob
from chemsmart.jobs.orca.settings import ORCAJobSettings


class ORCATDDFTJob(ORCAJob):
    """Electronic excitation calculation at the supplied geometry."""

    TYPE = "orcatd"

    @classmethod
    def settings_class(cls) -> Type[ORCAJobSettings]:
        return ORCAJobSettings
