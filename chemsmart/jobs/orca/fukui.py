"""
ORCA Fukui calculation job implementation.

Runs population calculations for the neutral, radical cation (N-1 electrons),
and radical anion (N+1 electrons) at the same geometry.

Post-processing is backend-independent via ``chemsmart run fukui``
(``chemsmart.analysis.fukui``).
"""

from chemsmart.analysis.fukui import ORCA_FUKUI_MODES
from chemsmart.jobs.fukui import FukuiChainMixin
from chemsmart.jobs.orca.job import ORCAJob
from chemsmart.jobs.orca.settings import ORCAJobSettings
from chemsmart.jobs.orca.singlepoint import ORCASinglePointJob


class ORCAFukuiJob(FukuiChainMixin, ORCAJob):
    """
    ORCA job class for Fukui charge-state calculations.

    Runs population calculations at the same geometry for the neutral,
    radical cation (N-1 electrons), and radical anion (N+1 electrons).

    Attributes:
        TYPE (str): Job type identifier ('orcafukui').
        molecule (Molecule): Neutral molecular structure.
        settings (ORCAJobSettings): Population-job configuration.
        label (str): Base job identifier used for file naming.
        jobrunner (JobRunner): Execution backend that runs the jobs.
        skip_completed (bool): If True, completed jobs are not rerun.
        mode (str): Charge partitioning mode for later analysis.
        radical_cation_charge (int, optional): Override radical-cation charge.
        radical_cation_multiplicity (int, optional): Override radical-cation
            multiplicity.
        radical_anion_charge (int, optional): Override radical-anion charge.
        radical_anion_multiplicity (int, optional): Override radical-anion
            multiplicity.
    """

    TYPE = "orcafukui"
    _supported_modes_label = "ORCA modes"

    def __init__(self, molecule, settings=None, **kwargs):
        if not isinstance(settings, ORCAJobSettings):
            raise ValueError(
                f"Settings must be instance of ORCAJobSettings for "
                f"{self.__class__.__name__}, but is {settings} instead!"
            )
        super().__init__(molecule=molecule, settings=settings, **kwargs)

    def _allowed_modes(self):
        return ORCA_FUKUI_MODES

    def _population_job_class(self):
        return ORCASinglePointJob

    def _population_settings(self, charge, multiplicity):
        settings = super()._population_settings(charge, multiplicity)
        if self.mode == "hirshfeld":
            self._append_route_parameter(settings, "Hirshfeld")
        return settings
