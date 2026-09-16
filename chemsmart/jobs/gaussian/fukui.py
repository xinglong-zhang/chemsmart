"""
Gaussian Fukui calculation job implementation.

Runs population calculations for the neutral, radical cation (N-1 electrons),
and radical anion (N+1 electrons) at the same geometry.

Post-processing is backend-independent via ``chemsmart run fukui``
(``chemsmart.analysis.fukui``).
"""

from chemsmart.analysis.fukui import FUKUI_MODES
from chemsmart.jobs.fukui import FukuiChainMixin
from chemsmart.jobs.gaussian.job import GaussianJob
from chemsmart.jobs.gaussian.settings import GaussianJobSettings
from chemsmart.jobs.gaussian.singlepoint import GaussianSinglePointJob
from chemsmart.jobs.gaussian.wbi import GaussianWBIJob


class GaussianFukuiJob(FukuiChainMixin, GaussianJob):
    """
    Gaussian job class for Fukui charge-state calculations.

    Runs population calculations at the same geometry for the neutral,
    radical cation (N-1 electrons), and radical anion (N+1 electrons).

    Attributes:
        TYPE (str): Job type identifier ('g16fukui').
        molecule (Molecule): Neutral molecular structure.
        settings (GaussianJobSettings): Population-job configuration.
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

    TYPE = "g16fukui"

    def __init__(self, molecule, settings=None, **kwargs):
        if not isinstance(settings, GaussianJobSettings):
            raise ValueError(
                f"Settings must be instance of GaussianJobSettings for "
                f"{self.__class__.__name__}, but is {settings} instead!"
            )
        super().__init__(molecule=molecule, settings=settings, **kwargs)

    def _allowed_modes(self):
        return FUKUI_MODES

    def _population_job_class(self):
        if self.mode == "nbo":
            return GaussianWBIJob
        return GaussianSinglePointJob

    def _population_settings(self, charge, multiplicity):
        settings = super()._population_settings(charge, multiplicity)
        if self.mode == "nbo":
            settings.jobtype = "wbi"
        elif self.mode in ("hirshfeld", "cm5"):
            self._append_route_parameter(settings, "pop=hirshfeld")
        return settings
