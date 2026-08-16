"""
Gaussian Fukui calculation job implementation.

Runs population calculations for the neutral, radical cation (N-1 electrons),
and radical anion (N+1 electrons) at the same geometry.

Post-processing is backend-independent via ``chemsmart run fukui``
(``chemsmart.analysis.fukui``).
"""

import logging

from chemsmart.analysis.fukui import (
    FUKUI_MODES,
    radical_ion_charge_and_multiplicity,
)
from chemsmart.jobs.chain import ChainMixin, JobPhase
from chemsmart.jobs.gaussian.job import GaussianJob
from chemsmart.jobs.gaussian.settings import GaussianJobSettings
from chemsmart.jobs.gaussian.singlepoint import GaussianSinglePointJob
from chemsmart.jobs.gaussian.wbi import GaussianWBIJob

logger = logging.getLogger(__name__)


class GaussianFukuiJob(ChainMixin, GaussianJob):
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

    def __init__(
        self,
        molecule,
        settings=None,
        label=None,
        jobrunner=None,
        skip_completed=True,
        mode="mulliken",
        radical_cation_charge=None,
        radical_cation_multiplicity=None,
        radical_anion_charge=None,
        radical_anion_multiplicity=None,
        **kwargs,
    ):
        if not isinstance(settings, GaussianJobSettings):
            raise ValueError(
                f"Settings must be instance of GaussianJobSettings for "
                f"{self.__class__.__name__}, but is {settings} instead!"
            )

        mode = mode.lower()
        if mode not in FUKUI_MODES:
            raise ValueError(
                f"Unknown Fukui mode {mode}. Supported modes are: "
                f"{', '.join(FUKUI_MODES)}."
            )

        super().__init__(
            molecule=molecule,
            settings=settings,
            label=label,
            jobrunner=jobrunner,
            skip_completed=skip_completed,
            **kwargs,
        )

        self.mode = mode
        self.radical_cation_charge = radical_cation_charge
        self.radical_cation_multiplicity = radical_cation_multiplicity
        self.radical_anion_charge = radical_anion_charge
        self.radical_anion_multiplicity = radical_anion_multiplicity

        self.neutral_job = None
        self.cation_job = None
        self.anion_job = None
        self.charge_jobs = []

        self._prepare_fukui_jobs()
        self.phases = [
            JobPhase(
                name="Fukui",
                jobs_factory=lambda: self.charge_jobs,
                stop_on_incomplete=True,
                require_complete=True,
                stop_message=(
                    "Fukui charge-state jobs incomplete, "
                    "halting serial execution."
                ),
            )
        ]

    def _prepare_fukui_jobs(self):
        """Prepare charge-state population jobs."""
        if self.settings is None:
            return
        self._create_charge_jobs(self.molecule)

    def _population_job_class(self):
        if self.mode == "nbo":
            return GaussianWBIJob
        return GaussianSinglePointJob

    def _population_settings(self, charge, multiplicity):
        settings = self.settings.copy()
        settings.charge = charge
        settings.multiplicity = multiplicity
        settings.freq = False
        if self.mode == "nbo":
            settings.jobtype = "wbi"
        else:
            settings.jobtype = "sp"
            if self.mode in ("hirshfeld", "cm5"):
                extra = "pop=hirshfeld"
                current = settings.additional_route_parameters
                if current:
                    settings.additional_route_parameters = f"{current} {extra}"
                else:
                    settings.additional_route_parameters = extra
        return settings

    def _create_charge_jobs(self, molecule):
        """Create neutral / cation / anion population jobs at ``molecule``."""
        charge = self.settings.charge
        multiplicity = self.settings.multiplicity
        job_cls = self._population_job_class()

        self.neutral_job = job_cls(
            molecule=molecule,
            settings=self._population_settings(charge, multiplicity),
            label=f"{self.label}_n",
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )

        cation_charge, cation_mult = radical_ion_charge_and_multiplicity(
            charge, multiplicity, +1
        )
        if self.radical_cation_charge is not None:
            cation_charge = self.radical_cation_charge
        if self.radical_cation_multiplicity is not None:
            cation_mult = self.radical_cation_multiplicity
        self.cation_job = job_cls(
            molecule=molecule,
            settings=self._population_settings(cation_charge, cation_mult),
            label=f"{self.label}_rc",
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )

        anion_charge, anion_mult = radical_ion_charge_and_multiplicity(
            charge, multiplicity, -1
        )
        if self.radical_anion_charge is not None:
            anion_charge = self.radical_anion_charge
        if self.radical_anion_multiplicity is not None:
            anion_mult = self.radical_anion_multiplicity
        self.anion_job = job_cls(
            molecule=molecule,
            settings=self._population_settings(anion_charge, anion_mult),
            label=f"{self.label}_ra",
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )

        self.charge_jobs = [
            self.neutral_job,
            self.cation_job,
            self.anion_job,
        ]
        return self.charge_jobs

    @property
    def neutral_output(self):
        """Parsed output for the neutral charge-state job."""
        if self.neutral_job is None:
            return None
        return self.neutral_job._output()

    @property
    def cation_output(self):
        """Parsed output for the radical-cation job."""
        if self.cation_job is None:
            return None
        return self.cation_job._output()

    @property
    def anion_output(self):
        """Parsed output for the radical-anion job."""
        if self.anion_job is None:
            return None
        return self.anion_job._output()
