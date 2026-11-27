"""Shared Fukui charge-state chain for Gaussian and ORCA jobs.

``FukuiChainMixin`` builds one Fukui phase with neutral, radical-cation,
and radical-anion population jobs at the same geometry.

Program jobs implement ``_allowed_modes``, ``_population_job_class``,
and ``_population_settings``.
"""

from chemsmart.analysis.fukui import radical_ion_charge_and_multiplicity
from chemsmart.jobs.chain import ChainMixin, JobPhase


class FukuiChainMixin(ChainMixin):
    """Neutral / radical-cation / radical-anion population jobs.

    Subclasses implement ``_allowed_modes``, ``_population_job_class``,
    and ``_population_settings``.
    """

    _supported_modes_label = "modes"

    def __init__(
        self,
        *args,
        mode="mulliken",
        radical_cation_charge=None,
        radical_cation_multiplicity=None,
        radical_anion_charge=None,
        radical_anion_multiplicity=None,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)

        mode = mode.lower()
        allowed = self._allowed_modes()
        if mode not in allowed:
            raise ValueError(
                f"Unknown Fukui mode {mode}. Supported "
                f"{self._supported_modes_label} are: "
                f"{', '.join(allowed)}."
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

        self._create_charge_jobs(self.molecule)
        self.phases = self._build_phases()

    def _allowed_modes(self):
        raise NotImplementedError

    def _population_job_class(self):
        raise NotImplementedError

    def _child_job(self, job_class, molecule, settings, label):
        return job_class(
            molecule=molecule,
            settings=settings,
            label=label,
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )

    def _append_route_parameter(self, settings, extra):
        current = settings.additional_route_parameters
        if current:
            settings.additional_route_parameters = f"{current} {extra}"
        else:
            settings.additional_route_parameters = extra
        return settings

    def _population_settings(self, charge, multiplicity):
        settings = self.settings.copy()
        settings.charge = charge
        settings.multiplicity = multiplicity
        settings.freq = False
        settings.jobtype = "sp"
        return settings

    def _create_charge_jobs(self, molecule):
        """Create neutral / cation / anion population jobs at ``molecule``."""
        charge = self.settings.charge
        multiplicity = self.settings.multiplicity
        job_cls = self._population_job_class()

        self.neutral_job = self._child_job(
            job_cls,
            molecule,
            self._population_settings(charge, multiplicity),
            f"{self.label}_n",
        )

        cation_charge, cation_mult = radical_ion_charge_and_multiplicity(
            charge, multiplicity, +1
        )
        if self.radical_cation_charge is not None:
            cation_charge = self.radical_cation_charge
        if self.radical_cation_multiplicity is not None:
            cation_mult = self.radical_cation_multiplicity
        self.cation_job = self._child_job(
            job_cls,
            molecule,
            self._population_settings(cation_charge, cation_mult),
            f"{self.label}_rc",
        )

        anion_charge, anion_mult = radical_ion_charge_and_multiplicity(
            charge, multiplicity, -1
        )
        if self.radical_anion_charge is not None:
            anion_charge = self.radical_anion_charge
        if self.radical_anion_multiplicity is not None:
            anion_mult = self.radical_anion_multiplicity
        self.anion_job = self._child_job(
            job_cls,
            molecule,
            self._population_settings(anion_charge, anion_mult),
            f"{self.label}_ra",
        )

        self.charge_jobs = [
            self.neutral_job,
            self.cation_job,
            self.anion_job,
        ]
        return self.charge_jobs

    def _build_phases(self):
        return [
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
