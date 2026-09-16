"""Shared redox exchange chain for Gaussian and ORCA jobs.

``RedoxChainMixin`` builds Opt (Ox, Red) → Ref Opt → SP → Ref SP.
Program jobs set ``_opt_job_class`` and ``_sp_job_class``.
"""

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.chain import ChainMixin, JobPhase


def reduced_charge_and_multiplicity(ox_charge, ox_multiplicity, n_electrons):
    """Return reduced-form charge and multiplicity after adding ``n`` electrons.

    Even ``n_electrons`` keeps multiplicity. Odd ``n_electrons`` turns a
    singlet into a doublet and otherwise decreases multiplicity by one.
    """
    red_charge = ox_charge - n_electrons
    if n_electrons % 2 == 0:
        red_multiplicity = ox_multiplicity
    elif ox_multiplicity == 1:
        red_multiplicity = 2
    else:
        red_multiplicity = max(ox_multiplicity - 1, 1)
    return red_charge, red_multiplicity


def _molecule_with_charge(molecule, charge, multiplicity):
    mol = molecule.copy()
    mol.charge = charge
    mol.multiplicity = multiplicity
    return mol


def _load_molecule(filepath, charge, multiplicity):
    mol = Molecule.from_filepath(filepath)
    return _molecule_with_charge(mol, charge, multiplicity)


class RedoxJobSettingsMixin:
    """Redox-specific fields mixed into Gaussian/ORCA job settings."""

    def __init__(
        self,
        n_electrons=None,
        red_file=None,
        red_charge=None,
        red_multiplicity=None,
        reference=None,
        ref_ox_file=None,
        ref_red_file=None,
        ref_ox_charge=None,
        ref_ox_multiplicity=None,
        ref_red_charge=None,
        ref_red_multiplicity=None,
        temperature=298.15,
        concentration=1.0,
        pressure=1.0,
        cutoff_entropy_grimme=100.0,
        cutoff_enthalpy=100.0,
        entropy_method=None,
        energy_units="hartree",
        **kwargs,
    ):
        from chemsmart.analysis.redox import resolve_redox_reference

        super().__init__(**kwargs)
        self.reference = resolve_redox_reference(reference)
        if n_electrons is None:
            n_electrons = self.reference.n_electrons
        elif n_electrons != self.reference.n_electrons:
            raise ValueError(
                f"n_electrons ({n_electrons}) must match the reference "
                f"couple {self.reference.name!r} "
                f"(n={self.reference.n_electrons})."
            )
        if n_electrons < 1:
            raise ValueError("n_electrons must be a positive integer.")
        self.n_electrons = n_electrons
        self.red_file = red_file
        self.red_charge = red_charge
        self.red_multiplicity = red_multiplicity
        self.ref_ox_file = ref_ox_file
        self.ref_red_file = ref_red_file
        self.ref_ox_charge = ref_ox_charge
        self.ref_ox_multiplicity = ref_ox_multiplicity
        self.ref_red_charge = ref_red_charge
        self.ref_red_multiplicity = ref_red_multiplicity
        self.temperature = temperature
        self.concentration = concentration
        self.pressure = pressure
        self.cutoff_entropy_grimme = cutoff_entropy_grimme
        self.cutoff_enthalpy = cutoff_enthalpy
        self.entropy_method = entropy_method
        self.energy_units = energy_units

    @property
    def resolved_ref_ox_file(self):
        if self.ref_ox_file is not None:
            return self.ref_ox_file
        return self.reference.ox_file

    @property
    def resolved_ref_red_file(self):
        if self.ref_red_file is not None:
            return self.ref_red_file
        return self.reference.red_file

    @property
    def has_reference_files(self):
        return self.resolved_ref_ox_file is not None

    def ox_charge_and_multiplicity(self, molecule):
        """Charge and multiplicity of the oxidized target."""
        if self.charge is not None:
            charge = self.charge
        elif molecule.charge is not None:
            charge = molecule.charge
        else:
            charge = 0
        if self.multiplicity is not None:
            multiplicity = self.multiplicity
        elif molecule.multiplicity is not None:
            multiplicity = molecule.multiplicity
        else:
            multiplicity = 1
        return charge, multiplicity

    def red_charge_and_multiplicity(self, ox_charge, ox_multiplicity):
        """Charge and multiplicity of the reduced target."""
        default_charge, default_mult = reduced_charge_and_multiplicity(
            ox_charge, ox_multiplicity, self.n_electrons
        )
        if self.red_charge is not None:
            charge = self.red_charge
        else:
            charge = default_charge
        if self.red_multiplicity is not None:
            multiplicity = self.red_multiplicity
        else:
            multiplicity = default_mult
        return charge, multiplicity

    def ox_molecule(self, molecule):
        """Oxidized target molecule with resolved charge and multiplicity."""
        charge, multiplicity = self.ox_charge_and_multiplicity(molecule)
        return _molecule_with_charge(molecule, charge, multiplicity)

    def red_molecule(self, molecule):
        """Reduced target from ``red_file`` or the oxidized geometry."""
        ox_charge, ox_mult = self.ox_charge_and_multiplicity(molecule)
        red_charge, red_mult = self.red_charge_and_multiplicity(
            ox_charge, ox_mult
        )
        if self.red_file is not None:
            return _load_molecule(self.red_file, red_charge, red_mult)
        return _molecule_with_charge(molecule, red_charge, red_mult)

    def _ref_ox_charge_and_multiplicity(self, molecule):
        if self.ref_ox_charge is not None:
            charge = self.ref_ox_charge
        elif self.reference.ox_charge is not None:
            charge = self.reference.ox_charge
        elif molecule.charge is not None:
            charge = molecule.charge
        else:
            charge = 0
        if self.ref_ox_multiplicity is not None:
            multiplicity = self.ref_ox_multiplicity
        elif self.reference.ox_multiplicity is not None:
            multiplicity = self.reference.ox_multiplicity
        elif molecule.multiplicity is not None:
            multiplicity = molecule.multiplicity
        else:
            multiplicity = 1
        return charge, multiplicity

    def _ref_red_charge_and_multiplicity(self, molecule, ox_charge, ox_mult):
        default_charge, default_mult = reduced_charge_and_multiplicity(
            ox_charge, ox_mult, self.n_electrons
        )
        if self.ref_red_charge is not None:
            charge = self.ref_red_charge
        elif self.reference.red_charge is not None:
            charge = self.reference.red_charge
        elif molecule.charge is not None:
            charge = molecule.charge
        else:
            charge = default_charge
        if self.ref_red_multiplicity is not None:
            multiplicity = self.ref_red_multiplicity
        elif self.reference.red_multiplicity is not None:
            multiplicity = self.reference.red_multiplicity
        elif molecule.multiplicity is not None:
            multiplicity = molecule.multiplicity
        else:
            multiplicity = default_mult
        return charge, multiplicity

    def _derived_ref_red_charge_and_multiplicity(self, ox_charge, ox_mult):
        default_charge, default_mult = reduced_charge_and_multiplicity(
            ox_charge, ox_mult, self.n_electrons
        )
        if self.ref_red_charge is not None:
            charge = self.ref_red_charge
        elif self.reference.red_charge is not None:
            charge = self.reference.red_charge
        else:
            charge = default_charge
        if self.ref_red_multiplicity is not None:
            multiplicity = self.ref_red_multiplicity
        elif self.reference.red_multiplicity is not None:
            multiplicity = self.reference.red_multiplicity
        else:
            multiplicity = default_mult
        return charge, multiplicity

    def ref_red_molecule(self, ox_mol, ox_charge, ox_mult):
        """Reduced reference from ``ref_red_file`` or the oxidized geometry."""
        if self.resolved_ref_red_file is not None:
            red_mol = Molecule.from_filepath(self.resolved_ref_red_file)
            red_charge, red_mult = self._ref_red_charge_and_multiplicity(
                red_mol, ox_charge, ox_mult
            )
            return _molecule_with_charge(red_mol, red_charge, red_mult)
        red_charge, red_mult = self._derived_ref_red_charge_and_multiplicity(
            ox_charge, ox_mult
        )
        return _molecule_with_charge(ox_mol, red_charge, red_mult)

    def reference_pair_molecules(self):
        """Load oxidized and reduced reference molecules.

        Raises:
            ValueError: If the oxidized reference geometry path is missing.
        """
        if not self.has_reference_files:
            raise ValueError(
                "Redox exchange calculations require an oxidized reference "
                "geometry (registry ox_file or ref_ox_file)."
            )
        ox_mol = Molecule.from_filepath(self.resolved_ref_ox_file)
        ox_charge, ox_mult = self._ref_ox_charge_and_multiplicity(ox_mol)
        ref_ox_mol = _molecule_with_charge(ox_mol, ox_charge, ox_mult)
        ref_red_mol = self.ref_red_molecule(ox_mol, ox_charge, ox_mult)
        return ref_ox_mol, ref_red_mol


class RedoxChainMixin(ChainMixin):
    """Ox/Red and reference Opt → SP workflow on top of ``ChainMixin``.

    Subclasses set ``_opt_job_class`` and ``_sp_job_class``.
    """

    _opt_job_class = None
    _sp_job_class = None

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if not self.settings.has_reference_files:
            raise ValueError(
                "Redox exchange calculations require an oxidized reference "
                "geometry (registry ox_file or ref_ox_file)."
            )
        self.ox_job = None
        self.red_job = None
        self.ref_ox_job = None
        self.ref_red_job = None
        self.ox_sp_job = None
        self.red_sp_job = None
        self.ref_ox_sp_job = None
        self.ref_red_sp_job = None
        self.opt_jobs = []
        self.ref_opt_jobs = []
        self.sp_jobs = None
        self.ref_sp_jobs = None
        self._prepare_opt_jobs()
        self.phases = self._build_phases()

    def _child_job(self, job_class, molecule, settings, label):
        return job_class(
            molecule=molecule,
            settings=settings,
            label=label,
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )

    def _apply_charge_mult(self, settings, charge, multiplicity):
        settings.charge = charge
        settings.multiplicity = multiplicity
        return settings

    def _opt_settings(self, charge, multiplicity):
        settings = self.settings.copy()
        settings.jobtype = "opt"
        settings.freq = True
        settings.remove_solvent()
        return self._apply_charge_mult(settings, charge, multiplicity)

    def _sp_settings(self, charge, multiplicity):
        settings = self.settings.copy()
        settings.jobtype = "sp"
        settings.freq = False
        return self._apply_charge_mult(settings, charge, multiplicity)

    def _prepare_opt_jobs(self):
        ox_mol = self.settings.ox_molecule(self.molecule)
        red_mol = self.settings.red_molecule(self.molecule)
        self.ox_job = self._child_job(
            self._opt_job_class,
            ox_mol,
            self._opt_settings(ox_mol.charge, ox_mol.multiplicity),
            f"{self.label}_ox_opt",
        )
        self.red_job = self._child_job(
            self._opt_job_class,
            red_mol,
            self._opt_settings(red_mol.charge, red_mol.multiplicity),
            f"{self.label}_red_opt",
        )
        self.opt_jobs = [self.ox_job, self.red_job]

        ref_ox_mol, ref_red_mol = self.settings.reference_pair_molecules()
        self.ref_ox_job = self._child_job(
            self._opt_job_class,
            ref_ox_mol,
            self._opt_settings(ref_ox_mol.charge, ref_ox_mol.multiplicity),
            f"{self.label}_RefOx_opt",
        )
        self.ref_red_job = self._child_job(
            self._opt_job_class,
            ref_red_mol,
            self._opt_settings(ref_red_mol.charge, ref_red_mol.multiplicity),
            f"{self.label}_RefRed_opt",
        )
        self.ref_opt_jobs = [self.ref_ox_job, self.ref_red_job]

    def _make_sp_job(self, opt_job, suffix):
        molecule = self._output_molecule(opt_job, opt_job.molecule)
        charge = opt_job.settings.charge
        multiplicity = opt_job.settings.multiplicity
        molecule = _molecule_with_charge(molecule, charge, multiplicity)
        return self._child_job(
            self._sp_job_class,
            molecule,
            self._sp_settings(charge, multiplicity),
            f"{self.label}_{suffix}",
        )

    def _sp_jobs_for_phase(self):
        if self.sp_jobs is None:
            self.ox_sp_job = self._make_sp_job(self.ox_job, "ox_sp")
            self.red_sp_job = self._make_sp_job(self.red_job, "red_sp")
            self.sp_jobs = [self.ox_sp_job, self.red_sp_job]
        return self.sp_jobs

    def _ref_sp_jobs_for_phase(self):
        if self.ref_sp_jobs is None:
            self.ref_ox_sp_job = self._make_sp_job(self.ref_ox_job, "RefOx_sp")
            self.ref_red_sp_job = self._make_sp_job(
                self.ref_red_job, "RefRed_sp"
            )
            self.ref_sp_jobs = [self.ref_ox_sp_job, self.ref_red_sp_job]
        return self.ref_sp_jobs

    def _build_phases(self):
        return [
            JobPhase(
                name="Opt",
                jobs_factory=lambda: self.opt_jobs,
                require_complete=True,
                stop_on_incomplete=True,
                stop_message="Opt jobs incomplete, halting serial execution.",
            ),
            JobPhase(
                name="Ref Opt",
                jobs_factory=lambda: self.ref_opt_jobs,
                require_complete=True,
                stop_on_incomplete=True,
                stop_message=(
                    "Ref Opt jobs incomplete, halting serial execution."
                ),
            ),
            JobPhase(
                name="SP",
                jobs_factory=self._sp_jobs_for_phase,
                require_complete=True,
                stop_on_incomplete=True,
                stop_message="SP jobs incomplete, halting serial execution.",
            ),
            JobPhase(
                name="Ref SP",
                jobs_factory=self._ref_sp_jobs_for_phase,
                require_complete=True,
                stop_on_incomplete=True,
                stop_message=(
                    "Ref SP jobs incomplete, halting serial execution."
                ),
            ),
        ]

    def _redox_output_files(self):
        return dict(
            ox_gas_file=self.ox_job.outputfile,
            red_gas_file=self.red_job.outputfile,
            ref_ox_gas_file=self.ref_ox_job.outputfile,
            ref_red_gas_file=self.ref_red_job.outputfile,
            ox_solv_file=(
                self.ox_sp_job.outputfile if self.ox_sp_job else None
            ),
            red_solv_file=(
                self.red_sp_job.outputfile if self.red_sp_job else None
            ),
            ref_ox_solv_file=(
                self.ref_ox_sp_job.outputfile if self.ref_ox_sp_job else None
            ),
            ref_red_solv_file=(
                self.ref_red_sp_job.outputfile if self.ref_red_sp_job else None
            ),
        )

    def compute_redox_potential(self):
        """Compute the exchange redox potential from completed child outputs."""
        from chemsmart.analysis.redox import compute_redox_potential

        thermo_kwargs = {
            "cutoff_entropy_grimme": self.settings.cutoff_entropy_grimme,
            "cutoff_enthalpy": self.settings.cutoff_enthalpy,
            "entropy_method": self.settings.entropy_method,
        }
        return compute_redox_potential(
            reference=self.settings.reference,
            n_electrons=self.settings.n_electrons,
            temperature=self.settings.temperature,
            concentration=self.settings.concentration,
            pressure=self.settings.pressure,
            **{
                key: value
                for key, value in thermo_kwargs.items()
                if value is not None
            },
            **self._redox_output_files(),
        )

    def format_redox_summary(self):
        """Return a formatted summary of the exchange redox calculation."""
        from chemsmart.analysis.redox import format_redox_summary

        return format_redox_summary(self.compute_redox_potential())
