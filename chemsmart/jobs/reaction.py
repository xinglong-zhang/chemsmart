"""Shared Endpoint Opt → Guess → Opt → SP chain for Gaussian and ORCA reaction jobs."""

import logging

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.chain import ChainMixin, JobPhase

logger = logging.getLogger(__name__)

PATH_SEARCH_SINGLE_STRUCTURE = (
    "Path search requires a single reactant geometry and a single product "
    "geometry. Combine fragments into one structure, or provide a TS guess "
    "with extra fragments to optimize without path search."
)


def validate_path_search_structures(reactant, product, ts_guess=None):
    """Require matching atom counts and atom order across path-search endpoints."""
    structures = [("reactant", reactant), ("product", product)]
    if ts_guess is not None:
        structures.append(("ts guess", ts_guess))

    for role, molecule in structures:
        if not isinstance(molecule, Molecule):
            raise ValueError(
                f"Path-search {role} must be a Molecule, got "
                f"{type(molecule).__name__}."
            )

    n_atoms = reactant.num_atoms
    ref_symbols = list(reactant.chemical_symbols)
    for role, molecule in structures[1:]:
        if molecule.num_atoms != n_atoms:
            raise ValueError(
                "Path-search structures must have the same number of atoms; "
                f"reactant has {n_atoms}, {role} has {molecule.num_atoms}."
            )
        if list(molecule.chemical_symbols) != ref_symbols:
            raise ValueError(
                "Path-search structures must have the same atom order; "
                f"reactant and {role} chemical symbols differ."
            )


class ReactionChainMixin(ChainMixin):
    """Reaction child-job construction on top of ``ChainMixin``.

    Subclasses set ``_opt_job_class``, ``_ts_job_class``, and
    ``_sp_job_class``, and implement ``_make_guess_job``.

    Path search (QST/NEB) runs only for exactly one reactant geometry and
    one product geometry. Extra fragments with a TS parent skip Guess.
    """

    _opt_job_class = None
    _ts_job_class = None
    _sp_job_class = None
    uses_neb = False

    def __init__(
        self,
        *args,
        reactants=None,
        products=None,
        ts_guess=None,
        no_path_search=False,
        no_sp=False,
        opt_settings=None,
        ts_settings=None,
        sp_settings=None,
        neb_settings=None,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)
        self.reactants = tuple(reactants or ())
        self.products = tuple(products or ())
        self.ts_guess = ts_guess
        self.no_path_search = bool(no_path_search)
        self.no_sp = bool(no_sp)
        self.opt_settings = opt_settings
        self.ts_settings = ts_settings
        self.sp_settings = sp_settings
        self.neb_settings = neb_settings
        self._configure_path_search()
        self._validate_path_search_endpoints()

        self.guess_job = None
        self.reactant_sp_jobs = None
        self.product_sp_jobs = None
        self.ts_sp_job = None
        self.sp_jobs = None
        self._create_opt_jobs()
        self.phases = self._build_reaction_phases()

    def _configure_path_search(self):
        """Skip or reject Guess when fragment counts are not a single R/P pair."""
        if self.no_path_search or not self.products:
            return
        n_reactants = len(self.reactants) if self.reactants else 1
        n_products = len(self.products)
        if n_reactants == 1 and n_products == 1:
            return
        if self.reactants:
            logger.info(
                "Skipping path search for %s: extra fragments with a TS "
                "geometry are optimized, not used as QST/NEB endpoints.",
                self.label,
            )
            self.no_path_search = True
            return
        raise ValueError(PATH_SEARCH_SINGLE_STRUCTURE)

    def _validate_path_search_endpoints(self):
        if not self.uses_path_search:
            return
        validate_path_search_structures(
            self.path_reactant,
            self.path_product,
            ts_guess=self.path_ts_guess,
        )

    @property
    def uses_path_search(self):
        return bool(self.products) and not self.no_path_search

    @property
    def path_reactant(self):
        return self.reactants[0] if self.reactants else self.molecule

    @property
    def path_product(self):
        return self.products[0] if self.products else None

    @property
    def path_ts_guess(self):
        if self.ts_guess is not None:
            return self.ts_guess
        if self.reactants:
            return self.molecule
        return None

    @property
    def reactant_molecules(self):
        if self.reactants:
            return self.reactants
        if self.uses_path_search:
            return (self.molecule,)
        return ()

    @property
    def endpoint_opt_jobs(self):
        if not self.uses_path_search:
            return ()
        return tuple(self.reactant_opt_jobs) + tuple(self.product_opt_jobs)

    def _make_guess_job(self, reactant, product, ts_guess=None):
        raise NotImplementedError

    def _prepare_guess_inputs(self):
        return None

    def _apply_charge_mult(self, settings, molecule):
        if molecule.charge is not None:
            settings.charge = molecule.charge
        if molecule.multiplicity is not None:
            settings.multiplicity = molecule.multiplicity
        return settings

    def _child_settings(self, template, jobtype, freq, molecule):
        base = template if template is not None else self.settings
        settings = base.copy()
        settings.jobtype = jobtype
        settings.freq = freq
        return self._apply_charge_mult(settings, molecule)

    def _ts_settings_for(self, molecule):
        return self._child_settings(self.ts_settings, "ts", True, molecule)

    def _optimized_endpoint_molecule(self, opt_job, fallback):
        if opt_job is None:
            return fallback
        return self._output_molecule(opt_job, fallback)

    def _optimized_reactant_for_guess(self):
        reactant_job = (
            self.reactant_opt_jobs[0] if self.reactant_opt_jobs else None
        )
        return self._optimized_endpoint_molecule(
            reactant_job, self.path_reactant
        )

    def _optimized_product_for_guess(self):
        product_job = (
            self.product_opt_jobs[0] if self.product_opt_jobs else None
        )
        return self._optimized_endpoint_molecule(
            product_job, self.path_product
        )

    def _ts_molecule(self):
        if self.uses_path_search and self.guess_job is not None:
            return self._output_molecule(self.guess_job, self.molecule)
        return self.molecule

    def _child_job(self, job_class, molecule, settings, label):
        return job_class(
            molecule=molecule,
            settings=settings,
            label=label,
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )

    def _opt_jobs_for_role(self, molecules, letter):
        n = len(molecules)
        suffixes = (
            (letter,)
            if n == 1
            else tuple(f"{letter}{i}" for i in range(1, n + 1))
        )
        return [
            self._child_job(
                self._opt_job_class,
                molecule,
                self._child_settings(self.opt_settings, "opt", True, molecule),
                f"{self.label}_{suffix}_opt",
            )
            for molecule, suffix in zip(molecules, suffixes)
        ]

    def _set_opt_jobs(self):
        if self.uses_path_search:
            self.opt_jobs = [self.ts_opt_job]
        else:
            self.opt_jobs = (
                list(self.reactant_opt_jobs)
                + [self.ts_opt_job]
                + list(self.product_opt_jobs)
            )

    def _create_opt_jobs(self):
        self.reactant_opt_jobs = self._opt_jobs_for_role(
            self.reactant_molecules, "R"
        )
        ts_molecule = self._ts_molecule()
        self.ts_opt_job = self._child_job(
            self._ts_job_class,
            ts_molecule,
            self._ts_settings_for(ts_molecule),
            f"{self.label}_TS_opt",
        )
        self.product_opt_jobs = self._opt_jobs_for_role(self.products, "P")
        self._set_opt_jobs()

    def _endpoint_opt_jobs_for_phase(self):
        return list(self.endpoint_opt_jobs)

    def _guess_jobs_for_phase(self):
        reactant = self._optimized_reactant_for_guess()
        product = self._optimized_product_for_guess()
        self.guess_job = self._make_guess_job(
            reactant,
            product,
            ts_guess=self.path_ts_guess,
        )
        self._prepare_guess_inputs()
        return [self.guess_job]

    def _opt_jobs_for_phase(self):
        if self.uses_path_search:
            ts_molecule = self._ts_molecule()
            self.ts_opt_job = self._child_job(
                self._ts_job_class,
                ts_molecule,
                self._ts_settings_for(ts_molecule),
                f"{self.label}_TS_opt",
            )
            self.opt_jobs = [self.ts_opt_job]
        return self.opt_jobs

    def _make_sp_job(self, opt_job):
        molecule = self._output_molecule(opt_job, opt_job.molecule)
        return self._child_job(
            self._sp_job_class,
            molecule,
            self._child_settings(self.sp_settings, "sp", False, molecule),
            f"{opt_job.label.removesuffix('_opt')}_sp",
        )

    def _sp_jobs_for_phase(self):
        if self.sp_jobs is None:
            self.reactant_sp_jobs = [
                self._make_sp_job(job) for job in self.reactant_opt_jobs
            ]
            self.ts_sp_job = self._make_sp_job(self.ts_opt_job)
            self.product_sp_jobs = [
                self._make_sp_job(job) for job in self.product_opt_jobs
            ]
            self.sp_jobs = (
                list(self.reactant_sp_jobs)
                + [self.ts_sp_job]
                + list(self.product_sp_jobs)
            )
        return self.sp_jobs

    def _build_reaction_phases(self):
        return [
            JobPhase(
                name="Endpoint Opt",
                jobs_factory=self._endpoint_opt_jobs_for_phase,
                skip_if=lambda: not self.uses_path_search,
                require_complete=True,
                stop_on_incomplete=True,
                stop_message=(
                    "Endpoint Opt jobs incomplete, halting serial execution."
                ),
            ),
            JobPhase(
                name="Guess",
                jobs_factory=self._guess_jobs_for_phase,
                skip_if=lambda: not self.uses_path_search,
                require_complete=True,
                stop_on_incomplete=True,
                stop_message=(
                    "Guess jobs incomplete, halting serial execution."
                ),
            ),
            JobPhase(
                name="Opt",
                jobs_factory=self._opt_jobs_for_phase,
                require_complete=True,
                stop_on_incomplete=True,
                stop_message="Opt jobs incomplete, halting serial execution.",
            ),
            JobPhase(
                name="SP",
                jobs_factory=self._sp_jobs_for_phase,
                skip_if=lambda: self.no_sp,
                require_complete=True,
                stop_on_incomplete=True,
                stop_message="SP jobs incomplete, halting serial execution.",
            ),
        ]
