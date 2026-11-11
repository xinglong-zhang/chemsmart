"""PySCF job settings.

``PySCFJobSettings`` subclasses ``MolecularJobSettings`` so that the project
YAML dialect, the ``keywords`` merge whitelist, and ``modify_solvent`` behave
exactly as they do for Gaussian and ORCA. It reuses the existing option
vocabulary (``functional``, ``ab_initio``, ``basis``, ``aux_basis``,
``defgrid``, ``scf_tol``, ``solvent_model``, ``solvent_id``) rather than
inventing PySCF-specific names, so one normalised vocabulary spans programs.

Fields the base class carries but PySCF cannot honour are listed in
``UNSUPPORTED`` and rejected by :meth:`validate`. A project YAML asking for a
genecp basis must not silently produce a plain def2-SVP calculation.
"""

import inspect
import logging
from functools import lru_cache

from ase.data import atomic_numbers

from chemsmart.jobs.settings import MolecularJobSettings

logger = logging.getLogger(__name__)

#: Settings ``MolecularJobSettings`` accepts. Its ``**kwargs`` discards any
#: other name, so a misspelt PySCF setting would vanish without a trace.
MOLECULAR_SETTINGS_FIELDS = frozenset(
    name
    for name, parameter in inspect.signature(
        MolecularJobSettings.__init__
    ).parameters.items()
    if parameter.kind is parameter.POSITIONAL_OR_KEYWORD and name != "self"
)

#: Solvent models accepted by ``solvent_model``, mapped to the PySCF call and
#: the ``with_solvent.method`` string. Verified against pyscf 2.14.0: the PCM
#: methods are exactly C-PCM, COSMO, IEF-PCM and SS(V)PE.
PYSCF_SOLVENT_MODELS = {
    "pcm": ("PCM", "IEF-PCM"),
    "iefpcm": ("PCM", "IEF-PCM"),
    "cpcm": ("PCM", "C-PCM"),
    "cosmo": ("PCM", "COSMO"),
    "ssvpe": ("PCM", "SS(V)PE"),
    "smd": ("SMD", None),
}

#: ``defgrid`` levels mapped to ``mf.grids.atom_grid`` (radial, angular).
#: Chosen to bracket ORCA's DEFGRID1/2/3 in cost, not to reproduce its exact
#: quadrature, which is not publicly specified.
PYSCF_DEFGRIDS = {
    "defgrid1": (50, 194),
    "defgrid2": (75, 302),
    "defgrid3": (99, 590),
}

#: Geometry-optimisation backends. ``geometric`` is the default because it is
#: what ``gpu4pyscf/drivers/opt_driver.py`` calls and what PySCF's own
#: ``[geomopt]`` extra installs first.
PYSCF_OPT_SOLVERS = ("geometric", "berny", "ase")

#: Execution engines. ``gpu`` routes through gpu4pyscf via ``.to_gpu()``.
PYSCF_ENGINES = ("cpu", "gpu")

#: Functionals for which the perturbative correlation term is not implemented
#: by this v1 mean-field-only backend.  Passing one of these names to
#: ``dft.KS`` can produce a converged energy that contains only the DFT part of
#: the advertised double hybrid, which is scientifically worse than a hard
#: failure.
PYSCF_DOUBLE_HYBRID_MARKERS = (
    "b2plyp",
    "b2gp-plyp",
    "b2gpplyp",
    "pwpb95",
    "pbe0-2",
    "wb97x-2",
    "qidh",
    "-dh",
    "revdsd",
    "dsd",
    "dhdft",
    "xyg3",
    "xygjos",
    "double-hybrid",
)

#: Characters that make PySCF read a functional string as an expression
#: (``pbe,pbe``, ``0.5*b3lyp``) instead of one named functional. CHEMSMART
#: requests functionals by name, so a functional is written one way only.
PYSCF_FUNCTIONAL_EXPRESSION_CHARACTERS = frozenset(",*+;() ")

#: Names PySCF accepts that do not denote what their spelling suggests,
#: mapped to the reason they are refused.
PYSCF_MISLEADING_FUNCTIONALS = {
    "b3lyp3": (
        "LibXC's B3LYP3 (functional 394) is not Gaussian's B3LYP, although "
        "ORCA calls Gaussian's variant 'VWN3'; for H2O/def2-SVP the two "
        "differ by 23 kcal/mol. Use 'b3lyp' for Gaussian's B3LYP (ORCA's "
        "B3LYP/G), 'b3lyp5' for ORCA's B3LYP, or 'hyb_gga_xc_b3lyp3' for "
        "LibXC functional 394 itself."
    ),
}

#: Functionals whose *name* means different things in different programs.
#: Maps a CHEMSMART functional name to (libxc name, note).
#:
#: Measured on H2O/def2-SVP at a fixed geometry (pyscf 2.14.0, ORCA 6.1.1):
#:
#:   PySCF 'b3lyp'  (VWN3, Gaussian's definition) = -76.35815198 Ha
#:   PySCF 'b3lyp5' (VWN5, ORCA's definition)     = -76.32100288 Ha
#:   ORCA  'B3LYP'                                = -76.32111960 Ha
#:
#: i.e. the two B3LYP definitions differ by **23.2 kcal/mol**, far larger
#: than any density-fitting or grid error (0.006 kcal/mol here). CHEMSMART
#: resolves ``b3lyp`` to the Gaussian/VWN3 definition, which is both libxc's
#: default and the dominant convention in the literature. Ask for
#: ``b3lyp5`` explicitly to reproduce an ORCA B3LYP number.
FUNCTIONAL_DIVERGENCES = {
    "b3lyp": (
        "b3lyp",
        "resolved to the Gaussian/VWN3 definition; ORCA's B3LYP is the "
        "VWN5 variant, which is 'b3lyp5' here and differs by ~23 kcal/mol "
        "for H2O/def2-SVP. Do not compare a PySCF 'b3lyp' energy with an "
        "ORCA 'B3LYP' energy.",
    ),
}


def resolve_functional(functional):
    """Return the libxc functional name for a CHEMSMART functional name.

    Warns when the requested name is known to denote different functionals
    in different programs, because such a mismatch is invisible in the
    output but shifts total energies by tens of kcal/mol.
    """
    if functional is None:
        return None
    key = str(functional).strip().lower()
    if key in FUNCTIONAL_DIVERGENCES:
        resolved, note = FUNCTIONAL_DIVERGENCES[key]
        _warn_once(f"Functional {functional!r}: {note}")
        return resolved
    return functional


@lru_cache(maxsize=None)
def _warn_once(message):
    """Log a message once per process; ``.xc`` is read many times per job."""
    logger.warning(message)


def is_double_hybrid_functional(functional):
    """Return whether *functional* names a known double-hybrid family."""
    if functional is None:
        return False
    normal = str(functional).strip().lower().replace("_", "-")
    return any(marker in normal for marker in PYSCF_DOUBLE_HYBRID_MARKERS)


class PySCFJobSettings(MolecularJobSettings):
    """Configuration for a PySCF calculation.

    Attributes beyond the shared molecular vocabulary:
        density_fit (bool): Apply ``.density_fit()`` to the mean-field object.
        opt_solver (str): One of ``PYSCF_OPT_SOLVERS``.
        opt_maxsteps (int): Geometry-optimisation step ceiling.
        engine (str): ``cpu``, ``gpu``, or None when not requested, in which
            case the CLI resolves it from the server. Recorded in the results
            provenance so every number states which engine produced it.
        scf_maxiter (int): ``mf.max_cycle``.
        point_group_tol (float): Geometry tolerance in Bohr for point-group
            detection, which sets the rotational symmetry number.
        forces (bool): Compute a gradient for a single point or Hessian job.
            An optimisation always stores the final forces.
    """

    #: Base-class fields PySCF cannot honour. Mapped to the value that means
    #: "unset", so ``validate`` can tell "absent" from "explicitly requested".
    UNSUPPORTED = {
        "semiempirical": None,
        "modred": None,
        "gen_genecp_file": None,
        "heavy_elements": None,
        "heavy_elements_basis": None,
        "light_elements_basis": None,
        "numfreq": False,
        "additional_route_parameters": None,
        "route_to_be_written": None,
        "input_string": None,
        "custom_solvent": None,
    }

    def __init__(
        self,
        ab_initio=None,
        functional=None,
        dispersion=None,
        basis=None,
        aux_basis=None,
        defgrid=None,
        scf_tol=None,
        scf_maxiter=None,
        charge=None,
        multiplicity=None,
        freq=False,
        density_fit=True,
        opt_solver="geometric",
        opt_maxsteps=100,
        engine=None,
        point_group_tol=1.0e-2,
        jobtype=None,
        title=None,
        solvent_model=None,
        solvent_id=None,
        **kwargs,
    ):
        # Passed under its private name by ``merge``. Not a project-YAML key:
        # the YAML loader derives it from the source file.
        project_yaml_digest = kwargs.pop("_project_yaml_digest", None)
        unknown = sorted(set(kwargs) - MOLECULAR_SETTINGS_FIELDS)
        if unknown:
            raise TypeError(
                f"Unknown PySCF settings: {', '.join(unknown)}. A setting "
                "that is not recognised would otherwise be ignored."
            )
        super().__init__(
            ab_initio=ab_initio,
            functional=functional,
            dispersion=dispersion,
            basis=basis,
            defgrid=defgrid,
            charge=charge,
            multiplicity=multiplicity,
            freq=freq,
            jobtype=jobtype,
            title=title,
            solvent_model=solvent_model,
            solvent_id=solvent_id,
            **kwargs,
        )
        self.aux_basis = aux_basis
        self.scf_tol = scf_tol
        self.scf_maxiter = scf_maxiter
        self.density_fit = density_fit
        self.opt_solver = opt_solver
        self.opt_maxsteps = opt_maxsteps
        self.engine = engine
        self.point_group_tol = point_group_tol
        self._project_yaml_digest = project_yaml_digest

    @classmethod
    def default(cls):
        """Return settings with every field at its unset default."""
        return cls()

    def copy(self):
        import copy as _copy

        return _copy.deepcopy(self)

    @property
    def project_yaml_digest(self):
        """SHA-256 of the project YAML that produced these settings."""
        return self._project_yaml_digest

    @project_yaml_digest.setter
    def project_yaml_digest(self, digest):
        self._project_yaml_digest = digest

    @property
    def resolved_engine(self):
        """Return the engine to run on; CPU unless GPU was requested."""
        return self.engine or "cpu"

    def merge(
        self, other, keywords=("charge", "multiplicity"), merge_all=False
    ):
        """Merge ``other`` into a copy of self.

        Mirrors ``ORCAJobSettings.merge``: without ``merge_all`` only the
        names in ``keywords`` cross over. A CLI override whose name was never
        appended to ``keywords`` is dropped here silently, which is the single
        most common bug in this codebase.
        """
        other_dict = other if isinstance(other, dict) else other.__dict__
        if not merge_all and keywords is not None:
            other_dict = {
                k: other_dict[k] for k in keywords if k in other_dict
            }
        merged_dict = self.__dict__.copy()
        merged_dict.update(other_dict)
        return type(self)(**merged_dict)

    @property
    def spin(self):
        """Return PySCF's ``mol.spin`` (= 2S = Nalpha - Nbeta).

        PySCF's ``spin`` is **not** the multiplicity. Passing a multiplicity
        straight into ``mol.spin`` silently computes a different electronic
        state, so this conversion lives in exactly one place.
        """
        if self.multiplicity is None:
            return None
        return self.multiplicity_to_spin(self.multiplicity)

    @staticmethod
    def multiplicity_to_spin(multiplicity):
        """Convert a multiplicity (2S+1) to PySCF's ``mol.spin`` (2S)."""
        return int(multiplicity) - 1

    @property
    def is_restricted(self):
        """Return whether an R (True) or U (False) reference is implied.

        Derived from multiplicity rather than exposed as a flag, so the
        reference cannot disagree with the requested electronic state.
        """
        return self.multiplicity in (None, 1)

    @property
    def method_name(self):
        """Return the method label used for provenance and the database."""
        if self.ab_initio:
            return str(self.ab_initio).lower()
        if self.functional:
            return str(self.functional).lower()
        return None

    @property
    def xc(self):
        """Return the libxc functional string, or None for a HF reference.

        Routed through :func:`resolve_functional`, which warns for names
        that denote different functionals in different programs.
        """
        if self.ab_initio and str(self.ab_initio).lower() == "hf":
            return None
        return resolve_functional(self.functional)

    def _check_functional_name(self):
        """Reject a functional string PySCF would not read as that name.

        Hyphenated names such as ``cam-b3lyp`` and ``b3lyp-vwn3`` (B3LYP
        minus VWN) can only be told apart with PySCF's own tables, so the
        driver checks those in the environment that runs the calculation.
        """
        functional = self.functional
        name = str(functional).strip().lower()
        if PYSCF_FUNCTIONAL_EXPRESSION_CHARACTERS.intersection(name):
            raise ValueError(
                f"Functional {functional!r} is a PySCF expression, not a "
                "named functional. Request the functional by name, for "
                "example 'pbe0' or 'b3lyp'."
            )
        if "-d3" in name or "-d4" in name:
            raise ValueError(
                f"Functional {functional!r} carries a dispersion suffix, "
                "which PySCF would apply on its own. Give the dispersion "
                "through the dispersion setting instead, for example "
                "functional 'b3lyp' with dispersion 'd3bj'."
            )
        if name in PYSCF_MISLEADING_FUNCTIONALS:
            raise ValueError(
                f"Functional {functional!r} is refused: "
                f"{PYSCF_MISLEADING_FUNCTIONALS[name]}"
            )

    def _check_solvent(self, solvent_model):
        """Reject a solvent model PySCF has no implementation for."""
        if str(solvent_model).lower() not in PYSCF_SOLVENT_MODELS:
            raise ValueError(
                f"Solvent model {solvent_model!r} is not available in PySCF. "
                f"Available models: {sorted(PYSCF_SOLVENT_MODELS)}"
            )

    def unsupported_requests(self):
        """Return ``[(field, value)]`` for every unsupported field that is set.

        Returns rather than raises so a caller can enumerate every problem at
        once instead of discovering them one exception at a time.
        """
        requested = []
        for field, unset in self.UNSUPPORTED.items():
            value = getattr(self, field, unset)
            if value != unset:
                requested.append((field, value))
        return requested

    def validate_electrons(self, molecule):
        """Raise if charge and multiplicity are impossible for a molecule.

        PySCF would otherwise stop at the start of the run, after the job has
        been launched. Charge and multiplicity must already be resolved.

        Raises:
            ValueError: If the electrons cannot form the requested spin state.
        """
        nuclear_charge = sum(
            atomic_numbers[symbol] for symbol in molecule.chemical_symbols
        )
        num_electrons = nuclear_charge - int(self.charge)
        unpaired = self.multiplicity_to_spin(self.multiplicity)
        problem = None
        if num_electrons < 0:
            problem = "leaves a negative number of electrons"
        elif (num_electrons - unpaired) % 2:
            parity = "even" if num_electrons % 2 else "odd"
            problem = f"needs an {parity} multiplicity"
        elif unpaired > num_electrons:
            problem = f"cannot have {unpaired} unpaired electrons"
        if problem:
            raise ValueError(
                f"Charge {self.charge} and multiplicity {self.multiplicity} "
                f"are impossible for {num_electrons} electrons (nuclear "
                f"charge {nuclear_charge}): this electron count {problem}. "
                f"Check -c/--charge and -m/--multiplicity."
            )

    def validate(self):
        """Raise if the settings request something PySCF cannot deliver.

        Raises:
            ValueError: If an unsupported field is set, or if a supported
                field holds a value outside its allowed set.
        """
        unsupported = self.unsupported_requests()
        if unsupported:
            details = ", ".join(f"{k}={v!r}" for k, v in unsupported)
            raise ValueError(
                f"These settings are not supported by the PySCF backend: "
                f"{details}.\nRemove them from the project YAML or the "
                f"command line; PySCF would otherwise run a different "
                f"calculation than the one requested."
            )

        if self.ab_initio is not None and str(self.ab_initio).lower() != "hf":
            raise ValueError(
                f"Unknown ab_initio method {self.ab_initio!r}; the PySCF "
                "v1 backend supports only 'hf' or a DFT functional."
            )

        if self.ab_initio is not None and self.functional is not None:
            raise ValueError(
                "Specify either 'ab_initio: hf' or 'functional', not both."
            )

        if self.functional is not None:
            self._check_functional_name()

        if is_double_hybrid_functional(self.functional):
            raise ValueError(
                f"Double-hybrid functional {self.functional!r} is not "
                "supported by the PySCF v1 backend because its "
                "perturbative correlation term would not be applied."
            )

        if self.solvent_model is None and self.solvent_id is not None:
            raise ValueError(
                "solvent_id cannot be applied without a solvent_model."
            )

        if self.solvent_model is not None:
            self._check_solvent(self.solvent_model)
            if not self.solvent_id:
                raise ValueError(
                    "A solvent_id is required for every PySCF solvent model "
                    "so the applied dielectric/environment is explicit."
                )

        if float(self.point_group_tol) <= 0:
            raise ValueError(
                f"point_group_tol must be positive, got "
                f"{self.point_group_tol!r}."
            )

        if self.opt_solver not in PYSCF_OPT_SOLVERS:
            raise ValueError(
                f"Unknown opt_solver {self.opt_solver!r}; "
                f"expected one of {PYSCF_OPT_SOLVERS}."
            )

        if self.engine is not None and self.engine not in PYSCF_ENGINES:
            raise ValueError(
                f"Unknown engine {self.engine!r}; "
                f"expected one of {PYSCF_ENGINES}."
            )

        if self.defgrid is not None:
            if str(self.defgrid).lower() not in PYSCF_DEFGRIDS:
                raise ValueError(
                    f"Unknown defgrid {self.defgrid!r}; "
                    f"expected one of {sorted(PYSCF_DEFGRIDS)}."
                )
            if self.ab_initio is not None:
                raise ValueError(
                    "defgrid is a DFT-only setting and cannot be applied to "
                    "an ab_initio HF calculation."
                )

        if not self.density_fit and self.aux_basis is not None:
            raise ValueError(
                "aux_basis cannot be applied when density_fit is disabled."
            )

        if self.multiplicity is not None and int(self.multiplicity) < 1:
            raise ValueError(
                f"Multiplicity must be >= 1, got {self.multiplicity!r}."
            )

        if self.xc is None and not self.ab_initio:
            raise ValueError(
                "No method specified: set 'functional' (DFT) or "
                "'ab_initio: hf' in the project YAML."
            )

        if not self.basis:
            raise ValueError(
                "No basis set specified; PySCF has no default basis."
            )
        return self
