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

import logging
import math
from numbers import Integral, Real

from chemsmart.jobs.settings import MolecularJobSettings

logger = logging.getLogger(__name__)

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
PYSCF_JOBTYPES = ("hess", "opt", "sp", "td")
PYSCF_RESPONSE_METHODS = ("tda", "tddft")
PYSCF_STATE_MANIFOLDS = ("singlet",)

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

#: Functionals whose *name* means different things in different programs or
#: target configurations. Maps the ChemSmart literal to an unambiguous PySCF
#: LibXC alias plus a public explanation. In PySCF 2.14, the bare ``B3LYP``
#: alias may be redirected by ``__config__.B3LYP_WITH_VWN5``; ``B3LYPG`` is
#: the explicit VWN3/Gaussian convention and is not changed by that switch.
FUNCTIONAL_DIVERGENCES = {
    "b3lyp": (
        "b3lypg",
        "resolved to PySCF's explicit B3LYPG (VWN3/Gaussian convention) so "
        "target-level B3LYP_WITH_VWN5 configuration cannot change the "
        "calculation. Request b3lyp5 explicitly for the VWN5 variant.",
    ),
}

# Public, deterministic semantics for aliases whose scientific interpretation
# can be stated without launching PySCF.  This registry intentionally stops at
# the alias/convention boundary: exact LibXC primitive IDs, factors, and hybrid
# coefficients remain target-environment observations produced by the child
# driver at execution time.
_REGISTERED_FUNCTIONAL_VARIANTS = {
    "b3lypg": {
        "functional_family": "b3lyp",
        "correlation_convention": "vwn3_gaussian",
        "rule_id": "pyscf.functional.b3lypg_vwn3_gaussian",
    },
    "b3lyp5": {
        "functional_family": "b3lyp",
        "correlation_convention": "vwn5",
        "rule_id": "pyscf.functional.b3lyp5_vwn5",
    },
}


def describe_functional_resolution(functional=None, *, ab_initio=None):
    """Describe the host-side functional literal applied by ChemSmart.

    The record is safe to expose during planning because it is derived by the
    same resolver used by :attr:`PySCFJobSettings.xc`.  It does not pretend to
    be the execution-time LibXC materialization: arbitrary literals and
    composite expressions remain explicitly unclassified until the target
    interpreter records ``libxc.parse_xc`` output in the result artifact.
    """

    method = str(ab_initio or "").strip().lower()
    if method == "hf":
        return {
            "schema_version": "chemsmart.pyscf-functional-resolution.v1",
            "status": "not_applicable",
            "requested_method_kind": "hf",
            "requested_literal": None,
            "normalized_requested_literal": "",
            "applied_xc": None,
            "normalized_applied_xc": "",
            "functional_family": "hartree_fock",
            "correlation_convention": "not_applicable",
            "source": "chemsmart.jobs.pyscf.settings.resolve_functional",
            "rule_id": "pyscf.functional.not_applicable_hf",
        }
    if functional is None or not str(functional).strip():
        return {
            "schema_version": "chemsmart.pyscf-functional-resolution.v1",
            "status": "missing",
            "requested_method_kind": "dft",
            "requested_literal": None,
            "normalized_requested_literal": "",
            "applied_xc": None,
            "normalized_applied_xc": "",
            "functional_family": "",
            "correlation_convention": "unresolved",
            "source": "chemsmart.jobs.pyscf.settings.resolve_functional",
            "rule_id": "pyscf.functional.missing",
        }

    requested = str(functional).strip()
    key = requested.lower()
    divergence = FUNCTIONAL_DIVERGENCES.get(key)
    # Call the executable resolver rather than reproducing its behaviour in
    # this evidence function.  Pass-through and composite strings therefore
    # retain the exact bytes that the writer will place in ``config["xc"]``.
    applied = str(resolve_functional(requested)).strip()
    applied_key = applied.lower()
    variant = _REGISTERED_FUNCTIONAL_VARIANTS.get(applied_key)
    if variant is None:
        return {
            "schema_version": "chemsmart.pyscf-functional-resolution.v1",
            "status": "literal_preserved",
            "requested_method_kind": "dft",
            "requested_literal": requested,
            "normalized_requested_literal": key,
            "applied_xc": applied,
            "normalized_applied_xc": applied_key,
            "functional_family": "unclassified_libxc_literal",
            "correlation_convention": "not_declared",
            "source": "chemsmart.jobs.pyscf.settings.resolve_functional",
            "rule_id": "pyscf.functional.literal_preserved",
        }
    return {
        "schema_version": "chemsmart.pyscf-functional-resolution.v1",
        "status": "registered_alias" if divergence else "explicit_variant",
        "requested_method_kind": "dft",
        "requested_literal": requested,
        "normalized_requested_literal": key,
        "applied_xc": applied,
        "normalized_applied_xc": applied_key,
        **variant,
        "source": "chemsmart.jobs.pyscf.settings.resolve_functional",
    }


def resolve_functional(functional):
    """Return the libxc functional name for a ChemSmart functional name.

    Warns when the requested name is known to denote different functionals
    in different programs, because such a mismatch is invisible in the
    output but shifts total energies by tens of kcal/mol.
    """
    if functional is None:
        return None
    key = str(functional).strip().lower()
    if key in FUNCTIONAL_DIVERGENCES:
        resolved, note = FUNCTIONAL_DIVERGENCES[key]
        logger.warning(f"Functional {functional!r}: {note}")
        return resolved
    return functional


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
        engine (str): ``cpu`` or ``gpu``; recorded in the results provenance
            so every number states which engine produced it.
        scf_maxiter (int): ``mf.max_cycle``.
    """

    #: Base-class fields PySCF cannot honour. Mapped to the value that means
    #: "unset", so ``validate`` can tell "absent" from "explicitly requested".
    UNSUPPORTED = {
        "semiempirical": None,
        "gen_genecp_file": None,
        "heavy_elements": None,
        "heavy_elements_basis": None,
        "light_elements_basis": None,
        "numfreq": False,
        "additional_route_parameters": None,
        "route_to_be_written": None,
        "modred": None,
        "input_string": None,
        "custom_solvent": None,
        "forces": False,
    }

    _INHERITED_SETTING_NAMES = frozenset(UNSUPPORTED)

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
        response_method=None,
        state_manifold=None,
        nstates=None,
        charge=None,
        multiplicity=None,
        freq=False,
        density_fit=False,
        opt_solver="geometric",
        opt_maxsteps=100,
        engine="cpu",
        jobtype=None,
        title=None,
        solvent_model=None,
        solvent_id=None,
        **kwargs,
    ):
        # Internal receipt metadata is deliberately private so it does not
        # become an accepted project-YAML setting. The YAML loader derives it
        # from the source artifact after validating the public settings keys.
        project_yaml_digest = kwargs.pop("_project_yaml_digest", None)
        unknown = sorted(set(kwargs).difference(self._INHERITED_SETTING_NAMES))
        if unknown:
            raise ValueError(
                "Unknown PySCF setting(s): "
                f"{', '.join(unknown)}. Refusing fields that would otherwise "
                "be silently ignored."
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
        self.response_method = response_method
        self.state_manifold = state_manifold
        self.nstates = nstates
        self.density_fit = density_fit
        self.opt_solver = opt_solver
        self.opt_maxsteps = opt_maxsteps
        self.engine = engine
        if project_yaml_digest is not None:
            self._project_yaml_digest = str(project_yaml_digest)

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
        return getattr(self, "_project_yaml_digest", None)

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
        return int(self.multiplicity) - 1

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

        if (self.charge is None) != (self.multiplicity is None):
            raise ValueError(
                "PySCF charge and multiplicity overrides must be supplied "
                "together, or both inherited from the molecular source."
            )

        if self.opt_solver not in PYSCF_OPT_SOLVERS:
            raise ValueError(
                f"Unknown opt_solver {self.opt_solver!r}; "
                f"expected one of {PYSCF_OPT_SOLVERS}."
            )

        if self.engine not in PYSCF_ENGINES:
            raise ValueError(
                f"Unknown engine {self.engine!r}; "
                f"expected one of {PYSCF_ENGINES}."
            )

        if self.jobtype not in PYSCF_JOBTYPES:
            raise ValueError(
                f"Unknown PySCF jobtype {self.jobtype!r}; expected exactly "
                f"one of {PYSCF_JOBTYPES}."
            )
        if type(self.density_fit) is not bool:
            raise ValueError(
                "density_fit must be a strict boolean, got "
                f"{self.density_fit!r}."
            )
        if type(self.freq) is not bool:
            raise ValueError(
                f"freq must be a strict boolean, got {self.freq!r}."
            )
        if self.jobtype != "hess" and self.freq:
            raise ValueError(
                f"PySCF {self.jobtype} requires freq=False; use a separate "
                "hess node so its input geometry is explicit."
            )
        if self.jobtype == "hess" and not self.freq:
            raise ValueError(
                "PySCF hess requires freq=True so the project setting "
                "matches the executed Hessian stage."
            )

        td_fields = {
            "response_method": self.response_method,
            "state_manifold": self.state_manifold,
            "nstates": self.nstates,
        }
        if self.jobtype == "td":
            response_method = str(self.response_method or "").strip().lower()
            if response_method not in PYSCF_RESPONSE_METHODS:
                raise ValueError(
                    "PySCF td requires response_method to be one of "
                    f"{PYSCF_RESPONSE_METHODS}, got {self.response_method!r}."
                )
            manifold = str(self.state_manifold or "").strip().lower()
            if manifold not in PYSCF_STATE_MANIFOLDS:
                raise ValueError(
                    "PySCF td currently supports only an explicit singlet "
                    f"state_manifold, got {self.state_manifold!r}."
                )
            if (
                isinstance(self.nstates, bool)
                or not isinstance(self.nstates, Integral)
                or int(self.nstates) <= 0
            ):
                raise ValueError(
                    "PySCF td requires nstates to be a positive integer, "
                    f"got {self.nstates!r}."
                )
            if self.ab_initio is not None or not self.functional:
                raise ValueError(
                    "PySCF td preview supports closed-shell DFT only; set a "
                    "functional and do not set ab_initio."
                )
            if self.multiplicity not in (None, 1):
                raise ValueError(
                    "PySCF td preview supports only a closed-shell singlet "
                    f"reference, got multiplicity={self.multiplicity!r}."
                )
            if self.solvent_model is not None or self.solvent_id is not None:
                raise ValueError(
                    "PySCF td preview is gas-phase only; solvent response is "
                    "not enabled."
                )
            if self.dispersion is not None:
                raise ValueError(
                    "PySCF td preview does not accept a ground-state "
                    "dispersion correction as an excited-state setting."
                )
            if str(self.engine or "cpu").strip().lower() != "cpu":
                raise ValueError(
                    "PySCF td is a CPU-only preview capability; "
                    "GPU4PySCF response calculations are not validated."
                )
        elif any(value is not None for value in td_fields.values()):
            populated = ", ".join(
                key for key, value in td_fields.items() if value is not None
            )
            raise ValueError(
                f"PySCF {populated} are valid only for the td jobtype."
            )
        if self.scf_tol is not None and (
            isinstance(self.scf_tol, bool)
            or not isinstance(self.scf_tol, Real)
            or not math.isfinite(float(self.scf_tol))
            or float(self.scf_tol) <= 0
        ):
            raise ValueError(
                f"scf_tol must be finite and > 0, got {self.scf_tol!r}."
            )
        for field in ("scf_maxiter", "opt_maxsteps"):
            value = getattr(self, field)
            if value is None:
                continue
            if (
                isinstance(value, bool)
                or not isinstance(value, Integral)
                or int(value) <= 0
            ):
                raise ValueError(
                    f"{field} must be a positive integer, got {value!r}."
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
