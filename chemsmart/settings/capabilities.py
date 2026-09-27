"""Immutable declarations for executable chemistry programs.

This module is deliberately independent of ``chemsmart.agent``.  It is the
single declaration that a future agent-harness merge can project its program
sets, project requirements, job-kind maps, and engine capabilities from.
"""

from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType
from typing import Mapping


def _validate_names(
    field: str, values: tuple[str, ...], *, allow_empty: bool = False
) -> None:
    """Require a deterministic tuple of unique, normalised identifiers."""

    if not values and not allow_empty:
        raise ValueError(f"{field} must not be empty")
    if values != tuple(sorted(set(values))):
        raise ValueError(f"{field} must be sorted and contain no duplicates")
    invalid = (
        not value.isidentifier() or value.lower() != value for value in values
    )
    if any(invalid):
        raise ValueError(f"{field} must contain lower-case identifiers")


@dataclass(frozen=True, order=True)
class EngineJobCapability:
    """One exact program-engine/job pairing and its execution boundary.

    ``preview_supported`` describes whether ChemSmart can compile and fake-run
    the pair.  ``execution_supported`` is deliberately narrower: it says the
    implementation is allowed to become executable after environment,
    approval, and validation gates pass.  Neither flag establishes that the
    current host is ready.
    """

    engine: str
    jobtype: str
    preview_supported: bool = True
    execution_supported: bool = True

    def __post_init__(self) -> None:
        _validate_names("engine", (self.engine,))
        _validate_names("jobtype", (self.jobtype,))
        if self.execution_supported and not self.preview_supported:
            raise ValueError("execution support requires preview support")


@dataclass(frozen=True)
class ProgramCapability:
    """CLI facts and the bounded agent surface for one chemistry program.

    ``jobtypes`` and ``engines`` describe the human CLI.  The exact
    ``engine_job_capabilities`` matrix is the narrower surface that Runtime V2
    may expose to a model.  ``agent_enabled`` keeps a human-only CLI program in
    this single canonical registry without silently advertising it to the
    agent.
    """

    program: str
    requires_project_configuration: bool
    supports_project_configuration: bool
    jobtypes: tuple[str, ...]
    project_owned_parameters: tuple[str, ...]
    engines: tuple[str, ...]
    project_parameter_domains: tuple[tuple[str, tuple[str, ...]], ...] = ()
    engine_job_capabilities: tuple[EngineJobCapability, ...] = ()
    #: How this program spells a driven coordinate on its own CLI.
    #: ``absolute_range`` takes the two endpoints and a point count;
    #: ``increment_steps`` takes a step size and an interval count and
    #: walks outward from the supplied geometry, so the increment
    #: carries the direction. A program that declares neither cannot be
    #: asked for a scan, and the refusal says so by name rather than by
    #: a branch on the program's name in the renderer.
    coordinate_idiom: str = ""
    agent_enabled: bool = True
    #: Top-level section names this program's project YAML accepts. The
    #: route-building programs group settings by phase (``gas`` for most job
    #: types, ``solv`` for ``sp``); PySCF keys sections by job type instead.
    #: Declaring it makes a wrong shape refusable at authoring time rather
    #: than surfacing deep inside a loader as an opaque AttributeError.
    project_section_names: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if (
            not self.program.isidentifier()
            or self.program.lower() != self.program
        ):
            raise ValueError("program must be a lower-case identifier")
        if (
            self.requires_project_configuration
            and not self.supports_project_configuration
        ):
            raise ValueError(
                "a program cannot require unsupported project configuration"
            )
        # A Click leaf has no child jobtype commands, so an empty tuple is a
        # meaningful inventory rather than a missing declaration.
        _validate_names("jobtypes", self.jobtypes, allow_empty=True)
        if self.project_owned_parameters:
            _validate_names(
                "project_owned_parameters", self.project_owned_parameters
            )
        _validate_names("engines", self.engines)
        _validate_names(
            "project_section_names",
            self.project_section_names,
            allow_empty=True,
        )
        domain_names = tuple(
            name for name, _values in self.project_parameter_domains
        )
        if domain_names != tuple(sorted(set(domain_names))):
            raise ValueError(
                "project_parameter_domains must be sorted and unique"
            )
        for name, values in self.project_parameter_domains:
            if name not in self.project_owned_parameters:
                raise ValueError(
                    "a parameter domain must target a project-owned parameter"
                )
            if not values or values != tuple(sorted(set(values))):
                raise ValueError(
                    "project parameter domain values must be sorted and unique"
                )
            if any(not value or value != value.lower() for value in values):
                raise ValueError(
                    "project parameter domain values must be lower-case"
                )
        if self.engine_job_capabilities:
            if self.engine_job_capabilities != tuple(
                sorted(set(self.engine_job_capabilities))
            ):
                raise ValueError(
                    "engine_job_capabilities must be sorted and unique"
                )
            for item in self.engine_job_capabilities:
                if item.engine not in self.engines:
                    raise ValueError(
                        "engine-job capability uses an undeclared engine"
                    )
                if item.jobtype not in self.jobtypes:
                    raise ValueError(
                        "engine-job capability uses an undeclared jobtype"
                    )
        if not isinstance(self.agent_enabled, bool):
            raise ValueError("agent_enabled must be boolean")

    @property
    def resolved_engine_job_capabilities(
        self,
    ) -> tuple[EngineJobCapability, ...]:
        """Return the exact matrix, deriving the legacy Cartesian declaration.

        Empty matrices remain valid migration input for the original registry,
        where engines and jobtypes were independent lists.  New non-Cartesian
        programs must declare their exact matrix.
        """

        if self.engine_job_capabilities:
            return self.engine_job_capabilities
        return tuple(
            EngineJobCapability(engine=engine, jobtype=jobtype)
            for engine in self.engines
            for jobtype in self.jobtypes
        )


# This is the project-owned option-name union used today by the agent harness
# for both Gaussian and ORCA.  Keeping it intact preserves existing synthesis
# and migration behaviour when those call sites become registry views.
_CURRENT_HARNESS_PROJECT_PARAMETERS = (
    "ab_initio",
    "additional_opt_options",
    "additional_route_parameters",
    "append_additional_info",
    "aux_basis",
    "basis",
    "custom_solvent",
    "defgrid",
    "dieze_tag",
    "dispersion",
    "extrapolation_basis",
    "functional",
    "geom_maxiter",
    "opt_convergence",
    "scf_algorithm",
    "scf_convergence",
    "scf_maxiter",
    "scf_tol",
    "semiempirical",
    "solvent_id",
    "solvent_model",
    "solvent_options",
    "solventfilename",
)


def _settable_parameters(names, module_name, class_names):
    """Keep only names some settings class for this program can accept.

    A project-owned parameter is advertised to the model as something it
    may set, and the model can only set one through project YAML -- where
    the loader refuses any key absent from the stage defaults. A name on
    this list that no settings class carries is therefore an instruction
    that cannot be followed: the model reads the capability, writes the
    key, and the loader rejects the project.

    Thirteen names were in that state -- four for ORCA and nine for
    Gaussian, the latter including four core SCF controls. Filtering here
    rather than by hand-editing the tuples keeps the guarantee true as
    the settings classes change, instead of true on the day someone last
    checked.

    Subclasses count, because the loader lifts a section into its own
    settings class for the jobtypes that have one.

    A field the module declares its Agent path refuses whenever it is set
    (``AGENT_REFUSED_FIELDS``: a whole input or route, free lines, a file
    path) is settable by a person and not offered here, for the same
    reason: the project tool would refuse the key the model was offered.
    """

    import importlib

    module = importlib.import_module(module_name)
    settable: set[str] = set()
    for class_name in class_names:
        cls = getattr(module, class_name, None)
        if cls is None:
            continue
        try:
            settable.update(cls.default().__dict__)
        except Exception:  # pragma: no cover - defensive
            settable.update(getattr(cls, "__dataclass_fields__", {}) or {})
    settable.difference_update(getattr(module, "AGENT_REFUSED_FIELDS", ()))
    return tuple(sorted(name for name in names if name in settable))


# ORCA has typed method controls that are scientifically stronger than the
# generic route-string escape hatch.  Advertising them through the canonical
# capability registry lets an agent choose explicit scalar-relativistic,
# reference, RI, frozen-core, and element-specific basis semantics in project
# YAML instead of hiding those choices in prose.
_ORCA_PROJECT_PARAMETERS = tuple(
    sorted(
        {
            *_CURRENT_HARNESS_PROJECT_PARAMETERS,
            "additional_solvent_options",
            # The broken-symmetry open-shell singlet, one typed request in
            # every program's words (HFTyp UHF + GuessMix here), instead of
            # FlipSpin on the simple-input line (R10 Q15 g1, ax41 ino2).
            "broken_symmetry",
            "dipole",
            "direction",
            "forces",
            "freq",
            "frozen_core",
            "frozen_core_electrons",
            "gbw",
            "heavy_elements",
            "heavy_elements_basis",
            "hessmode",
            "inithess",
            "light_elements_basis",
            "mdci_cutoff",
            "mdci_density",
            "joboption",
            "nimages",
            "preopt_ends",
            "nstates",
            "numfreq",
            "quadrupole",
            "reference",
            "relativistic",
            "response_method",
            "ri_approximation",
            # The other broken-symmetry mechanism: a flip of named centres
            # from the high-spin state (FlipSpin/FinalMs, written in ORCA's
            # 0-based numbering), for centres with more than one unpaired
            # electron, where the singlet mixing guess reached another
            # state (ino2's Ni(II)2, R11 truth-2, CUHK 2157086).
            "site_spin_flip",
            "state_manifold",
            # A saddle search's tuning controls. These are method
            # rationale and reusable across molecules, which is what
            # project YAML is for. Deliberately absent: the
            # hybrid-Hessian atom subset and the ScanTS coordinate,
            # which are facts about one molecule in one calculation and
            # belong on the workflow node, by the same reasoning that
            # keeps a scan's driven coordinate off the project; and the
            # initial-Hessian file pair, which arrives as a producer
            # edge and whose writer asserts the file exists.
            "full_scan",
            "numhess",
            "recalc_hess",
            "trust_radius",
            "tssearch_type",
            "vpt2",
            "vpt2_anharmonic_displacement",
            "vpt2_hessian_cutoff",
        }
    )
)
_ORCA_PROJECT_PARAMETERS = _settable_parameters(
    _ORCA_PROJECT_PARAMETERS,
    "chemsmart.jobs.orca.settings",
    (
        "ORCAJobSettings",
        "ORCAIRCJobSettings",
        "ORCANEBJobSettings",
        "ORCATSJobSettings",
    ),
)

# Gaussian exposes scientific controls the shared union omitted.
_GAUSSIAN_PROJECT_PARAMETERS = tuple(
    sorted(
        _CURRENT_HARNESS_PROJECT_PARAMETERS
        + (
            "additional_opt_options_in_route",
            "additional_solvent_options",
            # The broken-symmetry open-shell singlet: the unrestricted
            # method and its guess=mix written together, where a native
            # guess=mix on a restricted route ran restricted (R10 Q15 g1).
            "broken_symmetry",
            # What a reaction path is: which way it walks, how far, how
            # often the curvature is recomputed and how big a step is.
            # The CLI has taken all of these since the job existed and
            # the native writer writes them inside ``irc(...)``; the
            # project loader refused them, so the only channel a session
            # had was to append a bare ``maxpoints=50`` beside the IRC
            # keyword, which is not a Gaussian route keyword at all.
            "direction",
            "eqsolv",
            "flat_irc",
            "forces",
            "freq",
            "heavy_elements_basis",
            "guess",
            "jobtype",
            "link_route",
            "maxcycles",
            "maxpoints",
            "nstates",
            "numfreq",
            "predictor",
            "recalc_step",
            "recorrect",
            # The response and the manifold of a td stage, in the words
            # ORCA's and PySCF's settings take, so one td section is one
            # request in all three programs; ``states`` is Gaussian's own
            # older word for the manifold and is still read.
            "response_method",
            "root",
            "stable",
            "state_manifold",
            "states",
            "stepsize",
        )
    )
)
_GAUSSIAN_PROJECT_PARAMETERS = _settable_parameters(
    _GAUSSIAN_PROJECT_PARAMETERS,
    "chemsmart.jobs.gaussian.settings",
    (
        "GaussianJobSettings",
        "GaussianIRCJobSettings",
        "GaussianTDDFTJobSettings",
        "GaussianLinkJobSettings",
    ),
)

_PYSCF_PROJECT_PARAMETERS = (
    "ab_initio",
    "aux_basis",
    "basis",
    # A boolean like ``scf_stability``: the broken-symmetry open-shell
    # singlet, reached by following the restricted solution's own
    # RHF/RKS -> UHF/UKS instability (R10 Q18).
    "broken_symmetry",
    "cc_max_cycle",
    "defgrid",
    "density_fit",
    "dispersion",
    "excited_state_root",
    # A Hessian's derivative and its displacement: the loader applied
    # both from the day they existed, and neither was advertised, so the
    # validation receipt the review renders left them out and a human
    # approved a Hessian without seeing how it would be computed.
    "fd_step_angstrom",
    "freq",
    "frozen_core",
    "functional",
    "hessian_derivative",
    # The branch an irc node walks from its saddle; forward and backward
    # are a sign the host fixes on the transition vector, so two nodes on
    # one geometry walk opposite branches.
    "irc_direction",
    "nstates",
    "opt_maxsteps",
    "opt_solver",
    "response_method",
    "scf_maxiter",
    # A boolean with no enumerable domain, like ``density_fit``: it asks
    # PySCF whether the converged reference is a minimum in
    # orbital-rotation space, and the answer is recorded as an
    # observation that refuses nothing.
    "scf_stability",
    "scf_tol",
    "solvent_id",
    "solvent_model",
    "state_manifold",
    "td_max_cycle",
)


def _pyscf_parameter_domains() -> tuple[tuple[str, tuple[str, ...]], ...]:
    """PySCF's enumerable domains, imported from the settings vocabulary.

    The method, manifold and frozen-core words are the settings module's
    own tuples, so the declaration the model reads cannot drift from what
    ``validate()`` admits; ``frozen_core`` also takes an orbital count,
    which no enumeration can list.
    """

    from chemsmart.jobs.pyscf.settings import (
        PYSCF_AB_INITIO_METHODS,
        PYSCF_DEFGRIDS,
        PYSCF_FROZEN_CORE_AUTO,
        PYSCF_HESSIAN_DERIVATIVES,
        PYSCF_IRC_DIRECTIONS,
        PYSCF_OPT_SOLVERS,
        PYSCF_RESPONSE_METHODS,
        PYSCF_SOLVENT_MODELS,
        PYSCF_STATE_MANIFOLDS,
    )

    return tuple(
        sorted(
            (
                ("ab_initio", tuple(sorted(PYSCF_AB_INITIO_METHODS))),
                ("defgrid", tuple(sorted(PYSCF_DEFGRIDS))),
                ("frozen_core", (PYSCF_FROZEN_CORE_AUTO,)),
                (
                    "hessian_derivative",
                    tuple(sorted(PYSCF_HESSIAN_DERIVATIVES)),
                ),
                ("irc_direction", tuple(sorted(PYSCF_IRC_DIRECTIONS))),
                ("opt_solver", tuple(sorted(PYSCF_OPT_SOLVERS))),
                ("response_method", tuple(sorted(PYSCF_RESPONSE_METHODS))),
                ("solvent_model", tuple(sorted(PYSCF_SOLVENT_MODELS))),
                ("state_manifold", tuple(sorted(PYSCF_STATE_MANIFOLDS))),
            )
        )
    )


_XTB_PROJECT_PARAMETERS = (
    "charge",
    "gfn_version",
    "grad",
    "jobtype",
    "multiplicity",
    "optimization_level",
    "solvent_id",
    "solvent_model",
)


def _normalized_domain(values) -> tuple[str, ...]:
    """Lower-case, strip, dedupe and sort one vocabulary for a domain.

    Domains carry the invariant that values are sorted, unique and
    lower-case; the io tables carry mixed case and, historically, one entry
    with trailing whitespace, so the projection normalizes rather than
    trusting.
    """

    return tuple(
        sorted(
            {
                str(value).strip().lower()
                for value in values
                if str(value).strip()
            }
        )
    )


def _orca_opt_convergence_words() -> tuple[str, ...]:
    """The words ORCA's own convergence table accepts, read from it.

    A declared domain that is re-typed is a second source of truth; this
    reads the table the writer uses, so adding a preset there publishes
    it to the model in the same edit.
    """

    from chemsmart.jobs.orca.settings import (
        ORCA_OPT_CONVERGENCE_KEYWORDS,
    )

    return tuple(sorted(ORCA_OPT_CONVERGENCE_KEYWORDS))


def orca_response_domains() -> tuple[tuple[str, tuple[str, ...]], ...]:
    """What an ORCA td stage may be asked for, from the writer's tables.

    The response and manifold words are the ones every program's td stage
    takes; the tuples read here are the ones the writer checks and spells,
    so the declaration cannot name a word the writer does not translate
    (it had re-typed the manifold pair and would have kept doing so).
    """

    from chemsmart.jobs.orca.settings import (
        ORCA_TD_RESPONSE_METHODS,
        ORCA_TD_STATE_MANIFOLDS,
    )

    return (
        ("response_method", tuple(sorted(ORCA_TD_RESPONSE_METHODS))),
        ("state_manifold", tuple(sorted(ORCA_TD_STATE_MANIFOLDS))),
    )


def orca_method_domains() -> tuple[tuple[str, tuple[str, ...]], ...]:
    """Method vocabulary projected from the single-source ORCA io tables.

    Declaring these domains up front lets a planner learn that a keyword
    is outside ORCA's vocabulary before a whole stage is built, instead of
    at the safe preview's membership test. The tables are probe-verified
    against an installed ORCA binary, and the io module stays the single
    source, so the route parser, the preview verifier, and this receipt
    can never fork.
    """

    from chemsmart.io.orca import (
        ORCA_ALL_AUXILIARY_BASIS_SETS,
        ORCA_ALL_BASIS_SETS,
        ORCA_ALL_FUNCTIONALS,
        ORCA_ALL_SOLVENT_MODELS,
        ORCA_ALL_SOLVENTS,
    )
    from chemsmart.jobs.orca.settings import orca_functional_literal

    return (
        ("aux_basis", _normalized_domain(ORCA_ALL_AUXILIARY_BASIS_SETS)),
        ("basis", _normalized_domain(ORCA_ALL_BASIS_SETS)),
        # The literals ORCA's keywords apply (its B3LYP is b3lyp5, B3LYP/G
        # is b3lyp, BP86 is bp86-pw92), read through the writer's own table,
        # so the domain never offers a word the writer would spell as
        # another functional.
        (
            "functional",
            _normalized_domain(
                orca_functional_literal(word) for word in ORCA_ALL_FUNCTIONALS
            ),
        ),
        ("solvent_id", _normalized_domain(ORCA_ALL_SOLVENTS)),
        ("solvent_model", _normalized_domain(ORCA_ALL_SOLVENT_MODELS)),
    )


def gaussian_method_domains() -> tuple[tuple[str, tuple[str, ...]], ...]:
    """Gaussian vocabulary projected from the io tables.

    The functional list is manually curated from the Gaussian 16 Rev C.01
    keyword documentation -- the job-submission hold on this host forbids the
    binary probe ORCA received -- so its provenance is the documentation, not
    the binary. Dispersion covers both the native EmpiricalDispersion tokens
    and the aliases the loader normalizes.
    """

    from chemsmart.io.gaussian import (
        GAUSSIAN_ALL_FUNCTIONALS,
        GAUSSIAN_SOLVATION_MODELS,
    )
    from chemsmart.io.gaussian.route import GAUSSIAN_EMPIRICAL_DISPERSIONS

    dispersion = set(GAUSSIAN_EMPIRICAL_DISPERSIONS) | {
        "d2",
        "d3",
        "d3bj",
        "d3zero",
    }
    from chemsmart.jobs.gaussian.settings import gaussian_functional_literal

    return (
        ("dispersion", _normalized_domain(dispersion)),
        # The literals Gaussian's keywords apply (PBE1PBE is pbe0), read
        # through the writer's own table.
        (
            "functional",
            _normalized_domain(
                gaussian_functional_literal(word)
                for word in GAUSSIAN_ALL_FUNCTIONALS
            ),
        ),
        ("solvent_model", _normalized_domain(GAUSSIAN_SOLVATION_MODELS)),
    )


def gaussian_numerics_domains() -> tuple[tuple[str, tuple[str, ...]], ...]:
    """Gaussian's words for its SCF convergence and integration grid.

    Read from the tables the Gaussian writer spells them from and the
    route reader reads them back with, so the declaration offers no word
    the writer does not write. The grid words are Gaussian's own: ORCA's
    DEFGRID words and PySCF's are other quadratures and stay theirs.
    """

    from chemsmart.jobs.gaussian.settings import (
        GAUSSIAN_INTEGRATION_GRIDS,
        GAUSSIAN_SCF_CONVERGENCE,
    )

    return (
        ("defgrid", tuple(sorted(GAUSSIAN_INTEGRATION_GRIDS))),
        ("scf_convergence", tuple(sorted(GAUSSIAN_SCF_CONVERGENCE))),
    )


def gaussian_response_domains() -> tuple[tuple[str, tuple[str, ...]], ...]:
    """What a Gaussian td stage may be asked for, from the writer's tables.

    The words are ORCA's and PySCF's (``response_method``,
    ``state_manifold``); the table that spells each in Gaussian's grammar
    is the one read here, so the declaration cannot name a word the
    writer does not translate.
    """

    from chemsmart.jobs.gaussian.settings import (
        GAUSSIAN_TD_MANIFOLD_OPTIONS,
        GAUSSIAN_TD_RESPONSE_KEYWORDS,
    )

    return (
        ("response_method", tuple(sorted(GAUSSIAN_TD_RESPONSE_KEYWORDS))),
        ("state_manifold", tuple(sorted(GAUSSIAN_TD_MANIFOLD_OPTIONS))),
    )


def gaussian_path_domains() -> tuple[tuple[str, tuple[str, ...]], ...]:
    """What a Gaussian reaction path may be asked for, from the CLI itself.

    ``direction``, ``predictor`` and ``recorrect`` are ``click.Choice``
    options of ``run gaussian irc``, so the public command is already the
    authority on their vocabulary and this reads it rather than repeating
    it. Declaring them is what wires the settings: a setting advertised
    with no domain is a word the model may write and cannot check.

    The integer controls of the same job -- ``maxpoints``, ``maxcycles``,
    ``recalc_step``, ``stepsize`` -- and the boolean ``flat_irc`` stay
    undeclared here, because this table holds tuples of strings and a
    bound is not one.
    """

    from chemsmart.cli.gaussian.irc import irc

    domains = {}
    for parameter in irc.params:
        choices = getattr(getattr(parameter, "type", None), "choices", None)
        if choices:
            domains[parameter.name] = _normalized_domain(choices)
    return tuple(sorted(domains.items()))


def xtb_solvent_domains() -> tuple[tuple[str, tuple[str, ...]], ...]:
    """The xTB solvent vocabulary, previously invisible to the model."""

    from chemsmart.io.xtb import XTB_ALL_SOLVENT_IDS

    return (("solvent_id", _normalized_domain(XTB_ALL_SOLVENT_IDS)),)


def declared_jobtypes(program: str) -> tuple[str, ...]:
    """Job-type vocabulary projected from each program's own settings.

    The registry describes what a program can be asked for; the words
    themselves belong to the settings module that validates them, so a
    job type added there reaches the capability ladder, the conformance
    probe and the agent surface without a second list to remember.
    """

    if program == "pyscf":
        from chemsmart.jobs.pyscf.settings import PYSCF_JOBTYPES

        return tuple(sorted(PYSCF_JOBTYPES))
    if program == "xtb":
        from chemsmart.jobs.xtb.settings import XTBJobSettings

        return tuple(sorted(XTBJobSettings.JOBTYPES))
    return ()


def loader_project_section_names(program: str) -> tuple[str, ...]:
    """Project-section vocabulary projected from each concrete loader."""

    if program in {"gaussian", "orca"}:
        from chemsmart.jobs.settings import molecular_project_section_names

        return molecular_project_section_names(program)
    if program == "pyscf":
        from chemsmart.settings.pyscf import PYSCF_ALLOWED_SECTIONS

        return tuple(sorted(PYSCF_ALLOWED_SECTIONS))
    if program == "xtb":
        from chemsmart.settings.xtb import YamlXTBProjectSettingsBuilder

        return tuple(sorted(YamlXTBProjectSettingsBuilder.SECTIONS))
    return ()


PROGRAM_CAPABILITIES: Mapping[str, ProgramCapability] = MappingProxyType(
    {
        "gaussian": ProgramCapability(
            program="gaussian",
            coordinate_idiom="increment_steps",
            requires_project_configuration=True,
            supports_project_configuration=True,
            jobtypes=(
                "com",
                "crest",
                "dias",
                "irc",
                "link",
                "modred",
                "nci",
                "opt",
                "pka",
                "qrc",
                "resp",
                "scan",
                "sp",
                "td",
                "traj",
                "ts",
                "userjob",
                "wbi",
            ),
            project_owned_parameters=_GAUSSIAN_PROJECT_PARAMETERS,
            engines=("cpu",),
            engine_job_capabilities=(
                # A reaction path, walked from a saddle this Agent
                # found. Declared here because this flag is what admits a
                # node to approval at all, which is the order `ts` and
                # ORCA's `scan` were held to; the recorded runs keep it.
                EngineJobCapability(engine="cpu", jobtype="irc"),
                EngineJobCapability(
                    engine="cpu",
                    jobtype="link",
                    execution_supported=False,
                ),
                # A constrained optimisation, a relaxed scan and a
                # fixed-geometry response calculation are admitted to
                # approval for their first Agent runs (R10 Q7), the order
                # `ts`, `irc` and ORCA's `scan` were held to: this flag is
                # what admits a node to approval at all. Until those runs
                # are recorded in release.json they are admitted, not
                # qualified, and the flags are withdrawn if the runs do not
                # hold. Each had real engine runs through the human CLI on
                # CUHK (Slurm 2142374/2142393); what the Agent could not do
                # before was preview a scan at all (cb1acd9a, f412dac0).
                EngineJobCapability(engine="cpu", jobtype="modred"),
                # Qualified by real approved runs rather than by inspection.
                #
                # Gaussian 16 C.02 had never been driven through ChemSmart on
                # any target. Slurm jobs 2142374 and 2142393 on CUHK Charles
                # ran fifteen small jobs over every declared job type through
                # the human CLI: a distorted hydrogen peroxide relaxed to
                # O-O 1.4557 A (experiment 1.452) with six real modes, the
                # HCN/HNC saddle carried one imaginary mode at -1146.1 cm^-1
                # whose IRC branches reach HNC and HCN.
                #
                # The Agent then ran these two under the approval chain on two
                # chemically different molecules, each goal an `opt` feeding an
                # `sp` across a geometry edge, provider-free executor, both
                # nodes validated. See the release records for the run
                # evidence. The perturbation under which they reproduce is the
                # molecule: hydrogen peroxide and formaldehyde, distorted
                # differently, relaxed to different point groups.
                #
                # The declaration necessarily precedes the first approved
                # execution, because it is what admits a node to approval at
                # all -- the same order ORCA's scan was held to, and a goal
                # issued with these false settles `returned_to_human` with
                # "bounded execution has no executable jobs" before any
                # planning. It is withdrawn if the runs do not hold.
                #
                # `link` has real engine runs on this target through the
                # human CLI and no approved Agent execution. That is a
                # different fact and stays unclaimed.
                EngineJobCapability(engine="cpu", jobtype="opt"),
                EngineJobCapability(engine="cpu", jobtype="scan"),
                EngineJobCapability(engine="cpu", jobtype="sp"),
                EngineJobCapability(engine="cpu", jobtype="td"),
                EngineJobCapability(engine="cpu", jobtype="ts"),
            ),
            project_section_names=loader_project_section_names("gaussian"),
            project_parameter_domains=tuple(
                sorted(
                    (
                        ("states", ("50-50", "singlets", "triplets")),
                        *gaussian_method_domains(),
                        *gaussian_numerics_domains(),
                        *gaussian_path_domains(),
                        *gaussian_response_domains(),
                    )
                )
            ),
        ),
        "nciplot": ProgramCapability(
            program="nciplot",
            requires_project_configuration=False,
            supports_project_configuration=False,
            jobtypes=(),
            project_owned_parameters=(),
            engines=("cpu",),
            agent_enabled=False,
        ),
        "orca": ProgramCapability(
            program="orca",
            coordinate_idiom="absolute_range",
            requires_project_configuration=True,
            supports_project_configuration=True,
            jobtypes=(
                "inp",
                "irc",
                "modred",
                "neb",
                "opt",
                "pka",
                "qrc",
                "scan",
                "sp",
                "td",
                "ts",
            ),
            project_owned_parameters=_ORCA_PROJECT_PARAMETERS,
            engines=("cpu",),
            engine_job_capabilities=(
                # Qualified by a real run rather than by inspection: an
                # HCN -> HNC 1,2-hydrogen shift executed on this host through
                # the ChemSmart CLI at B3LYP/def2-SVP.  The transition state
                # carried exactly one imaginary mode (-1121.7 cm^-1) and its
                # analytic Hessian fed the IRC through `%irc InitHess read`;
                # the 41-frame path ORCA wrote to its own _IRC_Full_trj.xyz
                # sidecar starts at H-C 1.078 / H-N 2.202 angstrom and ends at
                # H-C 2.142 / H-N 1.005, so the saddle demonstrably connects
                # the two minima and the connectivity change is observable
                # rather than asserted.
                #
                # The approved Agent execution has since held: a TS -> IRC
                # workflow (ts feeding two irc nodes, each consuming the
                # converged saddle's geometry and Hessian as role-distinct
                # producer bindings) planned, reviewed, executed
                # provider-free, validated, and delivered claims
                # (qualification/irc-agent-path, 2026-08-23; -1121.05 cm^-1,
                # exactly one imaginary mode, wB97X-D3(BJ)/def2-TZVP).
                # `neb` stays preview-only: it has had no such run. `modred`
                # had one (below).
                EngineJobCapability(engine="cpu", jobtype="irc"),
                EngineJobCapability(
                    engine="cpu",
                    jobtype="modred",
                    # Constrained optimisation was declared for planning and
                    # preview first; execution is recorded from the run that
                    # exercised it under the approval chain (R9 r9o-g5, CUHK
                    # Slurm 2144929, 2026-09-23): one modred node held a
                    # BINOL biaryl torsion at 0.00 degrees, terminated
                    # normally at ORCA's cycle limit (3N = 108) without
                    # converging, was verified and extracted, and handed its
                    # structure across the producer edge to a saddle search
                    # that converged. The owner ruled that run a qualifying
                    # success; the record in release.json says exactly what
                    # ran, including that it did not converge.
                    execution_supported=True,
                ),
                EngineJobCapability(
                    engine="cpu",
                    jobtype="neb",
                    execution_supported=False,
                ),
                EngineJobCapability(engine="cpu", jobtype="opt"),
                # Qualified by a real run rather than by inspection: an ORCA
                # 6.1.1 relaxed scan executed on this host through the
                # ChemSmart CLI and its seven converged points match the
                # .relaxscanact.dat sidecar ORCA wrote beside them; the
                # Agent's own compiled invocation reproduces that native input
                # and its safe preview is green.
                #
                # The declaration necessarily precedes the first approved
                # Agent execution, because it is what admits a node to
                # approval at all -- leaving it false until an approved run
                # exists would keep the family unqualifiable for ever. It is
                # withdrawn if that run does not hold.
                EngineJobCapability(engine="cpu", jobtype="scan"),
                EngineJobCapability(engine="cpu", jobtype="sp"),
                EngineJobCapability(engine="cpu", jobtype="td"),
                EngineJobCapability(engine="cpu", jobtype="ts"),
            ),
            project_section_names=loader_project_section_names("orca"),
            # The literal domains below are merged with the probe-verified
            # method vocabulary (orca_method_domains) and sorted, because the
            # contract requires domain names in order and the method names --
            # functional, basis, aux_basis, solvent_* -- interleave with them.
            project_parameter_domains=tuple(
                sorted(
                    (
                        *orca_method_domains(),
                        # Derived from the writer's own keyword table
                        # rather than re-typed, so the declaration
                        # cannot drift from what the route accepts.
                        #
                        # This was the one code-level ORCA enum with no
                        # declared domain: `ORCA_OPT_CONVERGENCE_KEYWORDS`
                        # has existed in jobs/orca/settings.py while the
                        # model was told only that `opt_convergence` is
                        # settable, never what it may be set to. po3-r18
                        # declared `tight`, the validator reported
                        # `expected 'tight', observed None` because the
                        # input *reader* has no property for it, and the
                        # session then deleted the control the host's own
                        # rule had recommended. Publishing the domain is
                        # FUNDAMENTAL 1 directly: the model should not
                        # have to memorise a program's vocabulary.
                        (
                            "opt_convergence",
                            tuple(sorted(_orca_opt_convergence_words())),
                        ),
                        (
                            "ab_initio",
                            (
                                "ccsd(t)",
                                "dlpno-ccsd",
                                "dlpno-ccsd(t)",
                                "hf",
                                "mp2",
                                "rhf",
                                "uhf",
                            ),
                        ),
                        (
                            "defgrid",
                            (
                                "defgrid1",
                                "defgrid2",
                                "defgrid3",
                                "grid1",
                                "grid2",
                                "grid3",
                                "grid4",
                                "grid5",
                                "grid6",
                                "grid7",
                            ),
                        ),
                        ("direction", ("backward", "both", "down", "forward")),
                        ("dispersion", ("d2", "d3bj", "d3zero", "d4")),
                        (
                            "frozen_core",
                            ("fc_electrons", "fc_ewin", "fc_none"),
                        ),
                        ("mdci_cutoff", ("loose", "normal", "tight")),
                        ("reference", ("rhf", "rohf", "uhf")),
                        ("relativistic", ("dkh", "dkh2", "zora")),
                        (
                            "ri_approximation",
                            ("none", "ri", "rijcosx", "rijk"),
                        ),
                        *orca_response_domains(),
                    )
                )
            ),
        ),
        "pyscf": ProgramCapability(
            program="pyscf",
            requires_project_configuration=True,
            supports_project_configuration=True,
            jobtypes=declared_jobtypes("pyscf"),
            project_owned_parameters=_PYSCF_PROJECT_PARAMETERS,
            engines=("cpu", "gpu"),
            # PySCF keys sections by job type, and its loader also
            # accepts the legacy gas/solv pair and canonicalises it.
            project_section_names=loader_project_section_names("pyscf"),
            project_parameter_domains=_pyscf_parameter_domains(),
            engine_job_capabilities=(
                EngineJobCapability(engine="cpu", jobtype="hess"),
                # One branch of the intrinsic reaction coordinate from a
                # supplied saddle, walked by geomeTRIC 1.1.1 on PySCF's own
                # engine: the start's analytic Hessian on the walked
                # surface, the transition vector, the whole accepted path
                # and the endpoint every property belongs to are in the
                # artifact (result contract v8). Qualified first through
                # the human CLI on CUHK from ORCA OptTS saddles at matched
                # levels; no GPU row, because GPU4PySCF has run no IRC.
                EngineJobCapability(engine="cpu", jobtype="irc"),
                EngineJobCapability(engine="cpu", jobtype="opt"),
                EngineJobCapability(engine="cpu", jobtype="sp"),
                EngineJobCapability(engine="cpu", jobtype="td"),
                # A saddle search on PySCF's own surface (result contract
                # v9): geomeTRIC's partitioned rational-function step on
                # PySCF's engine, seeded with PySCF's analytic Hessian
                # there, so a PySCF IRC no longer needs a saddle located
                # by another program. Qualified through the human CLI on
                # CUHK: from a seed 0.24 A away it reproduced ORCA's own
                # HF/6-31G* H2CO/trans-HCOH saddle to 0.0001 A in every
                # interatomic distance, the Hessian at what it reached
                # carries the one imaginary mode the archived IRC
                # fixtures start from, and both branches walked from it
                # reach the two minima. No GPU row: GPU4PySCF has run no
                # saddle search.
                EngineJobCapability(engine="cpu", jobtype="ts"),
                EngineJobCapability(
                    engine="gpu",
                    jobtype="hess",
                    execution_supported=False,
                ),
                EngineJobCapability(
                    engine="gpu",
                    jobtype="opt",
                    execution_supported=False,
                ),
                EngineJobCapability(
                    engine="gpu",
                    jobtype="sp",
                    execution_supported=False,
                ),
            ),
        ),
        "xtb": ProgramCapability(
            program="xtb",
            requires_project_configuration=False,
            supports_project_configuration=True,
            jobtypes=declared_jobtypes("xtb"),
            project_owned_parameters=_XTB_PROJECT_PARAMETERS,
            engines=("cpu",),
            project_section_names=loader_project_section_names("xtb"),
            project_parameter_domains=(
                ("gfn_version", ("gfn0", "gfn1", "gfn2", "gfnff")),
                (
                    "optimization_level",
                    (
                        "crude",
                        "extreme",
                        "lax",
                        "loose",
                        "normal",
                        "sloppy",
                        "tight",
                        "vtight",
                    ),
                ),
                xtb_solvent_domains()[0],
                (
                    "solvent_model",
                    ("alpb", "cosmo", "cpcmx", "gbsa", "tmcosmo"),
                ),
            ),
            # xTB runs every declared job type on the one engine it has,
            # so its matrix is the product -- but it is declared rather
            # than left to the legacy Cartesian fallback, because the
            # conformance probe reads the matrix and fell back to a
            # hand-written core-stage set for any program that had none.
            engine_job_capabilities=tuple(
                EngineJobCapability(engine="cpu", jobtype=jobtype)
                for jobtype in declared_jobtypes("xtb")
            ),
        ),
    }
)

# Compatibility name for the current harness module.  This is an alias, not a
# second registry, so consumers cannot drift from PROGRAM_CAPABILITIES.
ENGINE_CAPABILITIES = PROGRAM_CAPABILITIES

# Immutable projections for the hard-coded sets and maps in the current
# harness. ``EXECUTABLE_PROGRAMS`` is the top-level executable Click inventory;
# ``COMPUTATIONAL_PROGRAMS`` is the primary molecular-method set and excludes
# the advanced NCIPLOT leaf.
KNOWN_PROGRAMS = frozenset(PROGRAM_CAPABILITIES)
EXECUTABLE_PROGRAMS = KNOWN_PROGRAMS
AGENT_PROGRAMS = frozenset(
    name
    for name, capability in PROGRAM_CAPABILITIES.items()
    if capability.agent_enabled and capability.resolved_engine_job_capabilities
)
PROJECT_PROGRAMS = frozenset(
    name
    for name, capability in PROGRAM_CAPABILITIES.items()
    if capability.supports_project_configuration
)
PROJECT_REQUIRED_PROGRAMS = frozenset(
    name
    for name, capability in PROGRAM_CAPABILITIES.items()
    if capability.requires_project_configuration
)
COMPUTATIONAL_PROGRAMS = PROJECT_PROGRAMS
PRIMARY_PROGRAMS = COMPUTATIONAL_PROGRAMS
PROJECT_PROGRAM_ORDER = tuple(sorted(PROJECT_PROGRAMS))

# These are direct Click child inventories, not a promise that an agent has a
# validated task model, concrete JobKind class, or renderer for every leaf.
PROGRAM_CLI_JOBTYPES: Mapping[str, tuple[str, ...]] = MappingProxyType(
    {
        name: capability.jobtypes
        for name, capability in PROGRAM_CAPABILITIES.items()
    }
)
PROGRAM_EXECUTION_ENGINES: Mapping[str, tuple[str, ...]] = MappingProxyType(
    {
        name: capability.engines
        for name, capability in PROGRAM_CAPABILITIES.items()
    }
)
AGENT_PROGRAM_PREVIEW_ENGINES: Mapping[str, tuple[str, ...]] = (
    MappingProxyType(
        {
            name: tuple(
                sorted(
                    {
                        item.engine
                        for item in capability.resolved_engine_job_capabilities
                        if item.preview_supported
                    }
                )
            )
            for name, capability in PROGRAM_CAPABILITIES.items()
            if name in AGENT_PROGRAMS
        }
    )
)
AGENT_PROGRAM_JOBTYPES: Mapping[str, tuple[str, ...]] = MappingProxyType(
    {
        name: tuple(
            sorted(
                {
                    item.jobtype
                    for item in capability.resolved_engine_job_capabilities
                    if item.preview_supported
                }
            )
        )
        for name, capability in PROGRAM_CAPABILITIES.items()
        if name in AGENT_PROGRAMS
    }
)
PROGRAM_PROJECT_OWNED_CLI_PARAMETERS: Mapping[str, tuple[str, ...]] = (
    MappingProxyType(
        {
            name: capability.project_owned_parameters
            for name, capability in PROGRAM_CAPABILITIES.items()
        }
    )
)

# Compatibility aliases for the names used by the current v2 harness. They
# remain references to the canonical immutable views above, not copies.
EngineCapability = ProgramCapability
PROGRAM_JOBTYPES = PROGRAM_CLI_JOBTYPES
PROGRAM_ENGINES = PROGRAM_EXECUTION_ENGINES
PROJECT_OWNED_PARAMETERS = PROGRAM_PROJECT_OWNED_CLI_PARAMETERS


def program_capability(program: str | None) -> ProgramCapability | None:
    """Return a capability after harmless case and whitespace normalisation."""

    return PROGRAM_CAPABILITIES.get(str(program or "").strip().lower())


def engine_capability(program: str | None) -> ProgramCapability | None:
    """Compatibility lookup for the current harness capability API."""

    return program_capability(program)


def requires_project_configuration(program: str | None) -> bool:
    """Return whether synthesis must resolve an approved project artifact."""

    capability = program_capability(program)
    return bool(
        capability is not None and capability.requires_project_configuration
    )


def supports_project_configuration(program: str | None) -> bool:
    """Return whether the program consumes project-YAML settings."""

    capability = program_capability(program)
    return bool(
        capability is not None and capability.supports_project_configuration
    )


def project_owns_parameter(program: str | None, parameter: str) -> bool:
    """Return whether an option value must come from approved project YAML."""

    capability = program_capability(program)
    normalised = str(parameter or "").strip().lower().replace("-", "_")
    return bool(
        capability is not None
        and normalised in capability.project_owned_parameters
    )


__all__ = [
    "AGENT_PROGRAM_PREVIEW_ENGINES",
    "AGENT_PROGRAM_JOBTYPES",
    "AGENT_PROGRAMS",
    "COMPUTATIONAL_PROGRAMS",
    "ENGINE_CAPABILITIES",
    "EXECUTABLE_PROGRAMS",
    "KNOWN_PROGRAMS",
    "PRIMARY_PROGRAMS",
    "PROGRAM_CAPABILITIES",
    "PROGRAM_CLI_JOBTYPES",
    "PROGRAM_ENGINES",
    "PROGRAM_EXECUTION_ENGINES",
    "PROGRAM_JOBTYPES",
    "PROGRAM_PROJECT_OWNED_CLI_PARAMETERS",
    "PROJECT_OWNED_PARAMETERS",
    "PROJECT_PROGRAMS",
    "PROJECT_PROGRAM_ORDER",
    "PROJECT_REQUIRED_PROGRAMS",
    "ProgramCapability",
    "EngineJobCapability",
    "EngineCapability",
    "engine_capability",
    "program_capability",
    "project_owns_parameter",
    "loader_project_section_names",
    "requires_project_configuration",
    "supports_project_configuration",
]
