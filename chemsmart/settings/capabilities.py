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
            "state_manifold",
            "vpt2",
            "vpt2_anharmonic_displacement",
            "vpt2_hessian_cutoff",
        }
    )
)

# Gaussian exposes scientific controls the shared union omitted.
_GAUSSIAN_PROJECT_PARAMETERS = tuple(
    sorted(
        _CURRENT_HARNESS_PROJECT_PARAMETERS
        + (
            "additional_opt_options_in_route",
            "additional_solvent_options",
            "eqsolv",
            "forces",
            "freq",
            "heavy_elements_basis",
            "guess",
            "jobtype",
            "link_route",
            "nstates",
            "numfreq",
            "root",
            "stable",
            "states",
        )
    )
)

_PYSCF_PROJECT_PARAMETERS = (
    "ab_initio",
    "aux_basis",
    "basis",
    "defgrid",
    "density_fit",
    "dispersion",
    "freq",
    "functional",
    "nstates",
    "opt_maxsteps",
    "opt_solver",
    "response_method",
    "scf_maxiter",
    "scf_tol",
    "solvent_id",
    "solvent_model",
    "state_manifold",
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

    return tuple(sorted({str(value).strip().lower() for value in values if str(value).strip()}))


def orca_method_domains() -> tuple[tuple[str, tuple[str, ...]], ...]:
    """Method vocabulary projected from the single-source ORCA io tables.

    A live session planned MN15 for ORCA; the YAML validated, the CLI
    compiled, the native input was generated, and only the safe preview's
    membership test against these very tables refused it -- after the whole
    stage was built. The tables existed all along; nothing declared them as
    domains, so the receipt the model reads first carried every enumerated
    vocabulary *except* the ones a chemist actually chooses by.

    The tables are probe-verified against the installed ORCA 6.1.1 binary
    (campaign keyword-probe report): all 94 prior entries were probed, one
    legacy keyword the binary rejects was removed, and fourteen keywords the
    binary accepts were added. The io module stays the single source, so the
    route parser, the preview verifier, and this receipt can never fork.
    """

    from chemsmart.io.orca import (
        ORCA_ALL_AUXILIARY_BASIS_SETS,
        ORCA_ALL_BASIS_SETS,
        ORCA_ALL_FUNCTIONALS,
        ORCA_ALL_SOLVENT_MODELS,
        ORCA_ALL_SOLVENTS,
    )

    return (
        ("aux_basis", _normalized_domain(ORCA_ALL_AUXILIARY_BASIS_SETS)),
        ("basis", _normalized_domain(ORCA_ALL_BASIS_SETS)),
        ("functional", _normalized_domain(ORCA_ALL_FUNCTIONALS)),
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
    return (
        ("dispersion", _normalized_domain(dispersion)),
        ("functional", _normalized_domain(GAUSSIAN_ALL_FUNCTIONALS)),
        ("solvent_model", _normalized_domain(GAUSSIAN_SOLVATION_MODELS)),
    )


def xtb_solvent_domains() -> tuple[tuple[str, tuple[str, ...]], ...]:
    """The xTB solvent vocabulary, previously invisible to the model."""

    from chemsmart.io.xtb import XTB_ALL_SOLVENT_IDS

    return (("solvent_id", _normalized_domain(XTB_ALL_SOLVENT_IDS)),)


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
                EngineJobCapability(
                    engine="cpu",
                    jobtype="irc",
                    execution_supported=False,
                ),
                EngineJobCapability(
                    engine="cpu",
                    jobtype="link",
                    execution_supported=False,
                ),
                EngineJobCapability(
                    engine="cpu",
                    jobtype="modred",
                    execution_supported=False,
                ),
                EngineJobCapability(
                    engine="cpu",
                    jobtype="opt",
                    execution_supported=False,
                ),
                EngineJobCapability(
                    engine="cpu",
                    jobtype="scan",
                    execution_supported=False,
                ),
                EngineJobCapability(
                    engine="cpu",
                    jobtype="sp",
                    execution_supported=False,
                ),
                EngineJobCapability(
                    engine="cpu",
                    jobtype="td",
                    execution_supported=False,
                ),
                EngineJobCapability(
                    engine="cpu",
                    jobtype="ts",
                    execution_supported=False,
                ),
            ),
            project_section_names=loader_project_section_names("gaussian"),
            project_parameter_domains=tuple(sorted((
                ("states", ("50-50", "singlets", "triplets")),
                *gaussian_method_domains(),
            ))),
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
                EngineJobCapability(
                    engine="cpu",
                    jobtype="irc",
                    execution_supported=False,
                ),
                EngineJobCapability(
                    engine="cpu",
                    jobtype="modred",
                    # Constrained optimisation is declared for planning and
                    # preview first.  Execution is claimed only once a real
                    # run on this target has qualified it, which is the same
                    # order irc and neb are held to above.
                    execution_supported=False,
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
            project_parameter_domains=tuple(sorted((
                *orca_method_domains(),
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
                ("frozen_core", ("fc_electrons", "fc_ewin", "fc_none")),
                ("mdci_cutoff", ("loose", "normal", "tight")),
                ("reference", ("rhf", "rohf", "uhf")),
                ("relativistic", ("dkh", "dkh2", "zora")),
                ("response_method", ("tda", "tddft")),
                ("ri_approximation", ("none", "ri", "rijcosx", "rijk")),
                (
                    "state_manifold",
                    ("singlet", "singlet_triplet"),
                ),
            ))),
        ),
        "pyscf": ProgramCapability(
            program="pyscf",
            requires_project_configuration=True,
            supports_project_configuration=True,
            jobtypes=("hess", "opt", "sp", "td"),
            project_owned_parameters=_PYSCF_PROJECT_PARAMETERS,
            engines=("cpu", "gpu"),
            # PySCF keys sections by job type, and its loader also
            # accepts the legacy gas/solv pair and canonicalises it.
            project_section_names=loader_project_section_names("pyscf"),
            project_parameter_domains=(
                ("ab_initio", ("hf",)),
                ("defgrid", ("defgrid1", "defgrid2", "defgrid3")),
                ("opt_solver", ("ase", "berny", "geometric")),
                ("response_method", ("tda", "tddft")),
                (
                    "solvent_model",
                    ("cosmo", "cpcm", "iefpcm", "pcm", "smd", "ssvpe"),
                ),
                ("state_manifold", ("singlet",)),
            ),
            engine_job_capabilities=(
                EngineJobCapability(engine="cpu", jobtype="hess"),
                EngineJobCapability(engine="cpu", jobtype="opt"),
                EngineJobCapability(engine="cpu", jobtype="sp"),
                EngineJobCapability(
                    engine="cpu",
                    jobtype="td",
                    execution_supported=False,
                ),
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
            jobtypes=("hess", "opt", "sp"),
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
    if capability.agent_enabled
    and capability.resolved_engine_job_capabilities
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
AGENT_PROGRAM_EXECUTION_ENGINES: Mapping[str, tuple[str, ...]] = (
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
    "AGENT_PROGRAM_EXECUTION_ENGINES",
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
