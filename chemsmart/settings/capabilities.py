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


@dataclass(frozen=True)
class ProgramCapability:
    """Agent-operability facts for one executable chemistry program."""

    program: str
    requires_project_configuration: bool
    supports_project_configuration: bool
    jobtypes: tuple[str, ...]
    project_owned_parameters: tuple[str, ...]
    engines: tuple[str, ...]

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

_PYSCF_PROJECT_PARAMETERS = (
    "ab_initio",
    "aux_basis",
    "basis",
    "defgrid",
    "density_fit",
    "dispersion",
    "freq",
    "functional",
    "opt_maxsteps",
    "opt_solver",
    "scf_maxiter",
    "scf_tol",
    "solvent_id",
    "solvent_model",
)


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
            project_owned_parameters=_CURRENT_HARNESS_PROJECT_PARAMETERS,
            engines=("cpu",),
        ),
        "nciplot": ProgramCapability(
            program="nciplot",
            requires_project_configuration=False,
            supports_project_configuration=False,
            jobtypes=(),
            project_owned_parameters=(),
            engines=("cpu",),
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
                "ts",
            ),
            project_owned_parameters=_CURRENT_HARNESS_PROJECT_PARAMETERS,
            engines=("cpu",),
        ),
        "pyscf": ProgramCapability(
            program="pyscf",
            requires_project_configuration=True,
            supports_project_configuration=True,
            jobtypes=("hess", "opt", "sp"),
            project_owned_parameters=_PYSCF_PROJECT_PARAMETERS,
            engines=("cpu", "gpu"),
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
    "EngineCapability",
    "engine_capability",
    "program_capability",
    "project_owns_parameter",
    "requires_project_configuration",
    "supports_project_configuration",
]
