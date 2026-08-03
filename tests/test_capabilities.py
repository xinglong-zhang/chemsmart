from dataclasses import FrozenInstanceError

import click
import pytest

from chemsmart.cli.gaussian import gaussian
from chemsmart.cli.nciplot import nciplot
from chemsmart.cli.orca import orca
from chemsmart.cli.pyscf import pyscf
from chemsmart.cli.subcommands import subcommands
from chemsmart.cli.xtb import xtb
from chemsmart.settings.capabilities import (
    COMPUTATIONAL_PROGRAMS,
    ENGINE_CAPABILITIES,
    EXECUTABLE_PROGRAMS,
    KNOWN_PROGRAMS,
    PRIMARY_PROGRAMS,
    PROGRAM_CAPABILITIES,
    PROGRAM_CLI_JOBTYPES,
    PROGRAM_ENGINES,
    PROGRAM_EXECUTION_ENGINES,
    PROGRAM_JOBTYPES,
    PROGRAM_PROJECT_OWNED_CLI_PARAMETERS,
    PROJECT_OWNED_PARAMETERS,
    PROJECT_PROGRAM_ORDER,
    PROJECT_PROGRAMS,
    PROJECT_REQUIRED_PROGRAMS,
    EngineCapability,
    ProgramCapability,
    engine_capability,
    program_capability,
    project_owns_parameter,
    requires_project_configuration,
    supports_project_configuration,
)

PROGRAM_COMMANDS = {
    "gaussian": gaussian,
    "nciplot": nciplot,
    "orca": orca,
    "pyscf": pyscf,
    "xtb": xtb,
}

CURRENT_HARNESS_PROJECT_PARAMETERS = (
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


def _click_jobtypes(command: click.Command) -> tuple[str, ...]:
    if isinstance(command, click.Group):
        return tuple(sorted(command.commands))
    return ()


def test_registry_agrees_with_live_click_leaf_inventory():
    assert frozenset(PROGRAM_COMMANDS) == KNOWN_PROGRAMS
    assert {command.name for command in subcommands}.intersection(
        KNOWN_PROGRAMS
    ) == KNOWN_PROGRAMS

    live_jobtypes = {
        name: _click_jobtypes(command)
        for name, command in PROGRAM_COMMANDS.items()
    }
    assert dict(PROGRAM_JOBTYPES) == live_jobtypes


def test_registry_and_all_derived_views_are_immutable():
    capability = PROGRAM_CAPABILITIES["pyscf"]

    with pytest.raises(FrozenInstanceError):
        setattr(capability, "program", "other")
    with pytest.raises(TypeError):
        PROGRAM_CAPABILITIES["other"] = capability
    with pytest.raises(TypeError):
        PROGRAM_JOBTYPES["pyscf"] = ("sp",)
    with pytest.raises(TypeError):
        PROGRAM_ENGINES["pyscf"] = ("cpu",)
    with pytest.raises(TypeError):
        PROJECT_OWNED_PARAMETERS["pyscf"] = ()


def test_derived_views_are_exact_registry_projections():
    assert ENGINE_CAPABILITIES is PROGRAM_CAPABILITIES
    assert EngineCapability is ProgramCapability
    assert EXECUTABLE_PROGRAMS is KNOWN_PROGRAMS
    assert PROJECT_PROGRAMS == frozenset(
        {"gaussian", "orca", "pyscf", "xtb"}
    )
    assert PROJECT_REQUIRED_PROGRAMS == frozenset(
        {"gaussian", "orca", "pyscf"}
    )
    assert COMPUTATIONAL_PROGRAMS is PROJECT_PROGRAMS
    assert PRIMARY_PROGRAMS is COMPUTATIONAL_PROGRAMS
    assert PROJECT_PROGRAM_ORDER == ("gaussian", "orca", "pyscf", "xtb")
    assert PROGRAM_JOBTYPES is PROGRAM_CLI_JOBTYPES
    assert PROGRAM_ENGINES is PROGRAM_EXECUTION_ENGINES
    assert PROJECT_OWNED_PARAMETERS is PROGRAM_PROJECT_OWNED_CLI_PARAMETERS
    assert dict(PROGRAM_JOBTYPES) == {
        name: capability.jobtypes
        for name, capability in PROGRAM_CAPABILITIES.items()
    }
    assert dict(PROGRAM_ENGINES) == {
        name: capability.engines
        for name, capability in PROGRAM_CAPABILITIES.items()
    }
    assert dict(PROJECT_OWNED_PARAMETERS) == {
        name: capability.project_owned_parameters
        for name, capability in PROGRAM_CAPABILITIES.items()
    }


def test_declared_capabilities_preserve_project_ownership_contract():
    assert PROJECT_OWNED_PARAMETERS["gaussian"] == (
        CURRENT_HARNESS_PROJECT_PARAMETERS
    )
    assert PROJECT_OWNED_PARAMETERS["orca"] == (
        CURRENT_HARNESS_PROJECT_PARAMETERS
    )
    assert PROJECT_OWNED_PARAMETERS["nciplot"] == ()
    assert PROJECT_OWNED_PARAMETERS["pyscf"] == (
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
    assert PROJECT_OWNED_PARAMETERS["xtb"] == (
        "charge",
        "gfn_version",
        "grad",
        "jobtype",
        "multiplicity",
        "optimization_level",
        "solvent_id",
        "solvent_model",
    )


def test_lookup_and_boolean_views_normalise_without_guessing():
    capability = PROGRAM_CAPABILITIES["pyscf"]
    assert program_capability(" PySCF ") is capability
    assert engine_capability("PYSCF") is capability
    assert program_capability(None) is None
    assert requires_project_configuration("pyscf")
    assert supports_project_configuration("orca")
    assert not requires_project_configuration("nciplot")
    assert not supports_project_configuration("unknown")
    assert project_owns_parameter("pyscf", "opt-solver")
    assert not project_owns_parameter("pyscf", "gpu")


@pytest.mark.parametrize("name, capability", PROGRAM_CAPABILITIES.items())
def test_registry_internal_invariants(name, capability):
    assert name == capability.program
    assert capability.jobtypes == tuple(sorted(set(capability.jobtypes)))
    assert capability.engines == tuple(sorted(set(capability.engines)))
    assert capability.project_owned_parameters == tuple(
        sorted(set(capability.project_owned_parameters))
    )
    assert not capability.requires_project_configuration or (
        capability.supports_project_configuration
    )


def test_program_capability_rejects_incoherent_declarations():
    with pytest.raises(ValueError, match="cannot require unsupported"):
        ProgramCapability(
            program="invalid",
            requires_project_configuration=True,
            supports_project_configuration=False,
            jobtypes=("sp",),
            project_owned_parameters=(),
            engines=("cpu",),
        )

    with pytest.raises(ValueError, match="sorted and contain no duplicates"):
        ProgramCapability(
            program="invalid",
            requires_project_configuration=False,
            supports_project_configuration=False,
            jobtypes=("sp", "sp"),
            project_owned_parameters=(),
            engines=("cpu",),
        )
