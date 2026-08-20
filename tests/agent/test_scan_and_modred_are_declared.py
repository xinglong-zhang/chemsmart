"""Coordinate scans and constrained optimisation reach the Agent's surface.

`scan` and `modred` are registered human CLI families on both route-building
programs and were declared for neither Agent, so every torsional profile and
every constrained optimisation in the recorded campaign needed a human to
hand-build the endpoint structures.

These tests exercise the live capability query the Agent actually issues, not
the declaration table, because the table is projected through a compatibility
view that derives the Agent's whole jobtype set from the narrower
engine-capability matrix — a family can be present in one and absent from the
other, which is precisely the state this change repairs.

Execution is deliberately *not* claimed here. A family becomes executable only
after a real run on this target qualifies it, which is the same order `irc` and
`neb` are already held to.
"""

from __future__ import annotations

import pytest

from chemsmart.agent.capabilities import load_program_capabilities
from chemsmart.agent.live_session import _conformance_jobtypes
from chemsmart.jobs.settings import molecular_project_section_names

_REGISTRY = load_program_capabilities()


def _pairs(program):
    item = next(
        entry for entry in _REGISTRY.programs if entry.program == program
    )
    return {
        capability.jobtype: capability
        for capability in item.engine_job_capabilities
    }


@pytest.mark.parametrize("program", ("orca", "gaussian"))
@pytest.mark.parametrize("jobtype", ("scan", "modred"))
def test_the_family_is_previewable_through_the_live_registry(program, jobtype):
    capability = _pairs(program).get(jobtype)

    assert capability is not None, (
        f"{program} {jobtype} is a registered human CLI family; the Agent "
        "cannot express it until the engine matrix declares it"
    )
    assert capability.preview_supported is True


#: The one pair a real run has qualified: an ORCA 6.1.1 relaxed scan executed
#: on this host through the ChemSmart CLI, whose converged points match the
#: .relaxscanact.dat sidecar ORCA wrote beside them, with the Agent's own
#: compiled invocation reproducing that native input and previewing green.
_QUALIFIED = {("orca", "scan")}


@pytest.mark.parametrize("program", ("orca", "gaussian"))
@pytest.mark.parametrize("jobtype", ("scan", "modred"))
def test_execution_is_claimed_only_where_a_run_qualified_it(program, jobtype):
    """The declaration must not quietly promise a run nobody has made.

    Adding a pair to `_QUALIFIED` is the deliberate act of claiming that
    evidence exists, so it cannot happen as a side effect of editing the
    matrix. Gaussian stays out of it regardless of evidence: this release
    does not claim Gaussian Agent execution at all.
    """

    expected = (program, jobtype) in _QUALIFIED

    assert _pairs(program)[jobtype].execution_supported is expected


@pytest.mark.parametrize("jobtype", ("scan", "modred"))
def test_gaussian_execution_is_never_claimed_here(jobtype):
    """A release-policy hold, not a missing-evidence one."""

    assert ("gaussian", jobtype) not in _QUALIFIED
    assert _pairs("gaussian")[jobtype].execution_supported is False


@pytest.mark.parametrize("program", ("orca", "gaussian"))
@pytest.mark.parametrize("jobtype", ("scan", "modred"))
def test_the_declaration_reaches_conformance(program, jobtype):
    """Conformance reads the same matrix, so the family must appear there."""

    assert jobtype in _conformance_jobtypes(program, "cpu")


@pytest.mark.parametrize("program", ("orca", "gaussian"))
@pytest.mark.parametrize("jobtype", ("scan", "modred"))
def test_the_project_layer_already_carries_the_family(program, jobtype):
    """Declaring a family the project loader cannot express would be a lie."""

    assert jobtype in molecular_project_section_names(program)


@pytest.mark.parametrize("program", ("orca", "gaussian"))
def test_nothing_previously_declared_was_lost(program):
    """Widening a surface must not narrow it somewhere else."""

    declared = set(_pairs(program))
    assert {"opt", "sp", "ts", "irc"} <= declared
