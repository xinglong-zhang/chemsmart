"""A declared project section must actually reach the job types that read it.

``ProgramCapability.project_section_names`` tells a caller which top-level keys
a program's project YAML accepts.  Acceptance is not the useful fact.  What a
chemist -- or a model -- needs to know is that settings written under an
accepted section *arrive* at the job they are running.

Those two came apart silently.  ``gas`` was declared and accepted for ORCA and
Gaussian, and ``sp`` read only ``solv``, so a project consisting of

    gas:
      functional: b3lyp
      basis: aug-cc-pVTZ

passed every gate the harness had, compiled a canonical command, round-tripped
through the CLI parser -- and then handed the writer settings with no level of
theory.  The run died at the bottom of the stack with a bare

    ValueError: Warning: neither ab initio nor DFT is specified!

after creating a zero-byte input file.  A live agent session read that as
"ORCA cannot materialise a single point in this environment" and switched
program, which is the hub silently narrowing its own program surface.

This gate derives the truth from the loader instead of restating a routing
table that would drift from it.  For every phase-keyed program the registry
declares, it writes a project containing exactly one phase section and asserts
that every job type fed by phase sections receives the level of theory in it.
Job types that own a dedicated section (``td``, ``qmmm``) are excluded, because
requiring their own section is the declared design rather than a gap.
"""

import os
import tempfile

import pytest

from chemsmart.jobs.settings import read_molecular_job_yaml
from chemsmart.settings.capabilities import PROGRAM_CAPABILITIES

#: Sections that group settings by phase rather than by job type.
PHASE_SECTIONS = ("gas", "solv")

#: Job types the loader deliberately routes to their own section, so a
#: phase-only project is not expected to feed them.
SECTION_OWNING_JOBTYPES = frozenset({"td", "qmmm"})

#: Job types with no settings of their own to merge.
NON_MOLECULAR_JOBTYPES = frozenset({"inp", "userjob"})


def _phase_keyed_programs():
    """Yield the programs whose project YAML groups settings by phase."""
    for name, capability in sorted(PROGRAM_CAPABILITIES.items()):
        sections = set(capability.project_section_names)
        if sections.issuperset(PHASE_SECTIONS):
            # PySCF declares the phase pair only as legacy migration input; it
            # keys its real sections by stage and is covered by its own tests.
            if sections - set(PHASE_SECTIONS) - set(capability.jobtypes):
                continue
            if set(capability.jobtypes) & sections:
                continue
            yield name, capability


def _configs_from_single_section(section, program):
    """Load a project whose only content is one phase section."""
    body = f"{section}:\n  functional: b3lyp\n  basis: aug-cc-pVTZ\n"
    with tempfile.TemporaryDirectory() as directory:
        path = os.path.join(directory, "project.yaml")
        with open(path, "w") as handle:
            handle.write(body)
        return read_molecular_job_yaml(path, program=program)


def test_the_registry_declares_at_least_one_phase_keyed_program():
    """Guard the gate itself: a silently empty sweep would prove nothing."""
    assert dict(_phase_keyed_programs())


@pytest.mark.parametrize("section", PHASE_SECTIONS)
def test_every_declared_phase_section_feeds_every_jobtype_that_reads_it(
    section,
):
    """The regression this gate exists for, swept over the whole registry."""
    starved = []
    for program, capability in _phase_keyed_programs():
        configs = _configs_from_single_section(section, program)
        for jobtype in sorted(capability.jobtypes):
            if jobtype in SECTION_OWNING_JOBTYPES:
                continue
            if jobtype in NON_MOLECULAR_JOBTYPES:
                continue
            settings = configs.get(jobtype)
            if settings is None:
                # A job type the loader does not build from project YAML at
                # all is a different question, answered by its own accessor.
                continue
            if settings.get("functional") is None and (
                settings.get("ab_initio") is None
            ):
                starved.append(f"{program}.{jobtype} from '{section}:'")

    assert not starved, (
        "these job types accept a declared section and then receive no level "
        f"of theory from it: {starved}. A project that passes every gate must "
        "not reach the input writer with nothing to write."
    )


@pytest.mark.parametrize("section", PHASE_SECTIONS)
def test_a_single_phase_section_never_leaves_a_single_point_doing_frequencies(
    section,
):
    """``sp`` means one energy; inheriting an optimisation's ``freq`` is wrong.

    Recorded rather than asserted for ``solv:``: that branch predates this gate
    and changing it would move existing projects, so the asymmetry is reported
    in the campaign notes instead of silently repaired here.
    """
    if section != "gas":
        pytest.skip("the solv-only branch is pre-existing behaviour")

    for program, _capability in _phase_keyed_programs():
        body = (
            "gas:\n  functional: b3lyp\n  basis: aug-cc-pVTZ\n  freq: true\n"
        )
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, "project.yaml")
            with open(path, "w") as handle:
                handle.write(body)
            configs = read_molecular_job_yaml(path, program=program)

        assert configs["sp"]["freq"] is False, program
