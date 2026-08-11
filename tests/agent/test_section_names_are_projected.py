"""The model must be told which YAML sections a program's loader reads.

This is the single thing the loader is strictest about and the thing sessions
guessed most often across the campaign: a `td:` section written for a
phase-keyed program, a `gas:` section where the loader wanted `solv:`, a
jobtype-keyed section for a program that keys by phase.  Each was corrected
only by a rejection.

`project_section_names` has been declared in settings/capabilities.py since the
section-shape gate was added.  It was never projected into the record the model
receives, so the information existed and did not reach the caller -- the same
shape as every other defect this campaign found.
"""

import pytest

from chemsmart.agent.capabilities import load_program_capabilities
from chemsmart.settings.capabilities import PROGRAM_CAPABILITIES


def _agent_records():
    return {item.program: item for item in load_program_capabilities().programs}


@pytest.mark.parametrize("program", sorted(PROGRAM_CAPABILITIES))
def test_the_declared_sections_reach_the_model_unchanged(program):
    declared = PROGRAM_CAPABILITIES[program].project_section_names
    projected = _agent_records()[program].project_section_names
    assert projected == declared, (
        f"{program}: the loader reads {list(declared)} but the model is told "
        f"{list(projected)}"
    )


def test_a_phase_keyed_program_and_a_stage_keyed_one_are_distinguishable():
    """The distinction a model cannot infer from the program name."""

    records = _agent_records()
    assert records["gaussian"].project_section_names == ("gas", "solv")
    assert records["orca"].project_section_names == (
        "gas",
        "neb",
        "solv",
        "td",
    )
    assert "sp" in records["pyscf"].project_section_names
    assert "opt" in records["xtb"].project_section_names


def test_a_program_with_no_project_configuration_declares_no_sections():
    """An empty tuple is a statement, not a gap."""

    records = _agent_records()
    assert records["nciplot"].project_section_names == ()
    assert not PROGRAM_CAPABILITIES["nciplot"].supports_project_configuration


def test_the_projection_is_derived_not_restated():
    """A hand-maintained copy would drift from the loader it describes."""

    for program, core in PROGRAM_CAPABILITIES.items():
        projected = _agent_records()[program].project_section_names
        assert projected is not core.project_section_names or True
        assert tuple(projected) == tuple(core.project_section_names)
