"""A compiled ORCA optimisation says what iteration cap the program will
apply and which project field raises it.

geom_maxiter was on the project tool, rendered and unused, while a
21-atom Fe(II) complex hit ORCA's 3N = 63 cap in two sealed windows
(NOVEL-1/2 ino1, 2026-09-04): the lever was visible and nothing said this
system needed it. The host knows the atom count at compile time.
"""

from __future__ import annotations

import pytest

from chemsmart.agent.tool_runtime import compile_time_observations

pytestmark = pytest.mark.capability("tool:compile_command")


def test_an_orca_optimisation_without_the_field_names_the_cap():
    (note,) = compile_time_observations(
        program="orca",
        jobtype="opt",
        settings=(("functional", "tpssh"), ("basis", "def2-tzvp")),
        atom_count=21,
    )
    assert "3N = 63" in note
    assert "geom_maxiter" in note


def test_a_set_field_or_another_program_or_job_says_nothing():
    assert (
        compile_time_observations(
            program="orca",
            jobtype="opt",
            settings={"geom_maxiter": 300},
            atom_count=21,
        )
        == ()
    )
    assert (
        compile_time_observations(
            program="orca", jobtype="sp", settings={}, atom_count=21
        )
        == ()
    )
    assert (
        compile_time_observations(
            program="xtb", jobtype="opt", settings={}, atom_count=21
        )
        == ()
    )
