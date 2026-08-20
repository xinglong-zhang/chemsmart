"""A scan with no driven coordinate is an optimisation under a false name.

The driven coordinate has to be threaded through every path that compiles a
node, and there are several: the prepare path, the amend rebuild, the approved
binding the executor rebuilds from, the repair path, and the review builder for
a node whose geometry comes from a producer. Each one that forgot it failed the
same way -- silently, by emitting a command that runs cleanly and computes
something other than what was planned.

The last of them was found on a real paper task. Two torsion scans were
specified correctly by the model, with the atom mapping translated properly out
of the paper's own labels, and reached the human review as a bare `scan` with no
range at all, because those nodes consumed a producer geometry and so compiled
through the review builder rather than the prepare path.

A guard at the compile boundary is worth more than fixing the fifth site,
because it makes the sixth fail loudly instead of quietly.
"""

from __future__ import annotations

import pytest

from chemsmart.agent.commands import native_coordinate_options

_SCAN = {
    "scan": {
        "kind": "dihedral",
        "atoms": [1, 2, 3, 4],
        "start": 0.0,
        "stop": 330.0,
        "points": 12,
    }
}
_HELD = {"constrained": [{"kind": "bond", "atoms": [1, 2]}]}


@pytest.mark.parametrize("jobtype", ("scan", "modred"))
def test_the_guard_names_the_substitution_it_prevents(jobtype):
    """Driven straight at the boundary, without the whole state machine."""

    from chemsmart.agent.commands import _COORDINATE_DRIVEN_JOBTYPES

    assert jobtype in _COORDINATE_DRIVEN_JOBTYPES


def test_an_ordinary_optimisation_is_not_caught_by_the_guard():
    from chemsmart.agent.commands import _COORDINATE_DRIVEN_JOBTYPES

    for jobtype in ("opt", "sp", "ts", "irc", "td", "neb", "hess"):
        assert jobtype not in _COORDINATE_DRIVEN_JOBTYPES


def test_a_driven_coordinate_still_renders_for_both_programs():
    """The guard must not fire on the case it is protecting."""

    assert native_coordinate_options("orca", _SCAN)
    assert native_coordinate_options("gaussian", _SCAN)
    assert native_coordinate_options("orca", _HELD)


def test_an_absent_specification_renders_nothing_to_pass():
    """Which is exactly the state the guard has to catch downstream."""

    assert native_coordinate_options("orca", None) == {}
    assert native_coordinate_options("orca", {}) == {}


def test_the_paper_task_mapping_is_expressible_as_written():
    """Atom order in the file need not match the labels in a paper.

    The real task described psi and theta by connectivity, and the model
    translated them onto this molecule's own numbering. Both must render.
    """

    psi = native_coordinate_options("orca", _SCAN)
    theta = native_coordinate_options(
        "orca",
        {
            "scan": {
                "kind": "dihedral",
                "atoms": [4, 3, 5, 6],
                "start": 0.0,
                "stop": 330.0,
                "points": 12,
            }
        },
    )

    assert psi["coordinates"] == "[[1, 2, 3, 4]]"
    assert theta["coordinates"] == "[[4, 3, 5, 6]]"
    assert psi["dist_start"] == "0.0" and psi["dist_end"] == "330.0"
    assert psi["num_steps"] == "12"
