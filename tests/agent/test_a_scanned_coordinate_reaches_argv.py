"""A driven coordinate is scientific intent and must survive to the command.

Which atoms are scanned, over what range, in how many points, is a fact about
*this molecule in this calculation* -- the same class of fact as charge and
multiplicity, which already live on the node. Freezing it into a reusable
project would be wrong, and it has no artifact to bind to, so before this
change there was no way for it to reach the compiler at all: the only per-node
channel carried hashed files.

These tests exercise the rendering the host actually performs, and check the
translation both programs need. They differ genuinely -- ORCA drives a
coordinate between absolute endpoints, Gaussian by a repeated increment -- and
one physical specification has to become both, which is the hub's job.
"""

from __future__ import annotations

import pytest

from chemsmart.agent._contracts import ContractError
from chemsmart.agent.commands import native_coordinate_options
from chemsmart.agent.workflows import (
    ArtifactInputIntentV1,
    ArtifactOutputIntentV1,
    CommandNodeIntentV1,
)

_TORSION = {
    "scan": {
        "kind": "dihedral",
        "atoms": [1, 2, 3, 4],
        "start": 0.0,
        "stop": 180.0,
        "points": 13,
    }
}


def _node(jobtype, coordinates):
    return CommandNodeIntentV1(
        node_id="scan-a",
        program="orca",
        jobtype=jobtype,
        project_role="torsion",
        dependencies=(),
        inputs=(
            ArtifactInputIntentV1(
                binding_id="filename",
                artifact_class="geometry_xyz",
                artifact_id="geom",
                producer_node_id="",
                producer_output_id="",
            ),
        ),
        expected_outputs=(
            ArtifactOutputIntentV1(
                output_id="profile", artifact_class="orca_output"
            ),
        ),
        unresolved_fields=(),
        internal_coordinates=coordinates,
    )


def test_one_specification_renders_into_each_program_s_own_idiom():
    """The same physical scan, two native renderings."""

    orca = native_coordinate_options("orca", _TORSION)
    gaussian = native_coordinate_options("gaussian", _TORSION)

    # ORCA drives between absolute endpoints.
    assert orca["dist_start"] == "0.0"
    assert orca["dist_end"] == "180.0"
    assert orca["num_steps"] == "13"

    # Gaussian repeats an increment, so 13 points is 12 steps of 15 degrees.
    assert gaussian["step_size"] == "15.0"
    assert gaussian["num_steps"] == "12"

    # Both take the coordinate list in the same literal form.
    assert orca["coordinates"] == "[[1, 2, 3, 4]]"
    assert gaussian["coordinates"] == orca["coordinates"]


def test_a_constraint_alone_renders_without_scan_options():
    values = native_coordinate_options(
        "orca", {"constrained": [{"kind": "bond", "atoms": [1, 19]}]}
    )

    assert values == {"coordinates": "[[1, 19]]"}


def test_a_scan_may_also_carry_constraints():
    values = native_coordinate_options(
        "orca",
        {**_TORSION, "constrained": [{"kind": "bond", "atoms": [5, 6]}]},
    )

    assert values["coordinates"] == "[[1, 2, 3, 4], [5, 6]]"


@pytest.mark.parametrize(
    "entry,message",
    (
        ({"kind": "bond", "atoms": [1, 2, 3]}, "defined by 2 atoms"),
        ({"kind": "dihedral", "atoms": [1, 2]}, "defined by 4 atoms"),
        ({"kind": "torsion", "atoms": [1, 2]}, "unsupported internal"),
        ({"kind": "bond", "atoms": [0, 2]}, "numbered from 1"),
        ({"kind": "bond", "atoms": [2, 2]}, "names an atom twice"),
    ),
)
def test_a_malformed_coordinate_is_refused_with_its_reason(entry, message):
    with pytest.raises(ContractError, match=message):
        native_coordinate_options("orca", {"constrained": [entry]})


def test_a_degenerate_range_is_refused():
    spec = {"scan": {**_TORSION["scan"], "stop": 0.0}}
    with pytest.raises(ContractError, match="start and stop are the same"):
        native_coordinate_options("orca", spec)


def test_a_scan_node_without_a_driven_coordinate_is_refused():
    """A scan is defined by what it drives; say so where it is written."""

    with pytest.raises(ContractError, match="names no driven coordinate"):
        _node("scan", {"constrained": [{"kind": "bond", "atoms": [1, 2]}]})


def test_a_constrained_optimisation_needs_a_constraint():
    with pytest.raises(ContractError, match="names no constrained"):
        _node("modred", _TORSION)


def test_a_driven_coordinate_on_a_plain_optimisation_is_refused():
    with pytest.raises(ContractError, match="does not drive one"):
        _node("opt", _TORSION)


def test_an_empty_specification_is_refused_rather_than_ignored():
    with pytest.raises(ContractError, match="omit the field instead"):
        _node("scan", {})


def test_a_node_without_coordinates_is_unaffected():
    node = _node("opt", None)
    assert node.internal_coordinates is None
    assert native_coordinate_options("orca", None) == {}
