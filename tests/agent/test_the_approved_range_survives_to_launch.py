"""What was displayed and approved has to be what gets launched.

The executor is a fresh process holding none of the planning session's state.
It rebuilds every approved node from its `ApprovedNodeBindingV1` alone, through
`synthesize_command`. A driven coordinate that lived only in the planning draft
was therefore not there to rebuild from.

The consequence was visible in the first approved scan execution. Two nodes
scanning the same molecule over different ranges -- 0-350 degrees in 36 points
and 90-150 degrees in 61 -- rebuilt to the *same* coordinate-free argv and so
reported the same `invocation_identity_sha256`, and both were refused with
"recompiled ChemSmart CLI operation differs from human review". The approval
chain behaved correctly; what it was comparing against had lost the range.

The coordinate now travels in the binding, which is where it belongs: the range
a human saw is part of what that human approved. It has to stay additive, so
approvals recorded before it existed keep their digests and stay replayable --
the same treatment `auxiliary_input_bindings` already gets.
"""

from __future__ import annotations

from dataclasses import replace

from chemsmart.agent.execution import (
    ApprovedNodeBindingV1,
    _approved_node_binding_body,
)

_COORDINATES = {
    "scan": {
        "kind": "dihedral",
        "atoms": [3, 1, 2, 4],
        "start": 0.0,
        "stop": 350.0,
        "points": 36,
    }
}


def _binding(**overrides):
    base = dict(
        node_id="hooh-scan-full",
        program="orca",
        engine="cpu",
        jobtype="scan",
        project_artifact_sha256="a" * 64,
        settings_sha256="b" * 64,
        charge=0,
        multiplicity=1,
        input_mode="initial",
        initial_artifact_id="geometry-hooh",
        initial_artifact_sha256="c" * 64,
        scientific_identity_sha256="d" * 64,
        producer_edge_sha256="",
    )
    base.update(overrides)
    return ApprovedNodeBindingV1(**base)


def test_a_node_that_drives_nothing_keeps_the_body_it_always_had():
    """Digest compatibility: an old approval must project unchanged."""

    body = _approved_node_binding_body(_binding())

    assert "internal_coordinates" not in body


def test_the_driven_range_reaches_the_canonical_body():
    body = _approved_node_binding_body(
        _binding(internal_coordinates=_COORDINATES)
    )

    assert body["internal_coordinates"]["scan"]["start"] == 0.0
    assert body["internal_coordinates"]["scan"]["stop"] == 350.0
    assert body["internal_coordinates"]["scan"]["points"] == 36


def test_two_ranges_of_one_molecule_stay_distinguishable():
    """The exact collision that let two different scans look identical."""

    wide = _binding(internal_coordinates=_COORDINATES)
    narrow = replace(
        wide,
        node_id="hooh-scan-gauche-window",
        internal_coordinates={
            "scan": {
                "kind": "dihedral",
                "atoms": [3, 1, 2, 4],
                "start": 90.0,
                "stop": 150.0,
                "points": 61,
            }
        },
    )

    assert _approved_node_binding_body(wide) != _approved_node_binding_body(
        narrow
    )


def test_the_binding_round_trips_through_its_own_mapping():
    """The bundle is JSON on disk, so the field has to survive the trip."""

    import json

    original = _binding(internal_coordinates=_COORDINATES)
    body = _approved_node_binding_body(original)
    revived = ApprovedNodeBindingV1(
        **{**body, "internal_coordinates": json.loads(
            json.dumps(body["internal_coordinates"])
        )}
    )

    assert revived.internal_coordinates["scan"]["atoms"] == [3, 1, 2, 4]


def test_the_revived_range_still_renders_the_same_options():
    """A JSON round trip must not change the argv this produces."""

    import json

    from chemsmart.agent.commands import native_coordinate_options

    direct = native_coordinate_options("orca", _COORDINATES)
    revived = native_coordinate_options(
        "orca", json.loads(json.dumps(_COORDINATES))
    )

    assert direct == revived
    assert direct["dist_start"] == "0.0"
    assert direct["dist_end"] == "350.0"
    assert direct["num_steps"] == "36"
