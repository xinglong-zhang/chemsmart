"""A source geometry carries its builder's symmetry, and the host says so.

Six live saddles in two goals came from exactly symmetric starts: an
idealised D4h iron complex whose triplet sat on a degenerate imaginary
pair, three PMe3 methyl rotors appended at torsions of exactly 60, 180
and 300 degrees, a planar CH2OH and a linear vinyl C-H (NOVEL-2 ino1,
NOVEL-3 ino3, H4). The host now states a point-group estimate at bind
and compile time, counts idealised appended coordinates, and offers a
seeded perturbation whose bytes are reproducible from the seed.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from chemsmart.agent._contracts import ContractError
from chemsmart.agent.guides import GUIDES_BY_ID
from chemsmart.agent.symmetry import (
    idealised_internal_coordinate_count,
    point_group_estimate,
    seeded_perturbation,
    symmetry_observation,
)
from tests.agent.test_a_geometry_edit_is_arithmetic_not_judgement import (
    _host_with,
)

_WATER = ("O", "H", "H"), (
    (0.0, 0.0, 0.1173),
    (0.0, 0.7572, -0.4692),
    (0.0, -0.7572, -0.4692),
)
_SQUARE_PLANAR = ("Fe", "Cl", "Cl", "Cl", "Cl"), (
    (0.0, 0.0, 0.0),
    (2.3, 0.0, 0.0),
    (-2.3, 0.0, 0.0),
    (0.0, 2.3, 0.0),
    (0.0, -2.3, 0.0),
)
_METHANE = ("C", "H", "H", "H", "H"), (
    (0.0, 0.0, 0.0),
    (0.63, 0.63, 0.63),
    (-0.63, -0.63, 0.63),
    (-0.63, 0.63, -0.63),
    (0.63, -0.63, -0.63),
)
_CHIRAL = ("C", "H", "F", "Cl", "Br"), (
    (0.0, 0.0, 0.0),
    (0.63, 0.63, 0.63),
    (-0.8, -0.8, 0.8),
    (-1.0, 1.0, -1.0),
    (1.2, -1.2, -1.2),
)


@pytest.mark.capability("tool:break_symmetry")
def test_the_estimate_finds_the_textbook_groups_and_says_its_tolerance():
    assert point_group_estimate(*_WATER)["point_group"] == "C2v"
    square = point_group_estimate(*_SQUARE_PLANAR)
    assert square["point_group"] == "D4h"
    assert "4 C2 ⊥" in square["elements"] and "σh" in square["elements"]
    assert point_group_estimate(*_METHANE)["point_group"].startswith("cubic")
    assert point_group_estimate(*_CHIRAL)["point_group"] == "C1"
    assert (
        point_group_estimate(
            ("C", "O", "O"), ((0, 0, 0), (1.16, 0, 0), (-1.16, 0, 0))
        )["point_group"]
        == "D∞h"
    )
    # A nearly symmetric start is C1 at 0.01 A and named at 0.1 A. A
    # triatomic with equal bond lengths is C2v about its bisector however
    # it is tilted, so the nudge stretches one O-H bond by 0.036 A.
    symbols, rows = _WATER
    nudged = np.asarray(rows).copy()
    nudged[1] = np.array([0.0, 0.0, 0.1173]) + 1.04 * (
        nudged[1] - np.array([0.0, 0.0, 0.1173])
    )
    # A planar triatomic keeps its molecular plane: Cs, never C1.
    assert point_group_estimate(symbols, nudged)["point_group"] == "Cs"
    assert (
        point_group_estimate(symbols, nudged, tolerance_angstrom=0.1)[
            "point_group"
        ]
        == "C2v"
    )
    sentence = symmetry_observation(symbols, nudged)
    assert sentence.startswith(
        "symmetry estimate: Cs within 0.01 Å (σ); C2v within 0.1 Å"
    )
    exact = symmetry_observation(*_SQUARE_PLANAR)
    assert exact.startswith("symmetry estimate: D4h within 0.01 Å")
    assert "break_symmetry" in exact
    assert "vibrational_mode_degeneracy_group" in exact


def test_idealised_appended_coordinates_are_counted_not_judged():
    receipts = [
        SimpleNamespace(dihedral_degrees=60.0, angle_degrees=109.4712206),
        SimpleNamespace(dihedral_degrees=180.0, angle_degrees=109.4712206),
        SimpleNamespace(dihedral_degrees=-60.0, angle_degrees=109.4712206),
        SimpleNamespace(dihedral_degrees=47.3, angle_degrees=111.2),
    ]
    counts = idealised_internal_coordinate_count(receipts)
    assert counts == {
        "appended_atoms": 4,
        "idealised_torsions": 3,
        "idealised_angles": 3,
    }


def test_a_seeded_perturbation_is_reproducible_and_bounded():
    symbols, rows = _SQUARE_PLANAR
    first, largest, rms = seeded_perturbation(
        rows, seed=7, amplitude_angstrom=0.05
    )
    second, _largest, _rms = seeded_perturbation(
        rows, seed=7, amplitude_angstrom=0.05
    )
    other, _l, _r = seeded_perturbation(rows, seed=8, amplitude_angstrom=0.05)
    assert np.array_equal(first, second)
    assert not np.array_equal(first, other)
    assert 0.0 < rms <= largest <= 0.05 + 1e-12
    assert np.allclose(first.mean(axis=0), np.asarray(rows).mean(axis=0))


_SQUARE_XYZ = (
    "5\nidealised square-planar start\n"
    "Fe 0.0 0.0 0.0\nCl 2.3 0.0 0.0\nCl -2.3 0.0 0.0\n"
    "Cl 0.0 2.3 0.0\nCl 0.0 -2.3 0.0\n"
)


@pytest.mark.capability("tool:break_symmetry")
def test_break_symmetry_through_the_host_leaves_a_hop_and_an_observation(
    tmp_path,
):
    host = _host_with(tmp_path, _SQUARE_XYZ, "fecl4")
    bound = host.dispatch(
        turn_id="t1",
        tool_name="bind_scientific_identity",
        arguments={
            "input_artifact_id": "fecl4",
            "charge": -2,
            "multiplicity": 5,
        },
    )["result"]
    (symmetry_line,) = bound["observations"]
    assert symmetry_line.startswith("symmetry estimate: D4h within 0.01 Å")

    reply = host.dispatch(
        turn_id="t2",
        tool_name="break_symmetry",
        arguments={
            "perturbed_artifact_id": "fecl4-broken",
            "input_artifact_id": "fecl4",
            "seed": 11,
            "amplitude_angstrom": 0.03,
        },
    )["result"]
    receipt = reply["symmetry_break"]
    assert receipt["seed"] == 11
    assert 0.0 < receipt["max_displacement_angstrom"] <= 0.03
    assert receipt["point_group_before"] == "D4h"
    assert receipt["point_group_after"] == "C1"
    assert receipt["atom_count"] == 5 and receipt["formula"] == "Cl4Fe"
    kept = host.symmetry_breaks[reply["artifact"]["sha256"]]
    assert kept.receipt_sha256 == receipt["receipt_sha256"]

    again = _host_with(tmp_path / "again", _SQUARE_XYZ, "fecl4")
    twin = again.dispatch(
        turn_id="t2",
        tool_name="break_symmetry",
        arguments={
            "perturbed_artifact_id": "fecl4-broken",
            "input_artifact_id": "fecl4",
            "seed": 11,
            "amplitude_angstrom": 0.03,
        },
    )["result"]
    assert twin["artifact"]["sha256"] == reply["artifact"]["sha256"]

    rebound = host.dispatch(
        turn_id="t3",
        tool_name="bind_scientific_identity",
        arguments={
            "input_artifact_id": "fecl4-broken",
            "charge": -2,
            "multiplicity": 5,
        },
    )["result"]
    assert rebound["observations"][0].startswith(
        "symmetry estimate: C1 within 0.01 Å, D4h within 0.1 Å"
    )

    # The schema gate refuses a non-integer seed before the transform.
    with pytest.raises(ContractError, match="seed must be integer"):
        host.dispatch(
            turn_id="t4",
            tool_name="break_symmetry",
            arguments={
                "perturbed_artifact_id": "fecl4-x",
                "input_artifact_id": "fecl4",
                "seed": 1.5,
                "amplitude_angstrom": 0.03,
            },
        )
    with pytest.raises(ContractError, match="amplitude"):
        host.dispatch(
            turn_id="t5",
            tool_name="break_symmetry",
            arguments={
                "perturbed_artifact_id": "fecl4-y",
                "input_artifact_id": "fecl4",
                "seed": 1,
                "amplitude_angstrom": 0.9,
            },
        )
    # The transform guards the same line behind the schema.
    from chemsmart.agent.execution import break_trusted_molecular_symmetry

    parent = host.artifacts["fecl4"]
    with pytest.raises(ContractError, match=r"\(0, 0.5\] angstrom"):
        break_trusted_molecular_symmetry(
            approved_workspace=host.approved_workspace,
            perturbed_artifact_id="fecl4-z",
            parent=parent,
            parent_identity_sha256="a" * 64,
            seed=1,
            amplitude_angstrom=0.9,
        )


def test_break_symmetry_is_a_structure_leaf_and_the_hop_renders():
    from chemsmart.agent.tool_specs import build_command_compiled_tool_surface
    from chemsmart.agent.tui.review import _geometry_lineage_panels

    assert "break_symmetry" in GUIDES_BY_ID["structure"].tools
    names = {
        item["function"]["name"]
        for item in build_command_compiled_tool_surface(
            guides=("structure",)
        ).tool_definitions
    }
    assert "break_symmetry" in names

    from tests.agent.test_a_geometry_edit_is_arithmetic_not_judgement import (
        _duck_review,
        _rendered,
    )

    review = _duck_review(
        (
            SimpleNamespace(
                node_id="opt",
                molecular_identity={
                    "geometry_lineage": [
                        {
                            "kind": "symmetry_break",
                            "parent_artifact_id": "fecl4",
                            "formula": "Cl4Fe",
                            "atom_count": 5,
                            "seed": 11,
                            "amplitude_angstrom": 0.03,
                            "max_displacement_angstrom": 0.0271,
                            "rms_displacement_angstrom": 0.0188,
                            "point_group_before": "D4h",
                            "point_group_after": "C1",
                            "atom_order_note": "parent atom i is perturbed atom i",
                            "min_interatomic_distance_angstrom": 2.28,
                            "close_contact_pairs": (),
                            "connectivity_changed": False,
                        }
                    ]
                },
            ),
        )
    )
    text = "\n".join(
        _rendered(panel) for panel in _geometry_lineage_panels(review)
    )
    assert "symmetry break" in text
    assert "seed 11" in text
    assert "D4h -> C1" in text
