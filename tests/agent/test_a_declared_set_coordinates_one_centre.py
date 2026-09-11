"""A model-declared set of atoms held at one distance from one centre.

The composition primitive placed one or two named pairs at a requested
distance and pushed **every** other interfragment pair outside
covalent-radii-plus-buffer. For multi-centre coordination the other
members of the set must sit at the *same* short distance, which that
inequality forbids: for Fe-C the generic floor is 2.380 A against a real
2.06 A.

The live consequence was not a refusal but a wrong molecule the host
accepted. ino3-r16 (2026-09-10) composed a formally correct C10H10Fe --
21 atoms, one connected piece, no clash -- whose ten Fe-C distances were
all 2.600 A, every bond 0.54 A too long, because the floor set them. One
attempt did refuse, with wording that invited a retry at a different
distance when the binding constraint was on the un-named pairs; the
session retried four times.

Nothing here teaches the host what a ring is. It is told which atoms and
which distance, it enforces geometry and the floor on everything it was
not told to exempt, and it reports what it achieved per atom.
"""

import math
import re

import pytest

from chemsmart.agent._contracts import (
    ContractError,
    TrustedArtifactRefV1,
    file_sha256,
)
from chemsmart.agent.execution import compose_trusted_molecular_arrangement
from chemsmart.utils.periodictable import covalent_radii

#: Cp with C-C 1.43 A; its circumradius is 1.2164 A.
_CC = 1.43
_R = _CC / (2 * math.sin(math.pi / 5))


def _ring(path, *, atoms=5, radius=_R, out_of_plane=0.0, radial_kick=0.0):
    rows = []
    for k in range(atoms):
        angle = 2 * math.pi * k / atoms
        rad = radius + (radial_kick if k == 0 else 0.0)
        z = out_of_plane if k == 0 else 0.0
        rows.append(
            f"C {rad*math.cos(angle):.6f} {rad*math.sin(angle):.6f} {z:.6f}"
        )
    for k in range(atoms):
        angle = 2 * math.pi * k / atoms
        rows.append(
            f"H {(radius+1.08)*math.cos(angle):.6f} "
            f"{(radius+1.08)*math.sin(angle):.6f} 0.000000"
        )
    path.write_text(f"{2*atoms}\nring\n" + "\n".join(rows) + "\n")
    return path


def _saddle(path, radius=1.4, rise=0.5, atoms=4):
    """Four points equidistant from their centroid and not coplanar.

    Needed to test the coplanarity refusal on its own: perturbing one
    atom out of plane also changes its distance from the centroid, so
    the equidistance check fires first and the test would pass for the
    wrong reason.
    """

    in_plane = math.sqrt(radius * radius - rise * rise)
    rows = []
    for k in range(atoms):
        angle = 2 * math.pi * k / atoms
        z = rise if k % 2 == 0 else -rise
        rows.append(
            f"C {in_plane*math.cos(angle):.6f} "
            f"{in_plane*math.sin(angle):.6f} {z:.6f}"
        )
    path.write_text(f"{atoms}\nsaddle\n" + "\n".join(rows) + "\n")
    return path


def _atom(path, symbol="Fe"):
    path.write_text(f"1\ncentre\n{symbol} 0.000000 0.000000 0.000000\n")
    return path


def _ref(path, artifact_id):
    return TrustedArtifactRefV1(
        artifact_id=artifact_id,
        kind="geometry_xyz",
        sha256=file_sha256(path),
        size_bytes=path.stat().st_size,
        path=str(path),
        cli_value=str(path),
    )


def _compose(tmp_path, **kwargs):
    workspace = tmp_path / "ws"
    (workspace / "artifacts").mkdir(parents=True, exist_ok=True)
    return compose_trusted_molecular_arrangement(
        approved_workspace=workspace,
        fragment_a_identity_sha256="a" * 64,
        fragment_b_identity_sha256="b" * 64,
        **kwargs,
    )


def test_a_declared_set_reaches_a_real_metal_ligand_distance(tmp_path):
    """The distance the generic floor forbids is now reachable."""

    centre = _atom(tmp_path / "fe.xyz")
    ring = _ring(tmp_path / "cp.xyz")
    floor = covalent_radii[26] + covalent_radii[6] + 0.3
    assert floor == pytest.approx(2.38, abs=1e-9)

    _artifact, receipt = _compose(
        tmp_path,
        composed_artifact_id="fe-cp",
        fragment_a=_ref(centre, "fe"),
        fragment_b=_ref(ring, "cp"),
        fragment_a_atom=1,
        fragment_b_atom=1,
        distance_angstrom=2.06,
        fragment_b_atoms=[1, 2, 3, 4, 5],
    )
    placement = receipt.placement
    assert placement["mode"] == "haptic_set"
    achieved = placement["achieved_set_distances_angstrom"]
    # Every named atom at the requested distance, not just the one the
    # old primitive would have held.
    assert len(achieved) == 5
    assert all(value == pytest.approx(2.06, abs=1e-4) for value in achieved)
    # And it is inside the floor, which is the whole point.
    assert max(achieved) < floor
    # The evidence a reader needs to see that the set really is one.
    assert placement["set_circumradius_angstrom"] == pytest.approx(
        _R, abs=1e-3
    )
    assert placement["set_out_of_plane_angstrom"] == pytest.approx(
        0.0, abs=1e-6
    )
    assert placement["clash_floor_exempt_pairs"] == [
        [1, 1],
        [1, 2],
        [1, 3],
        [1, 4],
        [1, 5],
    ]
    assert receipt.formula == "C5H5Fe"


def test_the_generic_floor_is_what_made_it_unreachable(tmp_path):
    """The defect, reproduced: a single contact cannot hold the others.

    This is the 2.600 A pathology from the live run, in miniature. The
    named pair is held at 2.06; every other ring carbon is pushed out by
    the floor, so the set is not a set.
    """

    centre = _atom(tmp_path / "fe.xyz")
    ring = _ring(tmp_path / "cp.xyz")
    _artifact, receipt = _compose(
        tmp_path,
        composed_artifact_id="fe-cp-one",
        fragment_a=_ref(centre, "fe"),
        fragment_b=_ref(ring, "cp"),
        fragment_a_atom=1,
        fragment_b_atom=1,
        distance_angstrom=2.06,
    )
    assert receipt.placement["mode"] == "contact"
    assert receipt.achieved_contact_distance_angstrom == pytest.approx(
        2.06, abs=1e-3
    )
    assert "achieved_set_distances_angstrom" not in receipt.placement
    # The other ring carbons are pushed out past the Fe-C floor, which
    # is exactly the 2.600 A pathology: only the named pair is held.
    rows = [
        row.split()
        for row in (tmp_path / "ws" / "artifacts" / "fe-cp-one.xyz")
        .read_text()
        .splitlines()[2:]
        if row.strip()
    ]
    fe = [float(x) for x in rows[0][1:4]]
    carbons = [
        math.dist(fe, [float(x) for x in row[1:4]])
        for row in rows[1:]
        if row[0] == "C"
    ]
    floor_c = covalent_radii[26] + covalent_radii[6] + 0.3
    assert min(carbons) == pytest.approx(2.06, abs=1e-3)
    assert sum(1 for value in carbons if value < floor_c) == 1, carbons


def test_one_and_two_contact_composition_are_unchanged(tmp_path):
    """The existing modes keep their receipts byte-for-byte in shape."""

    a = _ring(tmp_path / "a.xyz", atoms=3, radius=1.2)
    b = _ring(tmp_path / "b.xyz", atoms=3, radius=1.2)
    _art, single = _compose(
        tmp_path,
        composed_artifact_id="single",
        fragment_a=_ref(a, "a"),
        fragment_b=_ref(b, "b"),
        fragment_a_atom=1,
        fragment_b_atom=1,
        distance_angstrom=3.0,
    )
    assert single.placement["mode"] == "contact"
    assert set(single.placement) == {
        "schema_version",
        "mode",
        "fragment_a_atom",
        "fragment_b_atom",
        "distance_angstrom",
        "buffer_angstrom",
        "sphere_direction_samples",
        "axial_rotation_samples",
    }
    _art2, dual = _compose(
        tmp_path,
        composed_artifact_id="dual",
        fragment_a=_ref(a, "a"),
        fragment_b=_ref(b, "b"),
        fragment_a_atom=1,
        fragment_b_atom=1,
        distance_angstrom=3.0,
        fragment_a_atom_2=2,
        fragment_b_atom_2=2,
        distance_angstrom_2=3.2,
    )
    assert dual.placement["mode"] == "dual_contact"
    assert "fragment_b_atoms" not in dual.placement


def test_the_floor_still_holds_for_atoms_outside_the_declared_set(tmp_path):
    """Exemption is per named pair, never a blanket.

    The set's own hydrogens are not in the set, so they keep the generic
    floor -- and that is what makes the exemption safe to have.
    """

    centre = _atom(tmp_path / "fe.xyz")
    ring = _ring(tmp_path / "cp.xyz")
    _artifact, receipt = _compose(
        tmp_path,
        composed_artifact_id="fe-cp-floor",
        fragment_a=_ref(centre, "fe"),
        fragment_b=_ref(ring, "cp"),
        fragment_a_atom=1,
        fragment_b_atom=1,
        distance_angstrom=2.06,
        fragment_b_atoms=[1, 2, 3, 4, 5],
    )
    exempt = {
        tuple(pair) for pair in receipt.placement["clash_floor_exempt_pairs"]
    }
    # Exactly the five declared pairs, and no hydrogen among them.
    assert exempt == {(1, i) for i in range(1, 6)}
    floor_h = covalent_radii[26] + covalent_radii[1] + 0.3
    # Read the composed file and check every un-named atom is outside.
    lines = (tmp_path / "ws" / "artifacts" / "fe-cp-floor.xyz").read_text()
    rows = [row.split() for row in lines.splitlines()[2:] if row.strip()]
    fe = [float(x) for x in rows[0][1:4]]
    for index, row in enumerate(rows[1:], start=1):
        if index in (1, 2, 3, 4, 5):
            continue
        gap = math.dist(fe, [float(x) for x in row[1:4]])
        assert gap >= floor_h - 1e-6, (index, row[0], gap)


@pytest.mark.parametrize(
    "kwargs, expected",
    [
        ({"fragment_b_atoms": [1, 1, 2]}, "names the same atom twice"),
        ({"fragment_b_atoms": [1]}, "two or more atoms"),
        ({"fragment_b_atoms": [1, 99]}, "must be 1.."),
        (
            {"fragment_b_atoms": [2, 3, 4, 5]},
            "is not in fragment_b_atoms",
        ),
        (
            {
                "fragment_b_atoms": [1, 2, 3, 4, 5],
                "fragment_a_atom_2": 1,
                "fragment_b_atom_2": 2,
                "distance_angstrom_2": 3.0,
            },
            "two different placements",
        ),
    ],
)
def test_a_malformed_set_is_refused_mechanically(tmp_path, kwargs, expected):
    """Every refusal is arithmetic on indices, never chemistry."""

    centre = _atom(tmp_path / "fe.xyz")
    ring = _ring(tmp_path / "cp.xyz")
    with pytest.raises(ContractError, match=re.escape(expected)):
        _compose(
            tmp_path,
            composed_artifact_id="bad",
            fragment_a=_ref(centre, "fe"),
            fragment_b=_ref(ring, "cp"),
            fragment_a_atom=1,
            fragment_b_atom=1,
            distance_angstrom=2.06,
            **kwargs,
        )


def test_a_set_that_is_not_a_circle_is_refused_with_its_numbers(tmp_path):
    """One distance to one centre exists only for a set on a circle."""

    centre = _atom(tmp_path / "fe.xyz")
    skew = _ring(tmp_path / "skew.xyz", radial_kick=0.9)
    with pytest.raises(ContractError, match="not equidistant"):
        _compose(
            tmp_path,
            composed_artifact_id="skew",
            fragment_a=_ref(centre, "fe"),
            fragment_b=_ref(skew, "skew"),
            fragment_a_atom=1,
            fragment_b_atom=1,
            distance_angstrom=2.06,
            fragment_b_atoms=[1, 2, 3, 4, 5],
        )
    bent = _saddle(tmp_path / "bent.xyz")
    with pytest.raises(ContractError, match="not coplanar"):
        _compose(
            tmp_path,
            composed_artifact_id="bent",
            fragment_a=_ref(centre, "fe"),
            fragment_b=_ref(bent, "bent"),
            fragment_a_atom=1,
            fragment_b_atom=1,
            distance_angstrom=2.06,
            fragment_b_atoms=[1, 2, 3, 4],
        )


def test_a_distance_inside_the_rings_own_radius_is_refused(tmp_path):
    """Geometry, stated with both numbers and no chemistry."""

    centre = _atom(tmp_path / "fe.xyz")
    ring = _ring(tmp_path / "cp.xyz")
    with pytest.raises(ContractError, match="circumradius"):
        _compose(
            tmp_path,
            composed_artifact_id="tooclose",
            fragment_a=_ref(centre, "fe"),
            fragment_b=_ref(ring, "cp"),
            fragment_a_atom=1,
            fragment_b_atom=1,
            distance_angstrom=1.0,
            fragment_b_atoms=[1, 2, 3, 4, 5],
        )


def test_the_primitive_infers_no_chemistry(tmp_path):
    """It is element-blind and hapticity-blind.

    The same call with a different metal, a different ring size and a
    nonsense element behaves identically: the host is told the indices
    and the distance and never asks what they are. If any of these grew
    a special case, one of them would differ.
    """

    ring4 = _ring(tmp_path / "r4.xyz", atoms=4, radius=1.4)
    results = {}
    for symbol, element in (("Fe", 26), ("Ni", 28), ("Xe", 54)):
        centre = _atom(tmp_path / f"{symbol}.xyz", symbol=symbol)
        _artifact, receipt = _compose(
            tmp_path,
            composed_artifact_id=f"set-{symbol}",
            fragment_a=_ref(centre, symbol.lower()),
            fragment_b=_ref(ring4, "r4"),
            fragment_a_atom=1,
            fragment_b_atom=1,
            distance_angstrom=2.30,
            fragment_b_atoms=[1, 2, 3, 4],
        )
        results[symbol] = receipt.placement["achieved_set_distances_angstrom"]
        assert receipt.placement["mode"] == "haptic_set"
        del element
    # Identical geometry for every centre element: no per-element table,
    # no hapticity lookup, no name-based exception.
    assert results["Fe"] == results["Ni"] == results["Xe"]
    # And a four-membered set works exactly as a five-membered one does;
    # nothing privileges eta5.
    assert all(
        value == pytest.approx(2.30, abs=1e-4) for value in results["Fe"]
    )
