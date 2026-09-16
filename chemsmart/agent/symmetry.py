"""What the host can say about a starting structure's symmetry.

A source geometry carries its builder's symmetry: an exactly symmetric
start converges to the nearest stationary point of that symmetry, which
is a saddle when the true minimum is lower in symmetry. Observed live
six times before this module existed (an idealised D4h iron complex
whose triplet sat on a doubly degenerate imaginary pair; a PMe3 ligand
built with three exactly threefold methyl rotors, all three optimised
onto rotor saddles; a planar CH2OH and a linear vinyl C-H). Nothing here
judges: an estimate is stated with the tolerance that found it, an
idealised coordinate is counted, and a perturbation is applied only
when asked, by seed, so the same request gives the same bytes.
"""

from __future__ import annotations

import math
from typing import Any, Iterable, Sequence

import numpy as np

#: Torsions and angles a builder types exactly.  Multiples of 60 degrees
#: place threefold rotors on their saddles; 90 and 180 set planes and
#: lines; the tetrahedral angle is the idealised sp3 value.
_IDEALISED_TORSION_PERIOD_DEGREES = 60.0
_IDEALISED_ANGLES_DEGREES = (90.0, 109.4712206, 120.0, 180.0)
_IDEALISED_TOLERANCE_DEGREES = 1.0e-3


def _relative(positions: np.ndarray, masses: np.ndarray) -> np.ndarray:
    centre = (positions * masses[:, None]).sum(axis=0) / masses.sum()
    return positions - centre


def _matches(
    symbols: Sequence[str],
    positions: np.ndarray,
    transformed: np.ndarray,
    tolerance: float,
) -> bool:
    used: set[int] = set()
    for index, row in enumerate(transformed):
        found = False
        for candidate, other in enumerate(positions):
            if candidate in used or symbols[candidate] != symbols[index]:
                continue
            if float(np.linalg.norm(row - other)) <= tolerance:
                used.add(candidate)
                found = True
                break
        if not found:
            return False
    return True


def _rotation(axis: np.ndarray, angle: float) -> np.ndarray:
    unit = axis / np.linalg.norm(axis)
    x, y, z = unit
    c, s = math.cos(angle), math.sin(angle)
    t = 1.0 - c
    return np.array(
        [
            [t * x * x + c, t * x * y - s * z, t * x * z + s * y],
            [t * x * y + s * z, t * y * y + c, t * y * z - s * x],
            [t * x * z - s * y, t * y * z + s * x, t * z * z + c],
        ]
    )


def _reflection(normal: np.ndarray) -> np.ndarray:
    unit = normal / np.linalg.norm(normal)
    return np.eye(3) - 2.0 * np.outer(unit, unit)


def _candidate_directions(
    relative: np.ndarray, principal: np.ndarray
) -> list[np.ndarray]:
    """Directions worth testing as axes or plane normals: the principal
    axes, every atom's direction, every atom pair's midpoint and
    difference direction. Deterministic and finite."""

    found: list[np.ndarray] = [principal[:, k] for k in range(3)]
    rows = [row for row in relative if np.linalg.norm(row) > 1.0e-6]
    raw: list[np.ndarray] = list(rows)
    for i in range(len(rows)):
        for j in range(i + 1, len(rows)):
            mid = rows[i] + rows[j]
            diff = rows[i] - rows[j]
            if np.linalg.norm(mid) > 1.0e-6:
                raw.append(mid)
            if np.linalg.norm(diff) > 1.0e-6:
                raw.append(diff)
    found.extend(raw)
    # A perpendicular C2 or a vertical plane of a nearly symmetric
    # structure lies exactly in the plane normal to a principal axis;
    # a perturbed atom direction does not, so its projection is added.
    for k in range(3):
        axis = principal[:, k]
        for direction in raw:
            projected = direction - float(direction @ axis) * axis
            if np.linalg.norm(projected) > 1.0e-6:
                found.append(projected)
    unique: list[np.ndarray] = []
    for direction in found:
        unit = direction / np.linalg.norm(direction)
        if any(
            abs(float(abs(unit @ other)) - 1.0) < 1.0e-6 for other in unique
        ):
            continue
        unique.append(unit)
    return unique


def point_group_estimate(
    symbols: Sequence[str],
    positions: Iterable[Sequence[float]],
    *,
    tolerance_angstrom: float = 0.01,
) -> dict[str, Any]:
    """An estimate of the point group, with the elements that support it.

    Operations are tested against the molecule's own atoms within the
    stated tolerance; the label follows the ordinary flowchart over the
    elements found. Cubic groups are reported as the family they belong
    to when several C3 axes are found. An estimate, never a claim.
    """

    from chemsmart.utils.periodictable import PeriodicTable

    table = PeriodicTable()
    array = np.asarray(list(positions), dtype=float)
    count = len(array)
    if count < 2:
        return {
            "point_group": "K",
            "tolerance_angstrom": tolerance_angstrom,
            "elements": (),
            "note": "a single atom carries every symmetry element",
        }
    masses = np.array(
        [float(table.to_atomic_mass(symbol)) for symbol in symbols]
    )
    relative = _relative(array, masses)
    inertia = np.zeros((3, 3))
    for mass, row in zip(masses, relative):
        inertia += mass * (float(row @ row) * np.eye(3) - np.outer(row, row))
    moments, principal = np.linalg.eigh(inertia)
    tolerance = float(tolerance_angstrom)
    elements: list[str] = []

    def _has(matrix: np.ndarray) -> bool:
        return _matches(symbols, relative, relative @ matrix.T, tolerance)

    linear = bool(moments[0] < 1.0e-6 * max(moments[2], 1.0e-12))
    if linear:
        inversion = _has(-np.eye(3))
        return {
            "point_group": "D∞h" if inversion else "C∞v",
            "tolerance_angstrom": tolerance,
            "elements": ("C∞", "i") if inversion else ("C∞",),
        }
    inversion = _has(-np.eye(3))
    if inversion:
        elements.append("i")
    directions = _candidate_directions(relative, principal)
    axes: list[tuple[int, np.ndarray]] = []
    for direction in directions:
        for order in (6, 5, 4, 3, 2):
            if _has(_rotation(direction, 2.0 * math.pi / order)):
                axes.append((order, direction))
                break
    planes = [
        direction for direction in directions if _has(_reflection(direction))
    ]
    if not axes:
        if planes:
            label = "Cs"
            elements.append("σ")
        elif inversion:
            label = "Ci"
        else:
            label = "C1"
        return {
            "point_group": label,
            "tolerance_angstrom": tolerance,
            "elements": tuple(elements),
        }
    highest = max(order for order, _direction in axes)
    principal_axes = [d for order, d in axes if order == highest]
    if highest >= 3 and len(principal_axes) > 1:
        label = "cubic (Td, Oh or Ih family)"
        elements.append(f"{len(principal_axes)} C{highest} axes")
        return {
            "point_group": label,
            "tolerance_angstrom": tolerance,
            "elements": tuple(elements),
        }
    main = principal_axes[0]
    elements.append(f"C{highest}")
    # Classification uses an angular slack of about three degrees; the
    # geometric tolerance has already decided whether an element exists.
    slack = 0.05
    perpendicular = [
        d for order, d in axes if order == 2 and abs(float(d @ main)) < slack
    ]
    horizontal = any(abs(float(p @ main)) > math.cos(slack) for p in planes)
    vertical = [p for p in planes if abs(float(p @ main)) < slack]
    if len(perpendicular) >= highest:
        elements.append(f"{len(perpendicular)} C2 ⊥")
        if horizontal:
            elements.append("σh")
            label = f"D{highest}h"
        elif len(vertical) >= highest:
            elements.append(f"{len(vertical)} σd")
            label = f"D{highest}d"
        else:
            label = f"D{highest}"
    elif horizontal:
        elements.append("σh")
        label = f"C{highest}h"
    elif len(vertical) >= highest:
        elements.append(f"{len(vertical)} σv")
        label = f"C{highest}v"
    else:
        label = f"C{highest}"
    return {
        "point_group": label,
        "tolerance_angstrom": tolerance,
        "elements": tuple(elements),
    }


def symmetry_observation(
    symbols: Sequence[str], positions: Iterable[Sequence[float]]
) -> str:
    """One sentence a session and a reviewer can read beside a geometry."""

    array = np.asarray(list(positions), dtype=float)
    exact = point_group_estimate(symbols, array, tolerance_angstrom=0.01)
    if exact["point_group"] in {"C1", "K"}:
        loose = point_group_estimate(symbols, array, tolerance_angstrom=0.1)
        if loose["point_group"] == "C1":
            return (
                "symmetry estimate: C1 within 0.1 Å (no element found); an "
                "unsymmetric start carries no builder's saddle"
            )
        return (
            f"symmetry estimate: C1 within 0.01 Å, {loose['point_group']} "
            f"within 0.1 Å ({', '.join(loose['elements'])}); a nearly "
            "symmetric start relaxes freely"
        )
    loose = point_group_estimate(symbols, array, tolerance_angstrom=0.1)
    nearly = (
        f"; {loose['point_group']} within 0.1 Å"
        if loose["point_group"] != exact["point_group"]
        else ""
    )
    return (
        f"symmetry estimate: {exact['point_group']} within 0.01 Å "
        f"({', '.join(exact['elements']) or 'no element'}){nearly}; an "
        "exactly symmetric start converges to the nearest stationary point "
        "of that symmetry, which is a saddle when the minimum lies lower -- "
        "break_symmetry perturbs it by seed, and "
        "vibrational_mode_degeneracy_group shows whether an imaginary pair "
        "is degenerate"
    )


def _is_idealised_torsion(value: float) -> bool:
    remainder = abs(
        (
            value % _IDEALISED_TORSION_PERIOD_DEGREES
            + 0.5 * _IDEALISED_TORSION_PERIOD_DEGREES
        )
        % _IDEALISED_TORSION_PERIOD_DEGREES
        - 0.5 * _IDEALISED_TORSION_PERIOD_DEGREES
    )
    return remainder <= _IDEALISED_TOLERANCE_DEGREES


def _is_idealised_angle(value: float) -> bool:
    return any(
        abs(value - ideal) <= _IDEALISED_TOLERANCE_DEGREES
        for ideal in _IDEALISED_ANGLES_DEGREES
    )


def idealised_internal_coordinate_count(
    append_receipts: Iterable[Any],
    *,
    edit_receipts: Iterable[Any] = (),
) -> dict[str, Any]:
    """How many built coordinates were set to exactly idealised values.

    Counts requested torsions on the 60-degree lattice and requested
    angles at 90, 109.47, 120 or 180 degrees, to a millidegree. A
    builder that types 60/180/300 for three methyl hydrogens has placed
    a threefold rotor exactly on its saddle.

    Editing a coordinate onto the lattice is the same act as appending an
    atom onto it, and this counted only the appends -- so a live run that
    built its conformers by editing one torsion to exactly 0.00° and
    exactly 180.00° raised nothing, and the 0.00° structure was butane's
    syn-periplanar saddle, the top of the rotational barrier rather than
    a conformer. A bond length has no idealised value and is counted as
    a built coordinate and nothing more.
    """

    appended = 0
    built = 0
    idealised_torsions = 0
    idealised_angles = 0
    for receipt in append_receipts:
        appended += 1
        built += 1
        if _is_idealised_torsion(
            float(getattr(receipt, "dihedral_degrees", 0.0))
        ):
            idealised_torsions += 1
        if _is_idealised_angle(float(getattr(receipt, "angle_degrees", 0.0))):
            idealised_angles += 1
    for receipt in edit_receipts:
        built += 1
        operation = str(getattr(receipt, "operation", "") or "").lower()
        try:
            requested = float(getattr(receipt, "value_requested", 0.0))
        except (TypeError, ValueError):  # pragma: no cover
            continue
        if operation == "dihedral" and _is_idealised_torsion(requested):
            idealised_torsions += 1
        elif operation == "angle" and _is_idealised_angle(requested):
            idealised_angles += 1
    return {
        "appended_atoms": appended,
        "built_coordinates": built,
        "idealised_torsions": idealised_torsions,
        "idealised_angles": idealised_angles,
    }


def idealised_coordinate_observation(counts: dict[str, Any]) -> str:
    total = counts.get("built_coordinates", counts["appended_atoms"])
    return (
        f"{counts['idealised_torsions']} of {total} built coordinates "
        "sit at torsions on the exact 60° lattice and "
        f"{counts['idealised_angles']} at exactly idealised angles; an "
        "exactly idealised rotor starts on its own saddle"
    )


def seeded_perturbation(
    positions: Iterable[Sequence[float]],
    *,
    seed: int,
    amplitude_angstrom: float,
) -> tuple[np.ndarray, float, float]:
    """A deterministic displacement of every atom by at most the amplitude.

    Each atom moves by a vector drawn uniformly from the ball of that
    radius by the seeded generator; the mean displacement is removed so
    the centre stays put, and the steps are then rescaled together so
    that none exceeds the amplitude. Returns the new positions, the
    largest displacement achieved, and the RMS.
    """

    array = np.asarray(list(positions), dtype=float)
    generator = np.random.default_rng(int(seed))
    count = len(array)
    directions = generator.normal(size=(count, 3))
    norms = np.linalg.norm(directions, axis=1)
    norms[norms == 0.0] = 1.0
    directions = directions / norms[:, None]
    radii = float(amplitude_angstrom) * np.cbrt(generator.uniform(size=count))
    steps = directions * radii[:, None]
    steps = steps - steps.mean(axis=0)
    largest = float(np.linalg.norm(steps, axis=1).max())
    if largest > float(amplitude_angstrom) and largest > 0.0:
        steps = steps * (float(amplitude_angstrom) / largest)
    moved = array + steps
    lengths = np.linalg.norm(steps, axis=1)
    return moved, float(lengths.max()), float(np.sqrt((lengths**2).mean()))


__all__ = [
    "idealised_coordinate_observation",
    "idealised_internal_coordinate_count",
    "point_group_estimate",
    "seeded_perturbation",
    "symmetry_observation",
]
