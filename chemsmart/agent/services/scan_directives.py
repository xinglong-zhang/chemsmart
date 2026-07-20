"""Gaussian modred/scan directive parsing for agent job synthesis."""

from __future__ import annotations

import re
from typing import Any


def partition_modred_route_directives(
    value: Any,
) -> tuple[list[str], list[str]]:
    entries = value if isinstance(value, (list, tuple)) else [value]
    directives: list[str] = []
    route_options: list[str] = []
    for entry in entries:
        for chunk in re.split(r"[;\n]+", str(entry)):
            stripped = chunk.strip()
            if not stripped:
                continue
            directive = canonical_modred_directive(stripped)
            if directive is None:
                route_options.append(stripped)
            else:
                directives.append(directive)
    return directives, route_options


def canonical_modred_directive(value: str) -> str | None:
    tokens = [token for token in re.split(r"[,\s]+", value) if token]
    if not tokens:
        return None
    coordinate_type = tokens[0].upper()
    atom_count = {"B": 2, "A": 3, "D": 4}.get(coordinate_type)
    if atom_count is None or len(tokens) < atom_count + 2:
        return None
    try:
        atoms = [int(token) for token in tokens[1 : atom_count + 1]]
    except ValueError:
        return None
    if any(atom < 1 for atom in atoms):
        return None
    directive = tokens[atom_count + 1].upper()
    if directive == "F" and len(tokens) == atom_count + 2:
        return " ".join((coordinate_type, *(str(atom) for atom in atoms), "F"))
    if directive != "S" or len(tokens) != atom_count + 4:
        return None
    try:
        steps = int(tokens[atom_count + 2])
        step_size = float(tokens[atom_count + 3])
    except ValueError:
        return None
    if steps < 1:
        return None
    return " ".join(
        (
            coordinate_type,
            *(str(atom) for atom in atoms),
            "S",
            str(steps),
            str(step_size),
        )
    )


def parse_gaussian_scan_definition(scan_definition: str) -> dict[str, Any]:
    lines = [
        line.strip()
        for chunk in scan_definition.splitlines()
        for line in chunk.split(";")
        if line.strip()
    ]
    if not lines:
        raise ValueError("scan_definition cannot be empty.")

    atom_count_by_coordinate = {"B": 2, "A": 3, "D": 4}
    coords: list[list[int]] = []
    num_steps: list[int] = []
    step_sizes: list[float] = []
    constrained_coordinates: list[list[int]] = []

    for line in lines:
        tokens = line.split()
        coordinate_type = tokens[0].upper()
        expected_atoms = atom_count_by_coordinate.get(coordinate_type)
        if expected_atoms is None:
            raise ValueError(
                "scan_definition must start with B, A, or D followed by "
                "1-based atom indices."
            )
        if len(tokens) < expected_atoms + 2:
            raise ValueError(
                f"scan_definition {line!r} is incomplete. Expected a "
                "coordinate, indices, and a scan directive."
            )

        try:
            coordinate = [
                int(token) for token in tokens[1 : 1 + expected_atoms]
            ]
        except ValueError as exc:
            raise ValueError(
                f"scan_definition {line!r} has non-integer atom indices."
            ) from exc

        directive = tokens[1 + expected_atoms].upper()
        if directive == "F":
            constrained_coordinates.append(coordinate)
            continue

        if directive != "S" or len(tokens) != expected_atoms + 4:
            raise ValueError(
                "scan_definition entries must look like "
                "'D 1 2 3 4 S 10 36.0' or 'B 1 2 F'."
            )

        try:
            steps = int(tokens[2 + expected_atoms])
            step_size = float(tokens[3 + expected_atoms])
        except ValueError as exc:
            raise ValueError(
                f"scan_definition {line!r} has an invalid step count "
                "or step size."
            ) from exc

        coords.append(coordinate)
        num_steps.append(steps)
        step_sizes.append(step_size)

    if not coords:
        raise ValueError(
            "scan_definition must include at least one scan entry with "
            "an 'S <steps> <step_size>' directive."
        )

    modred: dict[str, Any] = {
        "coords": coords,
        "num_steps": num_steps,
        "step_size": step_sizes,
    }
    if constrained_coordinates:
        modred["constrained_coordinates"] = constrained_coordinates
    return modred


__all__ = [
    "canonical_modred_directive",
    "parse_gaussian_scan_definition",
    "partition_modred_route_directives",
]
