from __future__ import annotations

import re
from typing import Any


def extract_gaussian_route(content: str | None) -> str | None:
    """Return the first Gaussian route section, including continuation lines."""
    if not isinstance(content, str):
        return None
    lines = content.splitlines()
    for index, line in enumerate(lines):
        if not line.lstrip().startswith("#"):
            continue
        route_lines = [line.strip()]
        for continuation in lines[index + 1 :]:
            stripped = continuation.strip()
            if not stripped:
                break
            route_lines.append(stripped)
        return " ".join(route_lines)
    return None


def extract_orca_route(content: str | None) -> str | None:
    """Return ORCA route lines. Phase 1 does not attach hard gates to them."""
    if not isinstance(content, str):
        return None
    route_lines = [
        line.strip()
        for line in content.splitlines()
        if line.lstrip().startswith("!")
    ]
    if not route_lines:
        return None
    return " ".join(route_lines)


def extract_cartesian_state(
    content: str | None,
    *,
    software: str,
) -> dict[str, Any] | None:
    """Extract charge, multiplicity, and explicit Cartesian element symbols."""

    if not isinstance(content, str):
        return None
    if software == "orca":
        return _extract_orca_cartesian_state(content)
    if software == "gaussian":
        return _extract_gaussian_cartesian_state(content)
    raise ValueError(f"unsupported software: {software}")


def _extract_orca_cartesian_state(content: str) -> dict[str, Any] | None:
    lines = content.splitlines()
    for index, line in enumerate(lines):
        match = re.match(r"^\s*\*\s+xyz\s+(-?\d+)\s+(\d+)\s*$", line, re.I)
        if match is None:
            continue
        symbols: list[str] = []
        for coordinate in lines[index + 1 :]:
            if coordinate.strip().startswith("*"):
                break
            symbol = _coordinate_symbol(coordinate)
            if symbol:
                symbols.append(symbol)
        return {
            "charge": int(match.group(1)),
            "multiplicity": int(match.group(2)),
            "element_symbols": symbols,
        }
    return None


def _extract_gaussian_cartesian_state(content: str) -> dict[str, Any] | None:
    lines = content.splitlines()
    for index, line in enumerate(lines):
        state_values = _gaussian_state_values(line)
        if state_values is None:
            continue
        symbols: list[str] = []
        atom_layers: list[str | None] = []
        for coordinate in lines[index + 1 :]:
            if not coordinate.strip():
                break
            symbol = _coordinate_symbol(coordinate)
            if symbol is None:
                symbols = []
                atom_layers = []
                break
            symbols.append(symbol)
            atom_layers.append(_gaussian_coordinate_layer(coordinate))
        if symbols:
            result: dict[str, Any] = {
                "charge": state_values[0],
                "multiplicity": state_values[1],
                "element_symbols": symbols,
                "charge_multiplicity_pairs": [
                    state_values[offset : offset + 2]
                    for offset in range(0, len(state_values), 2)
                ],
            }
            if any(layer is not None for layer in atom_layers):
                result["atom_layers"] = atom_layers
                result["layer_atoms"] = {
                    layer: [
                        atom_index
                        for atom_index, observed in enumerate(
                            atom_layers, start=1
                        )
                        if observed == layer
                    ]
                    for layer in ("H", "M", "L")
                }
            return result
    return None


def _gaussian_state_values(line: str) -> list[int] | None:
    """Return Gaussian charge/multiplicity pairs, including ONIOM state lines."""

    tokens = line.split()
    if not tokens or len(tokens) % 2 or len(tokens) > 12:
        return None
    if any(re.fullmatch(r"-?\d+", token) is None for token in tokens):
        return None
    values = [int(token) for token in tokens]
    if any(multiplicity < 1 for multiplicity in values[1::2]):
        return None
    return values


def _gaussian_coordinate_layer(line: str) -> str | None:
    """Extract the primary ONIOM layer marker from one coordinate row."""

    tokens = line.split()
    # The first token is the element symbol. Any later H/M/L token before an
    # optional link-atom suffix is the atom's ONIOM partition marker.
    for token in tokens[1:]:
        if token in {"H", "M", "L"}:
            return token
    return None


def _coordinate_symbol(line: str) -> str | None:
    token = line.strip().split(maxsplit=1)[0] if line.strip() else ""
    match = re.match(r"^([A-Z][a-z]?)\b", token)
    if match is None:
        return None
    return match.group(1)
