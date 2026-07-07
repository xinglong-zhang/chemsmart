"""Shared helpers for TUI transcript cells."""

from __future__ import annotations

from pathlib import Path
from tempfile import NamedTemporaryFile

from chemsmart.io.molecules.structure import Molecule


def molecule_xyz(molecule: Molecule) -> str:
    lines = [str(len(molecule.symbols)), molecule.chemical_formula]
    for symbol, coords in zip(molecule.symbols, molecule.positions):
        x, y, z = [float(value) for value in coords]
        lines.append(f"{symbol:<2} {x: .6f} {y: .6f} {z: .6f}")
    return "\n".join(lines)


def write_temp_xyz(
    molecule: Molecule, directory: str | Path | None = None
) -> Path:
    target_dir = None if directory is None else str(directory)
    with NamedTemporaryFile(
        "w",
        suffix=".xyz",
        prefix="chemsmart-agent-",
        delete=False,
        dir=target_dir,
        encoding="utf-8",
    ) as handle:
        handle.write(molecule_xyz(molecule) + "\n")
        return Path(handle.name)
