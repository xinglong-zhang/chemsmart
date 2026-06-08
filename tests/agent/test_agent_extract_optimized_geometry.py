from __future__ import annotations

from pathlib import Path

import pytest

from chemsmart.agent.registry import ToolRegistry
from chemsmart.agent.tools import (
    build_gaussian_settings,
    build_job,
    extract_optimized_geometry,
)
from chemsmart.io.molecules.structure import Molecule


def _gaussian_opt_job(tmp_path: Path):
    molecule = Molecule(
        symbols=["O", "H", "H"],
        positions=[
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 0.0],
        ],
    )
    settings = build_gaussian_settings(
        "B3LYP",
        "6-31G*",
        charge=-1,
        multiplicity=2,
    )
    job = build_job(
        "gaussian.opt",
        molecule=molecule,
        settings=settings,
        label="water_opt",
    )
    job.set_folder(str(tmp_path))
    return job


@pytest.fixture()
def gaussian_optimized_log(tmp_path: Path) -> Path:
    logfile = tmp_path / "water_opt.log"
    logfile.write_text(
        """
 Random Gaussian output
                         Standard orientation:
 ---------------------------------------------------------------------
 Center     Atomic      Atomic             Coordinates (Angstroms)
 Number     Number       Type             X           Y           Z
 ---------------------------------------------------------------------
      1          8           0        0.000000    0.000000    0.100000
      2          1           0       -0.800000    0.000000   -0.500000
      3          1           0        0.800000    0.000000   -0.500000
 ---------------------------------------------------------------------
 More output
                         Standard orientation:
 ---------------------------------------------------------------------
 Center     Atomic      Atomic             Coordinates (Angstroms)
 Number     Number       Type             X           Y           Z
 ---------------------------------------------------------------------
      1          8           0        0.000000    0.000000    0.062600
      2          1           0       -0.792000    0.000000   -0.497300
      3          1           0        0.792000    0.000000   -0.497300
 ---------------------------------------------------------------------
 Normal termination
""".strip() + "\n",
        encoding="utf-8",
    )
    return logfile


def test_extract_optimized_geometry_reads_gaussian_log(
    tmp_path: Path,
    gaussian_optimized_log: Path,
):
    job = _gaussian_opt_job(tmp_path)

    molecule = extract_optimized_geometry(job)

    assert gaussian_optimized_log.exists()
    assert len(molecule) == 3
    assert molecule.symbols == ["O", "H", "H"]
    assert molecule.charge == -1
    assert molecule.multiplicity == 2
    assert molecule.positions[0][2] == pytest.approx(0.0626)


def test_extract_optimized_geometry_raises_if_no_log(tmp_path: Path):
    job = _gaussian_opt_job(tmp_path)

    with pytest.raises(FileNotFoundError):
        extract_optimized_geometry(job)


def test_registry_includes_extract_optimized_geometry():
    registry = ToolRegistry.default()

    tool_names = [tool.name for tool in registry.list_tools()]

    assert "extract_optimized_geometry" in tool_names
