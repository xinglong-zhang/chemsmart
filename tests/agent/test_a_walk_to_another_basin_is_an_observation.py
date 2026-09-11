"""An optimisation that walked far from its start, or made or broke a
bond, is an observation with the numbers, never a verdict.

A planar cyclohexane relaxed to the twist-boat in two live arms and
both delivered its energy as "the minimum": the stationary-point rule
types the order of a point, not which basin it is in. The host now
measures the walk against the structure it was given.
"""

import hashlib
from pathlib import Path

import numpy as np
import pytest

from chemsmart.agent._contracts import TrustedArtifactRefV1
from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1

_SN2 = (
    Path(__file__).resolve().parents[1]
    / "data"
    / "ORCATests"
    / "outputs"
    / "sn2_ts.out"
)


def _artifact(path: Path, kind: str, artifact_id: str) -> TrustedArtifactRefV1:
    return TrustedArtifactRefV1(
        artifact_id=artifact_id,
        kind=kind,
        sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        size_bytes=path.stat().st_size,
        path=str(path.resolve()),
        cli_value=str(path.resolve()),
    )


def _final_geometry():
    text = _SN2.read_text(errors="ignore")
    block = text.split("CARTESIAN COORDINATES (ANGSTROEM)")[-1].splitlines()[
        2:
    ]
    atoms = []
    for line in block:
        parts = line.split()
        if len(parts) != 4:
            break
        atoms.append((parts[0], np.array([float(v) for v in parts[1:4]])))
    return atoms


def _write_xyz(path: Path, atoms) -> Path:
    path.write_text(
        f"{len(atoms)}\nstart\n"
        + "".join(f"{s} {x:.6f} {y:.6f} {z:.6f}\n" for s, (x, y, z) in atoms)
    )
    return path


def _signals(tmp_path, atoms):
    start = _write_xyz(tmp_path / "start.xyz", atoms)
    evaluation = CommandCompiledToolHostV1._evaluate_execution_outputs(
        program="orca",
        jobtype="ts",
        charge=-1,
        multiplicity=1,
        expected_input_artifact=_artifact(start, "geometry_xyz", "start"),
        output_artifacts=(_artifact(_SN2, "orca_output", "result"),),
        exit_status=0,
    )
    return {item["signal_id"]: item for item in evaluation.anomalies}


@pytest.mark.capability("signal:geometry.heavy_atom_rmsd_ge_0.3")
def test_a_result_at_its_own_start_walks_nowhere(tmp_path):
    assert _signals(tmp_path, _final_geometry()) == {}


@pytest.mark.capability("signal:geometry.heavy_atom_rmsd_ge_0.3")
def test_a_heavy_atom_walk_is_measured(tmp_path):
    atoms = _final_geometry()
    heavy = [index for index, (symbol, _) in enumerate(atoms) if symbol != "H"]
    moved = list(atoms)
    symbol, position = moved[heavy[0]]
    moved[heavy[0]] = (symbol, position + np.array([1.2, 0.0, 0.0]))
    signals = _signals(tmp_path, moved)
    assert "geometry.heavy_atom_rmsd_ge_0.3" in signals
    assert (
        signals["geometry.heavy_atom_rmsd_ge_0.3"]["heavy_atom_rmsd_angstrom"]
        >= 0.3
    )


@pytest.mark.capability("signal:geometry.connectivity_changed")
def test_a_bond_made_or_broken_is_recorded_with_its_atoms(tmp_path):
    atoms = _final_geometry()
    hydrogens = [i for i, (symbol, _) in enumerate(atoms) if symbol == "H"]
    heavy = [i for i, (symbol, _) in enumerate(atoms) if symbol != "H"]
    moved = list(atoms)
    # Park one hydrogen far from everything: the start has one bond fewer.
    symbol, _ = moved[hydrogens[0]]
    moved[hydrogens[0]] = (symbol, np.array([30.0, 30.0, 30.0]))
    signals = _signals(tmp_path, moved)
    assert "geometry.connectivity_changed" in signals
    assert signals["geometry.connectivity_changed"]["bonds_made"]
    assert heavy  # the fixture carries heavy atoms to align on


@pytest.mark.capability("signal:geometry.connectivity_changed")
def test_an_equilibrium_bond_to_hydrogen_is_a_bond_to_the_perceiver():
    """The perceiver shrank every X-H buffer to 0.05 A, so the PH3
    minimum's 1.430 A P-H bond (cutoff 1.07+0.31+0.05) was "broken" on
    every repaired phosphine and a mode-displaced start's stretched N-H
    was "made" on every repaired ammonia (E4 window, 2026-09-03). The
    sensor asks whether topology changed; hydrogens get the same buffer
    as every other atom."""

    from chemsmart.agent.execution import _molecule_graph
    from chemsmart.io.molecules.structure import Molecule

    phosphine = Molecule(
        symbols=["P", "H", "H", "H"],
        positions=[
            [0.0, 0.0, 0.0],
            [1.430, 0.0, 0.0],
            [-0.715, 1.2384, 0.0],
            [-0.715, -1.2384, 0.0],
        ],
    )
    assert sorted(_molecule_graph(phosphine).edges) == [(1, 2), (1, 3), (1, 4)]
    # A start stepped 0.15 A along the umbrella mode keeps its bonds.
    stretched = Molecule(
        symbols=["N", "H", "H", "H"],
        positions=[
            [0.0, 0.0, 0.2],
            [1.15, 0.0, -0.1],
            [-0.575, 0.996, -0.1],
            [-0.575, -0.996, -0.1],
        ],
    )
    assert sorted(_molecule_graph(stretched).edges) == [(1, 2), (1, 3), (1, 4)]
