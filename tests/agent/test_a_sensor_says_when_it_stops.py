"""A floor a sensor stops at is declared, and reaching it is recorded.

The basin sensor and the same-structure sensor both measure a Kabsch
heavy-atom RMSD, and both stop below three heavy atoms because two
structures with fewer align exactly by construction. Neither said so.
That silence hid a real defect for a whole round: the same-structure
sensor joined a consumer to its producer by digest, which a handoff file
never matches, and every workspace small enough to be under the floor
reported nothing at all rather than reporting that it had not looked
(PySCF round 2, 2026-09-13).

The floor is a registered host policy now, and these are its boundary
probes: one heavy atom, two, and three, driving the real functions.
"""

from __future__ import annotations

from types import SimpleNamespace

import pytest

from chemsmart.agent._contracts import TrustedArtifactRefV1, file_sha256
from chemsmart.agent.tool_runtime import (
    SENSOR_HEAVY_ATOM_FLOOR,
    _basin_sensor_inputs,
    _same_structure_observations,
)

pytestmark = pytest.mark.capability("policy:sensor_heavy_atom_floor")

#: (name, symbols, positions) at one, two and three heavy atoms.
_MOLECULES = (
    (
        "water",
        ("O", "H", "H"),
        ((0.0, 0.0, 0.117), (0.0, 0.757, -0.469), (0.0, -0.757, -0.469)),
    ),
    (
        "hydrogen peroxide",
        ("O", "O", "H", "H"),
        (
            (0.0, 0.734, -0.053),
            (0.0, -0.734, -0.053),
            (0.839, 0.885, 0.423),
            (-0.839, -0.885, 0.423),
        ),
    ),
    (
        "hydrogen cyanide",
        ("H", "C", "N"),
        ((-1.064, 0.0, 0.0), (0.0, 0.0, 0.0), (1.156, 0.0, 0.0)),
    ),
)


def _heavy(symbols):
    return sum(1 for symbol in symbols if symbol != "H")


def _output(symbols, positions, *, shift=0.0):
    moved = tuple((x + shift, y, z) for x, y, z in positions)
    return SimpleNamespace(
        molecule=SimpleNamespace(
            chemical_symbols=list(symbols),
            positions=[list(row) for row in moved],
        ),
        final_energy=-76.0,
    )


def _artifact(tmp_path, name, symbols, positions):
    path = tmp_path / f"{name.replace(' ', '-')}.xyz"
    lines = [str(len(symbols)), name]
    for symbol, (x, y, z) in zip(symbols, positions):
        lines.append(f"{symbol} {x:.6f} {y:.6f} {z:.6f}")
    path.write_text("\n".join(lines) + "\n")
    return TrustedArtifactRefV1(
        artifact_id=f"input-{name}",
        kind="geometry_xyz",
        sha256=file_sha256(str(path)),
        size_bytes=path.stat().st_size,
        path=str(path),
        cli_value=str(path),
    )


@pytest.mark.parametrize("name,symbols,positions", _MOLECULES)
def test_the_basin_block_says_whether_it_measured(
    name, symbols, positions, tmp_path
):
    block = _basin_sensor_inputs(
        _artifact(tmp_path, name, symbols, positions),
        _output(symbols, positions, shift=0.05),
        "opt",
    )
    if _heavy(symbols) < SENSOR_HEAVY_ATOM_FLOOR:
        assert block.get("heavy_atom_rmsd_floor_applied") is True, name
        assert block.get("heavy_atom_count") == _heavy(symbols), name
        assert "heavy_atom_rmsd_angstrom" not in block, name
    else:
        assert "heavy_atom_rmsd_floor_applied" not in block, name
        assert block.get("heavy_atom_rmsd_angstrom") is not None, name


@pytest.mark.parametrize("name,symbols,positions", _MOLECULES)
def test_the_same_structure_sensor_reports_a_comparison_it_could_not_make(
    name, symbols, positions
):
    observations = _same_structure_observations(
        {}, "node-a", _output(symbols, positions)
    )
    below = _heavy(symbols) < SENSOR_HEAVY_ATOM_FLOOR
    if below:
        assert len(observations) == 1, name
        (observation,) = observations
        assert (
            observation["observation"] == "same_structure_comparison_not_made"
        )
        assert observation["heavy_atom_count"] == _heavy(symbols)
    else:
        # Above the floor with no receipts to compare against, the
        # sensor is silent because it looked and found no sibling.
        assert observations == (), name


def test_the_floor_is_the_declared_policy_and_not_a_second_number():
    from chemsmart.agent.rules import HOST_POLICIES

    owners = {policy_id: owner for policy_id, owner, _, _ in HOST_POLICIES}
    assert (
        owners["sensor_heavy_atom_floor"]
        == "chemsmart.agent.tool_runtime.SENSOR_HEAVY_ATOM_FLOOR"
    )
    assert SENSOR_HEAVY_ATOM_FLOOR == 3
