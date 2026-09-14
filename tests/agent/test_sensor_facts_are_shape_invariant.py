"""One quantity, two shapes, one set of facts.

Every program's reader answers the same selector vocabulary, and nothing
makes two readers agree on the *container*: PySCF's HDF5 path returns
numpy arrays where a log parser returns Python lists, and a sensor is
written against whichever one its author had in hand. Two idioms in this
host are shape-sensitive and neither fails loudly on the shape it was
not written for: ``if values:`` raises ``ValueError`` on an array of
more than one element, and ``isinstance(value, (list, tuple))`` is
simply false for an array, so the value takes a scalar branch that
raises or silently records the wrong thing.

This is a metamorphic test: the same numbers in two containers must
produce equal facts. It drives the host's own enumeration of its sensor
functions, so a sensor added without a driver here fails rather than
going unmeasured.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from chemsmart.agent._contracts import canonical_data
from chemsmart.agent.tool_runtime import (
    SENSOR_FACT_FUNCTIONS,
    _basin_sensor_inputs,
    _excited_state_sensor_facts,
    _gradient_anomaly,
    _imaginary_mode_sensor_inputs,
    _same_structure_observations,
    _scan_boundary_sensor,
)

pytestmark = pytest.mark.capability("gate:stationary_point_order")

_SYMBOLS = ("C", "C", "O", "H", "H", "H", "H")
_POSITIONS = (
    (0.000, 0.000, 0.000),
    (1.480, 0.000, 0.000),
    (2.100, 1.060, 0.000),
    (-0.540, 0.940, 0.000),
    (-0.540, -0.940, 0.000),
    (1.900, -0.960, 0.300),
    (3.060, 0.960, 0.100),
)
_FREQUENCIES = (-512.4, 118.7, 604.2, 1420.9, 3011.5)


def _molecule(positions, *, native):
    return SimpleNamespace(
        chemical_symbols=(
            list(_SYMBOLS) if native else np.array(_SYMBOLS, dtype=object)
        ),
        positions=(
            [list(row) for row in positions]
            if native
            else np.asarray(positions, dtype=float)
        ),
    )


def _output(*, native):
    frequencies = (
        list(_FREQUENCIES) if native else np.asarray(_FREQUENCIES, dtype=float)
    )
    excitations = (
        [0.1523, 0.2011, 0.2688]
        if native
        else np.asarray([0.1523, 0.2011, 0.2688], dtype=float)
    )
    counts = (
        {
            "nstates_requested": 3,
            "nstates_obtained": 3,
            "roots_filtered": 0,
            "unconverged_roots": [1],
        }
        if native
        else {
            "nstates_requested": np.int64(3),
            "nstates_obtained": np.int64(3),
            "roots_filtered": np.int64(0),
            "unconverged_roots": np.asarray([1], dtype=np.int64),
        }
    )
    return SimpleNamespace(
        molecule=_molecule(_POSITIONS, native=native),
        vibrational_frequencies=frequencies,
        final_energy=-192.114 if native else np.float64(-192.114),
        td_stage=counts,
        excitation_energies=excitations,
        excited_state_record={
            "root": 1 if native else np.int64(1),
            "root_gap_to_ground_end_ev": (
                3.052 if native else np.float64(3.052)
            ),
        },
    )


def _scan_profile(*, native):
    values = [(0.0, -152.10), (30.0, -152.07), (60.0, -152.04)]
    if native:
        return [
            {"coordinate": float(x), "energy": float(y)} for x, y in values
        ]
    return [
        {"coordinate": np.float64(x), "energy": np.float64(y)}
        for x, y in values
    ]


def _input_artifact(tmp_path, *, native):
    from chemsmart.agent._contracts import TrustedArtifactRefV1, file_sha256

    path = tmp_path / f"input-{'list' if native else 'array'}.xyz"
    lines = [str(len(_SYMBOLS)), "probe"]
    for symbol, (x, y, z) in zip(_SYMBOLS, _POSITIONS):
        lines.append(f"{symbol} {x:.6f} {y:.6f} {z:.6f}")
    path.write_text("\n".join(lines) + "\n")
    return TrustedArtifactRefV1(
        artifact_id="probe-input",
        kind="geometry_xyz",
        sha256=file_sha256(str(path)),
        size_bytes=path.stat().st_size,
        path=str(path),
        cli_value=str(path),
    )


def _drive(name, tmp_path, *, native):
    output = _output(native=native)
    if name == "_scan_boundary_sensor":
        return _scan_boundary_sensor(_scan_profile(native=native))
    if name == "_basin_sensor_inputs":
        return _basin_sensor_inputs(
            _input_artifact(tmp_path, native=native), output, "opt"
        )
    if name == "_same_structure_observations":
        return _same_structure_observations({}, "node-a", output)
    if name == "_imaginary_mode_sensor_inputs":
        return _imaginary_mode_sensor_inputs(
            output, output.vibrational_frequencies
        )
    if name == "_gradient_anomaly":
        return _gradient_anomaly(0.0042 if native else np.float64(0.0042))
    if name == "_excited_state_sensor_facts":
        return _excited_state_sensor_facts(output)
    raise AssertionError(f"no driver for {name}")


_DRIVEN = (
    "_scan_boundary_sensor",
    "_basin_sensor_inputs",
    "_same_structure_observations",
    "_imaginary_mode_sensor_inputs",
    "_gradient_anomaly",
    "_excited_state_sensor_facts",
)


def test_every_enumerated_sensor_has_a_driver_here():
    assert {item.__name__ for item in SENSOR_FACT_FUNCTIONS} == set(_DRIVEN)


@pytest.mark.parametrize("name", _DRIVEN)
def test_the_same_numbers_in_two_containers_give_the_same_facts(
    name, tmp_path
):
    from_lists = _drive(name, tmp_path, native=True)
    from_arrays = _drive(name, tmp_path, native=False)
    assert canonical_data(from_lists) == canonical_data(from_arrays), name
