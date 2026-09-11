"""A scan whose extremum sits on the grid's edge is observed with its
numbers, never refused.

Four of four scan extrema in NOVEL-3 sat on a boundary: two product-side
ring-opening scans rose monotonically to their last point (po3) and two
hydrogen-transfer scans reached their edge (po1); every "barrier
position" the extremum operation returned was the end of the grid, and
no host word said so (2026-09-05).
"""

from __future__ import annotations

import pytest

from chemsmart.agent.tool_runtime import _scan_boundary_sensor

pytestmark = [
    pytest.mark.capability("selector:orca:scan:scan_energies"),
    # The sensor this test is about: a scan whose extremum sits on the
    # grid's own edge has standing as an observation. It was pinned by
    # nothing, because the anomaly signals were string literals with no
    # declaration to name.
    pytest.mark.capability("signal:scan.extremum_at_grid_boundary"),
]


def _profile(energies, start=1.4, step=0.1):
    return tuple(
        {"coordinate": start + step * index, "energy": energy}
        for index, energy in enumerate(energies)
    )


def test_a_monotone_wall_is_observed_with_its_span_and_last_step():
    # po3's shape: strictly rising from the product, 12 points.
    energies = [-1077.036 + 0.0135 * i for i in range(12)]
    observation = _scan_boundary_sensor(_profile(energies))
    assert observation["signal_id"] == "scan.extremum_at_grid_boundary"
    assert observation["maximum_at_boundary"] is True
    assert observation["maximum_index"] == 11
    assert observation["maximum_coordinate"] == pytest.approx(2.5)
    assert observation["minimum_at_boundary"] is True
    assert observation["monotone"] is True
    assert observation["span_kcal_mol"] == pytest.approx(93.2, abs=0.5)
    assert observation["last_step_kcal_mol"] > 0


def test_an_interior_maximum_is_silent_and_a_short_scan_says_nothing():
    energies = [-9.90, -9.95, -9.80, -10.00, -9.90, -9.95]
    assert _scan_boundary_sensor(_profile(energies)) is None
    assert _scan_boundary_sensor(_profile([-10.0, -9.9])) is None


def test_a_boundary_minimum_alone_is_observed_but_not_called_a_wall():
    # A dissociation curve: minimum at the start, maximum in the interior
    # is impossible for a monotone rise, so make the interior the maximum.
    energies = [-10.00, -9.90, -9.80, -9.85, -9.90, -9.88]
    observation = _scan_boundary_sensor(_profile(energies))
    assert observation["minimum_at_boundary"] is True
    assert observation["maximum_at_boundary"] is False
    assert observation["monotone"] is False
