"""Two independent starts that converged on one structure is an observation.

A session that optimises twice and gets one answer reads the agreement
as convergence. Two live deliveries called two twist-boats "degenerate
chairs", and one dismissed a contradicting 1.85 kcal/mol as an artifact
(E4', 2026-09-03), while nothing in the host had ever compared two
sibling results to each other. The host measures it and says nothing
about what it means. It stays silent for the ordinary shape of a
workflow: a single point on an optimisation's own geometry, and two
nodes launched from the same geometry, are one structure by
construction.
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from chemsmart.agent.tool_runtime import _same_structure_observations
from chemsmart.analysis.result_readers import reader_for

_OUT = Path(__file__).resolve().parents[1] / "data" / "ORCATests" / "outputs"


def _receipt(name: str, node_id: str, *, inp: str = "", out: str = ""):
    path = _OUT / name
    return SimpleNamespace(
        node_id=node_id,
        state="valid",
        program="orca",
        input_artifact_sha256=inp,
        output_artifacts=(
            SimpleNamespace(path=str(path), sha256=out or f"sha.{name}"),
        ),
    )


def _output(name: str):
    return reader_for("orca").open_output(str(_OUT / name))


@pytest.mark.capability("signal:geometry.results_indistinguishable")
def test_two_independent_starts_on_one_structure_are_recorded():
    (found,) = _same_structure_observations(
        {"a": _receipt("phenol_pka_B.out", "opt-a", inp="start.a")},
        "opt-b",
        _output("phenol_pka_B_sp.out"),
        input_sha256="start.b",
        output_sha256s=("sha.b",),
    )
    assert found["signal_id"] == "geometry.results_indistinguishable"
    assert found["other_node_id"] == "opt-a"
    assert found["heavy_atom_rmsd_angstrom"] < 0.01
    assert "energy_difference_kcal_mol" in found


@pytest.mark.capability("signal:geometry.results_indistinguishable")
def test_two_different_structures_say_nothing():
    assert (
        _same_structure_observations(
            {"a": _receipt("fe3_doublet.out", "opt-a", inp="start.a")},
            "opt-b",
            _output("fe3_sextet.out"),
            input_sha256="start.b",
            output_sha256s=("sha.b",),
        )
        == ()
    )


@pytest.mark.capability("signal:geometry.results_indistinguishable")
def test_a_single_point_on_its_own_optimisation_is_no_surprise():
    """The commonest shape in the release: opt, then sp on that geometry."""

    parent = _receipt("phenol_pka_B.out", "opt-a", inp="start.a", out="geom.a")
    assert (
        _same_structure_observations(
            {"a": parent},
            "sp-b",
            _output("phenol_pka_B_sp.out"),
            input_sha256="geom.a",
            output_sha256s=("sha.b",),
        )
        == ()
    )


@pytest.mark.capability("signal:geometry.results_indistinguishable")
def test_two_siblings_from_one_geometry_are_no_surprise():
    """Three single points at three electronic states on one geometry."""

    assert (
        _same_structure_observations(
            {"a": _receipt("phenol_pka_B.out", "sp-a", inp="geom.shared")},
            "sp-b",
            _output("phenol_pka_B_sp.out"),
            input_sha256="geom.shared",
            output_sha256s=("sha.b",),
        )
        == ()
    )
