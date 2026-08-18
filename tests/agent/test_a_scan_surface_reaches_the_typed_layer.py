"""A relaxed scan is a surface, and no surface reached the typed layer.

ORCA had no scan parser at all -- zero occurrences of the word in its output
module -- so a torsional profile could only be read by a model looking at a
log, which is the hallucination channel the typed hub exists to close. The gap
was systemic rather than ORCA's alone: nothing anywhere carried a multi-point
series. `energies` is a flat concatenation of every SCF print in the file, with
no companion coordinate axis, and the trajectory selectors are gated to IRC and
expose only the first and last frame.

The fixture is an excerpt of a real ORCA 6.1.1 relaxed scan run through the
ChemSmart CLI on this host, kept verbatim; the complete 1 MB output is
preserved under the campaign's program-coverage directory. The excerpt parses
identically to the whole file, which is checked by the values below matching
the `.relaxscanact.dat` sidecar ORCA wrote alongside it.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

from chemsmart.agent._contracts import TrustedArtifactRefV1
from chemsmart.agent.postprocessing import extract_trusted_result_quantities
from chemsmart.analysis.result_quantities import QuantitySelectorV1
from chemsmart.analysis.result_readers import reader_for
from chemsmart.io.orca.output import ORCAOutput

_SCAN = "tests/data/ORCATests/outputs/hooh_relaxed_scan_excerpt.out"
_NOT_A_SCAN = "tests/data/ORCATests/outputs/CO2.out"

#: Exactly what ORCA wrote to the .dat sidecar for this run.
_SURFACE = (
    (0.0, -151.35504813),
    (30.0, -151.35532645),
    (60.0, -151.35587843),
    (90.0, -151.35525734),
    (120.0, -151.35112470),
    (150.0, -151.34532814),
    (180.0, -151.34272082),
)
_HARTREE_TO_KCAL = 627.5094740631


def _artifact(path):
    resolved = Path(path).resolve()
    return TrustedArtifactRefV1(
        artifact_id="scan-result",
        kind=reader_for("orca").artifact_kind,
        sha256=hashlib.sha256(resolved.read_bytes()).hexdigest(),
        size_bytes=resolved.stat().st_size,
        path=str(resolved),
        cli_value=str(resolved),
    )


def test_the_driven_coordinate_is_recovered_and_renumbered():
    """ORCA counts atoms from zero; every other surface here counts from one."""

    coordinate = ORCAOutput(_SCAN).scan_coordinate

    assert coordinate["kind"] == "dihedral"
    assert coordinate["atoms"] == (4, 3, 2, 1)
    assert coordinate["start"] == 0.0
    assert coordinate["stop"] == 180.0
    assert coordinate["points"] == 7


def test_the_surface_matches_what_orca_wrote():
    profile = ORCAOutput(_SCAN).scan_profile

    assert len(profile) == len(_SURFACE)
    for point, (coordinate, energy) in zip(profile, _SURFACE):
        assert point["coordinate"] == pytest.approx(coordinate)
        assert point["energy"] == pytest.approx(energy, abs=1e-8)


def test_the_selectors_deliver_the_surface_through_typed_extraction():
    """The production shape: the typed chain, not the parser directly."""

    receipt = extract_trusted_result_quantities(
        artifact=_artifact(_SCAN),
        program="orca",
        selectors=(
            QuantitySelectorV1(
                quantity_id="angles", selector="scan_coordinate_values"
            ),
            QuantitySelectorV1(quantity_id="e", selector="scan_energies"),
        ),
    )
    values = {item.quantity_id: item for item in receipt.quantities}

    assert list(values["angles"].value) == [point[0] for point in _SURFACE]
    assert values["e"].unit == "hartree"
    assert len(values["e"].value) == len(_SURFACE)


def test_the_barrier_is_reachable_from_the_existing_operations():
    """Height needs no new arithmetic once the surface is typed."""

    energies = [point["energy"] for point in ORCAOutput(_SCAN).scan_profile]
    barrier = (max(energies) - min(energies)) * _HARTREE_TO_KCAL

    # Hydrogen peroxide's torsional profile has one high and one low barrier;
    # this is the high one, and it is the number a session would report.
    assert barrier == pytest.approx(8.26, abs=0.05)


def test_the_step_count_is_independent_of_the_converged_surface():
    """How far a run got and how much converged are different questions."""

    assert ORCAOutput(_SCAN).scan_step_count == 7


def test_a_job_that_drove_nothing_reports_nothing():
    """No false surface on the ordinary optimisation everything else runs."""

    output = ORCAOutput(_NOT_A_SCAN)

    assert output.scan_coordinate is None
    assert output.scan_profile == ()
    assert output.scan_step_count == 0
