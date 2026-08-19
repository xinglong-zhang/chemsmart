"""What one session computes, the next session must be shown.

The workspace bootstrap registers normally terminated ORCA outputs so a fresh
session can analyse work that already exists. Its jobtype filter was a
hand-listed set that predated the scan family, so a completed relaxed scan was
silently excluded from registration.

The cost was observed, not imagined. A follow-on session opened a workspace
holding five validated results, among them two 72-point surfaces that had cost
about four engine-hours; it was shown the three optimisations and neither
scan, and rationally replanned both from scratch. The stop-point repair then
named the replanned scans as approval blockers -- working exactly as designed,
one layer above the actual defect.

There is an extra turn in how this one hid: before the reader learned to
classify a %geom Scan block, a scan output slipped through this filter
mislabelled as `opt`. Correcting the label is what armed the exclusion.

The fixture is the real ORCA 6.1.1 scan excerpt extended with the same run's
verbatim closing block, so it carries what the registration filter actually
inspects: normal termination, a final energy, the molecule, method and basis,
charge and multiplicity, and jobtype `scan`.
"""

from __future__ import annotations

import shutil
from pathlib import Path

from chemsmart.agent.live_session import _scan_orca_result_artifacts

_SCAN = "tests/data/ORCATests/outputs/hooh_relaxed_scan_excerpt.out"
_PLAIN = "tests/data/ORCATests/outputs/CO2.out"


def _workspace(tmp_path, *sources):
    workspace = tmp_path / "workspace"
    nodes = workspace / "nodes" / "some-node"
    nodes.mkdir(parents=True)
    for source in sources:
        shutil.copy(source, nodes / Path(source).name)
    return workspace


def test_a_completed_scan_is_registered_for_the_next_session(tmp_path):
    observations = _scan_orca_result_artifacts(
        _workspace(tmp_path, _SCAN)
    )

    assert [item.jobtype for item in observations] == ["scan"]
    scan = observations[0]
    assert scan.method == "b3lyp"
    assert scan.basis == "def2-svp"
    assert (scan.charge, scan.multiplicity) == (0, 1)


def test_the_scan_sits_beside_the_families_already_registered(tmp_path):
    observations = _scan_orca_result_artifacts(
        _workspace(tmp_path, _SCAN, _PLAIN)
    )

    assert sorted(item.jobtype for item in observations) == ["opt", "scan"]


def test_the_registered_artifact_reads_back_as_a_surface(tmp_path):
    """Registration is only worth having if the typed layer can use it."""

    from chemsmart.io.orca.output import ORCAOutput

    observations = _scan_orca_result_artifacts(_workspace(tmp_path, _SCAN))

    parsed = ORCAOutput(observations[0].artifact.path)
    assert len(parsed.scan_profile) == 7
    assert parsed.scan_point_records[0]["index"] == 1
