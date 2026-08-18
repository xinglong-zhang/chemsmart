"""A scan ends at a surface; carrying a structure off it is a choice.

ORCA writes the converged geometry of every relaxed-scan point beside its
output, 1-indexed as ``<stem>.001.xyz`` upward. Until now those files were
opaque ``program_output`` with nothing tying them to a coordinate or an energy,
and no tool could take one into a later calculation -- so a session could
compute a torsional profile, see exactly where the well was, and have no way to
optimise the structure sitting in it. Every observation of the cycle-038 paper
task died on that gap.

The repair is two-sided and the sides have different owners. The reader joins
each point to its geometry (`scan_point_records`) and the typed layer exposes
the indices (`scan_point_indices`), so the surface is legible; the model then
names the point it wants and `bind_scan_point_geometry` registers that
structure with its lineage. The host ranks nothing -- which point matters is
the scientist's judgement, and the mapping was verified against the real
25-point run before this was built: profile coordinate and file torsion agree
point for point.

Fixtures are the excerpt of a real ORCA 6.1.1 scan plus the per-point
geometries that run actually wrote, renamed only in stem to sit beside it.
"""

from __future__ import annotations

import hashlib
import shutil
from pathlib import Path

import pytest

from chemsmart.agent._contracts import ContractError, TrustedArtifactRefV1
from chemsmart.agent.runtime.event_store import RuntimeEventStore
from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1
from chemsmart.io.orca.output import ORCAOutput

_SCAN = "tests/data/ORCATests/outputs/hooh_relaxed_scan_excerpt.out"
_NOT_A_SCAN = "tests/data/ORCATests/outputs/CO2.out"


def _artifact(path, artifact_id="scan-result", kind="orca_output"):
    resolved = Path(path).resolve()
    return TrustedArtifactRefV1(
        artifact_id=artifact_id,
        kind=kind,
        sha256=hashlib.sha256(resolved.read_bytes()).hexdigest(),
        size_bytes=resolved.stat().st_size,
        path=str(resolved),
        cli_value=str(resolved),
    )


def _host(tmp_path, *artifacts):
    return CommandCompiledToolHostV1(
        event_store=RuntimeEventStore(
            tmp_path / "events.jsonl", session_id="s1"
        ),
        artifacts={item.artifact_id: item for item in artifacts},
        task_spec_sha256s=("a" * 64,),
        approved_workspace=tmp_path / "workspace",
    )


def _bind(host, **arguments):
    payload = host.dispatch(
        turn_id="t1",
        tool_name="bind_scan_point_geometry",
        arguments={"artifact_id": "scan-result", "program": "orca", **arguments},
    )
    assert payload["status"] == "ok"
    return payload["result"]


def test_every_point_reports_its_coordinate_energy_and_geometry():
    records = ORCAOutput(_SCAN).scan_point_records

    assert [record["index"] for record in records] == list(range(1, 8))
    assert records[0]["coordinate"] == 0.0
    assert records[-1]["coordinate"] == 180.0
    assert all(record["geometry_file"] for record in records)


def test_a_missing_point_file_is_reported_not_invented(tmp_path):
    """A truncated scan's converged points are data; absent files are gaps."""

    stem = tmp_path / "scan.out"
    shutil.copy(_SCAN, stem)
    for index in (1, 2, 3):  # deliberately omit 4..7
        shutil.copy(
            f"tests/data/ORCATests/outputs/hooh_relaxed_scan_excerpt.{index:03d}.xyz",
            tmp_path / f"scan.{index:03d}.xyz",
        )

    records = ORCAOutput(str(stem)).scan_point_records

    assert [r["geometry_file"] is not None for r in records] == (
        [True] * 3 + [False] * 4
    )


def test_a_run_that_scanned_nothing_has_no_points():
    assert ORCAOutput(_NOT_A_SCAN).scan_point_records == ()


def test_the_chosen_point_is_bound_with_its_lineage(tmp_path):
    host = _host(tmp_path, _artifact(_SCAN))

    result = _bind(host, point_index=3)

    bound = result["artifact"]
    assert bound["kind"] == "geometry_xyz"
    assert result["coordinate"] == 60.0
    assert result["point_count"] == 7
    assert result["selection_owner"] == "model"
    assert bound["artifact_id"] in host.artifacts


def test_the_bound_bytes_are_the_files_bytes(tmp_path):
    """Lineage means the exact structure ORCA wrote, digest and all."""

    host = _host(tmp_path, _artifact(_SCAN))

    result = _bind(host, point_index=5)

    source = Path(
        "tests/data/ORCATests/outputs/hooh_relaxed_scan_excerpt.005.xyz"
    ).read_bytes()
    assert result["artifact"]["sha256"] == hashlib.sha256(source).hexdigest()
    assert result["energy_hartree"] == pytest.approx(-151.35112470, abs=1e-8)


def test_a_point_that_does_not_exist_is_refused_by_range(tmp_path):
    host = _host(tmp_path, _artifact(_SCAN))

    with pytest.raises(ContractError, match="points 1 to 7"):
        host.dispatch(
            turn_id="t1",
            tool_name="bind_scan_point_geometry",
            arguments={
                "artifact_id": "scan-result",
                "point_index": 12,
                "program": "orca",
            },
        )


def test_a_result_with_no_surface_is_refused_by_name(tmp_path):
    host = _host(tmp_path, _artifact(_NOT_A_SCAN))

    with pytest.raises(ContractError, match="no scan surface"):
        host.dispatch(
            turn_id="t1",
            tool_name="bind_scan_point_geometry",
            arguments={
                "artifact_id": "scan-result",
                "point_index": 1,
                "program": "orca",
            },
        )


def test_the_indices_reach_the_typed_extraction_layer():
    """The model chooses from evidence, so the indices must be selectable."""

    from chemsmart.agent.postprocessing import extract_trusted_result_quantities
    from chemsmart.analysis.result_quantities import QuantitySelectorV1

    receipt = extract_trusted_result_quantities(
        artifact=_artifact(_SCAN),
        program="orca",
        selectors=(
            QuantitySelectorV1(
                quantity_id="idx", selector="scan_point_indices"
            ),
        ),
    )

    values = {item.quantity_id: item for item in receipt.quantities}
    assert [int(v) for v in values["idx"].value] == list(range(1, 8))
