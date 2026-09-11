"""The sensors also walk from the goal's own start.

NOVEL-3 ino3 (2026-09-05): a quartet restarted from a distorted geometry
relaxed a further 0.1 A from that input, so the walk sensor stayed
silent, while against the goal's original PH3 model the structure had
moved 0.6 A and closed an S...S contact to 2.81 A -- an incipient
disulfide that surfaced nowhere. A repair that restarts from a reached
or displaced geometry must not extinguish the anomaly it answers, so
the approved binding carries the root, the executor locates its bytes,
and the sensors compare against it when the input comparison is silent.
"""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from chemsmart.agent.execution import (
    ApprovedNodeBindingV1,
    _approved_node_binding_body,
)
from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1
from tests.agent.test_a_walk_to_another_basin_is_an_observation import (
    _SN2,
    _artifact,
    _final_geometry,
    _write_xyz,
)


def _shifted(atoms, dx):
    """The first heavy atom moved by dx along x, every other atom left in
    place: a walk of one atom, which no alignment removes (a shift of
    every heavy atom would be a rigid translation and no walk at all)."""

    moved = []
    done = False
    for symbol, (x, y, z) in atoms:
        if symbol != "H" and not done:
            moved.append((symbol, np.array([x + dx, y, z])))
            done = True
        else:
            moved.append((symbol, np.array([x, y, z])))
    return moved


@pytest.mark.capability("predicate:geometry.heavy_atom_rmsd_ge_0.3")
def test_a_silent_input_walk_still_fires_against_the_root(tmp_path):
    final = _final_geometry()
    # The node's input is the result's own geometry: no walk from it.
    start = _write_xyz(tmp_path / "start.xyz", final)
    # The goal's root is 0.5 A away on every heavy atom.
    root = _write_xyz(tmp_path / "root.xyz", _shifted(final, 1.5))
    evaluation = CommandCompiledToolHostV1._evaluate_execution_outputs(
        program="orca",
        jobtype="ts",
        charge=-1,
        multiplicity=1,
        expected_input_artifact=_artifact(start, "geometry_xyz", "start"),
        expected_root_artifact=_artifact(root, "geometry_xyz", "goal-root"),
        output_artifacts=(_artifact(_SN2, "orca_output", "result"),),
        exit_status=0,
    )
    signals = {item["signal_id"]: item for item in evaluation.anomalies}
    walk = signals["geometry.heavy_atom_rmsd_ge_0.3"]
    assert walk["reference"] == "goal_root"
    assert walk["reference_artifact_id"] == "goal-root"
    assert walk["heavy_atom_rmsd_angstrom"] >= 0.3

    # A root identical to the input adds nothing.
    same = CommandCompiledToolHostV1._evaluate_execution_outputs(
        program="orca",
        jobtype="ts",
        charge=-1,
        multiplicity=1,
        expected_input_artifact=_artifact(start, "geometry_xyz", "start"),
        expected_root_artifact=_artifact(start, "geometry_xyz", "start"),
        output_artifacts=(_artifact(_SN2, "orca_output", "result"),),
        exit_status=0,
    )
    assert "geometry.heavy_atom_rmsd_ge_0.3" not in {
        item["signal_id"] for item in same.anomalies
    }


def _binding(**extra):
    fields = {
        "node_id": "opt",
        "program": "orca",
        "engine": "cpu",
        "jobtype": "opt",
        "project_artifact_sha256": "e" * 64,
        "settings_sha256": "5" * 64,
        "charge": 0,
        "multiplicity": 1,
        "input_mode": "initial",
        "initial_artifact_id": "displaced",
        "initial_artifact_sha256": "d" * 64,
        "scientific_identity_sha256": "b" * 64,
        "producer_edge_sha256": "",
    }
    fields.update(extra)
    return ApprovedNodeBindingV1(**fields)


def test_the_binding_carries_the_root_and_keeps_old_digests_without_it():
    plain = _approved_node_binding_body(_binding())
    assert "root_artifact_sha256" not in plain
    rooted = _approved_node_binding_body(
        _binding(
            root_artifact_id="quartet-start", root_artifact_sha256="a" * 64
        )
    )
    assert rooted["root_artifact_id"] == "quartet-start"
    assert rooted["root_artifact_sha256"] == "a" * 64


def test_the_planning_host_walks_a_displaced_input_to_the_runs_input(
    tmp_path,
):
    from tests.agent.test_a_geometry_edit_is_arithmetic_not_judgement import (
        _host_with,
    )

    xyz = "3\nwater\nO 0.0 0.0 0.1173\nH 0.0 0.7572 -0.4692\nH 0.0 -0.7572 -0.4692\n"
    host = _host_with(tmp_path, xyz, "quartet-start")
    original = host.artifacts["quartet-start"]
    result_sha = "c" * 64
    displaced = _artifact(
        _write_xyz(tmp_path / "displaced.xyz", _final_geometry()),
        "geometry_xyz",
        "displaced",
    )
    host.artifacts["displaced"] = displaced
    host.mode_displacements[displaced.sha256] = SimpleNamespace(
        result_sha256=result_sha
    )
    evidence = tmp_path / "evidence"
    stream = evidence / ".chemsmart-agent" / "goals" / "g" / "runs" / "c1"
    stream.mkdir(parents=True)
    (stream / "events.jsonl").write_text(
        json.dumps(
            {
                "kind": "program_execution_observed",
                "payload": {
                    "record": {
                        "input_artifact_sha256": original.sha256,
                        "output_artifacts": [{"sha256": result_sha}],
                    }
                },
            }
        )
        + "\n"
    )
    host.run_evidence_root = Path(evidence)
    root = host._root_artifact_for("opt", displaced)
    assert root is not None and root.sha256 == original.sha256
    assert root.artifact_id == "quartet-start"
    # The goal's own start has no root.
    assert host._root_artifact_for("opt", original) is None


def test_the_executor_locates_the_root_beside_the_initial_artifacts(tmp_path):
    from chemsmart.agent.executor import _approved_initial_artifacts

    workspace = tmp_path / "ws"
    workspace.mkdir()
    start = _write_xyz(workspace / "displaced.xyz", _final_geometry())
    root = _write_xyz(workspace / "root.xyz", _shifted(_final_geometry(), 1.5))
    from chemsmart.agent.execution import file_sha256

    approval = SimpleNamespace(
        node_bindings=(
            _binding(
                initial_artifact_sha256=file_sha256(start),
                root_artifact_id="quartet-start",
                root_artifact_sha256=file_sha256(root),
            ),
        )
    )
    artifacts = _approved_initial_artifacts(workspace, approval)
    assert artifacts["quartet-start"].sha256 == file_sha256(root)
    assert artifacts["displaced"].sha256 == file_sha256(start)
