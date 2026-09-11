"""Host-owned evidence reaches any recorded result, verified by digest.

REACH-1 ino3 (2026-09-06): the seeded record and the run streams named
results from an earlier window by content digest; the session read the
digests off inspect_run and was refused ten times in 22 seconds while
the author's answer sat in one of the files. The bootstrap registered
only files under this workspace, outside its private root, and only
after parsing them. Owner ruling: any recorded result, by digest.
"""

from __future__ import annotations

import json

import pytest

from chemsmart.agent._contracts import file_sha256
from chemsmart.agent.live_session import (
    _scan_result_artifacts,
    discover_registered_result_artifacts,
)
from chemsmart.agent.workspace_record import (
    record_run,
    render_workspace_record,
)
from chemsmart.io.orca.output import ORCAOutput

_SP = (
    "|  1> ! sp wb97x-d3bj def2-svp\n"
    "|  2> * xyz 0 1\n|  3> O 0.0 0.0 0.1173\n|  4> H 0.0 0.7572 -0.4692\n"
    "|  5> H 0.0 -0.7572 -0.4692\n|  6> *\n****END OF INPUT****\n"
    "FINAL SINGLE POINT ENERGY     -76.400000000000\n"
    "                             ****ORCA TERMINATED NORMALLY****\n"
)


def _verified_stream(path, *, node_id, output, sha256, state="valid"):
    execution = {
        "kind": "program_execution_observed",
        "payload": {
            "node_id": node_id,
            "record": {
                "node_id": node_id,
                "input_artifact_sha256": "d" * 64,
                "started_at": "2026-09-06T12:00:00+00:00",
                "finished_at": "2026-09-06T12:00:05+00:00",
                "output_artifacts": [
                    {
                        "artifact_id": f"result.{node_id}.1",
                        "kind": "program_output",
                        "sha256": sha256,
                        "path": str(output),
                        "size_bytes": output.stat().st_size,
                        "cli_value": str(output),
                    }
                ],
            },
        },
    }
    verified = {
        "kind": "program_result_verified",
        "payload": {
            "node_id": node_id,
            "receipt_sha256": "5" * 64,
            "record": {
                "node_id": node_id,
                "program": "orca",
                "jobtype": "sp",
                "engine": "cpu",
                "state": state,
                "receipt_sha256": "5" * 64,
                "input_artifact_sha256": "d" * 64,
                "project_artifact_sha256": "e" * 64,
                "observations": {
                    "orca": {
                        "charge": 0,
                        "multiplicity": 1,
                        "energy_hartree": -76.4,
                        "normal_termination": True,
                    }
                },
                "output_artifacts": execution["payload"]["record"][
                    "output_artifacts"
                ],
            },
        },
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join(json.dumps(item) for item in (execution, verified)) + "\n"
    )


@pytest.mark.capability("tool:inspect_run")
def test_a_result_named_by_the_record_opens_by_id_without_a_parse(
    tmp_path, monkeypatch
):
    earlier = tmp_path / "earlier-window" / "nodes" / "cat-opt"
    earlier.mkdir(parents=True)
    output = earlier / "cat_opt.out"
    output.write_text(_SP)
    digest = file_sha256(output)
    workspace = tmp_path / "ws"
    stream = (
        workspace
        / ".chemsmart-agent"
        / "goals"
        / "novel3-goal-x"
        / "runs"
        / "cycle-1"
        / "events.jsonl"
    )
    _verified_stream(stream, node_id="cat-opt", output=output, sha256=digest)

    def no_parse(self):
        raise AssertionError("registration must not parse the structure")

    monkeypatch.setattr(ORCAOutput, "final_structure", property(no_parse))
    observations = _scan_result_artifacts(workspace)
    (found,) = [o for o in observations if o.artifact.sha256 == digest]
    assert found.artifact.artifact_id == f"orca-result-{digest[:16]}"
    assert found.artifact.kind == "orca_output"
    assert found.provenance_status == "recorded_run_result_by_digest"
    assert found.scientific_validation_state == "recorded_valid"
    assert found.charge == 0 and found.multiplicity == 1
    assert found.artifact.artifact_id in {
        a.artifact_id for a in discover_registered_result_artifacts(workspace)
    }, "the executor rebuilds the same id"


def test_bytes_that_no_longer_hash_to_the_record_are_not_registered(tmp_path):
    output = tmp_path / "elsewhere" / "x.out"
    output.parent.mkdir()
    output.write_text(_SP)
    workspace = tmp_path / "ws"
    stream = (
        workspace
        / ".chemsmart-agent"
        / "goals"
        / "g"
        / "runs"
        / "cycle-1"
        / "events.jsonl"
    )
    _verified_stream(stream, node_id="x", output=output, sha256="f" * 64)
    assert all(
        o.artifact.sha256 != "f" * 64
        for o in _scan_result_artifacts(workspace)
    )


def test_the_record_view_shows_the_id_and_the_run_outcome_carries_it(
    tmp_path,
):
    """The record view names the id beside each result; a run outcome
    names it beside each evidence digest. The outcome's ids are derived
    from the same digests its existing evidence field carries, so the
    rendered view is pinned here on a recorded stream and the outcome's
    comprehension on its own record type."""

    from chemsmart.agent.terminal_states import NodeTerminalStateV1

    output = tmp_path / "nodes" / "cat-opt" / "cat_opt.out"
    output.parent.mkdir(parents=True)
    output.write_text(_SP)
    digest = file_sha256(output)
    workspace = tmp_path / "ws"
    workspace.mkdir()
    stream = (
        workspace
        / ".chemsmart-agent"
        / "goals"
        / "g"
        / "runs"
        / "cycle-1"
        / "events.jsonl"
    )
    _verified_stream(stream, node_id="cat-opt", output=output, sha256=digest)
    record_run(
        workspace,
        goal_id="g",
        cycle=1,
        run_events_path=stream,
        run="goals/g/runs/cycle-1",
        review_file=None,
    )
    (row,) = render_workspace_record(workspace)["results"]
    assert row["artifact_id"] == f"orca-result-{digest[:16]}"
    node = NodeTerminalStateV1(
        node_id="cat-opt",
        program="orca",
        jobtype="opt",
        state="validated",
        evidence_artifact_sha256s=(digest,),
        evidence_artifact_ids=(f"orca-result-{digest[:16]}",),
    )
    assert node.public_record()["evidence_artifact_ids"] == (
        f"orca-result-{digest[:16]}",
    )
