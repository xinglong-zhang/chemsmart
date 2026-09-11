"""A rescan says what it knows, and no more.

The quarantine on a failed result is a property of the run, derived from
the sealed stream. The rescan that a later cycle performs is a property
of the file. The two never met: a run typed a failure in one cycle came
back in the next as an ordinary registered result whose record called
itself a validated native result, and nothing carried the failure
across. The bytes stay readable and claimable -- they are often the
finding -- and the record stops forgetting.
"""

from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path

import pytest

from chemsmart.agent.live_session import _scan_result_artifacts

_SADDLE = (
    Path(__file__).resolve().parents[1]
    / "data"
    / "ORCATests"
    / "outputs"
    / "sn2_ts.out"
)


def _workspace(tmp_path: Path, *, recorded_state: str | None) -> Path:
    workspace = tmp_path / "ws"
    (workspace / "nodes" / "opt-a").mkdir(parents=True)
    target = workspace / "nodes" / "opt-a" / "sn2_ts.out"
    shutil.copy(_SADDLE, target)
    if recorded_state is None:
        return workspace
    digest = hashlib.sha256(target.read_bytes()).hexdigest()
    stream = workspace / ".chemsmart-agent" / "goals" / "g" / "runs" / "c1"
    stream.mkdir(parents=True)
    (stream / "events.jsonl").write_text(
        json.dumps(
            {
                "kind": "program_result_verified",
                "payload": {
                    "receipt_sha256": "b" * 64,
                    "node_id": "opt-a",
                    "record": {
                        "node_id": "opt-a",
                        "state": recorded_state,
                        "output_artifacts": [{"sha256": digest}],
                    },
                },
            }
        )
        + "\n",
        encoding="utf-8",
    )
    return workspace


@pytest.mark.capability("tool:inspect_run")
def test_a_rescan_never_calls_itself_validated(tmp_path):
    (observation,) = _scan_result_artifacts(
        _workspace(tmp_path, recorded_state=None)
    )
    record = observation.public_record()
    assert "validated" not in record["provenance_status"]
    assert (
        record["provenance_status"] == "workspace_exact_parsed_native_result"
    )
    assert record["recorded_terminal_state"] == ""
    assert record["provenance_limitations"] == ()


@pytest.mark.capability("tool:inspect_run")
def test_a_recorded_failure_rides_the_rescan(tmp_path):
    (observation,) = _scan_result_artifacts(
        _workspace(tmp_path, recorded_state="invalid")
    )
    record = observation.public_record()
    assert record["recorded_terminal_state"] == "invalid"
    assert any(
        "invalid" in item and "opt-a" in item
        for item in record["provenance_limitations"]
    )
    # It stays readable: the artifact is registered exactly as before.
    assert record["artifact_class"] == "orca_output"
    assert record["jobtype"] == "ts"


@pytest.mark.capability("tool:inspect_run")
def test_a_valid_run_leaves_the_record_alone(tmp_path):
    (observation,) = _scan_result_artifacts(
        _workspace(tmp_path, recorded_state="valid")
    )
    record = observation.public_record()
    assert record["recorded_terminal_state"] == ""
    assert record["provenance_limitations"] == ()
