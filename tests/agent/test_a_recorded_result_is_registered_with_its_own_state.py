"""A result registered from a run record carries the state the record
holds, never a default.

``_recorded_result_artifacts`` read charge and multiplicity from the
program's observation block and defaulted to 0 and 1 when the block had
none -- and PySCF's block had none (its verification branch wrote only a
native-failure entry), so both PySCF results of an ORCA-opt-to-PySCF-sp
run registered as neutral singlets while their true state sat in the
same record's ``result_validation``.  "A default must never outrank a
binding" is the xTB audit's lesson arriving at registration.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from chemsmart.agent._contracts import file_sha256
from chemsmart.agent.live_session import _recorded_result_artifacts

FIXTURE = (
    Path(__file__).resolve().parents[1]
    / "data"
    / "PySCFTests"
    / "outputs"
    / "hydroxyl_sp"
    / "hydroxyl_sp_gas_phase.h5"
)


def _workspace_with_record(tmp_path: Path, observations: dict) -> Path:
    workspace = tmp_path / "workspace"
    stream = workspace / ".chemsmart-agent" / "runs" / "run-1"
    stream.mkdir(parents=True)
    event = {
        "kind": "program_result_verified",
        "payload": {
            "record": {
                "program": "pyscf",
                "jobtype": "sp",
                "state": "valid",
                "engine": "cpu",
                "output_artifacts": [
                    {
                        "kind": "pyscf_hdf5",
                        "path": str(FIXTURE),
                        "sha256": file_sha256(FIXTURE),
                    }
                ],
                "observations": observations,
            }
        },
    }
    (stream / "events.jsonl").write_text(json.dumps(event) + "\n")
    return workspace


@pytest.mark.capability("tool:inspect_run")
def test_the_validators_record_of_the_state_is_used_when_the_block_has_none(
    tmp_path,
):
    workspace = _workspace_with_record(
        tmp_path,
        {"pyscf": {}, "result_validation": {"charge": 0, "multiplicity": 2}},
    )
    (item,) = _recorded_result_artifacts(workspace)
    assert item.program == "pyscf"
    assert (item.charge, item.multiplicity) == (0, 2)
    record = item.public_record()
    assert record["multiplicity"] == 2
    assert not any(
        "no charge or multiplicity" in line
        for line in record["provenance_limitations"]
    )


@pytest.mark.capability("tool:inspect_run")
def test_a_record_with_no_state_registers_none_and_says_so(tmp_path):
    workspace = _workspace_with_record(tmp_path, {"pyscf": {}})
    (item,) = _recorded_result_artifacts(workspace)
    assert (
        item.charge is None and item.multiplicity is None
    ), "a default of 0 and 1 is a binding nobody made"
    assert any(
        "no charge or multiplicity" in line
        for line in item.public_record()["provenance_limitations"]
    )
