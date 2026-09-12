"""A two-node minimum is characterised by the Hessian that consumed it.

ORCA's opt+freq is one node, so a geometry read from it stands beside
the frequencies that characterise it. PySCF's and xTB's optimisation and
Hessian are two nodes joined by a validated handoff, and a number read
from the optimisation was worded "delivered from a result whose
stationary point is uncharacterised (no frequencies printed)" while the
Hessian on that exact geometry -- byte-identical through the handoff --
validated as a minimum in the same goal (PySCF round, g1 and g2,
2026-09-12). True per result, false per goal. The projection now joins
the Hessian's characterisation onto the optimisation whose reached
geometry it consumed, in the run stream and in the workspace record.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from chemsmart.agent.driver import _analysis_delivery
from chemsmart.agent.tool_runtime import _neutral_sensor_facts
from chemsmart.agent.workspace_record import (
    record_run,
    uncharacterised_artifacts,
)

from .test_the_sensors_read_pyscf_through_its_reader import _artifacts

pytestmark = pytest.mark.capability("program_jobtype:pyscf:cpu:hess")

OPT = "a" * 64
HESS = "b" * 64
EXTRACTION = "c" * 64


def _verified(node_id, jobtype, block, digest):
    return {
        "kind": "program_result_verified",
        "payload": {
            "node_id": node_id,
            "receipt_sha256": digest[:32] + "0" * 32,
            "status": "valid",
            "record": {
                "node_id": node_id,
                "program": "pyscf",
                "jobtype": jobtype,
                "state": "valid",
                "observations": {
                    "jobtype": jobtype,
                    "program": "pyscf",
                    "pyscf": block,
                },
                "output_artifacts": [{"kind": "pyscf_hdf5", "sha256": digest}],
            },
        },
    }


def _stream(path: Path, *, handoff: bool) -> Path:
    rows = [
        _verified("opt", "opt", {"vibrational_mode_count": 0}, OPT),
    ]
    if handoff:
        rows.append(
            {
                "kind": "optimized_geometry_handed_off",
                "payload": {
                    "producer_node_id": "opt",
                    "consumer_node_id": "hess",
                    "status": "validated_handoff",
                    "receipt_sha256": "d" * 64,
                },
            }
        )
    rows += [
        _verified(
            "hess",
            "hess",
            {
                "vibrational_mode_count": 6,
                "consequential_imaginary_mode_count": 0,
            },
            HESS,
        ),
        {
            "kind": "result_quantities_extracted",
            "payload": {
                "receipt_sha256": EXTRACTION,
                "artifact_sha256": OPT,
                "record": {"artifact_id": "result.opt", "program": "pyscf"},
            },
        },
        {
            "kind": "analysis_claims_recorded",
            "payload": {
                "receipt_sha256": "e" * 64,
                "record": {
                    "claims": [
                        {
                            "claim_id": "cn",
                            "quantity_id": "cn",
                            "display_value": 1.46,
                            "display_unit": "angstrom",
                            "dimension": [0, 1, 0, 0, 0, 0],
                            "source_receipt_sha256": EXTRACTION,
                        }
                    ]
                },
            },
        },
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join(json.dumps(row) for row in rows) + "\n", encoding="utf-8"
    )
    return path


def test_the_run_stream_joins_the_hessian_onto_the_optimisation(tmp_path):
    joined = _analysis_delivery(
        _stream(tmp_path / "joined.jsonl", handoff=True)
    )
    assert joined.uncharacterised_source_quantity_ids == ()
    alone = _analysis_delivery(
        _stream(tmp_path / "alone.jsonl", handoff=False)
    )
    assert alone.uncharacterised_source_quantity_ids == ("cn",)


def test_the_workspace_record_joins_the_hessian_onto_the_optimisation(
    tmp_path,
):
    workspace = tmp_path / "ws"
    record_run(
        workspace,
        goal_id="g",
        cycle=1,
        run_events_path=_stream(tmp_path / "joined.jsonl", handoff=True),
        run="goals/g/runs/cycle-1",
    )
    assert OPT not in uncharacterised_artifacts(workspace)
    other = tmp_path / "ws-alone"
    record_run(
        other,
        goal_id="g",
        cycle=1,
        run_events_path=_stream(tmp_path / "alone.jsonl", handoff=False),
        run="goals/g/runs/cycle-1",
    )
    assert OPT in uncharacterised_artifacts(other)


@pytest.mark.parametrize(
    ("case", "label", "jobtype", "count"),
    [
        ("water_hess", "water_hess", "hess", 3),
        ("water_opt", "water_opt", "opt", 0),
    ],
)
def test_the_neutral_block_says_how_many_modes_were_printed(
    case, label, jobtype, count
):
    block, _inputs = _neutral_sensor_facts(
        program="pyscf",
        jobtype=jobtype,
        multiplicity=1,
        output_artifacts=_artifacts(case, label),
        expected_input_artifact=None,
        expected_root_artifact=None,
    )
    assert block.get("vibrational_mode_count") == count
