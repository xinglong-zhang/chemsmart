"""A distance, angle or dihedral measured over indices that are not a
bonded chain is observed, beside the receipt, never refused.

A woken session measured the sulfone's F-C-C-S torsion with the indices
one atom off and reported F-C-C-O (78-82 deg) as a shifted well; the
true torsion is 62-71 deg, and the receipt that carried the positions
carried the perceived bonds too (NOVEL-2 po2, 2026-09-04).
"""

from __future__ import annotations

from pathlib import Path

import pytest

from chemsmart.agent._contracts import TrustedArtifactRefV1, file_sha256
from chemsmart.agent.runtime.event_store import RuntimeEventStore
from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1

pytestmark = pytest.mark.capability("operation:dihedral")

_RESULT = Path("tests/data/ORCATests/outputs/CO2.out")


def _host(tmp_path):
    resolved = _RESULT.resolve()
    artifact = TrustedArtifactRefV1(
        artifact_id="result.co2",
        kind="orca_output",
        sha256=file_sha256(resolved),
        size_bytes=resolved.stat().st_size,
        path=str(resolved),
        cli_value=str(resolved),
    )
    host = CommandCompiledToolHostV1(
        event_store=RuntimeEventStore(
            tmp_path / "events.jsonl", session_id="geom"
        ),
        task_spec_sha256s=("a" * 64,),
        approved_workspace=tmp_path,
    )
    host.artifacts[artifact.artifact_id] = artifact
    return host


def _measure(host, indices):
    extraction = host.dispatch(
        turn_id="t1",
        tool_name="extract_result_quantities",
        arguments={
            "artifact_id": "result.co2",
            "program": "orca",
            "selectors": [
                {"quantity_id": "positions", "selector": "positions"},
                {"quantity_id": "symbols", "selector": "symbols"},
            ],
        },
    )["result"]
    receipt = extraction["receipt_sha256"]
    refs = [
        {
            "node_id": f"a{k}",
            "operation": "ref",
            "reference": "xyz",
            "indices": [index],
        }
        for k, index in enumerate(indices)
    ]
    return host.dispatch(
        turn_id="t2",
        tool_name="evaluate_quantity_expression",
        arguments={
            "expression_id": "measure",
            "inputs": [
                {
                    "input_id": "xyz",
                    "receipt_sha256": receipt,
                    "quantity_id": "positions",
                }
            ],
            "nodes": refs
            + [
                {
                    "node_id": "value",
                    "operation": "angle" if len(indices) == 3 else "distance",
                    "input_ids": [ref["node_id"] for ref in refs],
                }
            ],
            "output_node_ids": ["value"],
        },
    )


def test_a_non_bonded_chain_is_observed_beside_the_receipt(tmp_path):
    # The fixture's atom order is O, O, C: the two oxygens are not bonded.
    reply = _measure(_host(tmp_path), [0, 1])
    assert reply["status"] == "ok"
    (observation,) = reply["observations"]
    assert observation["kind"] == "geometry_operation_over_non_bond"
    assert observation["unbonded_pairs"] == [[0, 1]]
    assert observation["elements"] == ["O", "O"]
    assert "bonded order" in observation["meaning"]


def test_a_bonded_chain_is_silent(tmp_path):
    # O(0)-C(2)-O(1) is the bonded chain in the fixture's atom order.
    reply = _measure(_host(tmp_path), [0, 2, 1])
    assert reply["status"] == "ok"
    assert "observations" not in reply
