"""A plan picks a normal mode by which atoms move, before the Hessian exists.

R9 xtb g3 (cephalexin) ran a second GFN2-xTB Hessian its session could not
assign. Naming the beta-lactam, amide and carboxyl C=O stretches means
reading each mode's per-atom participation, and a plan written before the
Hessian could select a row only by a literal index: "mode selection is a
post-execution data read; with no further agent cycle after this wave, the
adaptive assignment cannot be planned as a claim chain before execution".
The session was right not to carry indices over from the first Hessian:
the amide and carboxyl stretches swap rows between the two.

A null in ``ref``'s indices selects every element along that axis, so one
atom's participation in every mode is a column, and the mode in which a set
of atoms moves most is ``coordinate_at_maximum`` of their summed columns
against the frequencies -- planned with no index that names a mode.

An argmax alone does not say how far ahead of the rest the chosen mode was:
on both g3 Hessians the selected C=O shares were 0.96-0.98 and the runner-up
shares 0.36-0.63, low-frequency motions loading the same two atoms (R11
probe M, CUHK 2157069). So the host names the runner-up beside every
coordinate_at_maximum and coordinate_at_minimum, read from the same arrays.

The expected maximum is probe M's target step on this fixture, computed with
numpy outside the vocabulary before the change existed; the minimum is
computed the same way here, from the extraction reply's own values.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from chemsmart.agent._contracts import TrustedArtifactRefV1
from chemsmart.analysis.result_readers import reader_for

_RESULT = (
    Path(__file__).resolve().parents[1]
    / "data/XTBTests/outputs/acetaldehyde_hess/acetaldehyde_hess.out"
)
#: The carbonyl carbon and oxygen, zero-based in the result's atom order.
_CARBONYL = (2, 0)


def _host(tmp_path):
    from chemsmart.agent.runtime.event_store import RuntimeEventStore
    from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1

    path = _RESULT.resolve()
    artifact = TrustedArtifactRefV1(
        artifact_id="acetaldehyde-hess",
        kind=reader_for("xtb").artifact_kind,
        sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        size_bytes=path.stat().st_size,
        path=str(path),
        cli_value=str(path),
    )
    events = tmp_path / "events.jsonl"
    host = CommandCompiledToolHostV1(
        event_store=RuntimeEventStore(events, session_id="s"),
        artifacts={artifact.artifact_id: artifact},
        task_spec_sha256s=("a" * 64,),
        approved_workspace=tmp_path / "workspace",
    )
    return host, events


def _extract(host):
    reply = host.dispatch(
        turn_id="t-extract",
        tool_name="extract_result_quantities",
        arguments={
            "program": "xtb",
            "artifact_id": "acetaldehyde-hess",
            "selectors": [
                {
                    "quantity_id": "frequencies",
                    "selector": "vibrational_frequencies",
                },
                {
                    "quantity_id": "participation",
                    "selector": "vibrational_mode_atom_participation",
                },
            ],
        },
    )
    assert reply["status"] == "ok", reply
    result = reply["result"]
    values = {
        item["quantity_id"]: item["value"] for item in result["quantities"]
    }
    return result["receipt_sha256"], values


@pytest.mark.capability("operation:ref")
@pytest.mark.capability("operation:coordinate_at_maximum")
@pytest.mark.capability("operation:coordinate_at_minimum")
@pytest.mark.parametrize(
    "operation, reducer",
    (("coordinate_at_maximum", "max"), ("coordinate_at_minimum", "min")),
)
def test_a_plan_written_before_the_hessian_picks_a_mode_by_its_atoms(
    tmp_path, operation, reducer
):
    host, events = _host(tmp_path)
    receipt, values = _extract(host)

    # One column per atom: nothing in the plan names a mode.
    nodes = [
        {
            "node_id": f"atom{position}",
            "operation": "ref",
            "reference": "part",
            "indices": [None, atom],
        }
        for position, atom in enumerate(_CARBONYL)
    ]
    nodes += [
        {
            "node_id": "share",
            "operation": "add",
            "input_ids": ["atom0", "atom1"],
        },
        {
            "node_id": "nu",
            "operation": operation,
            "input_ids": ["share", "freq"],
        },
        {"node_id": "top", "operation": reducer, "input_ids": ["share"]},
    ]
    reply = host.dispatch(
        turn_id="t-select",
        tool_name="evaluate_quantity_expression",
        arguments={
            "expression_id": "carbonyl-stretch",
            "inputs": [
                {
                    "input_id": "part",
                    "receipt_sha256": receipt,
                    "quantity_id": "participation",
                },
                {
                    "input_id": "freq",
                    "receipt_sha256": receipt,
                    "quantity_id": "frequencies",
                },
            ],
            "nodes": nodes,
            "output_node_ids": ["nu", "top"],
        },
    )
    assert reply["status"] == "ok", reply
    outputs = {
        item["quantity_id"]: item["value"]
        for item in reply["result"]["outputs"]
    }

    # The oracle, computed outside the vocabulary: summed participation of
    # the two atoms in every mode, its extremum and the next in rank.
    table = np.asarray(values["participation"], dtype=float)
    frequencies = np.asarray(values["frequencies"], dtype=float)
    share = table[:, _CARBONYL[0]] + table[:, _CARBONYL[1]]
    order = np.argsort(share)
    if operation == "coordinate_at_maximum":
        order = order[::-1]
        # Probe M's target step, recorded before the change (EPISODE.md).
        assert (frequencies[order[0]], share[order[0]]) == (
            1798.58,
            0.9803607214428858,
        )
        assert (frequencies[order[1]], share[order[1]]) == (
            501.81,
            0.556488702259548,
        )
    chosen, runner_up = int(order[0]), int(order[1])

    assert outputs["nu"] == frequencies[chosen]
    assert outputs["top"] == share[chosen]

    # The runner-up travels with the selection: in the reply the session
    # reads, and in the record an approved chain leaves behind it.
    (seen,) = [
        item
        for item in reply.get("observations") or ()
        if item.get("kind") == "extremum_runner_up"
    ]
    assert seen["node_id"] == "nu"
    assert seen["operation"] == operation
    assert seen["selected"] == {
        "index": chosen,
        "coordinate": frequencies[chosen],
        "value": share[chosen],
    }
    assert seen["runner_up"] == {
        "index": runner_up,
        "coordinate": frequencies[runner_up],
        "value": share[runner_up],
    }
    recorded = [
        json.loads(line) for line in events.read_text().splitlines() if line
    ]
    (evaluated,) = [
        item
        for item in recorded
        if item.get("kind") == "quantity_expression_evaluated"
    ]
    assert evaluated["payload"]["extremum_observations"] == [seen]
