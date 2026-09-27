"""Probe: can a typed plan pick a vibrational mode by which atoms move,
before the Hessian that holds the modes exists?

The question comes from R9 xtb g3 (cephalexin), whose last wave ran a
GFN2-xTB Hessian that the session could not assign: "mode selection is a
post-execution data read; with no further agent cycle after this wave, the
adaptive assignment cannot be planned as a claim chain before execution"
(the graph's ``negative.an_analysis_chain_cannot_plan_a_mode_it_has_not_read``).
On the earlier Hessian the session had read the participation table and
chosen by eye: mode 99 (1866.33 cm-1) is the beta-lactam C=O, carried by
atoms 10 and 1 (0.649 + 0.320 = 0.969).

For each named result the probe drives the host's public tool surface,
``CommandCompiledToolHostV1.dispatch``, the way a session does, and records
every reply verbatim:

1. ``extract_result_quantities``: vibrational frequencies, the per-atom
   participation table (modes x atoms) and the atom symbols.
2. Blind constructions in today's vocabulary -- none may name a mode
   index, because a plan written before the Hessian cannot know it:
   ``coordinate_at_maximum`` over the table against the frequencies; a
   ``multiply`` by an atom mask; a column select in ``ref`` (as a null
   index and as an ``axis`` field); and ``sum``/``max`` over the table.
3. The target, computed from the extraction receipt's own values (numpy,
   outside the typed layer, as the witness's expected answer): for each
   atom set, the mode whose summed participation is largest, its
   frequency, its share, and the runner-up.
4. The static-index route a later wake cycle takes once the index has been
   read: ``ref`` [k, a] per atom, ``add``, ``ref`` frequencies [k], through
   ``evaluate_quantity_expression`` -- the typed receipt the goal needed.

It grades nothing and changes nothing: it reads each result where it lies.

    python probe_mode_selection.py OUT.jsonl LABEL:PROGRAM:PATH:SET[;SET...] ...

where SET is ``name=a+b`` with 0-based atom indices, for example
``route2:xtb:/path/run.out:lactam=10+1;amide=14+2``.
"""

from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

import numpy as np


def _host(workdir, artifacts):
    from chemsmart.agent._contracts import TrustedArtifactRefV1, file_sha256
    from chemsmart.agent.runtime.event_store import RuntimeEventStore
    from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1

    host = CommandCompiledToolHostV1(
        event_store=RuntimeEventStore(
            Path(workdir) / "events.jsonl", session_id="probe-mode-selection"
        ),
        task_spec_sha256s=("a" * 64,),
        approved_workspace=Path(workdir),
    )
    for artifact_id, program, path in artifacts:
        resolved = Path(path).resolve()
        host.artifacts[artifact_id] = TrustedArtifactRefV1(
            artifact_id=artifact_id,
            kind=f"{program}_output",
            sha256=file_sha256(resolved),
            size_bytes=resolved.stat().st_size,
            path=str(resolved),
            cli_value=str(resolved),
        )
    return host


def _call(host, turn, tool, arguments):
    try:
        reply = host.dispatch(
            turn_id=turn, tool_name=tool, arguments=arguments
        )
    except Exception as exc:  # noqa: BLE001 -- a raise is a reply too
        return {
            "status": "raised",
            "error_class": type(exc).__name__,
            "message": str(exc)[:1500],
        }
    if isinstance(reply, dict) and reply.get("status") != "ok":
        return {
            "status": reply.get("status"),
            "error_class": reply.get("error_class"),
            "message": str(reply.get("message") or reply.get("error") or "")[
                :1500
            ],
        }
    return reply


def _expression(host, turn, name, receipt, nodes, outputs, inputs=None):
    inputs = inputs or [
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
    ]
    return _call(
        host,
        turn,
        "evaluate_quantity_expression",
        {
            "expression_id": name,
            "inputs": inputs,
            "nodes": nodes,
            "output_node_ids": outputs,
        },
    )


def _summary(reply):
    if reply.get("status") != "ok":
        return reply
    result = reply.get("result") or {}
    values = {
        item.get("quantity_id"): item.get("value")
        for item in result.get("outputs") or ()
    }
    summary = {
        "status": "ok",
        "receipt_sha256": result.get("receipt_sha256"),
        "outputs": values,
    }
    # What the host says beside the receipt (E1: the runner-up of every
    # coordinate_at_maximum/minimum), recorded as the reply carried it.
    if reply.get("observations"):
        summary["observations"] = reply["observations"]
    return summary


def probe(label, program, path, sets, emit):
    workdir = tempfile.mkdtemp(prefix=f"probe-{label}-")
    artifact_id = f"result.{label}"
    host = _host(workdir, [(artifact_id, program, path)])
    extraction = _call(
        host,
        "t-extract",
        "extract_result_quantities",
        {
            "artifact_id": artifact_id,
            "program": program,
            "selectors": [
                {
                    "quantity_id": "frequencies",
                    "selector": "vibrational_frequencies",
                },
                {
                    "quantity_id": "participation",
                    "selector": "vibrational_mode_atom_participation",
                },
                {"quantity_id": "symbols", "selector": "symbols"},
            ],
        },
    )
    if extraction.get("status") != "ok":
        emit({"label": label, "step": "extract", "reply": extraction})
        return
    result = extraction["result"]
    receipt = result["receipt_sha256"]
    values = {
        item["quantity_id"]: item["value"] for item in result["quantities"]
    }
    table = np.asarray(values["participation"], dtype=float)
    frequencies = np.asarray(values["frequencies"], dtype=float)
    symbols = list(values["symbols"])
    emit(
        {
            "label": label,
            "step": "extract",
            "path": path,
            "artifact_sha256": result.get("artifact_sha256"),
            "receipt_sha256": receipt,
            "modes": int(frequencies.size),
            "atoms": len(symbols),
            "table_shape": list(table.shape),
            "lowest_frequencies": [
                round(float(x), 2) for x in frequencies[:3]
            ],
        }
    )

    # 2. Blind constructions: no mode index anywhere.
    first_set = next(iter(sets.values()))
    mask = [1.0 if i in first_set else 0.0 for i in range(len(symbols))]
    blind = {
        "B1_coordinate_at_maximum_over_table": [
            {
                "node_id": "pick",
                "operation": "coordinate_at_maximum",
                "input_ids": ["part", "freq"],
            },
        ],
        "B2_multiply_by_atom_mask": [
            {
                "node_id": "mask",
                "operation": "literal",
                "literal_value": mask,
                "literal_unit": "1",
            },
            {
                "node_id": "masked",
                "operation": "multiply",
                "input_ids": ["part", "mask"],
            },
            {
                "node_id": "pick",
                "operation": "coordinate_at_maximum",
                "input_ids": ["masked", "freq"],
            },
        ],
        "B3a_ref_column_as_null_index": [
            {
                "node_id": "col",
                "operation": "ref",
                "reference": "part",
                "indices": [None, first_set[0]],
            },
            {
                "node_id": "pick",
                "operation": "coordinate_at_maximum",
                "input_ids": ["col", "freq"],
            },
        ],
        "B3b_ref_column_as_axis_field": [
            {
                "node_id": "col",
                "operation": "ref",
                "reference": "part",
                "indices": [first_set[0]],
                "axis": 1,
            },
            {
                "node_id": "pick",
                "operation": "coordinate_at_maximum",
                "input_ids": ["col", "freq"],
            },
        ],
        "B4a_sum_over_table": [
            {"node_id": "pick", "operation": "sum", "input_ids": ["part"]},
        ],
        "B4b_max_over_table": [
            {"node_id": "pick", "operation": "max", "input_ids": ["part"]},
        ],
    }
    for name, nodes in blind.items():
        reply = _expression(
            host,
            f"t-{name}",
            name.lower().replace("_", "-"),
            receipt,
            nodes,
            ["pick"],
        )
        emit(
            {
                "label": label,
                "step": "blind",
                "construction": name,
                "reply": _summary(reply),
            }
        )

    # 3. The target, and 4. the static-index route a later cycle takes.
    for set_name, atoms in sets.items():
        # W: the candidate plan itself, written before the Hessian exists --
        # one column select per atom, their sum, the frequency at its
        # maximum and the maximum share. Today's tree refuses it at the
        # column select; the witness is green when it returns the target.
        plan = [
            {
                "node_id": f"col{i}",
                "operation": "ref",
                "reference": "part",
                "indices": [None, a],
            }
            for i, a in enumerate(atoms)
        ]
        running = "col0"
        for i in range(1, len(atoms)):
            plan.append(
                {
                    "node_id": f"sum{i}",
                    "operation": "add",
                    "input_ids": [running, f"col{i}"],
                }
            )
            running = f"sum{i}"
        plan += [
            {
                "node_id": "nu",
                "operation": "coordinate_at_maximum",
                "input_ids": [running, "freq"],
            },
            {"node_id": "top", "operation": "max", "input_ids": [running]},
        ]
        reply = _expression(
            host,
            f"t-plan-{set_name}",
            f"plan-{set_name}",
            receipt,
            plan,
            ["nu", "top"],
        )
        emit(
            {
                "label": label,
                "step": "candidate_plan",
                "set": set_name,
                "nodes": len(plan),
                "reply": _summary(reply),
            }
        )
        share = table[:, list(atoms)].sum(axis=1)
        order = np.argsort(share)[::-1]
        best, runner = int(order[0]), int(order[1])
        target = {
            "label": label,
            "step": "target",
            "set": set_name,
            "atoms": list(atoms),
            "atom_symbols": [symbols[i] for i in atoms],
            "mode_row": best,
            "mode_number": best + 1,
            "frequency_cm1": float(frequencies[best]),
            "share": float(share[best]),
            "per_atom": [float(table[best, i]) for i in atoms],
            "runner_up_row": runner,
            "runner_up_frequency_cm1": float(frequencies[runner]),
            "runner_up_share": float(share[runner]),
        }
        emit(target)
        nodes = [
            {
                "node_id": f"p{i}",
                "operation": "ref",
                "reference": "part",
                "indices": [best, a],
            }
            for i, a in enumerate(atoms)
        ]
        running = "p0"
        for i in range(1, len(atoms)):
            nodes.append(
                {
                    "node_id": f"s{i}",
                    "operation": "add",
                    "input_ids": [running, f"p{i}"],
                }
            )
            running = f"s{i}"
        nodes.append(
            {"node_id": "share", "operation": "ref", "reference": running}
        )
        nodes.append(
            {
                "node_id": "nu",
                "operation": "ref",
                "reference": "freq",
                "indices": [best],
            }
        )
        reply = _expression(
            host,
            f"t-static-{set_name}",
            f"static-{set_name}",
            receipt,
            nodes,
            ["share", "nu"],
        )
        emit(
            {
                "label": label,
                "step": "static_index_route",
                "set": set_name,
                "mode_row": best,
                "reply": _summary(reply),
            }
        )


def parse(spec):
    label, program, rest = spec.split(":", 2)
    path, sets_text = rest.rsplit(":", 1)
    sets = {}
    for item in sets_text.split(";"):
        name, atoms = item.split("=")
        sets[name] = tuple(int(a) for a in atoms.split("+"))
    return label, program, path, sets


def main(argv):
    import chemsmart

    out = Path(argv[1])
    with out.open("w", encoding="utf-8") as handle:

        def emit(row):
            text = json.dumps(row, sort_keys=True, default=str)
            handle.write(text + "\n")
            print(text)

        emit({"step": "import", "chemsmart": chemsmart.__file__})
        for spec in argv[2:]:
            probe(*parse(spec), emit)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
