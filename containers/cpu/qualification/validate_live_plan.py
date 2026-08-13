#!/usr/bin/env python3
"""Validate a provider-live, execution-free planning transaction.

The validator reads the public Runtime V2 stream produced by ``agent plan``.
It proves the event hash chain, reconstructs the canonical workflow records,
and requires one fully previewed materialization whose node set exactly matches
the latest scientific plan.  It never calls a provider or chemistry engine.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from chemsmart.agent._contracts import canonical_sha256
from chemsmart.agent.runtime.events import (
    OPTIMIZED_GEOMETRY_HANDED_OFF,
    PROGRAM_EXECUTED,
    RUNTIME_TERMINATED,
    SCIENTIFIC_WORKFLOW_MATERIALIZED,
    WORKFLOW_APPROVAL_CONSUMED,
    WORKFLOW_DATA_EDGE_BOUND,
    WORKFLOW_EXECUTION_STARTED,
    WORKFLOW_LAUNCH_RESERVED,
    WORKFLOW_NODE_STATE_CHANGED,
    WORKFLOW_PLANNED,
    RuntimeEvent,
)
from chemsmart.agent.runtime.records import (
    materialized_workflow_from_record,
    scientific_workflow_plan_from_record,
)
from chemsmart.agent.runtime.reducer import replay_events


_EXECUTION_EVENT_KINDS = frozenset(
    {
        PROGRAM_EXECUTED,
        OPTIMIZED_GEOMETRY_HANDED_OFF,
        WORKFLOW_APPROVAL_CONSUMED,
        WORKFLOW_DATA_EDGE_BOUND,
        WORKFLOW_EXECUTION_STARTED,
        WORKFLOW_LAUNCH_RESERVED,
        WORKFLOW_NODE_STATE_CHANGED,
    }
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def validate(workspace: Path) -> dict[str, object]:
    root = workspace.resolve(strict=True)
    streams = tuple(sorted(root.glob(".chemsmart-agent/runs/*/events.jsonl")))
    require(
        len(streams) == 1,
        "provider-live qualification requires exactly one Runtime event stream",
    )
    stream = streams[0]
    events = []
    with stream.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            require(
                line.endswith("\n"), f"unterminated Runtime line {line_number}"
            )
            value = json.loads(line)
            require(
                isinstance(value, dict),
                f"Runtime line {line_number} is not an object",
            )
            events.append(RuntimeEvent.from_dict(value))
    require(events, "provider-live Runtime stream is empty")
    replay_events(events, verify_hashes=True)
    require(
        len({event.session_id for event in events}) == 1,
        "provider-live Runtime stream contains multiple sessions",
    )
    require(
        not any(event.kind in _EXECUTION_EVENT_KINDS for event in events),
        "provider-live qualification observed an execution event",
    )

    planned = [event for event in events if event.kind == WORKFLOW_PLANNED]
    require(
        planned, "provider-live qualification has no scientific workflow plan"
    )
    plan_records = [
        event.payload.get("scientific_plan_record")
        for event in planned
        if event.payload.get("scientific_plan_record")
    ]
    require(
        plan_records, "scientific workflow event lacks its canonical record"
    )
    plan = scientific_workflow_plan_from_record(plan_records[-1])
    require(
        len(plan.nodes) == 2
        and {node.stage for node in plan.nodes} == {"sp", "opt"},
        "provider-live qualification requires exactly one SP and one optimization",
    )
    require(
        all(
            node.requested_program == node.program == "pyscf"
            and node.engine == "cpu"
            for node in plan.nodes
        ),
        "provider-live qualification requires PySCF CPU nodes",
    )
    require(
        not plan.edges
        and all(
            node.support_state == "resolvable" and not node.unresolved_fields
            for node in plan.nodes
        ),
        "provider-live SP and optimization must be independently resolvable",
    )

    observed = []
    for event in events:
        if event.kind != SCIENTIFIC_WORKFLOW_MATERIALIZED:
            continue
        record = event.payload["record"]
        require(
            event.payload["receipt_sha256"] == canonical_sha256(record),
            "materialized workflow event receipt differs from its record",
        )
        workflow = materialized_workflow_from_record(record)
        if workflow.plan_sha256 == plan.plan_sha256:
            observed.append(workflow)
    require(
        observed,
        "provider-live qualification has no materialization for its latest plan",
    )
    workflow = observed[-1]
    require(
        workflow.status == "previewed"
        and not workflow.unresolved_node_ids
        and all(node.state == "previewed" for node in workflow.nodes),
        "latest provider-live scientific workflow is not fully previewed",
    )
    planned_node_ids = tuple(sorted(node.node_id for node in plan.nodes))
    materialized_node_ids = tuple(node.node_id for node in workflow.nodes)
    require(
        materialized_node_ids == planned_node_ids,
        "fully previewed materialization differs from the planned node set",
    )
    require(
        len({node.input_artifact_sha256 for node in workflow.nodes}) == 1,
        "provider-live SP and optimization must consume one initial geometry",
    )

    terminal = [event for event in events if event.kind == RUNTIME_TERMINATED]
    require(terminal, "provider-live Runtime did not terminate durably")
    terminal_state = str(terminal[-1].payload.get("terminal_state") or "")
    require(
        terminal_state in {"complete", "planned"},
        f"provider-live Runtime terminated as {terminal_state or 'unknown'}",
    )
    return {
        "schema_version": "chemsmart.cpu-live-plan-qualification.v1",
        "status": "qualified",
        "event_count": len(events),
        "session_id": events[0].session_id,
        "terminal_state": terminal_state,
        "plan_sha256": plan.plan_sha256,
        "planned_node_ids": planned_node_ids,
        "materialized_workflow_sha256": workflow.materialized_sha256,
        "materialized_status": workflow.status,
        "unresolved_node_ids": workflow.unresolved_node_ids,
        "execution_event_count": 0,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = validate(args.workspace)
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        args.output.write_text(rendered, encoding="utf-8")
    print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
