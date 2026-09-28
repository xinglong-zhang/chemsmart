"""Replay approved routes at a tree: does its host still load, admit and parse them?

    HOME=<fence> PYTHONPATH=<tree> python replay_composition.py OUT.jsonl \\
        --list BUNDLES.txt [--label NAME]

``BUNDLES.txt`` holds one archived ``bundle.json`` path per line (``#``
comments allowed). Provider-free and engine-free: nothing here compiles an
argv, starts a program or reads a key. Every level calls the host functions
of the tree ``chemsmart`` resolves to (printed first, with a digest of its
Python sources), so the same list run on the producing tree and on the pin
separates what the pin changed from what the harness does.

Per bundle (c5 EPISODE.md, pre-registered levels):

R1 load -- ``live_session.load_workflow_execution_approval_bundle``: the
   executor's own loader (plan, toolchain, frozen producer rules, digests).
R2 edges -- on the approved plan parsed alone, the tree's
   ``producer_edge_selection_rule`` for every data edge beside the rule the
   frozen approval recorded, and ``admitted_producer_edge_rules`` over the
   frozen data edges (raises when the set is no longer admitted).
R3 analysis -- the toolchain parsed alone; operations it uses that the
   tree's vocabulary no longer holds; and each calculation node's
   (engine, stage) pair against the tree's execution pairs.

Each level is parsed independently of the others, so a refusal at R1 (a
digest or schema change anywhere in the bundle) does not hide whether the
route's composition itself is still admitted.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import traceback
from pathlib import Path


def _tree() -> dict:
    import chemsmart

    package = Path(chemsmart.__file__).resolve().parent
    digest = hashlib.sha256()
    for path in sorted(package.rglob("*.py")):
        digest.update(str(path.relative_to(package)).encode() + b"\0")
        digest.update(path.read_bytes())
    return {
        "chemsmart_file": str(Path(chemsmart.__file__).resolve()),
        "package_py_sha256": digest.hexdigest(),
    }


def _short(exc: BaseException) -> str:
    return f"{type(exc).__name__}: {str(exc)[:600]}"


def replay_one(path: Path) -> dict:
    from chemsmart.agent import live_session
    from chemsmart.agent.capabilities import load_program_capabilities
    from chemsmart.agent.execution import (
        admitted_producer_edge_rules,
        producer_edge_selection_rule,
    )
    from chemsmart.analysis.quantity_expressions import _OPERATIONS

    row: dict = {"bundle": str(path)}
    try:
        raw_document = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        row["error"] = _short(exc)
        return row
    raw = raw_document.get("workflow_execution_approval_bundle", raw_document)
    frozen_raw = raw.get("frozen_workflow_approval") or {}
    row["approval_id"] = frozen_raw.get("approval_id", "")

    # R1: the executor's own loader over the whole bundle.
    try:
        live_session.load_workflow_execution_approval_bundle(path)
        row["R1"] = "loaded"
    except Exception as exc:  # the refusal is the observation
        row["R1"] = "refused"
        row["R1_reason"] = _short(exc)

    # R2: the approved plan alone, then the tree's producer-edge rules.
    frozen_rules = {
        (
            item.get("source_node_id"),
            item.get("target_node_id"),
            item.get("artifact_class"),
        ): item.get("selection_rule", "")
        for item in frozen_raw.get("producer_edge_rules") or []
    }
    try:
        plan = live_session._parse_scientific_workflow_plan(
            raw["approved_scientific_plan"]
        )
        row["R2_plan"] = "parsed"
        row["R2_plan_sha256_equal"] = plan.plan_sha256 == raw[
            "approved_scientific_plan"
        ].get("plan_sha256")
        edges = []
        frozen_edges = []
        for edge in plan.edges:
            if edge.edge_kind != "data":
                continue
            key = (
                edge.source_node_id,
                edge.target_node_id,
                edge.artifact_class,
            )
            now = producer_edge_selection_rule(plan, edge)
            then = frozen_rules.get(key, "")
            edges.append(
                {
                    "edge": list(key),
                    "frozen": then,
                    "tree": now,
                    "equal": now == then,
                }
            )
            if then:
                frozen_edges.append(edge)
        row["R2_edges"] = edges
        row["R2_rules_equal"] = all(
            item["equal"] for item in edges if item["frozen"]
        )
        if frozen_edges:
            try:
                admitted_producer_edge_rules(
                    plan, tuple(frozen_edges), organ="c5 replay"
                )
                row["R2_admitted"] = "admitted"
            except Exception as exc:
                row["R2_admitted"] = "refused"
                row["R2_admitted_reason"] = _short(exc)
        else:
            row["R2_admitted"] = "no-frozen-edge"
        registry = load_program_capabilities()
        pairs = []
        for node in plan.nodes:
            capability = registry.get(node.program)
            executable = bool(
                capability is not None
                and (node.engine, node.stage)
                in capability.execution_engine_job_pairs
            )
            pairs.append([node.program, node.engine, node.stage, executable])
        row["R3_node_pairs"] = pairs
    except Exception as exc:
        row["R2_plan"] = "refused"
        row["R2_plan_reason"] = _short(exc)

    # R3: the toolchain alone, and the operations it names.
    toolchain_raw = raw.get("scientific_toolchain_plan")
    used = sorted(
        {
            expression.get("operation", "")
            for node in (toolchain_raw or {}).get("analysis_nodes") or []
            for expression in node.get("expression_nodes") or []
        }
    )
    row["R3_operations"] = used
    row["R3_missing_operations"] = [op for op in used if op not in _OPERATIONS]
    if toolchain_raw is None:
        row["R3"] = "no-toolchain"
    else:
        try:
            live_session._parse_scientific_toolchain_plan(toolchain_raw)
            row["R3"] = "parsed"
        except Exception as exc:
            row["R3"] = "refused"
            row["R3_reason"] = _short(exc)
    return row


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("out", type=Path)
    parser.add_argument("--list", required=True, type=Path)
    parser.add_argument("--label", default="")
    args = parser.parse_args(argv)
    tree = _tree()
    print(json.dumps({"label": args.label, **tree}))
    paths = [
        Path(line.strip())
        for line in args.list.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.startswith("#")
    ]
    with args.out.open("w", encoding="utf-8") as handle:
        for path in paths:
            try:
                row = replay_one(path)
            except Exception:
                row = {"bundle": str(path), "error": traceback.format_exc()}
            row.update({"label": args.label, **tree})
            handle.write(json.dumps(row, sort_keys=True) + "\n")
            print(
                f"{row.get('R1', '-'):8s} plan={row.get('R2_plan', '-'):8s} "
                f"rules_equal={row.get('R2_rules_equal', '-')!s:5s} "
                f"admitted={row.get('R2_admitted', '-'):14s} "
                f"toolchain={row.get('R3', '-'):12s} "
                f"missing_ops={row.get('R3_missing_operations', [])} {path}"
            )
    return 0


if __name__ == "__main__":
    sys.exit(main())
