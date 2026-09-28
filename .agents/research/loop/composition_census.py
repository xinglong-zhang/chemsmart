"""Composition census: every approved cycle's route, read from its bundle.

    python composition_census.py OUT.jsonl STRATUM=ROOT [STRATUM=ROOT ...]
    python composition_census.py --report OUT.jsonl [OUT2.jsonl ...]

Imports nothing from ``chemsmart``, so it runs on any Python, beside the
records it reads (the ax41 mirror on the Mac, CUHK records in place inside
a slot job). The unit is an **approved cycle**: one
``<agent-dir>/replays/<id>/bundle.json`` (the executor's own input),
de-duplicated by real path. Everything a row says is read from that bundle
or from a host stream or goal ledger in the same agent directory; nothing
is inferred from task text or narration.

Per cycle it records the calculation nodes (program, stage, engine, state,
heavy atoms, how the geometry arrived), the in-approval data edges with the
producer rule the frozen approval names, the geometry lineage of every node
(which operations built it, and whether it names an earlier result), the
analysis chain (kinds, operations, constants, predicates, and which
calculations each expression reaches), what executed (the last node states
the host recorded for this approval id, the handoffs, the analysis
settlements) and, for goals, the cycle's recorded workflow state and the
goal's settlement word.

Classes (pre-registered in c5's EPISODE.md): S one node; P several, no
data edge; M an in-approval data edge within one program; X one across
programs; L a node lifted from an earlier result (L-intra / L-cross); C a
side-by-side comparison (>= 2 programs, no crossing edge, one expression
reaching both); A an expression reaching >= 2 calculations.
"""

from __future__ import annotations

import hashlib
import json
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path

#: Job words that name a task rather than a stage, as the live CLI spells
#: them. A route whose argv carries one reaches task code.
TASK_JOB_WORDS = frozenset(
    {
        "pka",
        "dias",
        "nci",
        "resp",
        "wbi",
        "crest",
        "qrc",
        "traj",
        "neb",
        "link",
        "userjob",
        "com",
        "inp",
    }
)

#: Operations whose name is a task (the audit's hard cases).
TASK_NAMED_OPERATIONS = frozenset(
    {
        "gibbs_to_pka",
        "gibbs_to_redox_potential",
        "exponential_cbs_limit",
        "scf_exponential_cbs_limit",
        "scf_inverse_power_cbs_limit",
        "correlation_inverse_power_cbs_limit",
        "transition_state_crossover_temperature",
    }
)

#: Lineage kinds that start a node from an earlier calculation's result
#: rather than from a supplied or edited geometry.
RESULT_LINEAGE_HINTS = ("result_artifact_id", "result_sha256", "program")

MAX_STREAM_BYTES = 400 * 1024 * 1024


def _load(path: Path):
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError, UnicodeDecodeError):
        return None


def _bundles(root: Path):
    """Every bundle.json whose parent's parent is ``replays``, once."""

    seen_dirs: set[str] = set()
    seen_files: set[str] = set()
    for folder, dirs, files in os.walk(root, followlinks=True):
        real = os.path.realpath(folder)
        if real in seen_dirs:
            dirs[:] = []
            continue
        seen_dirs.add(real)
        dirs[:] = [
            d
            for d in dirs
            if d not in {"scratch", "__pycache__", ".git"}
            and not d.startswith("sealed")
        ]
        if "bundle.json" in files and Path(folder).parent.name == "replays":
            path = Path(folder) / "bundle.json"
            key = os.path.realpath(path)
            if key not in seen_files:
                seen_files.add(key)
                yield path


def _heavy(atom_order) -> int:
    return sum(1 for symbol in atom_order or () if str(symbol) != "H")


def _lineage(identity: dict) -> list[dict]:
    rows = []
    for entry in identity.get("geometry_lineage") or []:
        if not isinstance(entry, dict):
            continue
        rows.append(
            {
                "kind": entry.get("kind") or entry.get("operation") or "?",
                "operation": entry.get("operation", ""),
                "program": entry.get("program", ""),
                "names_result": any(
                    key in entry
                    for key in ("result_artifact_id", "result_sha256")
                ),
            }
        )
    for origin in (
        "composition",
        "derivation",
        "database_extraction",
        "geometry_edit",
    ):
        if identity.get(origin):
            rows.append(
                {
                    "kind": origin,
                    "operation": "",
                    "program": "",
                    "names_result": False,
                }
            )
    return rows


def _argv_job(argv, program: str) -> tuple[str, list[str]]:
    tokens = [str(item) for item in argv or ()]
    task_words = [t for t in tokens if t in TASK_JOB_WORDS]
    job = tokens[-1] if tokens else ""
    return job, task_words


def _calc_closure(analysis: dict, calc_ids: set[str]) -> dict[str, set[str]]:
    """For each analysis node, the calculation nodes it reaches."""

    by_id = {node["node_id"]: node for node in analysis}
    memo: dict[str, set[str]] = {}

    def reach(node_id: str, stack: tuple[str, ...] = ()) -> set[str]:
        if node_id in calc_ids:
            return {node_id}
        if node_id in memo:
            return memo[node_id]
        if node_id in stack or node_id not in by_id:
            return set()
        node = by_id[node_id]
        found: set[str] = set()
        if node.get("artifact_id"):
            # A root that reads a host-registered result of an earlier
            # workflow, not a calculation of this plan: a leaf of its own.
            found.add("reg:" + str(node["artifact_id"]))
        for item in node.get("inputs") or []:
            if (
                isinstance(item, dict)
                and item.get("source_kind") == "registered_result"
            ):
                found.add("reg:" + str(item.get("artifact_id", "")))
        refs = list(node.get("dependencies") or [])
        refs += [
            item.get("producer_node_id", "")
            for item in node.get("inputs") or []
            if isinstance(item, dict)
        ]
        for ref in refs:
            if ref:
                found |= reach(ref, stack + (node_id,))
        memo[node_id] = found
        return found

    return {node_id: reach(node_id) for node_id in by_id}


def _stream_index(agent_dir: Path) -> dict:
    """Host stream facts per approval id and per toolchain plan digest."""

    states: dict[str, dict[str, str]] = defaultdict(dict)
    workflow_state: dict[str, str] = {}
    handoffs: dict[str, list] = defaultdict(list)
    settled: dict[str, Counter] = defaultdict(Counter)
    reached: dict[str, dict] = {}
    for folder, dirs, files in os.walk(agent_dir, followlinks=False):
        dirs[:] = [d for d in dirs if d not in {"scratch", "__pycache__"}]
        if "events.jsonl" not in files:
            continue
        path = Path(folder) / "events.jsonl"
        try:
            if path.stat().st_size > MAX_STREAM_BYTES:
                continue
            handle = path.open(encoding="utf-8", errors="replace")
        except OSError:
            continue
        stream_approvals: set[str] = set()
        pending_handoffs: list = []
        with handle:
            for line in handle:
                if (
                    "workflow_node_state_changed" not in line
                    and "handed_off" not in line
                    and "workflow_analysis_node_settled" not in line
                    and "reached_geometry_bound" not in line
                ):
                    continue
                try:
                    event = json.loads(line)
                except ValueError:
                    continue
                kind = event.get("kind", "")
                payload = event.get("payload") or {}
                if kind == "workflow_node_state_changed":
                    record = payload.get("record") or {}
                    approval = record.get("approval_id", "")
                    if not approval:
                        continue
                    stream_approvals.add(approval)
                    for node in record.get("nodes") or []:
                        states[approval][node.get("node_id", "")] = node.get(
                            "state", ""
                        )
                    workflow_state[approval] = payload.get(
                        "workflow_state", record.get("state", "")
                    )
                elif "handed_off" in kind:
                    pending_handoffs.append(
                        (
                            payload.get("producer_node_id", ""),
                            payload.get("consumer_node_id", ""),
                            kind,
                            payload.get("status", ""),
                        )
                    )
                elif kind == "workflow_analysis_node_settled":
                    settled[payload.get("toolchain_plan_sha256", "")][
                        payload.get("state", "")
                    ] += 1
                elif kind == "reached_geometry_bound":
                    # bind_reached_geometry: a completed result's structure
                    # lifted into a new geometry artifact a later workflow
                    # starts from. The bundle's node review names only the
                    # artifact digest; this record names where it came from.
                    record = payload.get("record") or {}
                    digest = record.get("reached_artifact_sha256", "")
                    if digest:
                        reached[digest] = {
                            "program": record.get("program", ""),
                            "node": record.get("recorded_node_id", ""),
                            "state": record.get("recorded_terminal_state", ""),
                        }
        for approval in stream_approvals:
            handoffs[approval].extend(pending_handoffs)
    return {
        "states": states,
        "workflow_state": workflow_state,
        "handoffs": handoffs,
        "settled": settled,
        "reached": reached,
    }


def _ledgers(agent_dir: Path) -> dict[str, dict]:
    goals = {}
    base = agent_dir / "goals"
    if not base.is_dir():
        return goals
    for ledger in sorted(base.glob("*/ledger.jsonl")):
        cycles: dict[int, str] = {}
        settlement = ""
        goal_sha256 = ""
        try:
            lines = ledger.read_text(encoding="utf-8").splitlines()
        except OSError:
            continue
        for line in lines:
            try:
                row = json.loads(line)
            except ValueError:
                continue
            payload = row.get("payload") or {}
            if row.get("kind") == "run_recorded":
                try:
                    cycles[int(payload.get("cycle"))] = payload.get(
                        "workflow_state", ""
                    )
                except (TypeError, ValueError):
                    pass
            elif row.get("kind") == "goal_settled":
                settlement = payload.get("state", "")
            elif row.get("kind") == "goal_created":
                goal_sha256 = payload.get("goal_sha256", "")
        goals[ledger.parent.name] = {
            "cycles": cycles,
            "settlement": settlement,
            "goal_sha256": goal_sha256,
        }
    return goals


def census_row(
    stratum: str, root: Path, path: Path, index: dict, goals: dict
) -> dict | None:
    document = _load(path)
    if not isinstance(document, dict):
        return {"stratum": stratum, "bundle": str(path), "error": "unreadable"}
    # Content identity: the ax41 mirror holds byte-identical copies of
    # earlier goals' workspaces (R11 master, truth-4), so a route is counted
    # once per distinct bundle content, never once per path.
    content_sha256 = hashlib.sha256(path.read_bytes()).hexdigest()
    bundle = document.get("workflow_execution_approval_bundle", document)
    resolution = bundle.get("resolution") or {}
    frozen = bundle.get("frozen_workflow_approval") or {}
    plan = bundle.get("approved_scientific_plan") or {}
    approval = (
        resolution.get("approval_id")
        or frozen.get("approval_id")
        or (bundle.get("workflow_approval") or {}).get("approval_id", "")
    )
    rules = {
        (
            item.get("source_node_id"),
            item.get("target_node_id"),
            item.get("artifact_class"),
        ): item.get("selection_rule", "")
        for item in frozen.get("producer_edge_rules") or []
    }
    reviews = {
        item.get("node_id"): item for item in bundle.get("node_reviews") or []
    }
    nodes = []
    for node in plan.get("nodes") or []:
        node_id = node.get("node_id", "")
        review = reviews.get(node_id) or {}
        identity = review.get("molecular_identity") or {}
        coordinate = identity.get("coordinate_identity") or {}
        job, task_words = _argv_job(
            review.get("real_execution_argv"), node.get("program", "")
        )
        lineage = _lineage(identity)
        lifted = index["reached"].get(
            coordinate.get("geometry_artifact_sha256", "")
        )
        if lifted:
            lineage.append(
                {
                    "kind": "reached_geometry",
                    "operation": "",
                    "program": lifted["program"],
                    "names_result": True,
                }
            )
        nodes.append(
            {
                "id": node_id,
                "program": node.get("program", ""),
                "stage": node.get("stage", ""),
                "engine": node.get("engine", ""),
                "node_kind": node.get("node_kind", "program_call"),
                "charge": node.get("charge"),
                "multiplicity": node.get("multiplicity"),
                "formula": identity.get("formula", ""),
                "atoms": identity.get("atom_count"),
                "heavy": _heavy(identity.get("atom_order")),
                "coord_kind": coordinate.get("kind", ""),
                "identity_status": identity.get(
                    "identity_evidence_status", ""
                ),
                "lineage": lineage,
                "argv_job": job,
                "argv_task_words": task_words,
            }
        )
    by_id = {node["id"]: node for node in nodes}
    non_executable = set(bundle.get("non_executable_node_ids") or [])
    edges = []
    for edge in plan.get("edges") or []:
        if edge.get("edge_kind") != "data":
            continue
        src = by_id.get(edge.get("source_node_id"), {})
        tgt = by_id.get(edge.get("target_node_id"), {})
        rule = rules.get(
            (
                edge.get("source_node_id"),
                edge.get("target_node_id"),
                edge.get("artifact_class"),
            ),
            "",
        )
        if not rule and (
            edge.get("target_node_id") in non_executable
            or edge.get("source_node_id") in non_executable
        ):
            # An edge into a stage the review retained as non-executable
            # intent: displayed and approved as intent, never frozen.
            rule = "NON-EXECUTABLE"
        edges.append(
            {
                "src": edge.get("source_node_id"),
                "tgt": edge.get("target_node_id"),
                "class": edge.get("artifact_class", ""),
                "rule": rule,
                "src_program": src.get("program", ""),
                "src_stage": src.get("stage", ""),
                "tgt_program": tgt.get("program", ""),
                "tgt_stage": tgt.get("stage", ""),
                "cross": bool(
                    src and tgt and src.get("program") != tgt.get("program")
                ),
            }
        )
    toolchain = bundle.get("scientific_toolchain_plan") or {}
    analysis = [
        item
        for item in toolchain.get("analysis_nodes") or []
        if isinstance(item, dict)
    ]
    calc_ids = set(by_id)
    closure = _calc_closure(analysis, calc_ids)
    ops: Counter = Counter()
    constants: list[str] = []
    predicates: Counter = Counter()
    kinds: Counter = Counter()
    expressions = []
    registered_reads = sorted(
        {
            str(item["artifact_id"])
            for item in analysis
            if item.get("artifact_id")
        }
        | {
            str(entry.get("artifact_id", ""))
            for item in analysis
            for entry in item.get("inputs") or []
            if isinstance(entry, dict)
            and entry.get("source_kind") == "registered_result"
        }
    )
    for item in analysis:
        kind = item.get("analysis_kind") or item.get("kind") or ""
        kinds[kind] += 1
        for expression in item.get("expression_nodes") or []:
            ops[expression.get("operation", "")] += 1
            if expression.get("constant_name"):
                constants.append(expression["constant_name"])
        for rule in item.get("validation_rules") or []:
            if isinstance(rule, dict):
                predicates[rule.get("predicate", "")] += 1
        if kind == "quantity_expression":
            reached = sorted(closure.get(item.get("node_id", ""), set()))
            expressions.append(
                {
                    "node": item.get("node_id", ""),
                    "calc_nodes": [n for n in reached if n in by_id],
                    "registered": [n for n in reached if n.startswith("reg:")],
                    "programs": sorted(
                        {by_id[n]["program"] for n in reached if n in by_id}
                    ),
                    "task_ops": sorted(
                        {
                            e.get("operation", "")
                            for e in item.get("expression_nodes") or []
                            if e.get("operation") in TASK_NAMED_OPERATIONS
                        }
                    ),
                }
            )
    classes = []
    programs = sorted({node["program"] for node in nodes})
    if not edges:
        classes.append("S" if len(nodes) == 1 else ("P" if nodes else "none"))
    elif any(edge["cross"] for edge in edges):
        classes.append("X")
    else:
        classes.append("M")
    lifts = []
    for node in nodes:
        for entry in node["lineage"]:
            if entry["names_result"] or entry["program"]:
                lifts.append(
                    {
                        "node": node["id"],
                        "kind": entry["kind"],
                        "from_program": entry["program"],
                        "to_program": node["program"],
                        "to_stage": node["stage"],
                    }
                )
    if lifts:
        cross_lift = any(
            item["from_program"] and item["from_program"] != item["to_program"]
            for item in lifts
        )
        classes.append("L-cross" if cross_lift else "L-intra")
    crossing = any(edge["cross"] for edge in edges)
    if (
        len(programs) >= 2
        and not crossing
        and any(len(e["programs"]) >= 2 for e in expressions)
    ):
        classes.append("C")
    if any(len(e["calc_nodes"]) >= 2 for e in expressions):
        classes.append("A")
    if any(e["registered"] for e in expressions):
        # An expression over results of an earlier workflow: a cross-cycle
        # analysis composition (analysis-only revisions, reused results).
        classes.append("A-reg")
    shape_edges = sorted(
        {
            f"{e['src_program']}:{e['src_stage']}->{e['tgt_program']}:{e['tgt_stage']}[{e['rule'] or 'NO-RULE'}]"
            for e in edges
        }
    )
    shape_lifts = sorted(
        {
            f"lift:{i['kind']}:{i['from_program'] or '?'}->{i['to_program']}:{i['to_stage']}"
            for i in lifts
        }
    )
    shape = " ; ".join(shape_edges + shape_lifts) or (
        "nodes:"
        + ",".join(sorted({f"{n['program']}:{n['stage']}" for n in nodes}))
    )
    states = index["states"].get(approval, {})
    settled = index["settled"].get(toolchain.get("plan_sha256", ""), Counter())
    goal_id, cycle = "", None
    if approval.startswith("goal-") and "-cycle-" in approval:
        head, _, tail = approval.rpartition("-cycle-")
        goal_id = head[len("goal-") :]
        try:
            cycle = int(tail)
        except ValueError:
            cycle = None
    goal = goals.get(goal_id, {})
    return {
        "stratum": stratum,
        "bundle": os.path.relpath(path, root),
        "root": str(root),
        "approval_id": approval,
        "bundle_file_sha256": content_sha256,
        "goal_sha256": goal.get("goal_sha256", ""),
        "actor": resolution.get("actor", ""),
        "decision": resolution.get("decision", ""),
        "goal": goal_id,
        "cycle": cycle,
        "workflow_id": plan.get("workflow_id", ""),
        "nodes": nodes,
        "non_executable": list(bundle.get("non_executable_node_ids") or []),
        "edges": edges,
        "lifts": lifts,
        "analysis": {
            "registered_reads": registered_reads,
            "kinds": dict(kinds),
            "ops": dict(ops),
            "constants": constants,
            "predicates": dict(predicates),
            "expressions": expressions,
        },
        "classes": classes,
        "programs": programs,
        "shape": shape,
        "max_heavy": max((n["heavy"] for n in nodes), default=0),
        "exec": {
            "stream_found": bool(states),
            "node_states": states,
            "workflow_state": index["workflow_state"].get(approval, ""),
            "handoffs": index["handoffs"].get(approval, []),
            "analysis_settled": dict(settled),
        },
        "ledger": {
            "run_workflow_state": (
                (goal.get("cycles") or {}).get(cycle, "") if cycle else ""
            ),
            "settlement": goal.get("settlement", ""),
        },
    }


def run_census(out: Path, pairs: list[str]) -> None:
    rows = 0
    with out.open("w", encoding="utf-8") as handle:
        for pair in pairs:
            stratum, _, root_text = pair.partition("=")
            root = Path(root_text)
            cache: dict[str, tuple[dict, dict]] = {}
            for path in _bundles(root):
                agent_dir = path.parent.parent.parent
                key = str(agent_dir)
                if key not in cache:
                    cache[key] = (
                        _stream_index(agent_dir),
                        _ledgers(agent_dir),
                    )
                index, goals = cache[key]
                row = census_row(stratum, root, path, index, goals)
                if row is not None:
                    handle.write(json.dumps(row, sort_keys=True) + "\n")
                    rows += 1
    print(f"wrote {rows} rows to {out}")


VALIDATED = {"validated"}


def _executed(row: dict) -> bool:
    """A calculation validated, or an analysis-only cycle's chain ran."""

    states = row["exec"]["node_states"]
    if states and any(state in VALIDATED for state in states.values()):
        return True
    return not row["nodes"] and bool(
        row["exec"]["analysis_settled"].get("executed")
    )


def _distinct(rows: list[dict]) -> list[dict]:
    """One row per distinct bundle content.

    Among byte-identical copies the kept row is the one whose own agent
    directory recorded its execution (a stream naming the approval), then
    the one beside a goal record, then the first given.
    """

    best: dict[str, dict] = {}
    order: list[str] = []
    for row in rows:
        key = row.get("bundle_file_sha256") or row["root"] + row["bundle"]
        rank = (
            bool(row["exec"]["stream_found"]),
            bool(row["ledger"].get("settlement")),
        )
        if key not in best:
            order.append(key)
            best[key] = row
        elif rank > (
            bool(best[key]["exec"]["stream_found"]),
            bool(best[key]["ledger"].get("settlement")),
        ):
            best[key] = row
    return [best[key] for key in order]


def _goal_key(row: dict) -> tuple:
    return (row.get("goal_sha256") or row["root"], row["goal"])


def report(paths: list[str]) -> None:
    rows = []
    for path in paths:
        for line in Path(path).read_text(encoding="utf-8").splitlines():
            if line.strip():
                rows.append(json.loads(line))
    rows = [row for row in rows if "error" not in row]
    raw = len(rows)
    rows = _distinct(rows)
    strata = sorted({row["stratum"] for row in rows})
    print(
        f"approved cycles: {raw} bundle paths, {len(rows)} distinct by "
        f"bundle content, in strata {strata}"
    )
    for stratum in strata + ["ALL"]:
        chosen = [
            r for r in rows if stratum == "ALL" or r["stratum"] == stratum
        ]
        executed = [r for r in chosen if _executed(r)]
        classes = Counter(c for r in chosen for c in r["classes"])
        classes_exec = Counter(c for r in executed for c in r["classes"])
        goals = {_goal_key(r) for r in chosen if r["goal"]}
        print(
            f"\n[{stratum}] cycles {len(chosen)} (executed: >=1 node validated in a host stream: {len(executed)}); goals {len(goals)}"
        )
        print(
            "  classes (all / executed):",
            {k: (classes[k], classes_exec[k]) for k in sorted(classes)},
        )
    print(
        "\nroute shapes of executed M/X/L cycles (count, max heavy atoms, strata):"
    )
    shapes: dict[str, list] = defaultdict(list)
    for row in rows:
        if _executed(row) and any(
            c in {"M", "X", "L-intra", "L-cross"} for c in row["classes"]
        ):
            shapes[row["shape"]].append(row)
    for shape, members in sorted(
        shapes.items(), key=lambda kv: (-len(kv[1]), kv[0])
    ):
        print(
            f"  {len(members):4d}  heavy<={max(r['max_heavy'] for r in members):3d}  "
            f"{sorted({r['stratum'] for r in members})}  {shape}"
        )
    print(
        "\nin-approval data edges (all cycles / edge executed = consumer validated):"
    )
    edge_counts: Counter = Counter()
    edge_exec: Counter = Counter()
    for row in rows:
        states = row["exec"]["node_states"]
        for edge in row["edges"]:
            key = (
                edge["src_program"],
                edge["tgt_program"],
                edge["class"],
                edge["rule"] or "NO-RULE",
            )
            edge_counts[key] += 1
            if states.get(edge["tgt"]) in VALIDATED:
                edge_exec[key] += 1
    for key in sorted(edge_counts):
        print(f"  {edge_counts[key]:4d} / {edge_exec[key]:4d}  {key}")
    print("\ncross-cycle lifts (kind, from -> to; all / in executed cycles):")
    lift_counts: Counter = Counter()
    lift_exec: Counter = Counter()
    for row in rows:
        for lift in row["lifts"]:
            key = (
                lift["kind"],
                lift["from_program"] or "?",
                f"{lift['to_program']}:{lift['to_stage']}",
            )
            lift_counts[key] += 1
            if _executed(row):
                lift_exec[key] += 1
    for key in sorted(lift_counts):
        print(f"  {lift_counts[key]:4d} / {lift_exec[key]:4d}  {key}")
    print("\ngeometry origins on approved nodes (lineage kinds, all nodes):")
    origins = Counter(
        entry["kind"]
        for row in rows
        for node in row["nodes"]
        for entry in node["lineage"]
    )
    print("  ", dict(origins))
    print("\noperations (all cycles / executed cycles), task-named marked *:")
    op_all: Counter = Counter()
    op_exec: Counter = Counter()
    for row in rows:
        for op, count in row["analysis"]["ops"].items():
            op_all[op] += count
            if _executed(row):
                op_exec[op] += count
    for op in sorted(op_all, key=lambda o: -op_all[o]):
        mark = "*" if op in TASK_NAMED_OPERATIONS else " "
        print(f"  {mark} {op_all[op]:5d} / {op_exec[op]:5d}  {op}")
    print(
        "\ntask-named operations by cycle (stratum, bundle, calc nodes reached, programs):"
    )
    for row in rows:
        for expression in row["analysis"]["expressions"]:
            if expression["task_ops"]:
                print(
                    f"  {row['stratum']:10s} {'EXEC' if _executed(row) else 'plan'} "
                    f"{expression['task_ops']} n={len(expression['calc_nodes'])} {expression['programs']} {row['bundle']}"
                )
    print(
        "\nregistered-result reads (an earlier workflow's result read by "
        "this cycle's analysis; by the result's program prefix):"
    )
    reads = Counter(
        artifact.split("-result-")[0] if "-result-" in artifact else "?"
        for row in rows
        for artifact in row["analysis"].get("registered_reads", [])
    )
    print(
        "  ",
        dict(reads) or "none",
        "in",
        sum(1 for row in rows if row["analysis"].get("registered_reads")),
        "cycles",
    )
    print(
        "\nconstants selected:",
        dict(Counter(c for row in rows for c in row["analysis"]["constants"])),
    )
    print("\nargv task job words on approved nodes:")
    words = Counter(
        w
        for row in rows
        for node in row["nodes"]
        for w in node["argv_task_words"]
    )
    print("  ", dict(words) or "none")
    print(
        "\nnode kinds:",
        dict(
            Counter(node["node_kind"] for row in rows for node in row["nodes"])
        ),
    )
    print(
        "\nprogram:stage of approved nodes:",
        dict(
            Counter(
                f"{n['program']}:{n['stage']}"
                for row in rows
                for n in row["nodes"]
            )
        ),
    )
    print("\nlargest cross-program routes (X or L-cross, executed):")
    crossing = [
        r
        for r in rows
        if _executed(r) and ("X" in r["classes"] or "L-cross" in r["classes"])
    ]
    for row in sorted(crossing, key=lambda r: -r["max_heavy"])[:12]:
        print(
            f"  heavy {row['max_heavy']:3d}  {row['stratum']:10s} {row['shape'][:110]}  {row['bundle']}"
        )
    print("\nsettlements of goals with an executed X / L-cross cycle:")
    seen = set()
    for row in crossing:
        key = (row["root"], row["goal"])
        if row["goal"] and key not in seen:
            seen.add(key)
            print(
                f"  {row['stratum']:10s} {row['goal']:40s} {row['ledger']['settlement']}"
            )


#: The routes c5's EPISODE.md names for the replay whatever the shape rule
#: selects (matched against root + bundle path).
NAMED_ROUTES = (
    "pyscf-irc-20260920/goals/g2-hono",
    "xtb-ir-acetamide-pyscf-stability-r10",
    "r10/q11/goals/g1",
    "r10/q11/goals/g2",
    "standing-round/workspaces/e6-pcet-1/",
    "qualification/interop-fukui-path/",
    "qualification/pka-agent-path/",
)
SAMPLE_CAP = 40


def sample(paths: list[str]) -> None:
    """The pre-registered replay sample, one absolute bundle path per line.

    Executed cycles only: the first (by sorted root + bundle path) of every
    distinct route shape in M, X, L; the first of every distinct (program
    set, task-named operation set) in C, A, A-reg; the named routes first;
    above the cap, cross-program shapes before the rest, then by how many
    cycles share the shape.
    """

    rows = []
    for path in paths:
        for line in Path(path).read_text(encoding="utf-8").splitlines():
            if line.strip():
                rows.append(json.loads(line))
    rows = _distinct([row for row in rows if "error" not in row])
    executed = sorted(
        (row for row in rows if _executed(row)),
        key=lambda r: (r["root"], r["bundle"]),
    )
    representatives: dict[tuple, dict] = {}
    counts: Counter = Counter()
    for row in executed:
        keys = []
        if any(c in {"M", "X", "L-intra", "L-cross"} for c in row["classes"]):
            keys.append(("shape", row["shape"]))
        if any(c in {"C", "A", "A-reg"} for c in row["classes"]):
            task_ops = tuple(
                sorted(set(row["analysis"]["ops"]) & TASK_NAMED_OPERATIONS)
            )
            keys.append(("analysis", tuple(row["programs"]), task_ops))
        for key in keys:
            counts[key] += 1
            representatives.setdefault(key, row)

    def full(row: dict) -> str:
        return str(Path(row["root"]) / row["bundle"])

    named = [
        row
        for row in executed
        if any(marker in full(row) for marker in NAMED_ROUTES)
        and (
            any(c != "S" and c != "P" for c in row["classes"])
            or "qualification" in full(row)
        )
    ]

    def crossing(key: tuple, row: dict) -> bool:
        return (
            "X" in row["classes"]
            or "L-cross" in row["classes"]
            or (key[0] == "analysis" and len(key[1]) >= 2)
        )

    ranked = sorted(
        representatives.items(),
        key=lambda kv: (not crossing(*kv), -counts[kv[0]], full(kv[1])),
    )
    chosen: list[tuple[str, str]] = []
    seen: set[str] = set()
    for row in named:
        if full(row) not in seen and len(chosen) < SAMPLE_CAP:
            seen.add(full(row))
            chosen.append((full(row), "named"))
    for key, row in ranked:
        if full(row) not in seen and len(chosen) < SAMPLE_CAP:
            seen.add(full(row))
            chosen.append((full(row), f"{key[0]} x{counts[key]}"))
    print(
        f"# {len(representatives)} distinct keys among {len(executed)} "
        f"executed distinct cycles; {len(named)} named-route cycles; "
        f"cap {SAMPLE_CAP}"
    )
    for path, why in chosen:
        print(f"{path}\t{why}")


def main(argv: list[str]) -> int:
    if len(argv) >= 2 and argv[0] == "--report":
        report(argv[1:])
        return 0
    if len(argv) >= 2 and argv[0] == "--sample":
        sample(argv[1:])
        return 0
    if len(argv) < 2:
        print(__doc__)
        return 2
    run_census(Path(argv[0]), argv[1:])
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
