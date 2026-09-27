"""Re-sign archived goals' settlement words on the imported tree.

    HOME=<fence> PYTHONPATH=<tree> python .agents/research/loop/resign.py \
        OUT_DIR label=ROOT [label=ROOT ...]        # discover goal ledgers
    ... resign.py OUT_DIR @SPECS                   # lines "label|agent_dir|goal_id"

For every goal ledger under the roots (sealed, private and ``claude/``
paths excluded; deduplicated by goal id and goal digest), the goal's final
settle step is replayed through the imported tree's own driver: the ledger
is cut before its last ``goal_settled`` (or the ``reading_opened`` row that
held it), the cross-cycle state ``GoalDriver.resume`` rebuilds is restored,
and ``GoalDriver._settle`` -- or ``_delivery_settlement`` on the planning
path -- signs again. The executor's analysis word comes from the imported
tree's own walk when the approved chain dispatches nothing (an empty or
all-blocked chain), else from what the run's stream recorded (R10 Q24's
"walk" mode; its "stream" mode restated the signer and is not offered).

Each replayed word is compared with the archived one in the buckets R11
`truth` pre-registered: (a) the state differs; (b) the state is the same
and the reasons' content differs -- the sets of identifier-like tokens,
8-hex receipt prefixes and numbers they name; (c) only wording differs;
identical. Results stream to OUT_DIR/results.jsonl, one goal per line, so
a long census survives an interruption.

HOME must be fenced by the caller: a replayed achieved word writes
qualification rows, and the host store's path is bound at import. The
imported tree must be the one PYTHONPATH names first, or this refuses.
Adapted from R10 Q24's replay_settle.py (instruments/q24/tools), which was
adapted from R10 Q22's replay_final.py.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import sys
import traceback
from pathlib import Path
from types import SimpleNamespace

EXCLUDED_PARTS = ("sealed", "/private", "/claude/", "/grade-", "/grading/")


def rows(path: Path) -> list[dict]:
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return []
    out = []
    for line in text.splitlines():
        if line.strip():
            try:
                out.append(json.loads(line))
            except json.JSONDecodeError:
                continue
    return out


def write_rows(path: Path, items: list[dict]) -> None:
    path.write_text(
        "".join(json.dumps(item) + "\n" for item in items), encoding="utf-8"
    )


# -- discovery -------------------------------------------------------------


def discover(specs: list[str]) -> list[dict]:
    """Goal ledgers under label=ROOT specs, or listed in @file specs."""

    found: list[dict] = []
    for spec in specs:
        if spec.startswith("@"):
            for line in Path(spec[1:]).read_text().splitlines():
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                label, agent_dir, goal_id = line.split("|")
                found.append(
                    {"label": label, "agent": agent_dir, "goal_id": goal_id}
                )
            continue
        label, root = spec.split("=", 1)
        for folder, dirs, files in os.walk(root):
            dirs[:] = [
                d
                for d in dirs
                if not any(p.strip("/") in d for p in ("sealed", "private"))
                and d not in {"scratch", "claude"}
            ]
            if "ledger.jsonl" not in files:
                continue
            ledger = Path(folder) / "ledger.jsonl"
            text = str(ledger)
            if "/.chemsmart-agent/goals/" not in text:
                continue
            if any(part in text for part in EXCLUDED_PARTS):
                continue
            found.append(
                {
                    "label": label,
                    "agent": str(ledger.parents[2]),
                    "goal_id": ledger.parent.name,
                }
            )
    unique: dict[tuple[str, str], dict] = {}
    for item in found:
        goal_file = (
            Path(item["agent"]) / "goals" / item["goal_id"] / "goal.json"
        )
        try:
            digest = json.loads(goal_file.read_text()).get("goal_sha256", "")
        except (OSError, json.JSONDecodeError):
            digest = ""
        item["goal_sha256"] = digest
        key = (item["goal_id"], digest or item["agent"])
        if key in unique:
            unique[key].setdefault("copies", []).append(item["agent"])
            continue
        unique[key] = item
    return sorted(unique.values(), key=lambda i: (i["label"], i["agent"]))


# -- buckets ---------------------------------------------------------------

_HEX = re.compile(r"(?<![0-9a-f])([0-9a-f]{8})[0-9a-f]*(?![0-9a-f])")
_NUMBER = re.compile(r"(?<![\w.])-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?")
_TOKEN = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.:/+-]*[A-Za-z0-9]")


def content(reasons: list[str]) -> dict[str, list[str]]:
    """What a word's reasons name, apart from how they say it."""

    text = " ".join(str(item) for item in reasons or ())
    hexes = set(_HEX.findall(text))
    numbers = set()
    for raw in _NUMBER.findall(text):
        try:
            numbers.add(f"{float(raw):.6g}")
        except ValueError:
            continue
    tokens = {
        token
        for token in _TOKEN.findall(text)
        if any(mark in token for mark in "-_:.")
        and not _NUMBER.fullmatch(token)
    }
    return {
        "ids": sorted(tokens),
        "hex": sorted(hexes),
        "numbers": sorted(numbers),
    }


def bucket(archived: dict, replayed: dict | None) -> str:
    if not replayed:
        return "not_replayed"
    if replayed.get("state") != archived.get("state"):
        return "a_state"
    old = list(archived.get("reasons") or ())
    new = list(replayed.get("reasons") or ())
    if old == new:
        return "identical"
    if content(old) != content(new):
        return "b_content"
    return "c_wording"


# -- the replay (R10 Q24's, with its comments) -----------------------------


def stream_stamp(run_id: str) -> str:
    try:
        return run_id.split("-")[1][:15]
    except IndexError:
        return ""


def iso_stamp(text: str) -> str:
    return (
        text.replace("-", "").replace(":", "").split(".")[0].split("+")[0][:15]
    )


def terminal_of(stream: Path) -> str:
    terminal = ""
    for event in rows(stream):
        if event.get("kind") == "runtime_terminated":
            terminal = str(
                (event.get("payload") or {}).get("terminal_state") or ""
            )
    return terminal


def envelope_file(out: Path, goal_record: dict) -> Path:
    granted = goal_record["envelope"]
    lines = [
        "schema_version: chemsmart.bounded-execution-envelope.v1",
        "mode: bounded-local",
        "allowed_program_engines:",
    ]
    for program, engines in granted["allowed_program_engines"]:
        lines.append(f"  {program}:")
        lines.extend(f"  - {engine}" for engine in engines)
    lines += [
        "resources:",
        "  execution_target: run",
        f"  cores: {int(granted.get('cores') or 4)}",
        f"  memory_gb: {int(granted.get('memory_gb') or 16)}",
        "  gpu_count: 0",
        "  scratch_policy: server",
        "  node_timeout_seconds: 600",
        f"episode_wall_time_seconds: {int(granted['episode_wall_time_seconds'])}",
        "postprocess_reserve_seconds: 600",
        f"max_engine_calls: {int(granted['max_engine_calls'])}",
        f"scratch_root: {out / 'scratch'}",
    ]
    if int(granted.get("max_excursion_calls") or 0):
        lines.append(
            f"max_excursion_calls: {int(granted['max_excursion_calls'])}"
        )
    path = out / "envelope.yaml"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def bundle_path(
    agent: Path, kept: list[dict], goal_id: str, cycle: int
) -> Path | None:
    started = [
        e
        for e in kept
        if e["kind"] == "run_started"
        and int((e.get("payload") or {}).get("cycle") or 0) == cycle
    ]
    if started:
        recorded = str(
            (started[-1].get("payload") or {}).get("approval_file") or ""
        )
        if "/.chemsmart-agent/" in recorded:
            candidate = agent / recorded.split("/.chemsmart-agent/", 1)[1]
            if candidate.is_file():
                return candidate
    candidate = (
        agent / "replays" / f"goal-{goal_id}-cycle-{cycle}" / "bundle.json"
    )
    return candidate if candidate.is_file() else None


def executor_word_walk(
    bundle: Path | None, run_rows: list[dict], scratch: Path
) -> tuple[str, str]:
    """The imported tree's own walk over a chain that dispatches nothing;
    otherwise the recorded stream's word. Returns (word, how)."""

    recorded = any(
        r.get("kind")
        in {
            "workflow_analysis_node_settled",
            "workflow_analysis_report_rendered",
            "workflow_analysis_completion_refused",
            "analysis_completion_evaluated",
        }
        for r in run_rows
    )
    toolchain = None
    if not recorded and bundle is not None:
        from chemsmart.agent.live_session import (
            load_workflow_execution_approval_bundle,
        )

        loaded = load_workflow_execution_approval_bundle(bundle)
        toolchain = getattr(loaded, "scientific_toolchain_plan", None)
        if toolchain is None:
            return "", "bundle carries no toolchain"
    if toolchain is None and not recorded:
        return "", "no bundle and no analysis event"
    runnable = (
        [
            node
            for node in toolchain.analysis_nodes
            if getattr(node, "support_state", "") != "blocked_unsupported"
        ]
        if toolchain is not None
        else [None]
    )
    if not runnable:
        from chemsmart.agent.executor import ApprovedWorkflowExecutor
        from chemsmart.agent.runtime.event_store import RuntimeEventStore

        scratch.mkdir(parents=True, exist_ok=True)
        store = RuntimeEventStore(
            scratch / "walk-events.jsonl", session_id="replay-walk"
        )
        host = SimpleNamespace(
            event_store=store,
            artifacts={},
            execution_receipts={},
            activate_guides=lambda *args, **kwargs: None,
        )
        stub = SimpleNamespace(
            host=host,
            plan=SimpleNamespace(
                workflow_id=toolchain.workflow_id, plan_sha256="", nodes=()
            ),
            run_directory=scratch,
            approval_workspace=scratch,
        )
        _nodes, word, _receipts, _report = (
            ApprovedWorkflowExecutor._run_analysis_phase(stub, toolchain)
        )
        return word, (
            f"walked {len(toolchain.analysis_nodes)} analysis node(s), "
            "none runnable"
        )
    reports = [
        str((r.get("payload") or {}).get("report_kind") or "completed")
        for r in run_rows
        if r.get("kind") == "workflow_analysis_report_rendered"
    ]
    if reports:
        return (
            "partial" if reports[-1] == "partial" else "completed"
        ), "recorded report kind"
    if any(
        r.get("kind") == "workflow_analysis_completion_refused"
        for r in run_rows
    ):
        return "partial", "recorded completion refusal"
    states = {}
    for r in run_rows:
        if r.get("kind") == "workflow_analysis_node_settled":
            payload = r.get("payload") or {}
            states[str(payload.get("node_id"))] = str(payload.get("state"))
    if states and set(states.values()) <= {"executed", "blocked_unsupported"}:
        return "completed", "recorded node settlements, no report"
    return "partial", "recorded node settlements"


def replay(source: Path, goal_id: str, out: Path) -> dict:
    from chemsmart.agent import driver as drv
    from chemsmart.agent.goal import GoalLedger

    if out.exists():
        shutil.rmtree(out)
    workspace = out / "ws"
    agent = workspace / ".chemsmart-agent"
    shutil.copytree(
        source,
        agent,
        ignore=shutil.ignore_patterns(
            "*.gbw", "*.chk", "*.tmp", "scratch", "*.densities*"
        ),
        symlinks=True,
    )
    # The engine outputs the workspace record names live beside the agent
    # directory; link every sibling so result discovery reads them.
    for item in source.parent.iterdir():
        if item.name == ".chemsmart-agent":
            continue
        target = workspace / item.name
        if not target.exists():
            target.symlink_to(item.resolve())
    ledger_path = agent / "goals" / goal_id / "ledger.jsonl"
    entries = rows(ledger_path)
    at_cycle = int(os.environ.get("RESIGN_AT_CYCLE") or 0)
    if at_cycle:
        recorded = min(
            i
            for i, e in enumerate(entries)
            if e["kind"] == "run_recorded"
            and int((e.get("payload") or {}).get("cycle") or 0) == at_cycle
        )
        settled_at = min(
            i
            for i, e in enumerate(entries)
            if i > recorded
            and e["kind"]
            in {"goal_settled", "recovery_opened", "reading_opened"}
        )
        row = entries[settled_at]
        archived_row = {
            "state": (row.get("payload") or {}).get("state", row["kind"]),
            "reasons": list((row.get("payload") or {}).get("reasons") or ()),
        }
        entries = entries[: settled_at + 1]
        entries[settled_at] = {
            "kind": "goal_settled",
            "at": row.get("at"),
            "payload": archived_row,
        }
    settled_at = max(
        i for i, e in enumerate(entries) if e["kind"] == "goal_settled"
    )
    archived = dict(entries[settled_at]["payload"])
    archived.pop("evidence", None)
    cut = settled_at
    for index in range(settled_at - 1, -1, -1):
        kind = entries[index]["kind"]
        if kind == "reading_opened":
            payload = entries[index]["payload"]
            archived = {
                "state": payload.get("state"),
                "reasons": list(payload.get("reasons") or ()),
                "held_by_reading": True,
            }
            cut = index
            break
        if kind in {"run_recorded", "recovery_opened", "wake_composed"}:
            break
    settled_stamp = iso_stamp(str(entries[settled_at].get("at") or ""))
    kept = entries[:cut]
    write_rows(ledger_path, kept)
    archived_qualified = [
        str(e["payload"].get("id")) + "@" + str(e["payload"].get("node"))
        for e in entries[cut:]
        if e["kind"] == "qualified"
    ]
    cycles = max(
        [int((e.get("payload") or {}).get("cycle") or 0) for e in kept] + [1]
    )
    reasons0 = str((archived.get("reasons") or [""])[0])
    if ", planning session:" in reasons0 or ", goal re-wake:" in reasons0:
        return {"archived": archived, "replayed": None, "path": "typed-error"}
    run_rows_for_cycle = [
        e
        for e in kept
        if e["kind"] == "run_recorded"
        and int(e["payload"].get("cycle") or 0) == cycles
    ]
    named = [
        e["payload"]["run_id"]
        for e in kept
        if e["kind"] == "session_stream_recorded"
        and int(e["payload"].get("cycle") or 0) == cycles
    ]
    if named:
        session_stream = agent / "runs" / named[-1] / "events.jsonl"
    else:
        candidates = sorted(
            p
            for p in (agent / "runs").glob("live-*/events.jsonl")
            if not settled_stamp
            or stream_stamp(p.parent.name) <= settled_stamp
        )
        session_stream = candidates[-1] if candidates else None
    goal_record = json.loads(
        (agent / "goals" / goal_id / "goal.json").read_text(encoding="utf-8")
    )
    run_state = (
        run_rows_for_cycle[-1]["payload"].get("workflow_state")
        if run_rows_for_cycle
        else ""
    )
    if not run_rows_for_cycle or run_state == "analysis_only":
        stream = session_stream
        if run_state == "analysis_only":
            stream = (
                agent
                / "goals"
                / goal_id
                / "runs"
                / f"cycle-{cycles}"
                / "events.jsonl"
            )
        terminal = terminal_of(session_stream) if session_stream else ""
        if stream is None:
            return {
                "archived": archived,
                "replayed": None,
                "path": "planning-no-stream",
            }
        ledger = GoalLedger(agent / "goals" / goal_id)
        settled, reasons, _evidence = drv._delivery_settlement(
            ledger,
            goal_id=goal_id,
            events_path=stream,
            terminal=terminal,
            workspace=workspace,
        )
        return {
            "archived": archived,
            "archived_qualified": archived_qualified,
            "replayed": {"state": settled, "reasons": list(reasons)},
            "path": "planning",
            "stream": stream.parent.name,
            "terminal": terminal,
        }
    envelope = envelope_file(out, goal_record)
    driver_file = agent / "goals" / goal_id / "driver.json"
    if driver_file.exists():
        record = json.loads(driver_file.read_text(encoding="utf-8"))
        record["execution_envelope_file"] = str(envelope)
        record["stop_file"] = str(out / "STOP")
        driver_file.write_text(json.dumps(record), encoding="utf-8")
    empty = out / "empty-ws"
    empty.mkdir(parents=True, exist_ok=True)
    driver = drv.GoalDriver(
        task="replay",
        workspace=empty,
        execution_envelope_file=envelope,
        goal_id=goal_id,
        granted_by=str(goal_record.get("granted_by") or "replay"),
        max_revisions=int(goal_record.get("max_revisions") or 0),
        plan_session=lambda **kw: None,
        resolve_review=lambda **kw: None,
        execute_bundle=lambda **kw: None,
    )
    driver.workspace = workspace
    driver.goal_dir = agent / "goals" / goal_id
    driver.ledger = GoalLedger(driver.goal_dir)
    driver.goal = driver.ledger.load()
    driver.cycles = cycles
    driver.revisions_admitted = sum(
        1 for e in kept if e["kind"] == "revision_admitted"
    )
    earlier_runs = [
        item
        for item in kept
        if item["kind"] == "run_recorded"
        and int(item["payload"].get("cycle") or 0) < cycles
    ]
    if hasattr(driver, "_restore_standing_delivery"):
        # The tree's own restoration, the one its resume() calls, so the
        # replay rebuilds exactly the state the tree would.
        driver._restore_standing_delivery(earlier_runs)
    else:
        # Trees before that function (to c7102cc2): the loop their
        # resume() ran, restated.
        for item in earlier_runs:
            delivery = drv._analysis_delivery(
                agent
                / Path(*str(item["payload"].get("run") or "").split("/"))
                / "events.jsonl",
                inherited_rejected_artifacts=tuple(
                    sorted(driver.rejected_artifacts)
                ),
            )
            driver.rejected_artifacts.update(
                delivery.rejected_artifact_sha256s
            )
            if delivery.claims_rendered:
                driver.standing_stale = delivery.stale_quantity_ids
    run_directory = agent / "goals" / goal_id / "runs" / f"cycle-{cycles}"
    run_rows = rows(run_directory / "events.jsonl")
    word, how = executor_word_walk(
        bundle_path(agent, kept, goal_id, cycles), run_rows, out / "walk"
    )
    driver.run_directory = run_directory
    from chemsmart.agent.terminal_states import (
        derive_run_outcome,
        read_run_events,
    )

    try:
        driver.outcome = derive_run_outcome(
            read_run_events(run_directory / "events.jsonl")
        )
    except Exception:  # noqa: BLE001 - a stream with no workflow run
        driver.outcome = None
    driver.events_path = session_stream
    driver.execute_result = SimpleNamespace(
        status="completed" if run_state == "validated" else "partial",
        analysis_status=word,
    )
    driver.reading_turn = bool(archived.get("held_by_reading"))
    driver._settle()
    after = rows(ledger_path)[len(kept) :]
    replayed = None
    for e in after:
        if e["kind"] in ("goal_settled", "recovery_opened", "reading_opened"):
            p = e.get("payload") or {}
            replayed = {
                "kind": e["kind"],
                "state": p.get("state", e["kind"]),
                "reasons": list(p.get("reasons") or []),
            }
            break
    return {
        "archived": archived,
        "archived_qualified": archived_qualified,
        "replayed": replayed,
        "replayed_qualified": [
            str(e["payload"].get("id")) + "@" + str(e["payload"].get("node"))
            for e in after
            if e["kind"] == "qualified"
        ],
        "path": "run",
        "stream": run_directory.name,
        "executor_word": word,
        "executor_word_how": how,
    }


def main() -> None:
    if sys.argv[1] == "--list":
        # The census population, written before anything is replayed.
        goals = discover(sys.argv[3:])
        Path(sys.argv[2]).write_text(
            "".join(
                f"{g['label']}|{g['agent']}|{g['goal_id']}\n" for g in goals
            )
        )
        print(f"listed {len(goals)} goal(s) after dedupe")
        return
    import chemsmart

    imported = str(Path(chemsmart.__file__).resolve())
    expected = os.environ.get("PYTHONPATH", "").split(os.pathsep)[0]
    print("chemsmart imported from", chemsmart.__file__, flush=True)
    if not expected or not imported.startswith(str(Path(expected).resolve())):
        raise SystemExit(
            f"refusing: imported {imported}, PYTHONPATH names {expected!r}"
        )
    out_root = Path(sys.argv[1])
    out_root.mkdir(parents=True, exist_ok=True)
    goals = discover(sys.argv[2:])
    print(f"goals: {len(goals)} after dedupe", flush=True)
    results_file = out_root / "results.jsonl"
    with results_file.open("w", encoding="utf-8") as handle:
        for index, item in enumerate(goals):
            source = Path(item["agent"])
            goal_id = item["goal_id"]
            ledger = rows(source / "goals" / goal_id / "ledger.jsonl")
            result: dict = {}
            if not any(e.get("kind") == "goal_settled" for e in ledger):
                result = {"path": "unsettled"}
            else:
                work = (
                    out_root
                    / "work"
                    / hashlib.sha256(
                        f"{item['agent']}|{goal_id}".encode()
                    ).hexdigest()[:16]
                )
                try:
                    result = replay(source, goal_id, work)
                except Exception as exc:  # noqa: BLE001 - reported
                    result = {
                        "path": "error",
                        "error": f"{type(exc).__name__}: {str(exc)[:400]}",
                        "trace": traceback.format_exc()[-1500:],
                    }
                finally:
                    shutil.rmtree(work, ignore_errors=True)
            result.update(item)
            if "archived" in result:
                result["bucket"] = bucket(
                    result["archived"], result.get("replayed")
                )
            handle.write(json.dumps(result, default=str) + "\n")
            handle.flush()
            print(
                f"{index + 1:4d} {item['label']:10} {goal_id[:34]:34} "
                f"{result.get('path', '-'):9} {result.get('bucket', '-')}",
                flush=True,
            )


if __name__ == "__main__":
    main()
