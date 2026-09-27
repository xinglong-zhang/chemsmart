"""Which archived certificates named claims on a verdict the goal had answered.

    PYTHONPATH=<tree> python .agents/research/loop/certificate_census.py \
        OUT_DIR @SPECS [@SPECS ...]        # lines "label|agent_dir|goal_id"

R11 truth-3, item 1 (census 8). A completion certificate names every claim
that stands on a failed acceptance criterion no recorded decision has
answered, and the settlement trusts what it names. The certificate asked
that of its own host's records; the settlement asks it of every stream the
goal's ledger names. This reads, for every ``analysis_completion_evaluated``
row whose status is ``partial`` and whose findings name
``analysis.claim_on_failed_criterion.*``, in every stream a goal's ledger
names:

- host grain: the validations, decisions and lineage rows of the same
  stream before the completion row -- what the minting host held;
- faithful: the host-grain verdicts' ids equal the ``failed_criterion:``
  ids the completion itself carries (else "record insufficient", never
  counted either way);
- goal grain: the same joined with every other ledger-named stream's rows
  stamped before the completion;
- flip: every verdict unanswered at host grain is answered at goal grain
  (``flip_status`` when the completion's other findings are none, so the
  certificate itself would have passed);
- read_by: the ledger rows whose payload names the completion's 8-hex
  prefix, and whether it is its stream's last completion.

The verdict join and the stream set are the imported tree's own
(``chemsmart.agent.goal.failed_criteria``; ``chemsmart.agent.driver``'s
``_goal_streams``, ``_verdict_records``, ``_merge_verdict_records``), so
the census asks exactly what the settlement asks. Records only: nothing is
replayed, written or re-minted. Results stream to OUT_DIR/certificates.jsonl
and OUT_DIR/summary.json.
"""

from __future__ import annotations

import json
import sys
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from resign import discover  # noqa: E402

_PREFIX = "analysis.claim_on_failed_criterion."


def _stamp(value: str) -> datetime | None:
    try:
        return datetime.fromisoformat(str(value).replace("Z", "+00:00"))
    except ValueError:
        return None


def _lines(path: Path) -> list[str]:
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return []
    return [line for line in text.splitlines() if line.strip()]


def _parsed(line: str) -> dict:
    try:
        return json.loads(line)
    except json.JSONDecodeError:
        return {}


def _kind_of(agent: Path, stream: Path) -> str:
    relative = str(stream.relative_to(agent))
    return "run" if relative.startswith("goals/") else "session"


def census_goal(spec: dict, driver, failed_criteria) -> list[dict]:
    agent = Path(spec["agent"])
    goal_id = spec["goal_id"]
    goal_dir = agent / "goals" / goal_id
    ledger = driver.GoalLedger(goal_dir)
    entries = ledger.entries()
    streams = driver._goal_streams(ledger, agent.parent, goal_id)
    lines_of = {path: _lines(path) for path in streams}
    found: list[dict] = []
    for stream, lines in lines_of.items():
        completions = [
            index
            for index, line in enumerate(lines)
            if '"analysis_completion_evaluated"' in line
        ]
        for index in completions:
            event = _parsed(lines[index])
            if event.get("kind") != "analysis_completion_evaluated":
                continue
            payload = event.get("payload") or {}
            findings = list(
                (payload.get("record") or {}).get("findings") or ()
            )
            named = [item for item in findings if item.startswith(_PREFIX)]
            if payload.get("status") != "partial" or not named:
                continue
            at = _stamp(event.get("timestamp") or "")
            host = driver._verdict_records(lines[:index])
            host_verdicts = failed_criteria(
                host.validations,
                cited=host.cited,
                result_artifacts=host.result_artifacts,
                expression_sources=host.expression_sources,
            )
            carried = sorted(
                item
                for item in payload.get("anomaly_output_ids") or ()
                if str(item).startswith("failed_criterion:")
            )
            faithful = carried == sorted(
                verdict.observation_id for verdict in host_verdicts
            )
            own = {
                str(item.get("receipt_sha256") or "")
                for item in host.validations
            }
            others = []
            for other, other_lines in lines_of.items():
                if other == stream:
                    continue
                kept = [
                    line
                    for line in other_lines
                    if at is not None
                    and (_stamp(_parsed(line).get("timestamp") or "") or at)
                    < at
                ]
                others.append(driver._verdict_records(kept))
            goal = driver._merge_verdict_records(host, *others)
            goal_verdicts = [
                verdict
                for verdict in failed_criteria(
                    goal.validations,
                    cited=goal.cited,
                    result_artifacts=goal.result_artifacts,
                    expression_sources=goal.expression_sources,
                )
                if own.intersection(verdict.receipt_sha256s)
            ]
            host_open = sorted(
                verdict.label
                for verdict in host_verdicts
                if not verdict.answered
            )
            goal_open = sorted(
                verdict.label
                for verdict in goal_verdicts
                if not verdict.answered
            )
            answered_elsewhere = sorted(set(host_open) - set(goal_open))
            receipt = str(payload.get("receipt_sha256") or "")
            read_by = sorted(
                {
                    str(entry.get("kind") or "")
                    for entry in entries
                    if receipt[:8]
                    and receipt[:8] in json.dumps(entry.get("payload") or {})
                }
            )
            flip = bool(host_open) and not goal_open
            found.append(
                {
                    "label": spec["label"],
                    "goal_id": goal_id,
                    "agent": str(agent),
                    "stream": str(stream.relative_to(agent)),
                    "stream_kind": _kind_of(agent, stream),
                    "completion": receipt[:16],
                    "timestamp": event.get("timestamp") or "",
                    "faithful": faithful,
                    "carried_ids": carried,
                    "host_unanswered": host_open,
                    "goal_unanswered": goal_open,
                    "answered_elsewhere": answered_elsewhere,
                    "flip": flip,
                    "flip_status": flip
                    and all(item.startswith(_PREFIX) for item in findings),
                    "other_findings": [
                        item
                        for item in findings
                        if not item.startswith(_PREFIX)
                    ],
                    "last_in_stream": index == completions[-1],
                    "read_by": read_by,
                }
            )
    return found


def main(argv: list[str]) -> int:
    import chemsmart
    from chemsmart.agent import driver
    from chemsmart.agent.goal import failed_criteria

    out = Path(argv[1])
    out.mkdir(parents=True, exist_ok=True)
    specs = discover(argv[2:])
    print("chemsmart from", chemsmart.__file__, file=sys.stderr)
    results: list[dict] = []
    errors: list[dict] = []
    with (out / "certificates.jsonl").open("w", encoding="utf-8") as handle:
        for spec in specs:
            try:
                found = census_goal(spec, driver, failed_criteria)
            except (
                Exception
            ) as exc:  # noqa: BLE001 -- a census row, never a crash
                errors.append(
                    {
                        "goal_id": spec["goal_id"],
                        "agent": spec["agent"],
                        "error": f"{type(exc).__name__}: {exc}",
                    }
                )
                continue
            for row in found:
                handle.write(json.dumps(row, sort_keys=True) + "\n")
            results.extend(found)
    counted = [row for row in results if row["faithful"]]
    summary = {
        "goals": len(specs),
        "errors": errors,
        "partial_completions_on_a_criterion": len(results),
        "record_insufficient": len(results) - len(counted),
        "faithful": len(counted),
        "by_stream_kind": {
            kind: {
                "faithful": sum(
                    1 for r in counted if r["stream_kind"] == kind
                ),
                "answered_elsewhere": sum(
                    1
                    for r in counted
                    if r["stream_kind"] == kind and r["answered_elsewhere"]
                ),
                "flip": sum(
                    1
                    for r in counted
                    if r["stream_kind"] == kind and r["flip"]
                ),
                "flip_status": sum(
                    1
                    for r in counted
                    if r["stream_kind"] == kind and r["flip_status"]
                ),
            }
            for kind in ("run", "session")
        },
        "flips": [
            {
                key: row[key]
                for key in (
                    "label",
                    "goal_id",
                    "stream",
                    "completion",
                    "answered_elsewhere",
                    "flip_status",
                    "last_in_stream",
                    "read_by",
                )
            }
            for row in counted
            if row["answered_elsewhere"]
        ],
    }
    (out / "summary.json").write_text(
        json.dumps(summary, indent=1, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary["by_stream_kind"], sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
