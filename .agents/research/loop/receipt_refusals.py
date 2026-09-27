"""Was a receipt the host refused as "not one it minted" in fact minted by it?

    python .agents/research/loop/receipt_refusals.py OUT_DIR label=ROOT [...]

The decision gate (``decision.receipt_is_one_the_host_minted``) refuses a
citation of a receipt no registry of this session holds and no recorded
run stream of the workspace carries (``_recorded_run_receipt``). For every
such refusal the Agent met, this reads the digest prefix the refusal
names and searches every stream the goal wrote -- planning sessions and
runs, this cycle and earlier -- for an event that minted a receipt with
that prefix before the refusal. It classifies each refusal: minted in a
run stream; minted by an earlier planning session of the goal; minted
earlier in the same session; or minted nowhere the goal recorded. Only
the last is a refusal the records support. Imports nothing from
``chemsmart``. Written for R11 episode `truth`.
"""

from __future__ import annotations

import collections
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from word_reader import Goal, discover  # noqa: E402

GATE = "decision.receipt_is_one_the_host_minted"
_PREFIX = re.compile(r"\b([0-9a-f]{8})\b")


def gate_refusals(goal: Goal) -> list[dict]:
    """Every refusal at the gate in one goal's streams, classified by where
    the digest it named was minted before it (full digests included)."""

    out: list[dict] = []
    minted: list[tuple[str, str, str, str]] = []
    for stamp, role, path, event in goal.events:
        payload = event.get("payload") or {}
        digest = str(payload.get("receipt_sha256") or "")
        if digest and event.get("kind") not in {
            "tool_started",
            "tool_failed",
        }:
            minted.append((stamp, role, path, digest))
    for stamp, role, path, event in goal.events:
        if event.get("kind") != "tool_failed":
            continue
        payload = event.get("payload") or {}
        report = (
            payload.get("failure_report")
            or (payload.get("canonical_result") or {}).get("failure_report")
            or {}
        )
        if report.get("gate") != GATE:
            continue
        diagnosis = str(report.get("diagnosis") or "")
        match = _PREFIX.search(diagnosis)
        prefix = match.group(1) if match else ""
        before = [
            (s, r, p, d)
            for (s, r, p, d) in minted
            if prefix and d.startswith(prefix) and s <= stamp
        ]
        if not prefix:
            where = "no digest in the diagnosis"
        elif not before:
            where = "minted nowhere the goal recorded"
        elif any(p == path for (_s, _r, p, _d) in before):
            where = "minted earlier in the same session"
        elif any(r == "run" for (_s, r, _p, _d) in before):
            where = "minted in a run stream"
        else:
            where = "minted by an earlier planning session"
        kinds = sorted(
            {
                e.get("kind")
                for (_s, _r, p, d) in before
                for (_s2, _r2, p2, e) in goal.events
                if p2 == p
                and str((e.get("payload") or {}).get("receipt_sha256") or "")
                == d
            }
        )
        out.append(
            {
                "stamp": stamp,
                "session": Path(path).parent.name,
                "prefix": prefix,
                "where": where,
                "minted_by": kinds,
                "minted_in": sorted({str(p) for (_s, _r, p, _d) in before}),
                "digests": sorted({d for (_s, _r, _p, d) in before}),
                "diagnosis": diagnosis[:240],
                "tool": payload.get("tool"),
            }
        )
    return out


def main() -> None:
    out_dir = Path(sys.argv[1])
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for item in discover(sys.argv[2:]):
        goal = Goal(Path(item["agent"]), item["goal_id"])
        for refusal in gate_refusals(goal):
            rows.append(
                {
                    "label": item["label"],
                    "goal_id": item["goal_id"],
                    "session": refusal["session"],
                    "prefix": refusal["prefix"],
                    "where": refusal["where"],
                    "minted_by": refusal["minted_by"],
                    "diagnosis": refusal["diagnosis"],
                    "tool": refusal["tool"],
                }
            )
    (out_dir / "receipt_refusals.json").write_text(json.dumps(rows, indent=1))
    counts = collections.Counter(row["where"] for row in rows)
    print(f"{len(rows)} refusals at {GATE}")
    for where, n in counts.most_common():
        print(f"  {n:3d} {where}")
    for row in rows:
        print(
            f"  {row['goal_id'][:32]:32} {row['prefix']} {row['where'][:40]:40} {row['minted_by']}"
        )


if __name__ == "__main__":
    main()
