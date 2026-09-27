"""Re-read every refusal the decision gate made, through the imported tree.

    HOME=<fence> PYTHONPATH=<tree> python \
        .agents/research/loop/receipt_gate_replay.py OUT_DIR label=ROOT [...]

For every refusal at ``decision.receipt_is_one_the_host_minted`` the Agent
met (found and classified exactly as ``receipt_refusals.py`` finds them:
the digest prefix the refusal named, matched to the receipts the goal's
streams minted before it), the imported tree's own host is asked the two
questions the gate asks of each full digest: would it accept the digest
from a recorded stream (``_recorded_run_receipt``), and what does it say
the digest is (``_digest_names``). The replay host reads the archived
workspace in place as its ``run_evidence_root``, holds no session
registry, and is seeded with the anomalies the goal's ledger recorded
before the refusing session began, as a woken host is seeded. A digest
minted earlier in the same session is accepted by that session's own
registry, which no replay rebuilds, and is reported without a re-read.

For a digest no stream of the goal minted, the workspace's citable
streams and the seeded anomalies are searched for the prefix: anything
found there would change the diagnosis. Written for R11 episode `truth`
(Repair B's census).
"""

from __future__ import annotations

import collections
import json
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from receipt_refusals import gate_refusals  # noqa: E402
from word_reader import Goal, discover, rows, when  # noqa: E402

_CITABLE = (
    "goals/*/runs/*/events.jsonl",
    "runs/*/events.jsonl",
    "executions/*/events.jsonl",
)


def seeded_anomalies(goal: Goal, session: str) -> tuple[dict, ...]:
    """The anomalies the goal's ledger recorded before ``session`` began."""

    began = session.split("-")[1] if session.startswith("live-") else ""
    if len(began) >= 22:
        began = (
            f"{began[0:4]}-{began[4:6]}-{began[6:8]}T{began[9:11]}:"
            f"{began[11:13]}:{began[13:15]}.{began[15:21]}"
        )
    seen: dict[str, dict] = {}
    for entry in goal.ledger:
        if entry.get("kind") != "anomalies_observed":
            continue
        if began and when(str(entry.get("at") or "")) > began:
            continue
        for record in (entry.get("payload") or {}).get("anomalies") or ():
            digest = str(record.get("receipt_sha256") or "")
            if digest and digest not in seen:
                seen[digest] = dict(record)
    return tuple(seen.values())


def prefix_elsewhere(agent: Path, prefix: str) -> list[str]:
    """Citable streams of the workspace holding a receipt with ``prefix``."""

    found = []
    for pattern in _CITABLE:
        for stream in sorted(agent.glob(pattern)):
            text = stream.read_text(encoding="utf-8", errors="replace")
            if prefix not in text:
                continue
            for event in rows(stream):
                digest = str(
                    (event.get("payload") or {}).get("receipt_sha256") or ""
                )
                if digest.startswith(prefix) and not str(
                    event.get("kind") or ""
                ).startswith("tool_"):
                    found.append(f"{stream.parent.name}:{event.get('kind')}")
    return found


def main() -> None:
    import chemsmart
    from chemsmart.agent.runtime.event_store import RuntimeEventStore
    from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1

    out_dir = Path(sys.argv[1])
    out_dir.mkdir(parents=True, exist_ok=True)
    print("chemsmart imported from", chemsmart.__file__)
    scratch = Path(tempfile.mkdtemp(prefix="gate-replay-"))
    results = []
    for item in discover(sys.argv[2:]):
        agent = Path(item["agent"])
        goal = Goal(agent, item["goal_id"])
        for refusal in gate_refusals(goal):
            row = {
                "label": item["label"],
                "goal_id": item["goal_id"],
                **{key: refusal[key] for key in refusal if key != "stamp"},
            }
            seeded = seeded_anomalies(goal, refusal["session"])
            row["seeded_anomaly_prefixes"] = sorted(
                str(a.get("receipt_sha256") or "")[:8] for a in seeded
            )
            if refusal["where"] == "minted earlier in the same session":
                row["replayed"] = "not re-read: the session's own registry"
                results.append(row)
                continue
            host = CommandCompiledToolHostV1(
                event_store=RuntimeEventStore(
                    scratch / f"{len(results)}" / "events.jsonl",
                    session_id="gate-replay",
                ),
                artifacts={},
                task_spec_sha256s=("a" * 64,),
                approved_workspace=scratch / f"{len(results)}" / "ws",
                run_evidence_root=agent.parent,
                prior_anomaly_observations=seeded,
            )
            if refusal["digests"]:
                row["replayed"] = [
                    {
                        "digest": digest[:16],
                        "accepted": host._recorded_run_receipt(digest),
                        "diagnosis": host._digest_names(digest)
                        or "no digest this host minted",
                    }
                    for digest in refusal["digests"]
                ]
            else:
                row["elsewhere"] = prefix_elsewhere(agent, refusal["prefix"])
                row["seeded_with_prefix"] = [
                    p
                    for p in row["seeded_anomaly_prefixes"]
                    if p == refusal["prefix"]
                ]
                row["replayed"] = (
                    "no digest this host minted"
                    if not row["elsewhere"] and not row["seeded_with_prefix"]
                    else "FOUND ELSEWHERE"
                )
            results.append(row)
    (out_dir / "gate_replay.json").write_text(json.dumps(results, indent=1))
    counts: collections.Counter = collections.Counter()
    for row in results:
        replayed = row["replayed"]
        if isinstance(replayed, str):
            counts[(row["where"], replayed[:40])] += 1
            continue
        for item in replayed:
            text = item["diagnosis"]
            counts[
                (
                    row["where"],
                    "accepted" if item["accepted"] else "refused",
                    "names anomaly: route" if "anomaly:" in text else "",
                )
            ] += 1
    print(f"{len(results)} refusals re-read")
    for key, n in sorted(counts.items()):
        print(f"  {n:3d} {key}")
    for row in results:
        replayed = row["replayed"]
        if isinstance(replayed, str):
            continue
        for item in replayed:
            print(
                f"  {row['goal_id'][:24]:24} {row['prefix']} "
                f"{'ACCEPT' if item['accepted'] else 'refuse'} "
                f"{item['diagnosis'][:170]}"
            )


if __name__ == "__main__":
    main()
