"""Every refusal the Agent met, grouped by what refused it and how it ended.

    python .agents/research/loop/refusal_census.py OUT_DIR label=ROOT [...]

Reads every ``tool_failed`` event in the streams of every goal under the
roots (the same discovery and deduplication as ``word_reader.py``, whose
stream reading it reuses) and groups them by (tool, error class, gate) --
the gate a routed refusal names, else the refusal message with quoted
names, digests and numbers folded, so one rule is one group. For each group
it counts instances, goals and sessions; whether the session's next call of
the same tool succeeded (recovered in session); and it keeps the first
examples' invariant, diagnosis and route, or message. Classifying a group
-- protects an invariant, resolves a wrong reference, refuses a legitimate
choice, or names a route the host would refuse -- is a reading of those
examples against the code, and is written in the episode record, never
here. Written for R11 episode `truth`.
"""

from __future__ import annotations

import collections
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from word_reader import Goal, discover  # noqa: E402

_QUOTED = re.compile(r"'[^']*'|\"[^\"]*\"")
_HEX = re.compile(r"\b[0-9a-f]{8,64}\b")
_NUM = re.compile(r"-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?")


def template(text: str) -> str:
    text = _QUOTED.sub("<q>", str(text))
    text = _HEX.sub("<h>", text)
    text = _NUM.sub("<n>", text)
    return re.sub(r"\s+", " ", text).strip()[:110]


def main() -> None:
    out_dir = Path(sys.argv[1])
    out_dir.mkdir(parents=True, exist_ok=True)
    groups: dict[tuple[str, str, str], dict] = {}
    total = 0
    for item in discover(sys.argv[2:]):
        goal = Goal(Path(item["agent"]), item["goal_id"])
        by_stream: dict[str, list[dict]] = collections.defaultdict(list)
        for _stamp, role, path, event in goal.events:
            by_stream[path].append(event)
        for path, events in by_stream.items():
            for index, event in enumerate(events):
                if event.get("kind") != "tool_failed":
                    continue
                payload = event.get("payload") or {}
                result = payload.get("canonical_result") or {}
                report = payload.get("failure_report") or result.get(
                    "failure_report"
                ) or {}
                tool = str(payload.get("tool") or result.get("tool") or "")
                error = str(payload.get("error_class") or result.get("error_class") or "")
                gate = str(report.get("gate") or "")
                message = str(result.get("message") or report.get("diagnosis") or "")
                key = (tool, error, gate or template(message))
                recovered = None
                for later in events[index + 1 :]:
                    later_payload = later.get("payload") or {}
                    later_tool = str(
                        later_payload.get("tool")
                        or (later_payload.get("canonical_result") or {}).get("tool")
                        or ""
                    )
                    if later.get("kind") == "tool_succeeded" and later_tool == tool:
                        recovered = True
                        break
                    if later.get("kind") == "tool_failed" and later_tool == tool:
                        recovered = False
                        break
                group = groups.setdefault(
                    key,
                    {
                        "tool": tool,
                        "error_class": error,
                        "gate": gate,
                        "template": template(message) if not gate else "",
                        "instances": 0,
                        "goals": set(),
                        "streams": set(),
                        "recovered": 0,
                        "retried_and_refused": 0,
                        "not_retried": 0,
                        "examples": [],
                    },
                )
                total += 1
                group["instances"] += 1
                group["goals"].add(f"{item['label']}:{item['goal_id']}")
                group["streams"].add(path)
                if recovered is True:
                    group["recovered"] += 1
                elif recovered is False:
                    group["retried_and_refused"] += 1
                else:
                    group["not_retried"] += 1
                if len(group["examples"]) < 3:
                    group["examples"].append(
                        {
                            "goal": f"{item['label']}:{item['goal_id']}",
                            "invariant": str(report.get("invariant") or "")[:300],
                            "diagnosis": str(report.get("diagnosis") or "")[:400],
                            "route": str(report.get("route") or "")[:300],
                            "message": message[:500] if not gate else "",
                        }
                    )
    rows = []
    for group in groups.values():
        group["goals"] = sorted(group["goals"])
        group["streams"] = len(group["streams"])
        rows.append(group)
    rows.sort(key=lambda g: -g["instances"])
    (out_dir / "refusals.json").write_text(json.dumps(rows, indent=1))
    lines = [f"refusals {total} in {len(rows)} groups"]
    for g in rows:
        lines.append(
            f"{g['instances']:4d} goals {len(g['goals']):3d} rec {g['recovered']:3d} "
            f"again {g['retried_and_refused']:3d} none {g['not_retried']:3d}  "
            f"{g['tool'][:30]:30} {g['error_class'][:24]:24} {g['gate'] or g['template']}"
        )
    (out_dir / "refusals.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines[:60]))


if __name__ == "__main__":
    main()
