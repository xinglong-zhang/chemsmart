"""The literal census counted once per distinct report.

    python .agents/research/loop/literal_unique.py LITERAL_CLAIMS_JSON

``literal_claims.py`` counts every report file under its roots. An archive
that copies a workspace counts one rendered report several times: in the
ax41 mirror one goal's two reports sit, byte for byte, in 12 workspaces of
later rounds (under three goal ids, with no goal record beside them), and
the mirror's research directory holds copies of the two public cases the
repository also carries. The other census instruments deduplicate goals by
their ``goal_sha256``; these copies have no goal record, so the one identity
they share is their content. This groups the rows of a ``literal_claims.json``
by the sha256 of the report file each came from, keeps the rows of one file
per digest (the first path in sorted order), and prints the census's own
counters over those rows, beside the number of report files, distinct
reports and rows dropped as copies. A report that can no longer be read
keeps its rows under its own path. Imports nothing from ``chemsmart``.
Written for R11 truth-4's final-pin census (item 0), from a duplication found
in census 7b's local file before the census ran.
"""

from __future__ import annotations

import collections
import hashlib
import json
import sys
from pathlib import Path


def digest(path: str) -> str:
    try:
        return hashlib.sha256(Path(path).read_bytes()).hexdigest()
    except OSError:
        return "unreadable:" + path


def main() -> None:
    rows = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
    paths = sorted({row["report"] for row in rows})
    first: dict[str, str] = {}
    for path in paths:
        first.setdefault(digest(path), path)
    kept_paths = set(first.values())
    kept = [row for row in rows if row["report"] in kept_paths]
    counts: collections.Counter = collections.Counter()
    for row in kept:
        counts[row.get("class") or "(no class)"] += 1
        if row.get("pure"):
            counts["pure"] += 1
        for literal in row.get("literals") or ():
            counts[
                "physical literal "
                + ("carried" if literal.get("carried") else "foreign")
            ] += 1
            if literal.get("scaled"):
                counts[
                    "physical literal foreign, post hoc: " + literal["scaled"]
                ] += 1
    print(
        f"report files {len(paths)}; distinct reports {len(first)}; "
        f"rows {len(rows)}; rows in distinct reports {len(kept)}; "
        f"rows dropped as copies {len(rows) - len(kept)}"
    )
    for key, n in sorted(counts.items()):
        print(f"{n:6d} {key}")


if __name__ == "__main__":
    main()
