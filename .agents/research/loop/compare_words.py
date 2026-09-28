"""Compare two re-signings of one population, goal by goal.

    python .agents/research/loop/compare_words.py BASE NEW [--list]

BASE and NEW are ``resign.py`` results files (one goal per line). Goals are
joined on (agent directory, goal id), with a repository worktree prefix
(``.claude/worktrees/agent-<hex>``) normalised, so a census re-run from
another worktree joins the same public goals. For every goal both files
replayed, three classes of re-signed word are compared:

- the settlement or held word (W1, W2), in ``resign.py``'s buckets: a_state,
  b_content, c_wording (state and content equal, text differs), identical;
- the ``qualified`` rows the replayed word wrote (W3);
- the executor's analysis word (W4) and how it was obtained.

Goals present or replayed in one file only are counted and listed. The word
comparison is R11 truth-2's ``compare_trees.py`` (censuses 5 and 8),
unchanged; the W3 and W4 comparisons and the prefix normalisation were added
for truth-4's final-pin census (item 0).
"""

from __future__ import annotations

import collections
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from resign import bucket  # noqa: E402

WORKTREE = re.compile(r"/\.claude/worktrees/agent-[0-9a-f]+/")


def load(path: str) -> dict[tuple[str, str], dict]:
    out = {}
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if line.strip():
            row = json.loads(line)
            agent = WORKTREE.sub(
                "/.claude/worktrees/<worktree>/", row["agent"]
            )
            out[(agent, row["goal_id"])] = row
    return out


def main() -> None:
    base, new = load(sys.argv[1]), load(sys.argv[2])
    listing = "--list" in sys.argv
    counts: dict[str, collections.Counter] = {
        "word": collections.Counter(),
        "qualified": collections.Counter(),
        "executor_word": collections.Counter(),
        "population": collections.Counter(),
    }
    for key in sorted(set(base) | set(new)):
        if key not in base or key not in new:
            counts["population"]["in one file only"] += 1
            if listing:
                print("in one file only:", key[1], key in base, key in new)
            continue
        old_row, new_row = base[key], new[key]
        old, now = old_row.get("replayed"), new_row.get("replayed")
        if not old or not now:
            counts["population"]["not both replayed"] += 1
            if bool(old) != bool(now) and listing:
                print("replayed in one only:", key[1], bool(old), bool(now))
            continue
        counts["population"]["both replayed"] += 1
        verdict = bucket(
            {"state": old.get("state"), "reasons": old.get("reasons")}, now
        )
        if verdict == "identical" and old.get("reasons") != now.get("reasons"):
            verdict = "c_wording"
        counts["word"][verdict] += 1
        if verdict != "identical" and listing:
            print(
                verdict,
                key[0],
                key[1],
                old.get("state"),
                "->",
                now.get("state"),
            )
            for was, is_now in zip(
                old.get("reasons") or [], now.get("reasons") or []
            ):
                if was != is_now:
                    print("   -", was[:600])
                    print("   +", is_now[:600])
        same_rows = (old_row.get("replayed_qualified") or []) == (
            new_row.get("replayed_qualified") or []
        )
        counts["qualified"]["identical" if same_rows else "differs"] += 1
        if not same_rows and listing:
            print(
                "qualified differs:",
                key[1],
                old_row.get("replayed_qualified"),
                "->",
                new_row.get("replayed_qualified"),
            )
        executor = ("executor_word", "executor_word_how")
        same_word = all(old_row.get(k) == new_row.get(k) for k in executor)
        counts["executor_word"]["identical" if same_word else "differs"] += 1
        if not same_word and listing:
            print(
                "executor word differs:",
                key[1],
                [old_row.get(k) for k in executor],
                "->",
                [new_row.get(k) for k in executor],
            )
    for name, counter in counts.items():
        print(f"{name}: {dict(sorted(counter.items()))}")


if __name__ == "__main__":
    main()
