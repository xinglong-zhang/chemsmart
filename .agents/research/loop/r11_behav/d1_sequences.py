"""Exploratory: every D1 sample's call sequence after the cut, and whether a
ROUTE-class change appears anywhere within the probe turns (scratch)."""

import json
import sys

sys.path.insert(
    0,
    "/Users/hongjiseung/developer/chemsmart/.claude/worktrees/"
    "agent-a9e062387a514384b/.agents/research/loop",
)
import matched_outcomes as mo  # noqa: E402

S = "/private/tmp/claude-501/-Users-hongjiseung-developer-chemsmart/c946a1a6-a959-400c-ba03-fd5f4544cb43/scratchpad/behav/baseline"
for arm in ("d1-control-qwen", "d1-repaired-qwen", "d1-control-deepseek", "d1-repaired-deepseek"):
    eventually = 0
    for line in open(f"{S}/{arm}/samples.jsonl"):
        row = json.loads(line)
        acts = row.get("acts_after_cut") or []
        first = mo.classify_d1(row)
        routes_later = any(
            mo.classify_d1({"acts_after_cut": acts[i:]}) == "ROUTE" for i in range(len(acts))
        )
        eventually += routes_later
        seq = [f"{a['turn']}:{a['name']}({a['status'][0]})" for a in acts][:9]
        if arm.endswith("qwen"):
            print(arm, row["sample"], first, "route-within-probe" if routes_later else "-", " ".join(seq))
    print("ARM", arm, "route anywhere within the probe turns:", eventually, "/ 12")
