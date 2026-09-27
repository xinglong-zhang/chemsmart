"""The pre-registered baseline analysis (run only after the runner says DONE).

Counts per cell from matched_outcomes' rules; two-sided Fisher exact for
D1-M (repaired vs control, per model), D1-X (deepseek vs qwen, per tree)
and D2-X (deepseek vs qwen on LOOK, HYP, NOTICE). INFRA rows are counted
apart and never enter a denominator.
"""

import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

from scipy.stats import fisher_exact

sys.path.insert(
    0,
    "/Users/hongjiseung/developer/chemsmart/.claude/worktrees/"
    "agent-a9e062387a514384b/.agents/research/loop",
)
import matched_outcomes as mo  # noqa: E402

BASE = Path(
    "/private/tmp/claude-501/-Users-hongjiseung-developer-chemsmart/"
    "c946a1a6-a959-400c-ba03-fd5f4544cb43/scratchpad/behav/baseline"
)


def rows(arm):
    path = BASE / arm / "samples.jsonl"
    return [json.loads(line) for line in path.read_text().splitlines()] if path.exists() else []


def fisher(a, n1, b, n2):
    _odds, p = fisher_exact([[a, n1 - a], [b, n2 - b]])
    return p


d1 = {}
for tree in ("control", "repaired"):
    for model in ("deepseek", "qwen"):
        arm = f"d1-{tree}-{model}"
        labels = Counter(mo.classify_d1(row) for row in rows(arm))
        infra = labels.pop("INFRA", 0)
        n = sum(labels.values())
        d1[(tree, model)] = (labels["ROUTE"], n)
        print(f"D1 {tree:8} {model:8} ROUTE {labels['ROUTE']}/{n}  classes {dict(labels)}  INFRA {infra}")
for model in ("deepseek", "qwen"):
    (a, n1), (b, n2) = d1[("repaired", model)], d1[("control", model)]
    print(f"D1-M {model}: repaired {a}/{n1} vs control {b}/{n2}, diff {a - b}, Fisher p {fisher(a, n1, b, n2):.3g}")
for tree in ("control", "repaired"):
    (a, n1), (b, n2) = d1[(tree, "deepseek")], d1[(tree, "qwen")]
    print(f"D1-X {tree}: deepseek {a}/{n1} vs qwen {b}/{n2}, Fisher p {fisher(a, n1, b, n2):.3g}")

d2 = defaultdict(Counter)
counted = Counter()
for model in ("deepseek", "qwen"):
    for row in rows(f"d2-head-{model}"):
        outcome = mo.classify_d2(row)
        if outcome.get("INFRA"):
            d2[model]["INFRA"] += 1
            continue
        counted[model] += 1
        for key in ("LOOK", "HYP", "DECLARE", "NOTICE"):
            d2[model][key] += int(outcome[key])
for model in ("deepseek", "qwen"):
    n = counted[model]
    print(f"D2 {model:8} n {n}  " + "  ".join(f"{k} {d2[model][k]}/{n}" for k in ("LOOK", "HYP", "DECLARE", "NOTICE")) + f"  INFRA {d2[model]['INFRA']}")
for key in ("LOOK", "HYP", "NOTICE"):
    a, n1, b, n2 = d2["deepseek"][key], counted["deepseek"], d2["qwen"][key], counted["qwen"]
    print(f"D2-X {key}: deepseek {a}/{n1} vs qwen {b}/{n2}, Fisher p {fisher(a, n1, b, n2):.3g}")
