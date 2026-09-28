"""Arm M reading: R-DISSENT counts by cell (key opened), then the registered HYP rule etc."""

import json
import sys
from collections import Counter, defaultdict

from scipy.stats import fisher_exact

sys.path.insert(
    0,
    "/Users/hongjiseung/developer/chemsmart/.claude/worktrees/"
    "agent-a9e062387a514384b/.agents/research/loop",
)
import matched_outcomes as mo  # noqa: E402

S = "/private/tmp/claude-501/-Users-hongjiseung-developer-chemsmart/c946a1a6-a959-400c-ba03-fd5f4544cb43/scratchpad/behav"
key = json.load(open(f"{S}/dissent/dissent-key.json"))
verdicts = {}
for line in open(f"{S}/dissent/dissent-verdicts.txt"):
    parts = line.split()
    if parts and parts[0].startswith("D") and len(parts) >= 3 and parts[0][1:].isdigit():
        verdicts[parts[0]] = (parts[1], parts[2])
counts = defaultdict(lambda: [Counter(), Counter(), 0])
for code, (cell, _sample) in key.items():
    primary, secondary = verdicts[code]
    counts[cell][0][primary] += 1
    counts[cell][1][secondary] += 1
    counts[cell][2] += 1
for cell in sorted(counts):
    p, s, n = counts[cell]
    print(f"R-DISSENT {cell:18} n {n:2}  primary {dict(p)}  secondary {dict(s)}")


def fisher(a, n1, b, n2):
    return fisher_exact([[a, n1 - a], [b, n2 - b]])[1]


for reading, index in (("primary", 0), ("secondary", 1)):
    q = counts["d2-armM-qwen"][index]["OPPOSE"]
    d = counts["d2-armM-deepseek"][index]["OPPOSE"]
    nq, nd = counts["d2-armM-qwen"][2], counts["d2-armM-deepseek"][2]
    p = fisher(q, nq, d, nd)
    print(f"ARM M OPPOSE ({reading}): qwen {q}/{nq} vs deepseek {d}/{nd}, Fisher p {p:.3g}")

print("--- registered arm M rule (D2 outcomes)")
for model in ("deepseek", "qwen"):
    for root, cell in (("baseline", f"d2-head-{model}"), ("arm-m", f"d2-armM-{model}")):
        tally, n, infra = Counter(), 0, 0
        for line in open(f"{S}/{root}/{cell}/samples.jsonl"):
            outcome = mo.classify_d2(json.loads(line))
            if outcome.get("INFRA"):
                infra += 1
                continue
            n += 1
            for k in ("LOOK", "HYP", "DECLARE", "NOTICE"):
                tally[k] += int(outcome[k])
        print(f"{cell:18} n {n:2} " + " ".join(f"{k} {tally[k]}/{n}" for k in ("LOOK", "HYP", "DECLARE", "NOTICE")) + f" INFRA {infra}")
