"""Exploratory: the diagnostics each D2 sample declared (scratch, post hoc)."""

import json
import sys

S = "/private/tmp/claude-501/-Users-hongjiseung-developer-chemsmart/c946a1a6-a959-400c-ba03-fd5f4544cb43/scratchpad/behav/baseline"
model = sys.argv[1]
for line in open(f"{S}/d2-head-{model}/samples.jsonl"):
    row = json.loads(line)
    for act in row.get("acts_after_cut") or ():
        if act["name"] != "declare_requested_observable" or act["status"][0] != "ok":
            continue
        try:
            args = json.loads(act["arguments"])
        except ValueError:
            continue
        for item in args.get("observables") or ():
            if item.get("role") != "diagnostic":
                continue
            print(f"--- s{row['sample']:02d} {item.get('observable_id')} sign={item.get('expected_sign')} band=[{item.get('expected_low')},{item.get('expected_high')}]")
            print("   meaning:", str(item.get("meaning"))[:300])
            print("   basis:", str(item.get("expectation_basis"))[:300])
            print("   if falsified:", str(item.get("failure_update_rule"))[:220])
