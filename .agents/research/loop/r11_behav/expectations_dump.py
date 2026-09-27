"""Exploratory: declared expectations that bear on the regiochemical direction (scratch).

For each D2 sample, every declared observable (any role) whose meaning or
basis mentions a regioisomer direction, printed with its sign, band and
basis, for a hand classification: agrees with the chemist's expectation
(ester at C4 major), opposes it (ester at C5 / CF3-controlled major), or
no direction.
"""

import json
import re
import sys

S = "/private/tmp/claude-501/-Users-hongjiseung-developer-chemsmart/c946a1a6-a959-400c-ba03-fd5f4544cb43/scratchpad/behav"
root = sys.argv[1]
model = sys.argv[2]
DIRECTION = re.compile(r"regio|C4|C5|c4|c5|major|favou?r", re.IGNORECASE)
for line in open(f"{S}/{root}/samples.jsonl"):
    row = json.loads(line)
    seen = set()
    for act in row.get("acts_after_cut") or ():
        if act["name"] != "declare_requested_observable" or act["status"][0] != "ok":
            continue
        try:
            args = json.loads(act["arguments"])
        except ValueError:
            continue
        for item in args.get("observables") or ():
            if not item.get("expected_sign") and item.get("expected_low") is None:
                continue
            text = f"{item.get('meaning')} {item.get('expectation_basis')}"
            key = item.get("observable_id")
            if key in seen or not DIRECTION.search(text):
                continue
            seen.add(key)
            print(f"--- {model} s{row['sample']:02d} {key} role={item.get('role', 'requested')} sign={item.get('expected_sign')} band=[{item.get('expected_low')},{item.get('expected_high')}]")
            print("   meaning:", str(item.get("meaning"))[:260])
            print("   basis:", str(item.get("expectation_basis"))[:420])
