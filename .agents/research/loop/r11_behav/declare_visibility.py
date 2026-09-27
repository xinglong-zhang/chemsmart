"""Was the diagnostic sentence readable when each D2 declaration was composed? (scratch)

Per baseline D2 sample: the turn of the first declare_requested_observable
call; whether that tool was already callable (in view, its description --
which carries declare.diagnostic_has_standing -- readable) in the exposure
record in force at that turn; the first call's reply status; whether a
re-issue after schema_loaded carried byte-identical arguments; whether any
declare call carried role "diagnostic"; and search_capabilities calls.
"""

import json
from collections import Counter

S = "/private/tmp/claude-501/-Users-hongjiseung-developer-chemsmart/c946a1a6-a959-400c-ba03-fd5f4544cb43/scratchpad/behav/baseline"
TOOL = "declare_requested_observable"
for model in ("deepseek", "qwen"):
    tally = Counter()
    for line in open(f"{S}/d2-head-{model}/samples.jsonl"):
        row = json.loads(line)
        acts = row.get("acts_after_cut") or []
        declares = [a for a in acts if a["name"] == TOOL]
        tally["searches"] += sum(1 for a in acts if a["name"] == "search_capabilities")
        if not declares:
            tally["no_declare_call"] += 1
            continue
        first = declares[0]
        # exposure in force at the first declare's turn: last view recorded at
        # or before (turn index counts provider turns observed before it).
        views = [v for v in row.get("views") or () if v["turn_index"] <= first["turn"] - 1]
        in_view = bool(views) and TOOL in (views[-1].get("callable") or [])
        tally["in_view_when_first_composed"] += int(in_view)
        tally[f"first_status_{first['status'][0]}"] += 1
        if first["status"][0] == "schema_loaded" and len(declares) > 1:
            tally["reissued"] += 1
            tally["reissued_byte_identical"] += int(declares[1]["arguments"] == first["arguments"])
        tally["first_attempt_has_diagnostic"] += int('"diagnostic"' in first["arguments"])
        tally["any_declare_has_diagnostic"] += int(any('"diagnostic"' in d["arguments"] for d in declares))
    print(model, dict(sorted(tally.items())))
