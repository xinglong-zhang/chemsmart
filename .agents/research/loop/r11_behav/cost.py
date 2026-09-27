"""Sum provider token usage of baseline rows; prints nothing but counts (scratch)."""

import glob
import json
from collections import Counter

S = "/private/tmp/claude-501/-Users-hongjiseung-developer-chemsmart/c946a1a6-a959-400c-ba03-fd5f4544cb43/scratchpad/behav/baseline"
prompt, completion, samples = Counter(), Counter(), Counter()
for path in glob.glob(S + "/*/samples.jsonl"):
    for line in open(path):
        row = json.loads(line)
        arm = row["arm"]
        samples[arm] += 1
        for probe in row.get("probe") or ():
            usage = probe.get("usage") or {}
            prompt[arm] += int(usage.get("prompt_tokens") or 0)
            completion[arm] += int(usage.get("completion_tokens") or 0)
for arm in sorted(samples):
    print(f"{arm:26} samples {samples[arm]:2}  input {prompt[arm]:>9,}  output {completion[arm]:>8,}")
print(f"TOTAL input {sum(prompt.values()):,}  output {sum(completion.values()):,}")
