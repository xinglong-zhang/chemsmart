"""Print D2 samples' label-related text windows under random codes (scratch).

The key (code -> arm, sample) goes to blind-key.json, opened only after the
hand-read verdicts are written to blind-verdicts.txt. Windows: 400
characters around every mention of a label word, a masked file name,
IUPAC numbering or a swap word, in the assistant texts and call arguments
after the cut (the probe window).
"""

import json
import random
import re
import sys
from pathlib import Path

sys.path.insert(
    0,
    "/Users/hongjiseung/developer/chemsmart/.claude/worktrees/"
    "agent-a9e062387a514384b/.agents/research/loop",
)
import matched_outcomes as mo  # noqa: E402

S = Path(
    "/private/tmp/claude-501/-Users-hongjiseung-developer-chemsmart/"
    "c946a1a6-a959-400c-ba03-fd5f4544cb43/scratchpad/behav"
)
PATTERN = re.compile(
    r"label|FILE_C[45]|IUPAC|numbering|swapp?ed|transpos|mis-?label|"
    r"wrong way|revers|invert",
    re.IGNORECASE,
)
items = []
for model in ("deepseek", "qwen"):
    path = S / "baseline" / f"d2-head-{model}" / "samples.jsonl"
    for line in path.read_text().splitlines():
        row = json.loads(line)
        if row.get("infra"):
            continue
        texts = [str(t.get("content") or "") for t in row.get("texts_after_cut") or ()]
        texts += [str(a.get("arguments") or "") for a in row.get("acts_after_cut") or ()]
        windows = []
        for text in texts:
            for token, pattern in mo._FILES.items():
                text = pattern.sub(token, text)
            last = -1000
            for match in PATTERN.finditer(text):
                if match.start() < last + 400:
                    continue
                last = match.start()
                start = max(0, match.start() - 200)
                windows.append(text[start : start + 400].replace("\n", " "))
        items.append((f"d2-head-{model}", row["sample"], windows))
random.Random(20260928).shuffle(items)
key = {}
for number, (arm, sample, windows) in enumerate(items):
    code = f"S{number:02d}"
    key[code] = [arm, sample]
    print(f"===== {code} ({len(windows)} windows)")
    for window in windows[:8]:
        print("  -", window)
(S / "blind-key.json").write_text(json.dumps(key, indent=1))
