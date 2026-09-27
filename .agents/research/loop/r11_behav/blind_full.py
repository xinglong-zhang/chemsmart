"""Print full probe-window texts for blind codes (scratch; key not printed)."""

import json
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
key = json.loads((S / "blind-key.json").read_text())
for code in sys.argv[1:]:
    arm, sample = key[code]
    for line in (S / "baseline" / arm / "samples.jsonl").read_text().splitlines():
        row = json.loads(line)
        if row["sample"] != sample:
            continue
        print(f"===== {code}")
        for text in row.get("texts_after_cut") or ():
            content = str(text.get("content") or "")
            for token, pattern in mo._FILES.items():
                content = pattern.sub(token, content)
            print(f"--- turn {text['turn']}: {content[:2500]}")
