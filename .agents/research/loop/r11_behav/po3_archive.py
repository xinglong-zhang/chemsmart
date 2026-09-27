"""Score archived po3 sessions with the pre-registered D2 rules (scratch).

For every live session under the given roots whose workspace holds the
po3 product geometry (e57f5c9d...), read its public transcript and stream:
model (most common observed), first day, whether it is a goal's cycle-1
session (no previous_run in its goal context), and LOOK / HYP / NOTICE per
matched_outcomes, over the whole session and over its first 4 assistant
turns.
"""

import json
import os
import sys
from collections import Counter

sys.path.insert(
    0,
    "/Users/hongjiseung/developer/chemsmart/.claude/worktrees/"
    "agent-a9e062387a514384b/.agents/research/loop",
)
import matched_outcomes as mo  # noqa: E402


def rows_of(transcript, window=None):
    d = json.load(open(transcript))
    messages = d["transcript"] if isinstance(d, dict) else d
    replies = {m.get("tool_call_id"): str(m.get("content") or "") for m in messages if m.get("role") == "tool"}
    assistants = [m for m in messages if m.get("role") == "assistant"]
    if window:
        assistants = assistants[:window]
    acts, texts = [], []
    for n, turn in enumerate(assistants, 1):
        if turn.get("content"):
            texts.append({"turn": n, "content": str(turn["content"])})
        for call in turn.get("tool_calls") or ():
            f = call.get("function") or {}
            reply = replies.get(call.get("id"), "")
            try:
                parsed = json.loads(reply)
                status = [str(parsed.get("status") or ""), ""]
            except (ValueError, AttributeError):
                status = ["", ""]
            acts.append({"turn": n, "name": f.get("name"), "arguments": str(f.get("arguments") or ""), "status": status})
    context = ""
    for m in messages:
        if m.get("role") == "user":
            context = str(m.get("content") or "")
            break
    return {"acts_after_cut": acts, "texts_after_cut": texts}, context


def model_of(stream):
    models = Counter()
    first = ""
    for line in open(stream, encoding="utf-8", errors="replace"):
        try:
            e = json.loads(line)
        except ValueError:
            continue
        first = first or str(e.get("timestamp") or "")[:10]
        if e.get("kind") == "provider_turn_observed":
            m = str((e.get("payload") or {}).get("observed_model") or "")
            if m and m != "not_observed":
                models[m] += 1
    return (models.most_common(1)[0][0] if models else ""), first


for root in sys.argv[1:]:
    for base, _dirs, files in os.walk(root):
        if not base.split("/")[-1].startswith("live-") or "events.jsonl" not in files:
            continue
        workspace = base.split("/.chemsmart-agent/")[0]
        if not os.path.exists(os.path.join(workspace, "triazole-ester-at-c4.xyz")):
            continue
        transcripts = sorted(f for f in files if f.startswith("public-transcript-"))
        if not transcripts:
            continue
        model, day = model_of(os.path.join(base, "events.jsonl"))
        if not model:
            continue
        path = os.path.join(base, transcripts[-1])
        whole, context = rows_of(path)
        first4, _ = rows_of(path, window=4)
        cycle1 = '"previous_run":""' in context.replace(" ", "") or '"previous_run"' not in context
        a, b = mo.classify_d2(whole), mo.classify_d2(first4)
        label = os.path.relpath(workspace, root)
        print(
            f"{model.split('-')[0]:9} {day} c1={int(cycle1)} "
            f"whole L{int(a['LOOK'])} H{int(a['HYP'])} N{int(a['NOTICE'])} | "
            f"first4 L{int(b['LOOK'])} H{int(b['HYP'])} N{int(b['NOTICE'])} | {label}"
        )
