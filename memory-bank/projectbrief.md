# Project Brief — chemsmart agent v17 reasoning SFT

## One paragraph

chemsmart is a computational-chemistry CLI (Gaussian/ORCA job synthesis and
HPC submission) with an AI agent layer. We are building the SFT corpus for
**v17**: a local model that reasons (`<think>…</think>`) and then DIRECTLY
emits a valid `chemsmart` command. Data comes from REAL teacher-model runs
(DashScope reasoning models) through the actual unified agent loop; the
deterministic semantic gate acts as rejection sampling. This corpus feeds a
Qwen3-14B cloud validation train first, then a Qwen3-4B local (8GB MLX)
distillation.

## Why v17 (history in two lines)

- v16 (Qwen3-14B, tool_loop-only SFT): great eval_loss (0.862) but 6/26 on
  behavioral eval — it learned to CALL synthesize_command, never to DO the
  synthesis (that skill lived in tool-result messages, loss-masked).
- Fix: capture the teacher's synthesis reasoning + gate-approved command as
  the assistant target (`reasoning_synthesis` family) and drop delegation.

## Data flow (all paths repo-relative)

```
teacher run (scripts/training/reasoning_accum.py, DashScope)
  └─ episodes append → var/agent-training/runs/<model>/episodes/*.jsonl
       (APPEND-ONLY raw ledger; reasoning captured in synthesis.reasoning
        AND assistant reasoning_content)
  └─ scripts/training/audit_dataset.py  (quality/diversity/stability gates)
  └─ scripts/training/export_sft.py     (families below; fresh out-dir)
       reasoning_synthesis_sft.jsonl   ← THE v17 target family (trusted)
       reasoning_synthesis_review.jsonl← non-trainable candidates
       wrong_route_contrast.jsonl      ← router-boundary pairs
       tool_loop / compact_spec / command_answer / project_yaml / repair
```

## Key institutional knowledge

- Router confusion is the dominant failure: plain job requests get
  mis-routed into project-YAML authoring → no synthesis → gate-not-positive.
  The `--corpus command-hard` framing bypasses it (measured 3/4 vs 0/4).
- Diversity beats volume: template-generated data gates clean but transfers
  nothing (v13 proof). Skeleton ratio must stay >= 0.35.
- Only `synthesis.reasoning` qualifies a row as trusted; assistant-message
  reasoning is presentation logic (review-only).
- Quotas are per DashScope model id; dated snapshots have separate quotas.

## People/agents

The user (jiseung) directs; Claude/Codex/Cline agents operate under
`.clinerules/`. Any agent must keep this memory-bank current for the next.
