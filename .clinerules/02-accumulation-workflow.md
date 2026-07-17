# Accumulation Workflow (collection -> audit -> session export -> report)

Python interpreter: `/opt/anaconda3/envs/chemsmart/bin/python`. Run all
commands from the repo root `/Users/hongjiseung/developer/chemsmart`.

## 1. Teacher scan (start of session, or after any 403)

Only teachers used for the optional reasoning curriculum need an explicitly
returned `reasoning_content` field. Keep that trace isolated from tool-loop
execution truth and never request or infer hidden chain-of-thought. Known-good
reasoning teachers (verify liveness with a 1-token ping before batches):
- `qwen3-235b-a22b-thinking-2507` (strongest reasoning)
- `qwen3-next-80b-a3b-thinking`   (fast MoE)
- `qwen3-30b-a3b-thinking-2507`   (fastest; weakest routing)
- `qwen3.7-plus`                  (frontier; best routing, also reasons)
Endpoint: `https://dashscope-intl.aliyuncs.com/compatible-mode/v1`, key from
`api.env` field `alibaba-cloud-api-key` (via dotenv, never printed).
If all four are quota-dead, scan dated snapshot ids (separate quotas) or
stop collection and do audit/export/documentation work instead.

## 2. Collection batch (max 4 cases, then summarize)

```
/opt/anaconda3/envs/chemsmart/bin/python scripts/training/reasoning_accum.py \
    <model-id> --index <i> --stride <n> --limit 4 --corpus command-hard
```
- `--corpus command-hard` is the DEFAULT choice: prompts explicitly say
  "use the loaded project settings and synthesize ONLY the CLI command",
  which bypasses the known router confusion (job requests mis-routed into
  project-YAML authoring). Measured: command-hard 3/4 PASS vs 0/4 before.
- `--corpus workflow-hard` / `agentic-style`: use sparingly for trajectory
  diversity once trusted volume is healthy.
- Episodes append automatically to `var/agent-training/runs/<model>/`.
- A PASS needs: synthesis status `ready`, semantic gate `ok|warn`, and the
  expected job subcommand in the command. ASK (needs_clarification) is a
  safe non-count; investigate repeated WRONG/NOSYNTH per 03-quality-gates.

## 3. Audit (after every 2-3 batches)

```
/opt/anaconda3/envs/chemsmart/bin/python scripts/training/audit_dataset.py \
    --training-dir var/agent-training --include-runs \
    --strict --out var/agent-training/audit_<id>.json
```
- `--strict` must exit 0. If skeleton ratio drops below 0.35, STOP
  collecting repeats and vary request phrasing/kinds instead.
- Read `stability.multi_turn_session_share`, not the legacy per-record
  history count. Also inspect `quality.trainable_session_chains`,
  `quality.review_session_chains`, and `session_chain_skip_reasons`.
- `terminal_nosynth` is untrainable. A cross-turn reject that is later
  repaired is reported separately as `recovered_cross_turn_rejects`.

## 4. Export (fresh directory per session; never overwrite older exports)

```
/opt/anaconda3/envs/chemsmart/bin/python scripts/training/export_sft.py \
    --training-dir var/agent-training --include-runs \
    --out-dir var/agent-training/exports/v17_<date>_<label>
```
- `tool_loop_sft.jsonl` is reconstructed by `session_id`, not exported one
  turn at a time. It preserves correction and pause/resume context while
  removing repeated history prefixes.
- Failed or abandoned chains belong in `tool_loop_review.jsonl`; recovered
  reject-to-repair examples also belong in `repair_pairs.jsonl`.
- Do not pass a per-turn `--write-clean` output when the goal is an
  authoritative session-chain export. It can discard the rejected half of a
  later cross-turn repair.
- Never use `--include-rejected` for a production positive export.
- Read `manifest.json`: inspect all `written` and `skipped` counts plus
  `session_chains_seen`, `multi_turn_chains`, and
  `multi_turn_chain_share`. Reasoning counts are a separate curriculum, not
  the headline tool-loop count.
- Do not collect hidden chain-of-thought. Only explicitly surfaced reasoning
  may enter the isolated reasoning family after its own gates.

## 5. Contrast growth

When a batch produces a mis-routed episode (job request answered with
project-YAML authoring), pair it 1:1 with a corrected run of the SAME
request via command-hard framing. The exporter mints the contrast pair.
Never author these pairs by hand-editing episode files.

For ordinary command repair, keep the reject in the append-only ledger and
perform the corrected request in the same session. The exporter must emit the
successful full chain to `tool_loop_sft` and a separate reject-to-repair row
to `repair_pairs`; it must never promote the rejected command itself.

## 6. Close the loop

Update `memory-bank/progress.md` with the new audit + manifest numbers and
`memory-bank/activeContext.md` if focus changed. Report: batches run,
PASS/ASK/WRONG per teacher, positive/review chain counts, repair-pair delta,
multi-turn session share, coverage gaps, and quota status.
