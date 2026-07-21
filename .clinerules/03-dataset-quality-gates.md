# Dataset Quality Gates (what makes a row trainable)

## Trusted `tool_loop_sft` session chain - ALL must hold

1. Group deduplicated turn snapshots by `session_id`, sort by `turn`, and
   merge repeated history prefixes.
2. The terminal turn is not paused, denied, unresolved clarification, failed
   execution, failed tool result, or NOSYNTH/no-response.
3. A command workflow has a canonical `chemsmart ` command, semantic verdict
   `ok|warn`, and generated-input evidence.
4. The final assistant answer is non-empty and the message/tool sequence has
   no orphan tool result.
5. No secret or unmasked private path is present.

Anything failing these rules lands in `tool_loop_review.jsonl`. A later
successful correction can make the full chain positive and produce a
separate `repair_pairs.jsonl` row, but the rejected command is never positive
SFT by itself.

## Trusted `reasoning_synthesis` row - ALL must hold

1. Semantic gate verdict `ok` or `warn` (deterministic runtime gate, not
   your judgment).
2. `synthesis.reasoning` is an explicitly surfaced, user-auditable teacher
   trace. Never invent, request, or store hidden chain-of-thought.
3. Command starts with `chemsmart ` and violates no canonical rule
   (see 04-cli-ground-truth.md).
4. A real user request exists in the episode messages.
Anything failing these lands in `reasoning_synthesis_review.jsonl` — that is
correct behavior, not a bug to work around.

## Quality flags (audit_dataset.py) — untrainable episodes

`gate_reject_unrepaired` · `execute_failed` (rc != 0) · `user_denied` ·
`malformed_command` · `missing_user_message` · `empty_final_answer`
(tools ran but the agent never answered - trains silence) ·
`terminal_nosynth` (no useful terminal tool trajectory or response) ·
`secret_detected`. A reject with a same-turn or cross-turn positive repair is
reported as recovered and feeds `repair_pairs` contrast data.

Session-chain review reasons include `unresolved_pause`,
`unresolved_clarification`, `terminal_tool_error`, `terminal_no_response`,
`terminal_execute_failed`, `terminal_semantic_reject`,
`missing_generated_input_evidence`, and canonical route violations.

## Diversity discipline (the most expensive lesson of this project)

- v13 proof: 999 template-generated rows gated 100% clean, trained with
  healthy loss — and transferred ZERO behavior on held-out phrasing.
  Surface diversity is the lever, not volume.
- Gate: distinct query-skeleton ratio >= 0.35 (skeleton = lowercase,
  filenames→`<f>`, digits→`#`). Cap: max 20 rows or 5% share per skeleton.
- To raise a minority behavior, author MORE DISTINCT phrasings (structure,
  register, verbosity, EN/KO/ZH), never duplicate/oversample.
- Missing-slot "decline" style requests must contain NO atom/index cue
  words, or the harness treats the field as present.

## SOTA metrics — how to read them

- `multi_turn_session_share`: sessions with at least two user messages /
  deduplicated sessions. This is the headline multi-turn metric.
- `episode_share_in_multi_turn_sessions`: how much raw turn volume belongs to
  those sessions. Do not confuse it with session share.
- `record_level_multi_turn_episodes`: legacy evidence that one stored record
  already contained multiple user messages; never use it as the headline.
- `trainable_session_chains` / `review_session_chains`: strict terminal-gate
  split. Always report the dominant `session_chain_skip_reasons`.
- `trajectory_entropy` (Shannon over tool-transition bigrams): higher =
  more varied multi-step behavior. Current baseline ≈ 4.4.
- `distinct_tool_trajectories`: sha1 of the ORDERED tool sequence
  (deterministic). Baseline ≈ 181.
- `semantic_trajectory_duplicates`: counts duplicate MULTI-STEP (>=2 tools)
  trajectories only. Single-tool synthesize episodes are the expected norm
  and are NOT flooding — do not chase this number down for them.
- `average_schema_efficiency_score`: error-free fraction of tool calls,
  averaged over TOOL-BEARING episodes only (no-tool chat is excluded, not
  scored 0). Healthy baseline ≈ 0.90; a drop means teachers are emitting
  malformed tool calls — check the provider or the prompt.

## Report integrity

Numbers in any report must be paste-able from the tool output of THIS
session. If audit and export disagree (e.g. dedup counts), say so and show
both — never smooth over discrepancies.

Run authoritative exports directly from the append-only store with
`--training-dir var/agent-training --include-runs`. A per-turn
`--write-clean` file is not an authoritative input for session-chain,
pause/resume, or cross-turn-repair exports.
