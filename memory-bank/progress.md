# Progress Dashboard (numbers must come from a real audit/export run)

_Last measured: 2026-07-12, from `audit_dataset.py`, `export_sft.py`, and
`build_split_manifest.py`, all run with `--training-dir var/agent-training
--include-runs`. The authoritative export is
`var/agent-training/exports/production_readiness_20260712_codex/`. A
`--write-clean` per-turn file is NOT an authoritative input. Update this table
after every session's export — a stale dashboard silently misleads the next
agent._

_NOTE: the store includes Claude's isolated rare-kind run
`deepseek-r1-rk-20260712`; no raw ledger files were modified by this audit._

## Episode store (audit_dataset.py)

| metric | value | health |
|---|---|---|
| episodes raw / dedup / clean | **3311 / 2853 / 2499** | ok |
| clean after skeleton cap | **2466** | ok |
| distinct query skeletons / ratio | **1465 / 0.5135** | PASS (>=0.35) |
| health_failures | [] | PASS |
| job kinds covered (tag-based) | **25** | see caveat |
| multi-turn session share | **0.4596** | ok |
| trajectory_entropy | **5.2554** | ok |
| distinct_tool_trajectories | **410** | ok |
| semantic_trajectory_duplicates (multi-step) | **787** | watch |
| avg schema_efficiency (tool-bearing, n=2475) | **0.9458** | ok |

## Session chains

| bucket | count |
|---|---|
| trainable_session_chains | **633** |
| review_session_chains | **940** |
| session_count | **1573** |

## Export families (export_sft.py, latest full run)

| family | rows | note |
|---|---|---|
| **reasoning_synthesis (trusted)** | **347** | above the 150-300 target; curate/cap before training |
| reasoning_synthesis_review | **1636** | not trainable |
| tool_loop_sft | **623** | hidden-CoT scan: 0 `reasoning_content` fields |
| tool_loop_review | **950** | includes 1 legacy absolute home-path sample; redact before sharing |
| command_answer | **741** | |
| compact_spec | **703** | 3B-planner lineage (excluded from v17 mix) |
| project_yaml | 462 | |
| repair_pairs | **84** | |
| wrong_route_contrast | 43 | TARGET ~100 (was 14) |
| terminal_state_assertions | **83** | tau-bench-style end-state checks |

Export deterministic; split leakage is 0 session / 0 skeleton crossings, and
`reasoning_content` = 0 across all ten JSONL families. Positive families have
no home-path hits. Review-only data still contains one legacy
`/Users/bitplane/...` path and must be redacted before external distribution.
Residual tool-role `raw_response` is public synthesis JSON
(status/command/explanation/confidence), not hidden CoT; it is redundant bloat
and should be trimmed before a public release.

## Trajectory depth (tool_loop_sft)

avg user turns **2.073**, max 5, single-`synthesize_command`-only 46.5%,
positive tool-loop rows reaching `execute_chemsmart_command` **11**. 51
distinct providers contributed.

## Kind coverage & imbalance

The current audit reports **25 distinct job-kind tags** in the pipeline (the
extra tag is `gaussian.userjob`). This is neither the canonical 24-kind CLI
coverage nor a gated-positive count. Use the exporter-derived rare-kind table
below and a future per-kind gated audit for coverage claims.

**RARE-KIND CAMPAIGN COMPLETE (2026-07-12).** The authoritative positive export
now has >=20 gated-positive `tool_loop_sft` rows for all five rare kinds:
`gaussian.crest` **21**, `gaussian.dias` **20**, `orca.irc` **24**,
`orca.modred` **22**, `orca.qrc` **23**. The isolated run itself contains 112
hand-authored four-turn scenarios across 8 batches: 72 PASS, 34 ASK, and 6
WRONG. The 72 PASS rows are not the total kind coverage; the >=20 figures are
the cumulative gated-positive export counts. The run used
`agentic_workflow_accum.py --provider deepseek --models deepseek-reasoner
--run-label deepseek-r1-rk-20260712 --scenario-set rare4turn`, with no raw
ledger interference from the other collector. QRC wording must distinguish a
QRC from IRC, and CREST prompts must specify an explicit conformer count.

**qmmm correction:** QMMM remains emittable through the chained `opt qmmm`
sub-subcommand, e.g. `chemsmart run orca -p demo -f e.xyz opt qmmm -ha 1-12
-lm AMBER=...`. The earlier "qmmm has no CLI subcommand" claim was wrong.
The current audit sees 26 Gaussian-QMMM and 79 ORCA-QMMM tags, but these are
pipeline counts, not gated-positive counts.

## Quality / stability flags

quality: gate_reject_unrepaired **200**, terminal_nosynth 107, empty_final_answer 18,
missing_terminal_state **32**, execute_failed 6, terminal_state_assertion_failed **5**,
terminal_state_expected_failure 1.
stability: tool_error_unresolved 36, repeated_identical_tool_call 29,
execute_without_synthesize **8**. exact_duplicate_pairs **124**.
over-frequent skeletons: "explain the current project yaml..." 38, qm/mm
clarification pattern 56, "optimize <file>..." 27, "write the project yaml." 21.

## Session-chain skip reasons (why chains aren't trainable)

unresolved_clarification **355**, unresolved_pause **297**, terminal_tool_error **143**,
terminal_no_response 80, missing_generated_input_evidence **50**, canonical_freq_in_route **7**,
terminal_execute_failed 3, missing_terminal_state 3, canonical_empty_qmmm_layer_override 1.

## Teacher status

`gpt-5.4-mini` (OpenAI, paid) — hidden-CoT-free by construction; ~$8.10 of the
approved $9 left. Invoke with `--provider openai` (base_url pinned; never None).
DashScope thinking teachers alive on last scan (re-ping before use):
qwen3-235b-a22b-thinking-2507, qwen3-next-80b-a3b-thinking, qwen3-30b-a3b-thinking,
qwen3.7-plus. 51 providers total have contributed to the store.

## What's left (ordered)

1. ~~Strip hidden CoT from positives~~ **DONE** — recheck on every export.
2. Cap/curate `reasoning_synthesis` from 347 toward the 150-300 target.
3. ~~Rare kinds to >=20 gated positives~~ **DONE** — all five rare kinds pass.
4. Increase `wrong_route_contrast` from 43 toward ~100 mis-route/corrected twins.
5. Expand terminal receipts and long workflows; only 11 positive rows currently
   reach `execute_chemsmart_command`.
6. Fix `canonical_kind_coverage` to report gated positives directly; current
   audit kind counts remain tag/pipeline counts.
7. Redact legacy absolute paths from review exports and trim redundant public
   `raw_response` before external release.
8. Run a final per-kind held-out evaluation before v17 training; raw collection
   is improved but this is not yet a blanket production-readiness claim.

## Known issues

- **The semantic gate checks syntax and runtime, not intent.** `orca opt -f 1,3,2,4`
  gates ok (`-f` is overloaded) while freezing atoms when distances were asked
  (`modred -c`). Gate-ok != intent-correct.
- **`sub` cannot be gated on this machine** (no server config; discovery is not
  workspace-local). Keep collection corpora `run`-framed until a fixture is agreed.
- **Pruner synonym gaps** hide a subcommand so the model can't emit it; before
  blaming a teacher print `schema_variant_id(prune_schema_for_request(...))`.
- Router confusion (job request -> YAML authoring) is the dominant trusted killer;
  command-hard framing mitigates, contrast pairs cure.
