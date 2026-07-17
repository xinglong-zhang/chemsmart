# CLAUDE.md — ChemSmart agent SFT corpus

Working context for Claude Code in this repo. Cline reads `.clinerules/` +
`memory-bank/`; this file is the equivalent entry point for Claude.

**Read first:** `memory-bank/progress.md` (live numbers), then
`.clinerules/01-safety-invariants.md` and `.clinerules/04-cli-ground-truth.md`.

## Mission

Build a production-grade SFT corpus and fine-tuned local model for the ChemSmart
agent (Gaussian/ORCA CLI synthesis + HPC submission). Teachers run the real agent
loop; the deterministic semantic gate acts as rejection sampling. Only
gate-positive trajectories become trusted training rows.

## Hard rules (violating these has broken the machine or the data before)

- **Never** `pip install` / `conda install` / `brew install` on this Mac. It has
  broken the editable install twice. Use Colab or HPC.
- **Never** print, log, echo, or commit API keys. Secrets live only in `api.env`
  (gitignored) and `~/developer/chemsmart-finetune/v7/key.env`.
- **Never** modify or delete anything under `var/agent-training/runs/*/episodes/`.
  It is append-only evidence. Clean only via a new file, never in place.
- **Never** `git commit`, `git push`, `git stash`, or `git checkout -- <file>`
  unless explicitly asked. Push only to the fork (`Hongjiseung-ROK`), never to
  `origin` (`xinglong-zhang`). Commit messages must not mention "claude".
- **Never** edit `.clinerules/**` or `agent-orchestrator.yaml*` unless asked.
- **Never** fabricate training data with code-template loops. Diversity comes from
  hand-authored phrasings; a templated corpus gates 100% clean and transfers zero
  behavior (v13, measured).
- **Never** invent, request, or store hidden chain-of-thought. `synthesis.reasoning`
  must be an explicitly surfaced, user-auditable trace.
- **Never** report a number you did not just measure in the current session.
- Any executed `chemsmart` command must be a test run: `run` → `--fake --no-scratch`,
  `sub` → `--fake --test`. Never launch real Gaussian/ORCA compute from collection.
- Collection workspaces are throwaway temp dirs with exactly one project YAML.
  Do not write into the user's real `~/.chemsmart/`.

## Authoritative commands

Always export/audit from the append-only store. A `--write-clean` per-turn file is
**not** an authoritative input.

```bash
python scripts/training/audit_dataset.py --training-dir var/agent-training --include-runs
python scripts/training/export_sft.py   --training-dir var/agent-training --include-runs --out-dir <tmp>
```

Collection (`--corpus command-hard` framing defeats the router confusion that
sends plain job requests into project-YAML authoring):

```bash
# DashScope teachers (free quota; stop a model immediately on HTTP 403)
python scripts/training/reasoning_accum.py <model> --index 0 --stride 1 --limit 4 --corpus command-hard

# OpenAI teacher — pinned to api.openai.com, metered, budget-capped
python scripts/training/reasoning_accum.py gpt-5.4-mini --provider openai \
  --corpus gap-fill-wide --limit 21 --budget-usd 2.00
```

## State as of 2026-07-10 (measured this session)

| metric | value |
|---|---|
| episodes raw / dedup | 1,447 / 1,406 |
| session chains: trainable / review | 395 / 571 |
| trusted public-reasoning rows | 41 (target 150–300) |
| skeleton diversity ratio | 0.409 (gate ≥0.35) |
| schema efficiency | 0.9043 |
| kinds missing positives | 2 (`gaussian.qmmm`, `orca.qmmm` — both impossible) |
| health_failures | [] |

`gpt-5.4-mini` costs ~$0.008/case (~$0.033 when the repair loop fires and falls
back to the full schema). ~$8.10 of the approved $9 remains.

## Current tasks (ordered)

1. **Strip hidden CoT from positive export families.** `tool_loop_sft.jsonl` still
   carries **846** assistant `reasoning_content` fields (789,199 chars) and **358**
   tool messages with `raw_response` (1,550,986 chars). Both violate the no-hidden-CoT
   policy and must be removed before any fine-tune. Add a regression test.
2. **Grow trusted reasoning rows 41 → 150–300.** Hand-author more distinct phrasings;
   do not oversample. `gpt-5.4-mini` is the safest teacher: it emits no
   `reasoning_content`, so its rows are hidden-CoT-free by construction.
3. **Lift the 8 rare kinds to ≥20 rows** (`gaussian.crest`, `gaussian.dias`,
   `gaussian.nci`, `orca.irc`, `orca.modred`, `orca.scan`, `orca.sp`, `orca.ts`).
4. **Fix the audit's kind taxonomy.** `COMMAND_CANONICAL_AGENT_KINDS` lists 24
   CLI-emittable kinds; the real count is **22**. `qmmm` has no CLI subcommand for
   either program and can never gate-pass — move it out of `command_emittable`.
   (`gaussian.tddft` is a naming mismatch only; the subcommand is `td`.)
5. **Decide the `sub` fixture.** Submit-framed requests cannot be gated here: the
   gate's `sub --test --fake` fails with "Could not find any submitters for scheduler
   type: None", and server discovery is not workspace-local. Needs a `SLURM.yaml` in
   the user's real `~/.chemsmart/server/` — **ask before writing it.**
6. **Build a fixed held-out agentic benchmark.** Separate fake-semantic-pass from
   real runner end-state; measure `pass^k` before claiming production reliability.
7. **Grow long workflows** (`project YAML → command → dry-run → repair/execute`).
   Only 5 positive chains reach `execute_chemsmart_command`; avg 1.55 user turns.

## Traps that have each cost a failed run

- `OpenAIProvider(base_url=None)` silently falls back to a **third-party gateway**
  (`providers.py` does `base_url or self.gateway_url`). Always pin an explicit URL.
- **The semantic gate checks syntax and runtime, not intent.** `orca opt -f 1,3,2,4`
  gates `ok` (`-f` is overloaded: `--filename` at program level, `--freeze-atoms` on
  `opt`) while freezing atoms when the request wanted distances frozen (`modred -c`).
  Gate-ok ≠ correct.
- **A pruner synonym gap makes the teacher look wrong when it is not.** If
  `schema_prune._JOBKIND_PATTERNS` can't map the request, the subcommand is absent
  from the pruned schema and the model *cannot* emit it. Before blaming a teacher for
  a wrong jobkind, print `schema_variant_id(prune_schema_for_request(...))`.
- Qwen3 thinking models need `enable_thinking=False` when served, or they ramble in
  `<think>` and emit zero tool calls.
