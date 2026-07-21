# Safety Invariants (these override every other instruction)

## NEVER — hard prohibitions

1. NEVER run `pip install`, `conda install`, `brew install`, or any package
   manager on this machine. It has broken the user's editable install twice.
   All heavy compute/deps belong on Colab/HPC, not here.
2. NEVER print, log, echo, or commit API keys or tokens. Secrets live ONLY in
   `api.env` (repo root, gitignored) and `~/developer/chemsmart-finetune/v7/
   key.env`. Read them with dotenv; never paste their values anywhere.
3. NEVER modify, rewrite, or delete files under
   `var/agent-training/runs/*/episodes/` — the raw episode ledger is
   APPEND-ONLY evidence. Cleaning happens exclusively via
   `audit_dataset.py --write-clean` into a NEW file.
4. NEVER `git commit`, `git push`, `git stash`, or `git checkout -- <file>`.
   Leave all changes in the working tree; the user handles version control.
5. NEVER edit `.clinerules/**` or `agent-orchestrator.yaml*` unless the user
   explicitly asks in the current conversation.
6. NEVER fabricate training data with code-template loops (for-loop over a
   few sentence skeletons). Proven failure (v13): the model memorizes the
   template and the behavior does not transfer. Only REAL teacher-model runs
   through the agent loop produce training rows.
7. NEVER report a number you did not just measure. Every count in a report
   or in memory-bank/progress.md must come from an actual
   `audit_dataset.py` / `export_sft.py` run in this session.

## API quota discipline (DashScope free tier)

- Run AT MOST 4 collection cases per batch, then stop and summarize before
  the next batch.
- On HTTP 403 containing "free quota": STOP that model immediately (no
  retries), switch to another live teacher, and note it in progress.md.
  Quotas are per-model-id; dated snapshot ids (e.g. `qwen-plus-2025-12-01`)
  carry separate quotas from the base id.
- Qwen3 OPEN-WEIGHT models reject non-streaming calls unless
  `extra_body={"enable_thinking": false}` is sent — but for REASONING
  teachers we WANT thinking: do NOT send that flag to `*-thinking-*` models.
- Before a matrix/batch over a new model id, ping it first: catalog entries
  can 404 on chat (`qwen2-7b-instruct` does).

## Execution safety

- Any `chemsmart` command you execute must be a test run: `run` gets
  `--fake --no-scratch`, `sub` gets `--fake --test`. Never launch real
  Gaussian/ORCA compute from a collection session.
- Workspaces for collection are throwaway temp dirs with exactly ONE project
  YAML (two programs' YAMLs make the yaml_check gate refuse everything).
- Do not write anything into the user's real `~/.chemsmart/` config.
