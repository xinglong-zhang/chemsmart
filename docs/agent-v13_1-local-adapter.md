# ChemSmart Agent v13.1 Local Adapter

This document records the dataset state, validation evidence, and runtime
adapter contract used to deploy the v13.1 local chemsmart command-synthesis
model.

## Model

- Deployable merged model:
  `Smilesjs/chemsmart-qwen2.5-coder-3b-instruct-v13_1`
- LoRA/checkpoint repos kept for audit:
  `Smilesjs/chemsmart-qwen2.5-coder-3b-instruct-lora-v13_1` and
  `Smilesjs/chemsmart-qwen2.5-coder-3b-instruct-lora-v13_1-mixed-checkpoints`
- Training base noted in the v13.1 report:
  `unsloth/Qwen2.5-Coder-3B-Instruct-bnb-4bit`
- Final recipe: LoRA, learning rate `5e-6`, 3 epochs.

## Dataset Layout

Source folder:
`/Users/hongjiseung/developer/chemsmart-finetune/v13_1`

Top-level files:

- `Result_v13_1.md`: training and eval report.
- `config/`: training and evaluation config.
- `scripts/`: generation, semantic gate, Colab train/eval/merge scripts.
- `reports/build_v13_1_summary.json`: dataset build summary.
- `artifacts/train/`: training logs and checkpoints.
- `artifacts/eval/`: clean eval and checkpoint-eval artifacts.

Data files:

- `data/train-00000-of-00002.jsonl`: 3,415 records.
- `data/train-00001-of-00002.jsonl`: 3,415 records.
- `data/eval_representative_v131.jsonl`: 733 records.
- `data/eval_boundary_v131.jsonl`: 179 records.
- `data/eval_regression_ts_route_v131.jsonl`: 70 records.

Combined checked records: 7,812.

Build summary:

- training records: 6,830.
- base records: 5,369.
- authored augment records kept: 1,461 from 1,480 raw authored records.
- intents: 5,458 workflow, 912 decline, 260 advisory, 200 chitchat.
- decline proportion: 13.4%.
- augmentation buckets: 519 workflow twins, 622 decline, 150 ORCA parse,
  50 TS extra, 120 TS clean.
- query diversity: distinct query ratio 1.0.

## What v13.1 Fixed

v13.1 was built after v13 failed to improve generalization. The effective
changes were:

- Diverse, non-templated augmentation instead of mild template rewrites.
- Restoring the empirically stable `5e-6` learning rate after a `1e-4` run
  overfit early.
- Cleaning CLI-semantic targets before training, especially atom-index
  settings and route parameter structures.

The critical dataset bug was that JSON parsing and click resilient parsing
passed while the real runtime parser failed. v13.1 corrected:

- `gaussian.dias.fragment_indices`: flat fragment-1 integer list only.
- `high_level_atoms`, `low_level_atoms`, `freeze_atoms`: integer arrays that
  render to comma/range strings.
- `additional_route_parameters`: string only, never a dict.

The semantic gate rendered every SPEC to CLI and called
`chemsmart.utils.utils.get_list_from_string_range` on atom-index flags:

- 7,812 / 7,812 pass.

## Clean Eval Evidence

Clean full26 representative eval:
`v13_1/artifacts/eval/chemsmart-v131-clean-eval-bg-20260628-232719`

Summary:

- cases: 733.
- command correctness: 698/733 = 0.952.
- acceptance correctness: 698/733 = 0.952.
- adapter valid: 710/733.
- strict CLI valid: 710/733.
- required fields present: 726/733.
- settings allowed: 733/733.
- harness verdicts: 733 ok, 0 reject.
- average latency: 1.154 seconds.
- maximum GPU memory observed: 11,212 MB.

The eval used the corrected v13.1 representative split, not the earlier v12
full26 file that produced misleading harness rejects.

## Runtime Smoke Evidence

Colab session:
`chemsmart-v131-agent-e2e-fixed-20260629-150553`

Local artifact:
`/private/tmp/chemsmart-v131-agent-e2e-fixed-20260629-150553/result.json`

Key outcomes:

- Gaussian one-turn local run with project config passed the runtime semantic
  gate and generated `water_opt_fake.com`.
- Missing server and missing project cases were rejected with explicit runtime
  evidence instead of producing silently invalid commands.
- `sub --test --fake` with server and project reached `status=ready` with
  semantic verdict `warn`; the warning records that submit dry-run did not
  leave a generated input artifact in the current working directory.
- DIAS with explicit fragment atoms rendered `--fragment-indices` after the
  `dias` subcommand and generated the expected fragment and full-system
  Gaussian inputs.
- DIAS/QM/MM requests without required atom selections still exposed model-side
  hallucinated atom indices. Runtime missing-field gates must remain mandatory.

The current detailed status report is maintained in
`docs/agent-current-status.md`.

## Runtime Adapter Contract

The model emits compact SPEC:

```json
{"intent":"workflow","jobs":[{"id":1,"kind":"gaussian.opt","file":"m.xyz","charge":0,"mult":1}]}
```

The runtime owns:

- project defaults.
- functional, basis, solvent, dispersion, grid, SCF settings.
- labels unless explicitly safe and useful.
- canonical TS route tokens that chemsmart derives itself.

The adapter path is:

1. parse JSON.
2. apply `postprocess_v8.py`.
3. apply conservative kind disambiguation in `kind_disambiguator.py`.
4. render with `v8_adapter.py`.
5. validate with click parser.
6. validate atom-index flags with `get_list_from_string_range`.

This keeps known failures out of the generated input path instead of expecting
the model to solve runtime-owned formatting perfectly.

## Local Use

The repo ships a ready-to-merge example config at
`chemsmart/agent/local/agent.v13_1.yaml.example`. The active provider should be
`local_chemsmart_v13_1` and should point `base_model_id` at the merged HF repo.

## Method References Used In Project Notes

- QLoRA-style LoRA fine-tuning was used as an implementation technique, but
  the project report records that the common high learning-rate heuristic did
  not fit this dataset; the empirical v11/v12 rate of `5e-6` was retained.
- Abstention/decline behavior is treated as a harder SFT target, consistent
  with the R-Tuning-style lesson recorded in the v13.1 report. Runtime harness
  enforcement remains mandatory for missing structural inputs.
- Harness engineering patterns from the `oh-my-openagent` / `oh-my-codex`
  investigation are applied as deterministic runtime contracts:
  small context, compact model output, deterministic adapter, real parser
  validation, and artifact evidence.

## Residual Risks

- Model-side decline recall is still incomplete. Missing atom indices,
  missing fragments, and unavailable prior geometry must remain deterministic
  harness failures.
- Some ORCA TS/NEB/QRC misses remain model kind-selection errors, not adapter
  formatting failures.
- Two-step chains passed in clean eval for the available
  `gaussian.opt -> orca.sp` cases, but the clean representative split did not
  contain 3-step chains. Production registration should treat 3-step chains as
  requiring fresh validation.
