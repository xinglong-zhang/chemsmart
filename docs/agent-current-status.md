# ChemSmart Agent Current Status

Date: 2026-06-29

Latest documentation addendum: 2026-07-16

This document summarizes the current local chemsmart CLI agent state after the
v13.1 model integration, runtime semantic gate work, and Colab smoke testing.

## 2026-07-16 TUI And Workspace Addendum

The current TUI keeps model behavior subordinate to CHEMSMART CLI and
workspace truth while reducing transcript clutter:

- project YAML is discovered only from the launch workspace under
  `.chemsmart/{gaussian,orca}/`; one candidate auto-loads and multiple
  candidates require an explicit `Shift+Tab` selection;
- `/write-project` uses the latest provider-call-matched validation receipt and
  asks whether to overwrite the active YAML or create a new project name;
- generated commands remain the final visible artifact, while deterministic
  parsing, semantic/intent evidence, repair attempts, and intermediate answers
  collapse after completion into a reversible Tool chain row;
- response cells open a mouse-selectable copy view without replacing Rich
  rendering in the transcript;
- local calculations remain asynchronous and visible through the status strip
  and `Ctrl+B` calculation monitor, with chemistry receipts separated from raw
  logs.

Verification on the integrated branch: `150 passed` for the TUI plus
project-YAML slice and `1062 passed, 5 warnings` for the full `tests/agent`
suite.

Historical model-performance numbers below are retained as dated evidence;
they are not a claim about a newer model release.

## 2026-07-04 Addendum

The newest changes since the v13.1 local-agent report are mostly runtime,
adapter, and TUI grounding improvements rather than a new model-performance
benchmark:

- job-synthesis output remains CHEMSMART CLI first: user-facing job setup must
  produce a `chemsmart run ...` or `chemsmart sub ...` command before generated
  input evidence;
- API/frontier providers can explain, critique, or repair a generated command,
  but their answer is grounded by the deterministic command parser and runtime
  semantic gate;
- the TUI now renders a deterministic command interpretation cell and a
  collapsible public decision trace for API-routed turns;
- `/` opens a slash-command palette with prefix filtering, including `/init`
  for project-YAML build mode;
- project-YAML tools can extract a literature protocol, render CHEMSMART
  `gas:`/`solv:` YAML, validate it, critique it, and write it only after
  approval;
- basis-set phrase resolution is now a read-only `search_basis_sets` tool backed
  by a Basis Set Exchange catalog, returning compact top-k candidates instead
  of injecting all basis names into the prompt;
- Gaussian TS route invariants now reject duplicate `opt=(...)` blocks and
  runtime-owned TS tokens leaked outside the TS opt block.

Focused local verification for the latest doc update:

```text
74 passed
```

from:

```text
tests/agent/harness/test_basis_catalog.py
tests/agent/test_command_answerer.py
tests/agent/test_model_command_parser.py
tests/agent/test_project_yaml_tools.py
tests/agent/test_synthesis.py
tests/agent/tui/test_synthesis_mode.py
tests/agent/tui/test_slash_commands.py
tests/agent/tui/test_track_b_ux.py
```

## Executive Summary

The v13.1 local agent is usable for local review-and-approve command synthesis.
It is not yet ready for fully autonomous production execution.

The deployable model is:

`Smilesjs/chemsmart-qwen2.5-coder-3b-instruct-v13_1`

The production stack is intentionally hybrid:

1. the model emits a compact workflow SPEC;
2. deterministic postprocessing normalizes known model drift;
3. the v8 adapter renders the real chemsmart CLI;
4. the runtime semantic gate executes the generated command in safe fake/dry
   mode;
5. failed runtime evidence is surfaced as reject/warn metadata instead of
   relying on JSON validity alone.

## Current Runtime Evidence

Latest integration commits on the main-integration branch:

- `f8b7f4c` fixed local generator postprocessor wiring.
- `b958c1f` added runtime semantic validation for synthesized commands.
- `125642d` fixed fake-runner selection and structural option rendering.
- `1c0c004` treats successful submit dry-runs without observed generated input
  as a warning rather than a hard rejection.

Focused local validation after these changes:

- `tests/test_jobrunner_fake_selection.py`
- `tests/agent/test_v8_adapter.py`
- `tests/agent/harness/test_command_semantics.py`
- `tests/agent/test_synthesis.py`

Result: 35 passed.

The broader `tests/agent` suite was previously at 612 passed / 1 failed before
the final submit-warning adjustment. The remaining failure was an unrelated TUI
wrapping assertion where `crest_best.xyz` was rendered across a line break.

## Model Performance

The most important measured result is the corrected v13.1 representative full26
eval, not the earlier contaminated v13/v12 evals.

| Evaluation | Scope | Result |
|---|---:|---:|
| v13.1 representative full26 clean eval | 733 cases | 698/733 command correctness = 0.952 |
| v13.1 acceptance correctness | 733 cases | 698/733 = 0.952 |
| v13.1 adapter validity | 733 cases | 710/733 |
| v13.1 strict CLI validity | 733 cases | 710/733 |
| v13.1 settings allowed | 733 cases | 733/733 |
| v13.1 harness verdicts | 733 cases | 733 ok, 0 reject |
| v13.1 behavioral score in finetune report | corrected full26 | 0.947 |
| v13.1 semantic dataset gate | 7,812 records | 7,812/7,812 pass |
| historical v8 easy45 baseline | 45 cases | 33/45 = 0.733 |
| historical v8 native multi-turn baseline | 8 follow-ups | 7/8 correct |

The v13.1 numbers and the historical v8 numbers should not be mixed as the
same benchmark. v13.1 is the current deployable model; v8 remains useful only
as a baseline showing the magnitude of improvement after the dataset and
harness refactor.

## One-Turn Behavior

One-turn fully specified commands are now the strongest path.

In the 2026-06-29 Colab smoke test, this one-turn request:

```text
Optimize /content/water.xyz with Gaussian; charge 0 multiplicity 1.
```

validated as:

```text
chemsmart run gaussian -p test -f /content/water.xyz -c 0 -m 1 opt
```

The runtime semantic verdict was `ok`, and fake execution generated
`water_opt_fake.com` with route:

```text
# opt freq m062x def2svp
```

The corrected v13.1 eval also shows that the former `opt+freq` route drift was
removed in the model/harness path. The v13.1 result report records `opt+freq`
drift as 0 after previously observing 22 cases in v11 and 26 cases in v13.

## Multi-Turn And Missing-Information Behavior

The current multi-turn story has two layers.

At the runtime layer, missing environment information is now caught reliably by
the semantic gate. The Colab smoke test showed this chain:

1. submit request without server config: rejected with missing valid server
   evidence;
2. same request after server setup but without project config: rejected with
   missing valid project evidence;
3. same command after adding `-p test`: accepted as `ready` with semantic
   verdict `warn`.

The warning on step 3 was:

```text
cmd.semantic.submit_generated_input_not_observed
```

That means `sub --test --fake` succeeded far enough to validate the command
shape, but the dry-run did not leave a generated input artifact in the current
working directory. This is acceptable as a warning for submit dry-runs, but it
is weaker evidence than `run --fake --no-scratch`, which directly produces and
inspects input files.

At the model layer, broad v13.1 multi-turn benchmarking is still pending. The
historical v8 native multi-turn benchmark was 7/8, and the current runtime can
carry conversation/session state, but v13.1 has only been smoke-tested for
environment fallback and command correction rather than evaluated over a full
multi-turn suite.

## Structural Job Behavior

Structural jobs improved at the adapter/runtime level.

The Colab smoke test verified that DIAS with explicit fragments now renders the
option after the `dias` subcommand:

```text
chemsmart run gaussian -p test -f /content/water.xyz -c 0 -m 1 -l water_12_01 dias --fragment-indices 1,2
```

The semantic verdict was `ok`, and fake execution generated the expected three
Gaussian inputs for fragment 1, fragment 2, and the combined system.

The remaining problem is model-side abstention. When the user asks for DIAS or
QM/MM without required atom selections, the model can still hallucinate atom
indices instead of asking for missing information. This is why the harness must
remain the source of truth for missing structural fields.

## Production Readiness

Ready now:

- local v13.1 provider wiring;
- compact SPEC to real chemsmart CLI rendering;
- real parser validation for atom-index options;
- safe runtime semantic validation for generated commands;
- Gaussian `run --fake --no-scratch` generated-input inspection;
- session evidence through semantic verdicts, failed rule ids, and generated
  input metadata.

Not production-autonomous yet:

- model-side decline recall is still weak for missing structural inputs;
- v13.1 needs a dedicated multi-turn benchmark rather than only historical v8
  multi-turn evidence;
- submit dry-run evidence is currently weaker than run dry-run evidence;
- ORCA advanced workflows need the same per-kind semantic stress testing as the
  Gaussian paths already exercised.

Operational recommendation:

Use v13.1 locally in review-and-approve mode. Treat `semantic_verdict=ok` as
usable evidence, `warn` as user-review required, and `reject` as a hard block.
Do not enable unattended submit execution until missing structural input
declines and submit artifact evidence are fully covered by runtime gates.
