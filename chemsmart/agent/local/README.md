# Local v13.1 model adapter for chemsmart agent

Serves `Smilesjs/chemsmart-qwen2.5-coder-3b-instruct-v13_1`, the deployable
merged v13.1 command-synthesis model. The model emits compact SPEC JSON; the
runtime adapter expands it into real `chemsmart run|sub` commands and validates
against the chemsmart CLI plus runtime atom-index parser.

## Architecture

```
user -> chemsmart agent ask
      -> SynthesisSession
      -> LocalProvider
      -> local/generator.py
      -> postprocess_v8.py
      -> v8_adapter.py
      -> chemsmart CLI validation
```

- `loader.py` loads the merged v13.1 HF model by default. If
  `adapter_repo_id` is set, it falls back to the older base+PEFT LoRA path.
- `mlx_loader.py` loads the converted 4-bit MLX model when `runtime: mlx`.
  This is the Apple Silicon path; it keeps the same compact SPEC,
  postprocessor, adapter, and harness.
- `generator.py` greedy-decodes compact JSON SPEC.
- `postprocess_v8.py` strips runtime-owned leaks, canonicalizes opt+freq, and
  removes TS route tokens that chemsmart derives itself.
- `v8_adapter.py` renders commands and runs both click parsing and semantic
  range parsing for `fragment_indices`, `high_level_atoms`, `low_level_atoms`,
  and `freeze_atoms`.
- `adapter.py` converts compact SPEC results into the public synthesis schema.

## Local configuration

`~/.chemsmart/agent/agent.yaml` should contain:

```yaml
active: local_chemsmart_v13_1
providers:
  local_chemsmart_v13_1:
    type: local
    model: chemsmart-qwen2.5-coder-3b-instruct-v13_1
    base_model_id: Smilesjs/chemsmart-qwen2.5-coder-3b-instruct-v13_1
    adapter_repo_id: ""
    hf_token_env: HF_TOKEN
    hf_token: ""
    runtime: ""
    project: ""  # optional selector for a same-named workspace project
  local_chemsmart_v13_1_mlx4:
    type: local
    model: chemsmart-qwen2.5-coder-3b-instruct-v13_1-mlx-4bit
    base_model_id: Smilesjs/chemsmart-qwen2.5-coder-3b-instruct-v13_1-mlx-4bit
    adapter_repo_id: ""
    hf_token_env: HF_TOKEN
    hf_token: ""
    runtime: mlx
    project: ""  # optional selector for a same-named workspace project
```

Install runtime dependencies once:

```bash
pip install "huggingface_hub>=0.34.0,<1.0" "transformers==4.56.2" \
    "accelerate==1.10.0"
```

For PEFT adapter experiments, also install:

```bash
pip install "peft==0.16.0" "bitsandbytes==0.47.0"
```

For Apple Silicon MLX 4-bit inference, use a separate environment or accept
the newer MLX stack pins:

```bash
pip install "mlx-lm==0.31.3"
```

`mlx-lm==0.31.3` currently depends on `transformers>=5` and
`huggingface_hub>=1.5`, which intentionally differs from the PyTorch local
provider stack. Select `active: local_chemsmart_v13_1_mlx4` to use the MLX
runtime.

If runtime semantic validation reports missing project settings, create or
select a workspace-local YAML such as
`<workspace>/.chemsmart/gaussian/test.yaml`. The optional `project: test`
provider value can select that same-named workspace candidate, but it cannot
replace an absent file. Project/method/solvent resolution remains
runtime-owned; the model target must not emit those fields.

Then verify:

```bash
chemsmart agent doctor
chemsmart agent ask "run a Gaussian opt+freq on mols/oxetane.xyz"
```

## v13.1 validation snapshot

Clean full26 representative eval, 2026-06-28:

- model: `Smilesjs/chemsmart-qwen2.5-coder-3b-instruct-v13_1`
- cases: 733
- command correctness: 698/733 = 0.952
- adapter validity: 710/733
- strict CLI validity: 710/733
- harness verdicts: 733 ok, 0 reject
- opt+freq drift: eliminated in v13.1 report

Dataset-level semantic gate:

- records checked: 7,812
- pass: 7,812
- gate: render compact SPEC to command, then run the real
  `get_list_from_string_range` parser on atom-index arguments.

Known residual limitation: model-side declines remain incomplete. Runtime
harness/validation must continue to enforce missing structural inputs rather
than relying only on the model to abstain.
