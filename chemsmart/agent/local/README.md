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
