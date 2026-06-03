# Local LoRA model adapter for chemsmart agent

Serves [`Smilesjs/chemsmart-qwen2.5-7b-lora`](https://huggingface.co/Smilesjs/chemsmart-qwen2.5-7b-lora)
(a Qwen2.5-7B-Instruct + LoRA fine-tune on the chemsmart v4 dataset) behind an
OpenAI-compatible HTTP endpoint so that `chemsmart agent ask` can call it with
no changes to the provider plumbing.

## Architecture

```
        +---------------------+        HTTP            +-------------------+
user -> | chemsmart agent ask | <--------------------> | local/server.py   |
        | (SynthesisSession)  |  /v1/chat/completions  | FastAPI shim      |
        +---------------------+                        +---------+---------+
                                                                 |
                                                                 v
                                                 +-------------------------------+
                                                 | local/generator.py            |
                                                 |  + local/loader.py (Qwen+LoRA)|
                                                 |  + local/postprocessor.py     |
                                                 |  + local/adapter.py           |
                                                 +-------------------------------+
```

* `loader.py` — loads Qwen2.5-7B-Instruct in 4-bit NF4 + applies the LoRA adapter.
* `generator.py` — greedy decode → planner JSON (with the V4 postprocessor applied).
* `postprocessor.py` — token-level repairs (filepath fidelity, label rule,
  step refs, kind whitelist guard). Moved from `chemsmart-finetune/v4/`.
* `adapter.py` — converts the atomic `{steps, rationale, intent}` planner JSON
  to the `{status, command, …}` schema expected by
  `chemsmart.agent.synthesis.SynthesisSession`.
* `server.py` — OpenAI-compatible FastAPI shim. Designed to run on Colab
  (T4/L4/A100) and be tunnelled with ngrok or Cloudflare Tunnel.
* `planner_compact.md` — the exact 4.9 KB system prompt the LoRA was trained on.

## Quickstart (one command on Mac/Linux)

```bash
chemsmart config agent
# → choose `local` at the provider prompt
# → paste a HF read token (or leave blank to use $HF_TOKEN at runtime)
# → confirm pre-fetch — base model (~15 GB) + adapter (~177 MB) are pulled to
#   ~/.cache/huggingface/hub/ in one shot.
chemsmart agent doctor   # verifies tokenizer loads + reports model id
chemsmart agent ask "opt of h2o at b3lyp/6-31g*"
```

The `chemsmart config agent` flow writes an `agent.yaml` entry like:

```yaml
active: local_chemsmart_v4
providers:
  local_chemsmart_v4:
    type: local
    model: chemsmart-qwen2.5-7b-lora
    base_model_id: Qwen/Qwen2.5-7B-Instruct
    adapter_repo_id: Smilesjs/chemsmart-qwen2.5-7b-lora
    hf_token_env: HF_TOKEN
    hf_token: ""
    runtime: ""       # auto: cuda | mps | cpu
```

The `LocalProvider` is lazy — `chemsmart agent doctor` only touches the
tokenizer (a few MB), so config-time validation is fast. The 15 GB safetensor
weights are only paged in when the first `chemsmart agent ask` call lands.

## Quickstart (Colab via OpenAI-compatible HTTP)

```python
# Cell 1 — pin the V4 inference stack (huggingface_hub<1.0!)
%%capture
!pip install -q --upgrade "huggingface_hub>=0.34.0,<1.0"
!pip install -q --upgrade "transformers==4.56.2" "peft==0.16.0" \
    "accelerate==1.10.0" "bitsandbytes==0.47.0"
!pip install -q sentencepiece protobuf fastapi uvicorn pyngrok

# Cell 2 — install chemsmart from your fork (read-only path is fine)
!pip install -q -e /content/drive/MyDrive/chemsmart  # or git+https://...

# Cell 3 — set the HF read token (rotate after first run)
import os
os.environ["HF_TOKEN"] = "hf_REPLACE_ME"

# Cell 4 — load model
from chemsmart.agent.local.loader import load_lora_model
bundle = load_lora_model()

# Cell 5 — start server + ngrok tunnel
from chemsmart.agent.local.server import build_app
from pyngrok import ngrok
import nest_asyncio, uvicorn, threading
nest_asyncio.apply()
app = build_app(bundle)
public = ngrok.connect(8000, "http")
print("Public base_url:", f"{public.public_url}/v1")
threading.Thread(
    target=uvicorn.run,
    args=(app,),
    kwargs={"host": "0.0.0.0", "port": 8000, "log_level": "warning"},
    daemon=True,
).start()
```

Then copy the printed `Public base_url` into your local
`~/.chemsmart/agent/agent.yaml`:

```yaml
active: local_chemsmart_v4
providers:
  local_chemsmart_v4:
    type: openai
    model: chemsmart-qwen2.5-7b-lora
    base_url: https://<ngrok-tunnel>/v1
    api_key: local-not-used
```

Finally:

```bash
chemsmart agent doctor                     # should resolve local_chemsmart_v4
chemsmart agent ask "opt of h2o at b3lyp/6-31g*"
```

## What the adapter can and cannot synthesize

| build_job.kind | Mapped CLI | Status |
|---|---|---|
| `gaussian.{sp,opt,ts,irc,scan,modred,tddft}` | `chemsmart {sub,run} gaussian {sp,opt,ts,irc,modred,td}` | `ready` |
| `orca.{sp,opt,ts,irc,scan}` | `chemsmart {sub,run} orca {sp,opt,ts,irc,scan}` | `ready` |
| `gaussian.{resp,nci,dias,crest,traj,qmmm,qrc,wbi,freq}` | — | `needs_clarification` |
| `orca.{modred,neb,qmmm,qrc}` | — | `needs_clarification` |

`adapter.py` fails-closed on the 13 non-canonical kinds so the user always sees
"this kind has no direct chemsmart CLI subcommand" rather than an invalid
command. Extend `_KIND_TO_CLI` when the corresponding CLI subcommands ship.

## Runtime dependencies (must be installed once)

The local provider keeps torch / transformers / peft out of chemsmart's core
install footprint. Install them in your active env before running
`chemsmart agent ask` with `type: local`:

```bash
pip install "huggingface_hub>=0.34.0,<1.0" "transformers==4.56.2" \
    "peft==0.16.0" "accelerate==1.10.0" "bitsandbytes==0.47.0"
# bitsandbytes is CUDA-only — skip on Mac/CPU-only boxes
```

`chemsmart agent doctor` will surface a precise install line if any of these
are missing.

## Base model variant

The LoRA was trained on `unsloth/Qwen2.5-7B-Instruct-bnb-4bit` (a pre-quantized
clone of Qwen2.5-7B-Instruct). The loader defaults to the upstream
`Qwen/Qwen2.5-7B-Instruct` because PEFT applies adapters by layer name
(`q_proj`, `k_proj`, …), which are identical between the two repos. If you
hit layer-name mismatches, override `base_model_id` in `agent.yaml` to the
Unsloth repo.

## Token handling

Never hardcode `HF_TOKEN` in source. The loader reads it from the environment
or an explicit `hf_token=` kwarg. Tokens copied into Colab notebooks should be
**read-only** tokens and rotated after each session — the HF token panel at
<https://huggingface.co/settings/tokens> shows recent usage.

## Audit reference

* Compatibility verdict and three forward paths are in
  `chemsmart-finetune/v4/compatibility_audit.md`.
* The atomic-plan ↔ single-CLI mismatch is also documented in the auto-memory
  entry `chemsmart-v4-agent-ask-mismatch`.
