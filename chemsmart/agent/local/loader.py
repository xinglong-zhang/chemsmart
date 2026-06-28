"""Load the local v13.1 chemsmart command-synthesis model.

Designed for Colab-T4/L4/A100 or a dedicated GPU box. The HF token must be
supplied via the ``HF_TOKEN`` environment variable; never hardcode it.

Example:
    >>> from chemsmart.agent.local.loader import load_lora_model
    >>> bundle = load_lora_model()
    >>> bundle.model, bundle.tokenizer
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Any

MODEL_REPO_ID = "Smilesjs/chemsmart-qwen2.5-coder-3b-instruct-v13_1"
BASE_MODEL_ID = MODEL_REPO_ID
ADAPTER_REPO_ID = ""
DEFAULT_MAX_SEQ_LENGTH = 4096


@dataclass
class LoadedModel:
    """Loaded model + tokenizer + decoding config."""

    model: Any
    tokenizer: Any
    max_seq_length: int
    base_model_id: str
    adapter_repo_id: str
    model_repo_id: str


def detect_runtime() -> str:
    """Return the best available runtime: ``cuda``, ``mps``, or ``cpu``."""
    try:
        import torch
    except ImportError:
        return "cpu"
    if torch.cuda.is_available():
        return "cuda"
    if getattr(torch.backends, "mps", None) and torch.backends.mps.is_available():
        return "mps"
    return "cpu"


def load_lora_model(
    base_model_id: str = BASE_MODEL_ID,
    adapter_repo_id: str = ADAPTER_REPO_ID,
    hf_token: str | None = None,
    max_seq_length: int = DEFAULT_MAX_SEQ_LENGTH,
    device_map: str | None = None,
    runtime: str | None = None,
) -> LoadedModel:
    """Load a merged model, or base + LoRA when ``adapter_repo_id`` is set.

    Loading strategy is chosen from the available hardware:

    * **cuda** — bfloat16 merged model, or 4-bit NF4 for PEFT LoRA mode.
    * **mps** — Apple Silicon, fp16 weights (no bitsandbytes).
    * **cpu** — fp32 weights, very slow but functional for smoke tests.

    Args:
        base_model_id: HF repo for the merged model, or base model in PEFT
            mode.
        adapter_repo_id: Optional HF repo for the LoRA adapter (PEFT format).
        hf_token: Optional explicit token. Falls back to ``HF_TOKEN`` env var.
            Only required for an initial download; ignored when the model is
            already in ``~/.cache/huggingface/hub/``.
        max_seq_length: Max context length for tokenization defaults.
        device_map: Optional ``accelerate`` device_map override.
        runtime: Force ``cuda``/``mps``/``cpu``; defaults to :func:`detect_runtime`.

    Returns:
        :class:`LoadedModel` with ``model.eval()`` already called.

    Raises:
        RuntimeError: If the active ``huggingface_hub`` is incompatible with the
            transformers pin requires ``huggingface_hub<1.0``.
    """
    _enforce_huggingface_hub_compat()

    token = hf_token or os.environ.get("HF_TOKEN") or None
    if token:
        from huggingface_hub import login

        login(token=token, add_to_git_credential=False)

    import torch
    from transformers import AutoModelForCausalLM, AutoTokenizer

    chosen_runtime = (runtime or detect_runtime()).lower()
    tokenizer = AutoTokenizer.from_pretrained(base_model_id, token=token)

    if chosen_runtime == "cuda":
        if adapter_repo_id:
            from transformers import BitsAndBytesConfig

            bnb = BitsAndBytesConfig(
                load_in_4bit=True,
                bnb_4bit_quant_type="nf4",
                bnb_4bit_use_double_quant=True,
                bnb_4bit_compute_dtype=torch.float16,
            )
            base = AutoModelForCausalLM.from_pretrained(
                base_model_id,
                quantization_config=bnb,
                device_map=device_map or "auto",
                token=token,
            )
        else:
            base = AutoModelForCausalLM.from_pretrained(
                base_model_id,
                torch_dtype=torch.bfloat16,
                device_map=device_map or "auto",
                token=token,
            )
    elif chosen_runtime == "mps":
        base = AutoModelForCausalLM.from_pretrained(
            base_model_id,
            torch_dtype=torch.float16,
            device_map=device_map,
            token=token,
        )
        base = base.to("mps")
    else:
        base = AutoModelForCausalLM.from_pretrained(
            base_model_id,
            torch_dtype=torch.float32,
            device_map=device_map,
            token=token,
        )

    model = base
    if adapter_repo_id:
        from peft import PeftModel

        model = PeftModel.from_pretrained(
            base,
            adapter_repo_id,
            is_trainable=False,
            token=token,
        )
    model.eval()
    return LoadedModel(
        model=model,
        tokenizer=tokenizer,
        max_seq_length=max_seq_length,
        base_model_id=base_model_id,
        adapter_repo_id=adapter_repo_id,
        model_repo_id=base_model_id if not adapter_repo_id else adapter_repo_id,
    )


def _enforce_huggingface_hub_compat() -> None:
    """Block the well-known ``huggingface_hub>=1.0`` ↔ ``transformers==4.56`` clash."""
    try:
        import importlib.metadata as _m

        from packaging.version import Version
    except ImportError:
        return

    try:
        hub_version = Version(_m.version("huggingface_hub"))
    except Exception:
        return

    if hub_version >= Version("1.0.0"):
        raise RuntimeError(
            f"huggingface_hub=={hub_version} conflicts with the local stack "
            "(transformers==4.56.2 requires huggingface_hub>=0.34.0,<1.0). "
            "Reinstall: pip install 'huggingface_hub>=0.34.0,<1.0'"
        )
