"""Load a 4-bit MLX chemsmart command-synthesis model.

This runtime is for Apple Silicon/Metal. It intentionally stays separate from
``loader.py`` because current MLX-LM releases depend on the Transformers 5 /
Hugging Face Hub 1.x stack, while the PyTorch local provider is pinned to the
Transformers 4.x stack.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Any

MODEL_REPO_ID = (
    "Smilesjs/chemsmart-qwen2.5-coder-3b-instruct-v13_1-mlx-4bit"
)
DEFAULT_MAX_SEQ_LENGTH = 4096


@dataclass
class LoadedMLXModel:
    """Loaded MLX model + tokenizer + decoding config."""

    model: Any
    tokenizer: Any
    max_seq_length: int
    base_model_id: str
    adapter_repo_id: str
    model_repo_id: str
    backend: str = "mlx"


def load_mlx_model(
    model_id: str = MODEL_REPO_ID,
    hf_token: str | None = None,
    max_seq_length: int = DEFAULT_MAX_SEQ_LENGTH,
) -> LoadedMLXModel:
    """Load a converted MLX model from a local path or Hugging Face repo.

    Args:
        model_id: Local MLX model directory or Hugging Face model repo.
        hf_token: Optional token for private/gated repos. Falls back to
            ``HF_TOKEN``.
        max_seq_length: Prompt budget advertised to the generator.
    """
    token = hf_token or os.environ.get("HF_TOKEN") or None
    if token:
        from huggingface_hub import login

        login(token=token, add_to_git_credential=False)

    try:
        from mlx_lm import load
    except Exception as exc:  # pragma: no cover - depends on optional runtime
        raise RuntimeError(
            "MLX runtime requires Apple Silicon/Metal and mlx-lm. Install in "
            "an MLX-compatible environment with: pip install 'mlx-lm==0.31.3'"
        ) from exc

    model, tokenizer = load(model_id)
    return LoadedMLXModel(
        model=model,
        tokenizer=tokenizer,
        max_seq_length=max_seq_length,
        base_model_id=model_id,
        adapter_repo_id="",
        model_repo_id=model_id,
    )
