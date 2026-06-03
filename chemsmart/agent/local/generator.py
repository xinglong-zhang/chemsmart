"""Greedy decode wrapper around the loaded V4 LoRA model.

Returns the parsed planner JSON (with the V4 postprocessor already applied).
The synthesis layer can then translate the atomic plan into a single CLI
command via :mod:`chemsmart.agent.local.adapter`.
"""

from __future__ import annotations

import json
import re
from importlib import resources
from typing import Any

from chemsmart.agent.local.postprocessor import postprocess

_DEFAULT_MAX_NEW_TOKENS = 1024
_DEFAULT_TEMPERATURE = 0.0
_SYSTEM_PROMPT_FILE = "planner_compact.md"


def load_system_prompt() -> str:
    """Return the planner_compact system prompt shipped next to this module."""
    return resources.files("chemsmart.agent.local").joinpath(
        _SYSTEM_PROMPT_FILE
    ).read_text(encoding="utf-8")


def render_chat(
    tokenizer: Any,
    user_query: str,
    history: list[dict[str, str]] | None = None,
    system_prompt: str | None = None,
) -> str:
    """Render a Qwen2.5 chat-format prompt with optional history."""
    system_prompt = system_prompt or load_system_prompt()
    messages: list[dict[str, str]] = [
        {"role": "system", "content": system_prompt}
    ]
    if history:
        messages.extend(history)
    messages.append({"role": "user", "content": user_query})
    return tokenizer.apply_chat_template(
        messages,
        tokenize=False,
        add_generation_prompt=True,
    )


def generate_plan(
    bundle: Any,
    user_query: str,
    history: list[dict[str, str]] | None = None,
    max_new_tokens: int = _DEFAULT_MAX_NEW_TOKENS,
    temperature: float = _DEFAULT_TEMPERATURE,
    apply_postprocessor: bool = True,
) -> dict[str, Any]:
    """Greedy-decode a planner JSON for ``user_query``.

    Args:
        bundle: A :class:`~chemsmart.agent.local.loader.LoadedModel` instance.
        user_query: Natural-language request.
        history: Optional prior turns as Qwen chat messages.
        max_new_tokens: Cap on generated tokens.
        temperature: Decoding temperature; ``0.0`` triggers greedy.
        apply_postprocessor: Whether to apply V4 token-level repairs.

    Returns:
        A planner dict ``{"steps": [...], "rationale": str, "intent": str}``,
        optionally repaired in place by :func:`postprocess`.

    Raises:
        ValueError: If the model output cannot be parsed as JSON after a
            single recovery attempt.
    """
    import torch

    tokenizer = bundle.tokenizer
    model = bundle.model
    prompt = render_chat(tokenizer, user_query, history=history)
    inputs = tokenizer(prompt, return_tensors="pt").to(model.device)

    with torch.no_grad():
        outputs = model.generate(
            **inputs,
            max_new_tokens=max_new_tokens,
            do_sample=temperature > 0.0,
            temperature=max(temperature, 1e-5),
            pad_token_id=tokenizer.eos_token_id,
        )

    generated = outputs[0][inputs["input_ids"].shape[1]:]
    raw = tokenizer.decode(generated, skip_special_tokens=True)
    plan = _parse_planner_json(raw)
    if apply_postprocessor:
        plan = postprocess(plan, user_query)
    return plan


_JSON_BLOCK_RE = re.compile(r"\{.*\}", re.DOTALL)


def _parse_planner_json(raw: str) -> dict[str, Any]:
    """Extract the outermost JSON object from a model's raw decode output."""
    try:
        return json.loads(raw)
    except json.JSONDecodeError:
        pass

    match = _JSON_BLOCK_RE.search(raw)
    if not match:
        raise ValueError(f"no JSON object in model output: {raw[:200]!r}")
    return json.loads(match.group(0))
