from __future__ import annotations

from functools import lru_cache
from pathlib import Path


@lru_cache(maxsize=None)
def load_prompt(name: str) -> str:
    prompt_path = Path(__file__).with_name(name)
    return prompt_path.read_text(encoding="utf-8")
