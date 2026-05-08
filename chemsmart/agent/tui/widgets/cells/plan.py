"""Plan transcript cell."""

from __future__ import annotations

from rich.text import Text

from .base import BaseCell


class PlanCell(BaseCell):
    def __init__(self, text: str) -> None:
        super().__init__(Text(text), title="Plan", classes="plan-cell")
