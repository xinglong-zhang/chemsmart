"""Plan transcript cell."""

from __future__ import annotations

from rich.text import Text

from .base import BaseCell, no_data_text


class PlanCell(BaseCell):
    def __init__(self, text: str) -> None:
        renderable = Text(text) if text.strip() else no_data_text()
        super().__init__(renderable, title="Plan", classes="plan-cell")
