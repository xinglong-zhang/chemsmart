"""Error transcript cell."""

from __future__ import annotations

from rich.console import Group
from rich.text import Text

from .base import BaseCell


class ErrorCell(BaseCell):
    def __init__(
        self, title: str, message: str, details: dict | None = None
    ) -> None:
        lines = [Text(title, style="bold red"), Text(message)]
        if details:
            stage = details.get("stage") or details.get("tool")
            if stage:
                lines.append(Text(f"Context: {stage}"))
        super().__init__(
            Group(*lines),
            title="Error",
            classes="error-cell",
        )
