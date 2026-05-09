"""User message transcript cell."""

from __future__ import annotations

from rich.text import Text

from .base import BaseCell


class UserMessageCell(BaseCell):
    def __init__(self, text: str) -> None:
        classes = "user-cell"
        if not text.startswith("/"):
            classes += " user-request-cell"
        super().__init__(Text(text), title="You", classes=classes)
