"""User message transcript cell."""

from __future__ import annotations

from rich.text import Text
from textual.app import ComposeResult
from textual.containers import Right
from textual.widget import Widget

from .base import BaseCell


class UserMessageCell(Widget):
    can_focus = True

    DEFAULT_CSS = """
    UserMessageCell {
        width: 100%;
        height: auto;
    }

    UserMessageCell > Right {
        width: 1fr;
        height: auto;
    }

    UserMessageCell:focus BaseCell {
        border-left: outer $accent;
        background-tint: $accent 8%;
    }
    """

    def __init__(self, text: str) -> None:
        super().__init__()
        self.source_text = text
        self._is_slash_command = text.startswith("/")
        bubble_classes = "user-cell"
        if not self._is_slash_command:
            bubble_classes += " user-request-cell"
        self._bubble = BaseCell(
            Text(text),
            title="You",
            classes=bubble_classes,
            shrink=not self._is_slash_command,
        )

    @property
    def bubble(self) -> BaseCell:
        return self._bubble

    @property
    def renderable(self):
        return self._bubble.renderable

    @renderable.setter
    def renderable(self, value) -> None:
        self._bubble.update(value)

    def compose(self) -> ComposeResult:
        if self._is_slash_command:
            yield self._bubble
            return

        with Right():
            yield self._bubble
