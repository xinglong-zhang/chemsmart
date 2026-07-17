"""Agent/system message transcript cell."""

from __future__ import annotations

from rich.markdown import Markdown
from textual import events

from .base import BaseCell


class AgentMessageCell(BaseCell):
    def __init__(self, text, *, title: str = "Agent") -> None:
        self.source_text = text if isinstance(text, str) else None
        renderable = (
            Markdown(text or "_no data_") if isinstance(text, str) else text
        )
        super().__init__(renderable, title=title, classes="agent-cell")
        if self.source_text:
            self.tooltip = "Click to select and copy this response"

    def on_click(self, event: events.Click) -> None:
        if not self.source_text:
            return
        event.stop()
        self.post_message(
            self.CopyRequested(
                self.source_text, self.border_title or "Response"
            )
        )
