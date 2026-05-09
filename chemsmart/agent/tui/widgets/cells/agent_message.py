"""Agent/system message transcript cell."""

from __future__ import annotations

from rich.markdown import Markdown

from .base import BaseCell


class AgentMessageCell(BaseCell):
    def __init__(self, text, *, title: str = "Agent") -> None:
        self.source_text = text if isinstance(text, str) else None
        renderable = (
            Markdown(text or "_no data_") if isinstance(text, str) else text
        )
        super().__init__(renderable, title=title, classes="agent-cell")
