"""Agent/system message transcript cell."""

from __future__ import annotations

from rich.markdown import Markdown

from .base import BaseCell


class AgentMessageCell(BaseCell):
    def __init__(self, text: str, *, title: str = "Agent") -> None:
        super().__init__(Markdown(text), title=title, classes="agent-cell")
