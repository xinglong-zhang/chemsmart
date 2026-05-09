"""Agent/system message transcript cell."""

from __future__ import annotations

from rich.markdown import Markdown

from .base import BaseCell


class AgentMessageCell(BaseCell):
    def __init__(self, text, *, title: str = "Agent") -> None:
        if text == "Transcript cleared.":
            text = (
                "대화가 비워졌습니다. 새 요청을 입력하거나 /sessions 로 "
                "이전 기록을 열어보세요."
            )
        self.source_text = text if isinstance(text, str) else None
        renderable = (
            Markdown(text or "_no data_") if isinstance(text, str) else text
        )
        super().__init__(renderable, title=title, classes="agent-cell")
