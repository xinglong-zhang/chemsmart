"""One-line controller for a completed turn's tool evidence."""

from __future__ import annotations

from textual import events
from textual.binding import Binding
from textual.message import Message

from .base import BaseCell


class ToolChainToggleCell(BaseCell):
    BINDINGS = [
        Binding("enter", "toggle", "Toggle tools", show=False),
        Binding("space", "toggle", "Toggle tools", show=False),
    ]

    class Toggled(Message):
        def __init__(self, turn_id: str) -> None:
            super().__init__()
            self.turn_id = turn_id

    def __init__(
        self,
        *,
        turn_id: str,
        tool_count: int,
        expanded: bool = False,
    ) -> None:
        self.turn_id = turn_id
        self.tool_count = tool_count
        self.expanded = expanded
        super().__init__(
            self._summary(),
            title=self._title(),
            classes="agent-cell tool-chain-toggle-cell",
        )

    def action_toggle(self) -> None:
        self.post_message(self.Toggled(self.turn_id))

    def on_click(self, event: events.Click) -> None:
        event.stop()
        self.action_toggle()

    def set_expanded(self, expanded: bool) -> None:
        self.expanded = bool(expanded)
        self.border_title = self._title()
        self.update(self._summary())

    def set_tool_count(self, tool_count: int) -> None:
        self.tool_count = max(0, int(tool_count))
        self.update(self._summary())

    def _title(self) -> str:
        return f"Tool chain {'▾' if self.expanded else '▸'}"

    def _summary(self) -> str:
        state = "shown" if self.expanded else "hidden after completion"
        return (
            f"{self.tool_count} tool/validation step(s) · {state} · "
            "Enter, Space, or click to toggle"
        )
