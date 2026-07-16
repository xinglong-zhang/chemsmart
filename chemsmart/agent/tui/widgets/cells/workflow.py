"""Collapsible workflow summary cell."""

from __future__ import annotations

from rich.console import Group
from rich.text import Text
from textual.binding import Binding
from textual.events import Click

from .base import BaseCell

_STATUS_ORDER = {
    "failed": 0,
    "running": 1,
    "queued": 2,
    "cancelled": 3,
    "done": 4,
}


class WorkflowCell(BaseCell):
    BINDINGS = [
        Binding("enter", "toggle_expand", "Expand", show=False),
        Binding("space", "toggle_expand", "Expand", show=False),
    ]

    def __init__(
        self,
        title: str,
        children: list[dict],
        expanded: bool = False,
    ) -> None:
        self.workflow_title = title
        self.workflow_children = list(children)
        self.expanded = expanded
        super().__init__(
            self._build_renderable(),
            title="Workflow",
            classes="workflow-cell",
        )

    def action_toggle_expand(self) -> None:
        self.set_expanded(not self.expanded)

    def on_click(self, _event: Click) -> None:
        self.action_toggle_expand()

    def set_expanded(self, expanded: bool) -> None:
        self.expanded = bool(expanded)
        self.update(self._build_renderable())

    def _build_renderable(self):
        aggregate = min(
            (
                child.get("status") or "done"
                for child in self.workflow_children
            ),
            key=lambda status: _STATUS_ORDER.get(str(status), 99),
            default="done",
        )
        header = Text(f"{self.workflow_title} · {aggregate}", style="bold")
        header.append("  [enter expand]", style="dim")
        if not self.expanded:
            summary = ", ".join(
                f"{child.get('name')}={child.get('status')}"
                for child in self.workflow_children
            )
            return Group(header, Text(summary))
        lines = [header]
        for child in self.workflow_children:
            lines.append(
                Text(
                    f"• {child.get('name')} · {child.get('status')} · {child.get('detail') or ''}".rstrip(),
                )
            )
        return Group(*lines)
