"""Transcript cell for tool-use requests and outcomes."""

from __future__ import annotations

from rich.console import Group
from rich.text import Text

from chemsmart.agent.tui.tool_meta import (
    pretty_tool_args,
    tool_risk_badge,
    tool_side_effect_summary,
    tool_status_style,
)

from .base import BaseCell


class ToolCallCell(BaseCell):
    def __init__(
        self,
        *,
        tool: str,
        status: str,
        description: str,
        arguments: dict | None = None,
        note: str | None = None,
        queue_index: int | None = None,
        queue_total: int | None = None,
        session_rule_active: bool = False,
    ) -> None:
        self.tool = tool
        self.status = status
        self.description = description or tool
        self.arguments = dict(arguments or {})
        self.note = note or ""
        self.queue_index = queue_index
        self.queue_total = queue_total
        self.session_rule_active = session_rule_active
        super().__init__(
            self._build_renderable(),
            title="Tool use",
            classes="tool-call-cell",
        )

    def _build_renderable(self):
        risk, risk_style = tool_risk_badge(self.tool)
        status_style = tool_status_style(self.status)
        label = self.status.replace("_", " ")
        heading = Text.assemble(
            ("● ", status_style),
            (label.upper(), f"bold {status_style}"),
            (" ", "dim"),
            (self.tool, "bold"),
            (" ", "dim"),
            (f"[{risk}]", risk_style),
        )
        if (
            self.queue_index is not None
            and self.queue_total is not None
            and self.queue_total > 0
        ):
            heading.append(
                f"  ({self.queue_index}/{self.queue_total})", style="dim"
            )

        details = [
            Text(self.description),
            Text(
                f"effect: {tool_side_effect_summary(self.tool)}",
                style="dim",
            ),
        ]
        if self.session_rule_active:
            details.append(Text("session rule: active", style="dim"))
        if self.note:
            details.append(Text(self.note, style=status_style))
        if self.arguments:
            details.append(Text(""))
            details.append(Text(pretty_tool_args(self.arguments), style="dim"))
        return Group(heading, *details)
