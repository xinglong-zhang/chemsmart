"""Compact, mutable transcript cell for one provider tool call."""

from __future__ import annotations

from time import monotonic

from rich.console import Group
from rich.text import Text
from textual.binding import Binding
from textual.events import Click

from chemsmart.agent.tui.tool_meta import (
    pretty_tool_args,
    render_tool_result_detail,
    render_tool_result_summary,
    tool_risk_badge,
    tool_side_effect_summary,
    tool_status_style,
)

from .base import BaseCell

_AUTO_EXPAND_STATUSES = frozenset(
    {"warn", "reject", "partial", "error", "denied", "skipped", "interrupted"}
)
_TERMINAL_STATUSES = frozenset(
    {
        "ok",
        "warn",
        "reject",
        "partial",
        "error",
        "denied",
        "skipped",
        "interrupted",
    }
)


class ToolCallCell(BaseCell):
    BINDINGS = [
        Binding("enter", "toggle_expand", "Expand", show=False),
        Binding("space", "toggle_expand", "Expand", show=False),
    ]

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
        result: dict | None = None,
        provider_call_id: str = "",
        expanded: bool = False,
    ) -> None:
        self.tool = tool
        self.status = status
        self.description = description or tool
        self.arguments = dict(arguments or {})
        self.note = note or ""
        self.queue_index = queue_index
        self.queue_total = queue_total
        self.session_rule_active = session_rule_active
        self.result = dict(result) if isinstance(result, dict) else None
        self.provider_call_id = provider_call_id
        self.expanded = expanded or status in _AUTO_EXPAND_STATUSES
        self._started_at = monotonic()
        self._finished_at: float | None = (
            self._started_at if status in _TERMINAL_STATUSES else None
        )
        super().__init__(
            self._build_renderable(),
            title=self._title(),
            classes="tool-call-cell",
        )

    @property
    def elapsed_seconds(self) -> float:
        end = (
            self._finished_at if self._finished_at is not None else monotonic()
        )
        return max(0.0, end - self._started_at)

    @property
    def elapsed_text(self) -> str:
        return f"{self.elapsed_seconds:.1f}s"

    def update_lifecycle(
        self,
        *,
        status: str,
        description: str = "",
        arguments: dict | None = None,
        note: str | None = None,
        queue_index: int | None = None,
        queue_total: int | None = None,
        session_rule_active: bool | None = None,
        result: dict | None = None,
    ) -> None:
        self.status = status
        if description:
            self.description = description
        if arguments:
            self.arguments = dict(arguments)
        if note is not None:
            self.note = note
        if queue_index is not None:
            self.queue_index = queue_index
        if queue_total is not None:
            self.queue_total = queue_total
        if session_rule_active is not None:
            self.session_rule_active = session_rule_active
        if result is not None:
            self.result = dict(result)
        if status in _TERMINAL_STATUSES:
            self._finished_at = monotonic()
        if status in _AUTO_EXPAND_STATUSES:
            self.expanded = True
        self.border_title = self._title()
        self.update(self._build_renderable())

    def on_click(self, _event: Click) -> None:
        self.action_toggle_expand()

    def action_toggle_expand(self) -> None:
        self.set_expanded(not self.expanded)

    def set_expanded(self, expanded: bool) -> None:
        self.expanded = bool(expanded)
        self.border_title = self._title()
        self.update(self._build_renderable())

    def _title(self) -> str:
        return f"Tool · {'▾' if self.expanded else '▸'}"

    def activity_snapshot(self) -> dict[str, str]:
        detail = [
            self.description,
            f"effect: {tool_side_effect_summary(self.tool)}",
        ]
        if self.arguments:
            detail.extend(["", pretty_tool_args(self.arguments)])
        if self.note:
            detail.extend(["", self.note])
        return {
            "tool": self.tool,
            "status": self.status,
            "elapsed": self.elapsed_text,
            "detail": "\n".join(detail),
        }

    @staticmethod
    def _compact_text(value: str, *, limit: int = 140) -> str:
        compact = " ".join(value.split())
        if len(compact) <= limit:
            return compact
        return f"{compact[: limit - 1].rstrip()}…"

    def _heading(self, compact: str = "") -> Text:
        risk, risk_style = tool_risk_badge(self.tool)
        status_style = tool_status_style(self.status)
        heading = Text.assemble(
            ("● ", status_style),
            (self.status.replace("_", " ").upper(), f"bold {status_style}"),
            (" · ", "dim"),
            (self.tool, "bold"),
            (" · ", "dim"),
            (f"{risk}", risk_style),
            (" · ", "dim"),
            (self.elapsed_text, "dim"),
        )
        if self.queue_index is not None and self.queue_total:
            heading.append(
                f" · {self.queue_index}/{self.queue_total}", style="dim"
            )
        if compact:
            heading.append(" · ", style="dim")
            heading.append(self._compact_text(compact), style="dim")
        heading.append("  [enter expand]", style="dim")
        return heading

    def _build_renderable(self):
        result_summary = (
            render_tool_result_summary(self.tool, self.result)
            if self.result is not None
            else None
        )
        compact = self.note or result_summary or self.description
        if not self.expanded:
            return Group(self._heading(compact))

        status_style = tool_status_style(self.status)
        details = [
            self._heading(),
            Text(self.description),
            Text(
                f"effect: {tool_side_effect_summary(self.tool)}", style="dim"
            ),
        ]
        if self.session_rule_active:
            details.append(Text("session rule: active", style="dim"))
        if self.note:
            details.append(Text(self.note, style=status_style))
        if self.arguments:
            details.extend(
                [Text(""), Text(pretty_tool_args(self.arguments), style="dim")]
            )
        if self.result is not None:
            result_detail = render_tool_result_detail(self.tool, self.result)
            if result_summary is not None or result_detail is not None:
                details.append(Text(""))
                result_style = (
                    "error" if self.result.get("error") is not None else "dim"
                )
                details.append(
                    Text.assemble(
                        ("Result", "bold"),
                        (" · ", "dim"),
                        (result_summary or "details", result_style),
                    )
                )
                if result_detail is not None:
                    details.extend(result_detail)
        return Group(*details)
