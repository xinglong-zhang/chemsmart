"""Error transcript cell."""

from __future__ import annotations

import json

from rich.console import Group
from rich.text import Text
from textual.binding import Binding

from .base import BaseCell


class ErrorCell(BaseCell):
    BINDINGS = [Binding("e", "toggle_details", "Expand", show=False)]

    def __init__(
        self,
        title: str,
        message: str,
        details: dict | str | None = None,
    ) -> None:
        self.error_title = title or "Error"
        self.message = message or "no data"
        self.details = details
        self.show_details = False
        super().__init__(
            self._build_renderable(),
            title="Error",
            classes="error-cell",
        )

    def action_toggle_details(self) -> None:
        if not self._details_text():
            return
        self.show_details = not self.show_details
        self.update(self._build_renderable())

    def _details_text(self) -> str:
        if self.details is None:
            return ""
        if isinstance(self.details, str):
            return self.details.strip()
        payload = dict(self.details)
        stage = payload.pop("stage", None) or payload.pop("tool", None)
        parts = []
        if stage:
            parts.append(f"Context: {stage}")
        stack = payload.pop("traceback", None) or payload.pop("stack", None)
        if payload:
            parts.append(json.dumps(payload, indent=2, sort_keys=True))
        if stack:
            parts.append(str(stack))
        return "\n\n".join(part for part in parts if part).strip()

    def _build_renderable(self):
        summary = Text.assemble(
            ("✖ ", "bold error"),
            (f"{self.error_title}: ", "bold error"),
            self.message,
        )
        details = self._details_text()
        if not details:
            return Group(summary)
        hint = "  [e expand]" if not self.show_details else "  [e collapse]"
        summary.append(hint, style="dim")
        if not self.show_details:
            return Group(summary, Text(details.splitlines()[0], style="dim"))
        return Group(summary, Text(details, style="dim"))
