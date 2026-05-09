"""Runtime validation transcript cell."""

from __future__ import annotations

from rich.console import Group
from rich.text import Text

from .base import BaseCell


class RuntimeValidationCell(BaseCell):
    def __init__(self, validation: dict | None) -> None:
        self.validation = validation or {}
        super().__init__(
            self._build_renderable(),
            title="Runtime validation",
            classes="runtime-cell",
        )

    def _build_renderable(self):
        status = str(self.validation.get("ok") or "unknown")
        lines = [Text(f"Status: {status}", style="bold")]
        local_issues = self.validation.get("local_issues") or []
        remote_unknown = self.validation.get("remote_unknown") or []
        if local_issues:
            lines.append(Text("Local issues:", style="bold red"))
            lines.extend(Text(f"- {issue}") for issue in local_issues)
        if remote_unknown:
            lines.append(Text("Remote unknowns:", style="bold yellow"))
            lines.extend(Text(f"- {issue}") for issue in remote_unknown)
        if not local_issues and not remote_unknown:
            lines.append(
                Text(
                    (
                        "All runtime checks passed."
                        if self.validation
                        else "No runtime validation data."
                    ),
                    style="dim" if not self.validation else "",
                )
            )
        return Group(*lines)
