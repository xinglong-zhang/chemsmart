"""Submission preview transcript cell."""

from __future__ import annotations

from pathlib import Path

from rich.console import Group
from rich.syntax import Syntax
from rich.text import Text

from .base import BaseCell


class SubmissionPreviewCell(BaseCell):
    def __init__(self, preview: dict | None) -> None:
        self.preview = preview or {}
        super().__init__(
            self._build_renderable(),
            title="Submission preview",
            classes="submission-cell",
        )

    def _build_renderable(self):
        script_path = self.preview.get("script_path")
        transport = self.preview.get("transport") or "unknown"
        duplicate = (self.preview.get("duplicate_check") or {}).get(
            "duplicate"
        )
        command = self.preview.get("command_executed")
        lines = [
            Text("⇡ Remote submission has not been sent yet", style="bold"),
            Text(""),
            Text(f"Transport: {transport}", style="bold"),
            Text(f"Duplicate check: {'yes' if duplicate else 'no'}"),
            Text(f"Command: {command or 'no data'}"),
        ]
        if script_path:
            lines.append(Text(f"Script: {script_path}"))
            path = Path(script_path)
            if path.exists():
                lines.extend(
                    [
                        Text(""),
                        Syntax(
                            path.read_text(encoding="utf-8", errors="ignore")
                            or "# no data",
                            "bash",
                            theme="ansi_dark",
                            word_wrap=True,
                        ),
                    ]
                )
        elif not self.preview:
            lines.append(Text("No submission preview data.", style="dim"))
        return Group(*lines)
