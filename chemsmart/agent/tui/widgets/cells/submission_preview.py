"""Submission preview transcript cell."""

from __future__ import annotations

from pathlib import Path

from rich.console import Group
from rich.syntax import Syntax
from rich.text import Text

from .base import BaseCell


class SubmissionPreviewCell(BaseCell):
    def __init__(self, preview: dict) -> None:
        self.preview = preview
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
            Text(f"Transport: {transport}", style="bold"),
            Text(f"Duplicate check: {'yes' if duplicate else 'no'}"),
        ]
        if command:
            lines.append(Text(f"Command: {command}"))
        if script_path:
            lines.append(Text(f"Script: {script_path}"))
            path = Path(script_path)
            if path.exists():
                lines.extend(
                    [
                        Text(""),
                        Syntax(
                            path.read_text(encoding="utf-8", errors="ignore"),
                            "bash",
                            theme="ansi_dark",
                            word_wrap=True,
                        ),
                    ]
                )
        return Group(*lines)
