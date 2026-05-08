"""Dry-run input transcript cell."""

from __future__ import annotations

from difflib import unified_diff

from rich.console import Group
from rich.syntax import Syntax
from rich.text import Text
from textual.binding import Binding

from .base import BaseCell


class DryRunInputCell(BaseCell):
    BINDINGS = [Binding("d", "toggle_diff", "Diff", show=False)]

    def __init__(
        self,
        content: str,
        inputfile: str | None = None,
        previous_content: str | None = None,
    ) -> None:
        self.content = content
        self.inputfile = inputfile
        self.previous_content = previous_content
        self.show_diff = False
        super().__init__(
            self._build_renderable(),
            title="Dry run input",
            classes="dry-run-cell",
        )

    def action_toggle_diff(self) -> None:
        if self.previous_content is None:
            return
        self.show_diff = not self.show_diff
        self.update(self._build_renderable())

    def _build_renderable(self):
        header = Text()
        if self.inputfile:
            header.append("Input file: ", style="bold")
            header.append(self.inputfile)
        else:
            header.append("Dry-run input generated", style="bold")
        if self.previous_content is not None:
            header.append("  [d diff]", style="dim")

        if self.show_diff and self.previous_content is not None:
            diff = (
                "\n".join(
                    unified_diff(
                        self.previous_content.splitlines(),
                        self.content.splitlines(),
                        fromfile="previous",
                        tofile="current",
                        lineterm="",
                    )
                )
                or "(no changes)"
            )
            body = Syntax(diff, "diff", theme="ansi_dark", word_wrap=True)
        else:
            body = Syntax(
                self.content,
                "text",
                theme="ansi_dark",
                word_wrap=True,
            )
        return Group(header, Text(""), body)
