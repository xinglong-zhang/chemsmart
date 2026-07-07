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
        command: str | None = None,
        cli_grounded: bool = False,
        cli_grounding_issue: str | None = None,
    ) -> None:
        self.content = content or "# no data"
        self.inputfile = inputfile
        self.previous_content = previous_content
        self.command = command
        self.cli_grounded = cli_grounded
        self.cli_grounding_issue = cli_grounding_issue
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
        filename = (
            self.inputfile.rsplit("/", maxsplit=1)[-1]
            if self.inputfile
            else "dry-run input"
        )
        header = Text()
        header.append("✓ ", style="success")
        if self.inputfile:
            header.append(f"{filename} ready", style="bold")
        else:
            header.append("Dry-run input ready", style="bold")
        if self.previous_content is not None:
            header.append("  [d diff]", style="dim")
        path_line = Text()
        if self.inputfile:
            path_line.append("path: ", style="dim")
            path_line.append(self.inputfile, style="dim")

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
        lines = [header]
        if self.inputfile:
            lines.append(path_line)
        lines.append(Text(""))
        lines.append(Text("Generated chemsmart CLI command", style="bold"))
        if self.command:
            lines.append(
                Syntax(
                    self.command,
                    "bash",
                    theme="ansi_dark",
                    word_wrap=True,
                )
            )
        else:
            warning = Text()
            warning.append("CLI-grounding missing", style="bold red")
            issue = self.cli_grounding_issue or (
                "dry-run input did not include a generated chemsmart command"
            )
            warning.append(f": {issue}", style="yellow")
            lines.append(warning)
        if self.command and not self.cli_grounded:
            warning = Text()
            warning.append("CLI-grounding warning", style="bold yellow")
            issue = self.cli_grounding_issue or (
                "command evidence was present but not marked grounded"
            )
            warning.append(f": {issue}", style="yellow")
            lines.append(warning)
        lines.extend([Text(""), body])
        return Group(*lines)
