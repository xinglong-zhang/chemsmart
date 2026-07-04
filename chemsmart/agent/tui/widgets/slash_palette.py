"""Inline slash-command suggestion palette."""

from __future__ import annotations

from rich.table import Table
from rich.text import Text
from textual.widgets import Static


class SlashCommandPalette(Static):
    DEFAULT_CSS = """
    SlashCommandPalette {
        height: auto;
        max-height: 10;
        margin: 0 1;
        padding: 0 1;
        border: round $accent;
        background: $panel;
    }

    SlashCommandPalette.hidden {
        display: none;
    }
    """

    def __init__(self, *args, **kwargs) -> None:
        super().__init__("", *args, **kwargs)
        self.matches: list[tuple[str, str]] = []
        self.add_class("hidden")

    def show_matches(
        self,
        *,
        query: str,
        matches: list[tuple[str, str]],
    ) -> None:
        self.matches = list(matches)
        if not matches:
            self.update(Text("No slash commands match.", style="dim"))
            self.remove_class("hidden")
            return
        table = Table.grid(expand=True)
        table.add_column(no_wrap=True, style="bold")
        table.add_column(ratio=1, style="dim")
        for command, description in matches:
            table.add_row(command, description)
        title = "/" if not query else f"/{query}"
        self.border_title = f"Commands: {title}"
        self.update(table)
        self.remove_class("hidden")

    def hide(self) -> None:
        self.matches = []
        self.update("")
        self.add_class("hidden")
