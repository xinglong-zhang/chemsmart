"""Focusable slash-command palette with keyboard and mouse parity."""

from __future__ import annotations

from dataclasses import dataclass

from rich.text import Text
from textual.widgets import OptionList
from textual.widgets.option_list import Option


@dataclass(slots=True, frozen=True)
class SlashPaletteItem:
    command: str
    description: str
    enabled: bool = True
    unavailable_reason: str = ""
    shortcut: str = ""
    aliases: tuple[str, ...] = ()


class SlashCommandPalette(OptionList):
    DEFAULT_CSS = """
    SlashCommandPalette {
        height: auto;
        max-height: 12;
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
        super().__init__(*args, **kwargs)
        self.matches: list[tuple[str, str]] = []
        self.items: list[SlashPaletteItem] = []
        self.add_class("hidden")

    @property
    def is_open(self) -> bool:
        return "hidden" not in self.classes

    def show_matches(
        self,
        *,
        query: str,
        matches: list[tuple[str, str]] | list[SlashPaletteItem],
    ) -> None:
        items = [
            item
            if isinstance(item, SlashPaletteItem)
            else SlashPaletteItem(item[0], item[1])
            for item in matches
        ]
        self.items = items
        self.matches = [(item.command, item.description) for item in items]
        self.clear_options()
        if not items:
            self.add_option(Option(Text("No slash commands match.", style="dim"), disabled=True))
        else:
            self.add_options([self._option(item) for item in items])
            self.highlighted = self._first_enabled_index()
        title = "/" if not query else f"/{query}"
        self.border_title = f"Commands: {title} · ↑↓ select · Tab complete · Enter run"
        self.remove_class("hidden")

    def _option(self, item: SlashPaletteItem) -> Option:
        prompt = Text.assemble(
            (item.command, "bold" if item.enabled else "dim"),
            ("  ", "dim"),
            (item.description, "dim"),
        )
        if item.shortcut:
            prompt.append(f"  [{item.shortcut}]", style="accent")
        if item.unavailable_reason:
            prompt.append(f"  unavailable: {item.unavailable_reason}", style="warning")
        return Option(prompt, id=item.command, disabled=not item.enabled)

    def _first_enabled_index(self) -> int | None:
        for index, item in enumerate(self.items):
            if item.enabled:
                return index
        return None

    def move_highlight(self, delta: int) -> None:
        if not self.items:
            return
        start = self.highlighted if self.highlighted is not None else -1
        for offset in range(1, len(self.items) + 1):
            index = (start + delta * offset) % len(self.items)
            if self.items[index].enabled:
                self.highlighted = index
                self.scroll_to_highlight()
                return

    def selected_item(self) -> SlashPaletteItem | None:
        index = self.highlighted
        if index is None or index >= len(self.items):
            return None
        item = self.items[index]
        return item if item.enabled else None

    def hide(self) -> None:
        self.matches = []
        self.items = []
        self.clear_options()
        self.add_class("hidden")
