"""Searchable request-history overlay."""

from __future__ import annotations

from rich.text import Text
from textual.app import ComposeResult
from textual.binding import Binding
from textual.containers import Vertical
from textual.screen import ModalScreen
from textual.widgets import Input, OptionList
from textual.widgets.option_list import Option


class HistorySearchOverlay(ModalScreen[str | None]):
    BINDINGS = [
        Binding("up", "move(-1)", "Previous", show=False, priority=True),
        Binding("down", "move(1)", "Next", show=False, priority=True),
        Binding("enter", "choose", "Use request", show=False, priority=True),
        Binding("escape", "cancel", "Cancel", show=False, priority=True),
    ]

    DEFAULT_CSS = """
    HistorySearchOverlay { align: center middle; }
    #history-modal {
        width: 88;
        height: 22;
        max-width: 96%;
        max-height: 90%;
        border: round $accent;
        padding: 1 2;
        background: $surface;
    }
    #history-query { margin-bottom: 1; }
    #history-options { height: 1fr; }
    """

    def __init__(self, history: list[str]) -> None:
        super().__init__()
        self.history = list(dict.fromkeys(reversed(history)))
        self.filtered: list[str] = []

    def compose(self) -> ComposeResult:
        with Vertical(id="history-modal"):
            yield Input(placeholder="Search request history", id="history-query")
            yield OptionList(id="history-options")

    def on_mount(self) -> None:
        self.query_one("#history-modal", Vertical).border_title = "Request history"
        self._filter("")
        self.query_one("#history-query", Input).focus()

    def on_input_changed(self, event: Input.Changed) -> None:
        if event.input.id == "history-query":
            self._filter(event.value)

    def on_option_list_option_selected(
        self, event: OptionList.OptionSelected
    ) -> None:
        if event.option_list.id != "history-options":
            return
        self._dismiss_index(event.option_index)

    def _filter(self, query: str) -> None:
        lowered = query.strip().lower()
        self.filtered = [
            request
            for request in self.history
            if not lowered or lowered in request.lower()
        ]
        options = self.query_one("#history-options", OptionList)
        options.clear_options()
        if not self.filtered:
            options.add_option(
                Option(Text("No matching requests.", style="dim"), disabled=True)
            )
            return
        options.add_options(
            [
                Option(Text(request), id=f"history-{index}")
                for index, request in enumerate(self.filtered)
            ]
        )
        options.highlighted = 0

    def action_move(self, delta: int) -> None:
        options = self.query_one("#history-options", OptionList)
        if not self.filtered:
            return
        current = options.highlighted if options.highlighted is not None else 0
        options.highlighted = (current + delta) % len(self.filtered)
        options.scroll_to_highlight()

    def action_choose(self) -> None:
        options = self.query_one("#history-options", OptionList)
        index = options.highlighted if options.highlighted is not None else 0
        self._dismiss_index(index)

    def _dismiss_index(self, index: int) -> None:
        if 0 <= index < len(self.filtered):
            self.dismiss(self.filtered[index])

    def action_cancel(self) -> None:
        self.dismiss(None)
