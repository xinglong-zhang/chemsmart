"""Keyboard-navigable snapshot of the current tool loop."""

from __future__ import annotations

from textual.app import ComposeResult
from textual.binding import Binding
from textual.containers import Vertical
from textual.screen import ModalScreen
from textual.widgets import OptionList, Static
from textual.widgets.option_list import Option


class ToolActivityOverlay(ModalScreen[None]):
    BINDINGS = [Binding("escape", "close", "Close", show=False, priority=True)]

    DEFAULT_CSS = """
    ToolActivityOverlay { align: center middle; }
    #activity-modal {
        width: 92;
        height: 28;
        border: round $accent;
        padding: 1 2;
        background: $surface;
    }
    #activity-context { height: auto; margin-bottom: 1; }
    #activity-list { height: 10; }
    #activity-detail { height: 1fr; margin-top: 1; }
    """

    def __init__(
        self,
        entries: list[dict[str, str]],
        *,
        context: dict[str, str] | None = None,
    ) -> None:
        super().__init__()
        self.entries = entries
        self.context = dict(context or {})

    def compose(self) -> ComposeResult:
        options = [
            Option(
                f"{entry.get('status', 'unknown').upper():<12} "
                f"{entry.get('tool', 'tool')}  {entry.get('elapsed', '')}",
                id=str(index),
            )
            for index, entry in enumerate(self.entries)
        ]
        if not options:
            options = [Option("No tool activity in this turn.", disabled=True)]
        with Vertical(id="activity-modal"):
            yield Static(self._context_text(), id="activity-context")
            yield OptionList(*options, id="activity-list")
            yield Static(
                "Select a tool to inspect public evidence.",
                id="activity-detail",
            )

    def _context_text(self) -> str:
        phase = self.context.get("phase", "idle")
        operation = self.context.get("operation", "Ready")
        progress = self.context.get("progress", "no active tool")
        return f"phase: {phase}\noperation: {operation}\nworkflow: {progress}"

    def on_mount(self) -> None:
        self.query_one("#activity-context", Static).border_title = "Workflow"
        listing = self.query_one("#activity-list", OptionList)
        listing.border_title = "Tool activity · ↑↓ inspect · Esc close"
        if self.entries:
            listing.highlighted = len(self.entries) - 1
            self._show_detail(len(self.entries) - 1)
        listing.focus()

    def on_option_list_option_highlighted(
        self, event: OptionList.OptionHighlighted
    ) -> None:
        self._show_detail(event.option_index)

    def _show_detail(self, index: int) -> None:
        if index < 0 or index >= len(self.entries):
            return
        entry = self.entries[index]
        detail = self.query_one("#activity-detail", Static)
        detail.border_title = entry.get("tool", "Tool")
        detail.update(entry.get("detail", "No public detail recorded."))

    def action_close(self) -> None:
        self.dismiss(None)
