"""Mouse-selectable plain-text view of a rendered response."""

from __future__ import annotations

from textual.app import ComposeResult
from textual.binding import Binding
from textual.containers import Vertical
from textual.screen import ModalScreen
from textual.widgets import Static, TextArea


class ResponseCopyOverlay(ModalScreen[None]):
    """Keep rich transcript rendering while offering precise text selection."""

    BINDINGS = [
        Binding("a", "select_all", "Select all", show=False, priority=True),
        Binding("c", "copy_selection", "Copy", show=False, priority=True),
        Binding("escape", "close", "Close", show=False, priority=True),
    ]

    DEFAULT_CSS = """
    ResponseCopyOverlay { align: center middle; }
    #response-copy-modal {
        width: 96%;
        height: 82%;
        border: round $accent;
        padding: 1 2;
        background: $surface;
    }
    #response-copy-help {
        height: 2;
        color: $text-muted;
    }
    #response-copy-text {
        width: 100%;
        height: 1fr;
        border: round $panel;
    }
    """

    def __init__(self, text: str, *, title: str = "Response") -> None:
        super().__init__()
        self.text = text
        self.title = title

    def compose(self) -> ComposeResult:
        with Vertical(id="response-copy-modal"):
            yield Static(
                "Drag to select text. The current selection is copied "
                "automatically · A select all · C copy · Esc close",
                id="response-copy-help",
            )
            yield TextArea(
                self.text,
                read_only=True,
                soft_wrap=True,
                show_line_numbers=False,
                id="response-copy-text",
            )

    def on_mount(self) -> None:
        area = self.query_one("#response-copy-text", TextArea)
        area.border_title = self.title
        area.focus()

    def on_text_area_selection_changed(
        self, event: TextArea.SelectionChanged
    ) -> None:
        if event.text_area.id != "response-copy-text":
            return
        selected = event.text_area.selected_text
        if selected:
            self.app.copy_to_clipboard(selected)

    def action_select_all(self) -> None:
        area = self.query_one("#response-copy-text", TextArea)
        area.select_all()
        self.action_copy_selection()

    def action_copy_selection(self) -> None:
        selected = self.query_one("#response-copy-text", TextArea).selected_text
        if selected:
            self.app.copy_to_clipboard(selected)
            self.notify("Selection copied", timeout=1.5)

    def action_close(self) -> None:
        self.dismiss(None)
