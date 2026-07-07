"""Simple text-entry overlay."""

from __future__ import annotations

from textual.app import ComposeResult
from textual.binding import Binding
from textual.containers import Vertical
from textual.screen import ModalScreen
from textual.widgets import Input, Static


class TextPromptOverlay(ModalScreen[str | None]):
    BINDINGS = [
        Binding("enter", "submit", "Submit", show=False),
        Binding("escape", "cancel", "Cancel", show=False),
    ]

    DEFAULT_CSS = """
    TextPromptOverlay {
        align: center middle;
    }

    #text-prompt-modal {
        width: 88;
        height: auto;
        border: round $primary;
        padding: 1 2;
        background: $surface;
    }

    #text-prompt-input {
        margin-top: 1;
    }
    """

    def __init__(self, *, title: str, prompt: str, initial: str = "") -> None:
        super().__init__()
        self.title = title
        self.prompt = prompt
        self.initial = initial

    def compose(self) -> ComposeResult:
        with Vertical(id="text-prompt-modal"):
            prompt = Static(self.prompt, id="text-prompt-summary")
            prompt.border_title = self.title
            yield prompt
            yield Input(value=self.initial, id="text-prompt-input")

    def on_mount(self) -> None:
        self.query_one("#text-prompt-input", Input).focus()

    def action_submit(self) -> None:
        value = self.query_one("#text-prompt-input", Input).value.strip()
        if not value:
            return
        self.dismiss(value)

    def action_cancel(self) -> None:
        self.dismiss(None)
