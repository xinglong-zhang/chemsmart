"""Bottom composer widget."""

from __future__ import annotations

from textual.binding import Binding
from textual.message import Message
from textual.widgets import TextArea


class Composer(TextArea):
    DEFAULT_CSS = """
    Composer {
        height: 5;
        min-height: 5;
        border: round $primary;
    }
    """

    BINDINGS = [
        Binding("enter", "submit", show=False, priority=True),
        Binding("ctrl+j", "insert_newline", show=False, priority=True),
        Binding("shift+enter", "insert_newline", show=False, priority=True),
    ]

    class Submitted(Message):
        def __init__(self, composer: "Composer", text: str) -> None:
            super().__init__()
            self.composer = composer
            self.text = text

    def __init__(self) -> None:
        super().__init__(
            text="",
            soft_wrap=True,
            show_line_numbers=False,
            id="composer",
        )
        self.border_title = "Request"

    def action_submit(self) -> None:
        text = self.text.strip()
        if not text:
            return
        self.post_message(self.Submitted(self, self.text))

    def action_insert_newline(self) -> None:
        self.insert("\n")

    def clear_text(self) -> None:
        self.load_text("")
