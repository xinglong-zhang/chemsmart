"""Base class for transcript cells."""

from __future__ import annotations

from rich.text import Text
from textual.message import Message
from textual.widgets import Static


class BaseCell(Static):
    can_focus = True

    DEFAULT_CSS = """
    BaseCell {
        width: 100%;
        height: auto;
        margin: 0 0 1 0;
        padding: 0 1 0 2;
        border: none;
        border-left: outer $surface;
        background: $panel;
        background-tint: $foreground 3%;
    }

    BaseCell:focus {
        border-left: outer $accent;
        background-tint: $accent 8%;
    }
    """

    class CopyRequested(Message):
        def __init__(self, text: str, title: str) -> None:
            super().__init__()
            self.text = text
            self.title = title

    def __init__(
        self,
        renderable,
        *,
        title: str,
        classes: str | None = None,
        expand: bool = False,
        shrink: bool = False,
    ) -> None:
        super().__init__(
            renderable,
            classes=classes,
            expand=expand,
            shrink=shrink,
        )
        self.border_title = title


def no_data_text(message: str = "no data") -> Text:
    return Text(message, style="dim")
