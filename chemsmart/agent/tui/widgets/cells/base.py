"""Base class for transcript cells."""

from __future__ import annotations

from textual.widgets import Static


class BaseCell(Static):
    DEFAULT_CSS = """
    BaseCell {
        width: 100%;
        height: auto;
        margin: 0 0 1 0;
        padding: 0 1;
        border: round $surface;
        background: $panel;
    }
    """

    def __init__(
        self,
        renderable,
        *,
        title: str,
        classes: str | None = None,
    ) -> None:
        super().__init__(renderable, classes=classes)
        self.border_title = title
