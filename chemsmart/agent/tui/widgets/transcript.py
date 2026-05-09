"""Scrollable transcript widget."""

from __future__ import annotations

from textual.containers import Vertical, VerticalScroll
from textual.widget import Widget


class Transcript(VerticalScroll):
    DEFAULT_CSS = """
    Transcript {
        border: round $panel;
        padding: 0 1;
        overflow-x: hidden;
        overflow-y: auto;
    }

    #cells {
        width: 100%;
        height: auto;
    }
    """

    def compose(self):
        yield Vertical(id="cells")

    def add_cell(self, cell: Widget) -> None:
        self.query_one("#cells", Vertical).mount(cell)
        self.call_after_refresh(self.scroll_end, animate=False)

    def clear_cells(self) -> None:
        cells = self.query_one("#cells", Vertical)
        for child in list(cells.children):
            child.remove()
