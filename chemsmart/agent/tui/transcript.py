"""The transcript: one widget per block so rows can change state in place.

A `RichLog` is append-only, which is why the previous interface could show
nothing until a whole session finished. Here every block is its own widget:
a running tool row can later mute itself, an error can turn red, and the
view stays pinned to the bottom through Textual's anchor until the human
scrolls away. An `ExportRecorder` keeps every renderable independently of
widget trimming, so `/export` always has the full transcript.
"""

from __future__ import annotations

import io
from dataclasses import dataclass, field
from typing import Any

from rich.console import Console
from rich.text import Text
from textual.containers import VerticalScroll
from textual.widgets import Static

#: Keep at most this many blocks mounted; ExportRecorder retains everything.
MAX_MOUNTED_BLOCKS = 400
#: How many oldest blocks are trimmed together when the cap is exceeded.
TRIM_BATCH = 100


@dataclass
class ExportRecorder:
    """Everything ever written, for /export -- independent of trimming."""

    entries: list[Any] = field(default_factory=list)

    def record(self, renderable: Any) -> None:
        self.entries.append(renderable)

    def export_text(self, *, width: int = 100) -> str:
        console = Console(record=True, width=width, file=io.StringIO())
        for renderable in self.entries:
            console.print(renderable)
        return console.export_text()


class TranscriptBlock(Static):
    """One transcript entry; subclasses/states restyle it in place."""


class ToolRow(Static):
    """A one-line humanized tool row that mutates as its call settles.

    States follow the theme's single rule: running = text, finished =
    muted, failed = error.
    """

    def __init__(self, text: Text, *, state: str = "running") -> None:
        super().__init__(text)
        self.state = state

    def settle(self, text: Text, *, state: str) -> None:
        self.state = state
        self.update(text)


class TranscriptView(VerticalScroll):
    """Scrolling column of blocks, pinned to the bottom until released."""

    def __init__(self, *, id: str | None = None) -> None:  # noqa: A002
        super().__init__(id=id)
        self.recorder = ExportRecorder()
        self._trimmed = 0
        self._trim_marker: Static | None = None

    def on_mount(self) -> None:
        self.anchor()

    def add_block(self, renderable: Any) -> TranscriptBlock:
        self.recorder.record(renderable)
        block = TranscriptBlock(renderable)
        self.mount(block)
        self._enforce_cap()
        return block

    def add_row(self, text: Text, *, state: str = "running") -> ToolRow:
        self.recorder.record(text)
        row = ToolRow(text, state=state)
        self.mount(row)
        self._enforce_cap()
        return row

    def settle_row(self, row: ToolRow, text: Text, *, state: str) -> None:
        """Mutate a mounted row in place and keep the export record true."""

        self.recorder.record(text)
        row.settle(text, state=state)

    def _enforce_cap(self) -> None:
        blocks = [
            child
            for child in self.children
            if isinstance(child, (TranscriptBlock, ToolRow))
        ]
        if len(blocks) <= MAX_MOUNTED_BLOCKS:
            return
        for child in blocks[:TRIM_BATCH]:
            child.remove()
        self._trimmed += TRIM_BATCH
        marker = Text(
            f"· · · {self._trimmed} earlier blocks trimmed from view; "
            "/export saves the complete transcript",
            style="dim",
        )
        if self._trim_marker is None:
            self._trim_marker = Static(marker)
            self.mount(self._trim_marker, before=0)
        else:
            self._trim_marker.update(marker)


__all__ = [
    "ExportRecorder",
    "MAX_MOUNTED_BLOCKS",
    "ToolRow",
    "TranscriptBlock",
    "TranscriptView",
]
