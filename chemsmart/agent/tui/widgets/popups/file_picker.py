"""File picker overlay for @-insertion."""

from __future__ import annotations

from pathlib import Path

from textual.app import ComposeResult
from textual.binding import Binding
from textual.containers import Vertical
from textual.screen import ModalScreen
from textual.widgets import ListItem, ListView, Static

from chemsmart.agent.tui.services.file_index import iter_candidate_files


class FilePickerOverlay(ModalScreen[str | None]):
    BINDINGS = [
        Binding("escape", "cancel", "Cancel", show=False),
        Binding("enter", "select_current", "Select", show=False),
    ]

    DEFAULT_CSS = """
    FilePickerOverlay {
        align: center middle;
    }

    #file-picker-modal {
        width: 88;
        height: 24;
        border: round $secondary;
        padding: 1 2;
        background: $surface;
    }

    #file-picker-list {
        height: 1fr;
        margin-top: 1;
    }
    """

    def __init__(self, cwd: str | Path) -> None:
        super().__init__()
        self.cwd = Path(cwd)
        self._candidates = iter_candidate_files(self.cwd)

    def compose(self) -> ComposeResult:
        with Vertical(id="file-picker-modal"):
            summary = Static(
                "Insert a nearby file path. Preferred: .xyz/.log/.com/.inp/.gjf/.out",
                id="file-picker-summary",
            )
            summary.border_title = "Files"
            yield summary
            items = [
                ListItem(Static(str(path.relative_to(self.cwd))))
                for path in self._candidates[:20]
            ] or [ListItem(Static("No candidate files found."))]
            yield ListView(*items, id="file-picker-list")

    def on_mount(self) -> None:
        self.query_one("#file-picker-list", ListView).focus()

    def action_select_current(self) -> None:
        if not self._candidates:
            self.dismiss(None)
            return
        list_view = self.query_one("#file-picker-list", ListView)
        index = list_view.index or 0
        index = max(0, min(index, len(self._candidates) - 1))
        self.dismiss(str(self._candidates[index].relative_to(self.cwd)))

    def action_cancel(self) -> None:
        self.dismiss(None)
