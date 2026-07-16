"""File picker overlay for @-insertion."""

from __future__ import annotations

from concurrent.futures import Future
from pathlib import Path

from textual import work
from textual.app import ComposeResult
from textual.binding import Binding
from textual.containers import Vertical
from textual.screen import ModalScreen
from textual.widgets import ListItem, ListView, Static

from chemsmart.agent.tui.services.file_index import (
    iter_candidate_files,
    request_candidate_file_refresh,
)


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
        self.cwd = Path(cwd).resolve()
        self._candidates = iter_candidate_files(self.cwd)
        self._refresh_future = request_candidate_file_refresh(self.cwd)

    def compose(self) -> ComposeResult:
        with Vertical(id="file-picker-modal"):
            summary = Static(
                "Insert a nearby file path. Preferred: .xyz/.log/.com/.inp/.gjf/.out",
                id="file-picker-summary",
            )
            summary.border_title = "Files"
            yield summary
            yield ListView(*self._candidate_items(), id="file-picker-list")

    def on_mount(self) -> None:
        self.query_one("#file-picker-list", ListView).focus()
        if self._refresh_future is not None:
            self._await_candidate_refresh(self._refresh_future)

    def _candidate_items(self) -> list[ListItem]:
        if self._candidates:
            return [
                ListItem(Static(str(path.relative_to(self.cwd))))
                for path in self._candidates[:20]
            ]
        message = (
            "Indexing workspace files..."
            if self._refresh_future is not None
            else "No candidate files found."
        )
        return [ListItem(Static(message))]

    @work(
        thread=True,
        exclusive=True,
        exit_on_error=False,
        group="file-picker-index",
        name="file-picker-index",
    )
    def _await_candidate_refresh(self, future: Future) -> None:
        try:
            candidates = list(future.result())
        except Exception:
            candidates = []
        self.app.call_from_thread(self._apply_candidates, candidates)

    def _apply_candidates(self, candidates: list[Path]) -> None:
        if not self.is_mounted:
            return
        self._candidates = candidates
        self._refresh_future = None
        self._render_candidate_items()

    @work(
        exclusive=True,
        exit_on_error=False,
        group="file-picker-render",
        name="file-picker-render",
    )
    async def _render_candidate_items(self) -> None:
        list_view = self.query_one("#file-picker-list", ListView)
        await list_view.clear()
        await list_view.extend(self._candidate_items())
        list_view.index = 0
        list_view.focus()

    def action_select_current(self) -> None:
        if not self._candidates:
            self.dismiss(None)
            return
        list_view = self.query_one("#file-picker-list", ListView)
        index = list_view.index or 0
        index = max(0, min(index, len(self._candidates) - 1))
        self.dismiss(str(self._candidates[index].relative_to(self.cwd)))

    def on_list_view_selected(self, event: ListView.Selected) -> None:
        event.stop()
        self.action_select_current()

    def action_cancel(self) -> None:
        self.dismiss(None)
