"""Jobs overlay for active and recent jobs."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from rich.table import Table
from textual.app import ComposeResult
from textual.binding import Binding
from textual.containers import Vertical
from textual.screen import ModalScreen
from textual.widgets import Static

from chemsmart.agent.tui.services.job_poller import JobStatusUpdated


@dataclass(slots=True)
class JobsPanelAction:
    action: str
    job_id: str | None = None


class JobsPanel(ModalScreen[JobsPanelAction | None]):
    BINDINGS = [
        Binding("escape", "close", "Close", show=False, priority=True),
        Binding("up", "cursor_up", "Up", show=False),
        Binding("down", "cursor_down", "Down", show=False),
        Binding("enter", "jump_to_session", "Resume", show=False),
        Binding("c", "cancel_selected", "Cancel", show=False),
        Binding("e", "extract_selected", "Extract", show=False),
    ]

    DEFAULT_CSS = """
    JobsPanel {
        align: center middle;
    }

    #jobs-panel-modal {
        width: 112;
        height: auto;
        max-height: 28;
        border: round $primary;
        padding: 1 2;
        background: $surface;
    }
    """

    def __init__(self, jobs: dict[str, dict[str, Any]]) -> None:
        super().__init__()
        self.jobs = {
            job_id: dict(snapshot) for job_id, snapshot in jobs.items()
        }
        self._ordered_ids: list[str] = []
        self._selected_index = 0

    def compose(self) -> ComposeResult:
        with Vertical(id="jobs-panel-modal"):
            summary = Static(id="jobs-panel-summary")
            summary.border_title = "Jobs"
            yield summary

    def on_mount(self) -> None:
        self._refresh_table()

    def on_job_status_updated(self, message: JobStatusUpdated) -> None:
        if message.job_id not in self.jobs:
            self.jobs[message.job_id] = {"job_id": message.job_id}
        self.jobs[message.job_id].update(message.fields)
        self._refresh_table()

    def action_close(self) -> None:
        self.dismiss(None)

    def action_cursor_up(self) -> None:
        if self._ordered_ids:
            self._selected_index = max(0, self._selected_index - 1)
            self._refresh_table()

    def action_cursor_down(self) -> None:
        if self._ordered_ids:
            self._selected_index = min(
                len(self._ordered_ids) - 1,
                self._selected_index + 1,
            )
            self._refresh_table()

    def action_cancel_selected(self) -> None:
        if job_id := self._selected_job_id():
            self.dismiss(JobsPanelAction("cancel", job_id))

    def action_extract_selected(self) -> None:
        if job_id := self._selected_job_id():
            self.dismiss(JobsPanelAction("extract", job_id))

    def action_jump_to_session(self) -> None:
        if job_id := self._selected_job_id():
            self.dismiss(JobsPanelAction("resume", job_id))

    def _selected_job_id(self) -> str | None:
        if not self._ordered_ids:
            return None
        self._selected_index = max(
            0, min(self._selected_index, len(self._ordered_ids) - 1)
        )
        return self._ordered_ids[self._selected_index]

    def _refresh_table(self) -> None:
        rows = sorted(self.jobs.values(), key=_job_sort_key)
        self._ordered_ids = [str(row.get("job_id")) for row in rows]
        if self._selected_index >= len(self._ordered_ids):
            self._selected_index = max(0, len(self._ordered_ids) - 1)
        table = Table(show_header=True, header_style="bold")
        for column in [
            "job_id",
            "name",
            "scheduler",
            "status",
            "started",
            "runtime",
            "host",
        ]:
            table.add_column(column)
        if not rows:
            table.add_row("—", "No jobs found", "", "", "", "", "")
        for index, row in enumerate(rows):
            style = "reverse" if index == self._selected_index else ""
            table.add_row(
                str(row.get("job_id") or ""),
                str(row.get("name") or ""),
                str(row.get("scheduler") or ""),
                str(row.get("status") or ""),
                str(row.get("started") or ""),
                str(row.get("runtime") or ""),
                str(row.get("host") or ""),
                style=style,
            )
        self.query_one("#jobs-panel-summary", Static).update(table)


def _job_sort_key(job: dict[str, Any]) -> tuple[int, str, str]:
    status_order = {
        "running": 0,
        "queued": 1,
        "failed": 2,
        "cancelled": 3,
        "done": 4,
    }
    return (
        status_order.get(str(job.get("status")), 9),
        str(job.get("raw_started") or ""),
        str(job.get("job_id") or ""),
    )
