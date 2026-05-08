"""Inline job status cell."""

from __future__ import annotations

from rich.console import Group
from rich.table import Table
from rich.text import Text

from .base import BaseCell

_PROGRESS = {
    "queued": 0.2,
    "running": 0.6,
    "done": 1.0,
    "failed": 1.0,
    "cancelled": 1.0,
}
_STATUS_STYLE = {
    "queued": "warning",
    "running": "accent",
    "done": "success",
    "failed": "error",
    "cancelled": "dim",
}


class JobStatusCell(BaseCell):
    def __init__(self, job_id: str, snapshot: dict) -> None:
        self.job_id = job_id
        self.snapshot = dict(snapshot)
        super().__init__(
            self._build_renderable(),
            title="Job status",
            classes="job-status-cell",
        )

    def apply_update(self, fields: dict) -> None:
        self.snapshot.update(fields)
        self.update(self._build_renderable())

    def _build_renderable(self):
        status = str(self.snapshot.get("status") or "queued")
        name = str(self.snapshot.get("name") or self.job_id)
        progress = _render_progress(_PROGRESS.get(status, 0.0), status)
        heading = Text(f"{name}  [{self.job_id}]", style="bold")
        heading.append(f"  {status}", style=_STATUS_STYLE.get(status, "text"))
        table = Table.grid(expand=True)
        table.add_column(ratio=2)
        table.add_column(ratio=3)
        table.add_row("scheduler", str(self.snapshot.get("scheduler") or ""))
        table.add_row("started", str(self.snapshot.get("started") or ""))
        table.add_row("runtime", str(self.snapshot.get("runtime") or ""))
        table.add_row("host", str(self.snapshot.get("host") or ""))
        return Group(heading, progress, table)


def _render_progress(fraction: float, status: str, width: int = 24) -> Text:
    filled = max(0, min(width, round(width * fraction)))
    empty = width - filled
    bar = Text("█" * filled, style=_STATUS_STYLE.get(status, "accent"))
    bar.append("░" * empty, style="dim")
    return bar
