"""Keyboard-first monitor for local calculations and scheduler jobs."""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from rich.console import Group
from rich.table import Table
from rich.text import Text
from textual.app import ComposeResult
from textual.binding import Binding
from textual.containers import Vertical
from textual.screen import ModalScreen
from textual.widgets import Input, OptionList, Static
from textual.widgets.option_list import Option


@dataclass(slots=True, frozen=True)
class CalculationMonitorAction:
    action: str
    run_id: str


class CalculationMonitor(ModalScreen[CalculationMonitorAction | None]):
    BINDINGS = [
        Binding("escape", "close", "Close", show=False, priority=True),
        Binding("enter", "toggle_detail", "Details", show=False),
        Binding("l", "toggle_log", "Log", show=False),
        Binding("/", "search_log", "Search log", show=False),
        Binding("e", "extract", "Extract", show=False),
        Binding("c", "cancel", "Cancel", show=False),
    ]

    DEFAULT_CSS = """
    CalculationMonitor { align: center middle; }
    #calculation-monitor {
        width: 118;
        height: 34;
        max-width: 96%;
        max-height: 92%;
        border: round $accent;
        padding: 1 2;
        background: $surface;
    }
    #calculation-list { height: 9; }
    #calculation-detail { height: 1fr; margin-top: 1; overflow-y: auto; }
    #calculation-search { display: none; height: 3; margin-top: 1; }
    #calculation-search.visible { display: block; }
    """

    def __init__(
        self,
        runs: list[dict[str, Any]],
        jobs: dict[str, dict[str, Any]] | None = None,
    ) -> None:
        super().__init__()
        self.runs = {
            str(run.get("run_id")): dict(run)
            for run in runs
            if run.get("run_id")
        }
        self.jobs = {key: dict(value) for key, value in (jobs or {}).items()}
        self._ordered_ids: list[str] = []
        self._show_log = False
        self._log_query = ""

    def compose(self) -> ComposeResult:
        with Vertical(id="calculation-monitor"):
            yield OptionList(id="calculation-list")
            yield Input(
                placeholder="Search the selected calculation log",
                id="calculation-search",
            )
            yield Static(id="calculation-detail")

    def on_mount(self) -> None:
        listing = self.query_one("#calculation-list", OptionList)
        listing.border_title = "Calculations · Enter details · L log · / search · E extract · C cancel"
        self.query_one("#calculation-detail", Static).border_title = "Receipt"
        self._refresh_list()
        listing.focus()

    def update_run(self, run: dict[str, Any]) -> None:
        run_id = str(run.get("run_id") or "")
        if not run_id:
            return
        self.runs[run_id] = dict(run)
        self._refresh_list(preserve_selection=run_id)

    def update_job(self, job_id: str, fields: dict[str, Any]) -> None:
        self.jobs.setdefault(job_id, {"job_id": job_id}).update(fields)
        self._refresh_list(preserve_selection=f"job:{job_id}")

    def on_option_list_option_highlighted(
        self, event: OptionList.OptionHighlighted
    ) -> None:
        if event.option_id:
            self._render_detail(str(event.option_id))

    def on_input_submitted(self, event: Input.Submitted) -> None:
        if event.input.id != "calculation-search":
            return
        self._log_query = event.value.strip().lower()
        event.input.remove_class("visible")
        self.query_one("#calculation-list", OptionList).focus()
        if run_id := self._selected_id():
            self._show_log = True
            self._render_detail(run_id)

    def action_close(self) -> None:
        self.dismiss(None)

    def action_toggle_detail(self) -> None:
        self._show_log = False
        if run_id := self._selected_id():
            self._render_detail(run_id)

    def action_toggle_log(self) -> None:
        self._show_log = not self._show_log
        if run_id := self._selected_id():
            self._render_detail(run_id)

    def action_search_log(self) -> None:
        search = self.query_one("#calculation-search", Input)
        search.add_class("visible")
        search.value = self._log_query
        search.focus()

    def action_extract(self) -> None:
        if run_id := self._selected_id():
            self.dismiss(CalculationMonitorAction("extract", run_id))

    def action_cancel(self) -> None:
        if run_id := self._selected_id():
            self.dismiss(CalculationMonitorAction("cancel", run_id))

    def _refresh_list(self, preserve_selection: str | None = None) -> None:
        listing = self.query_one("#calculation-list", OptionList)
        selected = preserve_selection or self._selected_id()
        rows = sorted(
            self.runs.values(),
            key=lambda run: (
                str(run.get("status"))
                not in {"validating", "starting", "running"},
                str(run.get("updated_at") or run.get("started_at") or ""),
            ),
        )
        options: list[Option] = []
        self._ordered_ids = []
        for run in rows:
            run_id = str(run.get("run_id"))
            self._ordered_ids.append(run_id)
            options.append(
                Option(
                    _run_option_text(run),
                    id=run_id,
                )
            )
        for job_id, job in sorted(self.jobs.items()):
            synthetic_id = f"job:{job_id}"
            self._ordered_ids.append(synthetic_id)
            options.append(Option(_job_option_text(job), id=synthetic_id))
        listing.clear_options()
        if options:
            listing.add_options(options)
            target = (
                self._ordered_ids.index(selected)
                if selected in self._ordered_ids
                else 0
            )
            listing.highlighted = target
            self._render_detail(self._ordered_ids[target])
        else:
            listing.add_option(
                Option("No calculations recorded.", disabled=True)
            )
            self.query_one("#calculation-detail", Static).update(
                "Run a validated command with /run or submit with /submit."
            )

    def _selected_id(self) -> str | None:
        listing = self.query_one("#calculation-list", OptionList)
        index = listing.highlighted
        if index is None or index < 0 or index >= len(self._ordered_ids):
            return None
        return self._ordered_ids[index]

    def _render_detail(self, selected_id: str) -> None:
        detail = self.query_one("#calculation-detail", Static)
        if selected_id.startswith("job:"):
            job = self.jobs.get(selected_id.removeprefix("job:"), {})
            detail.border_title = "Scheduler receipt"
            detail.update(_job_detail(job))
            return
        run = self.runs.get(selected_id, {})
        detail.border_title = (
            "Live log" if self._show_log else "Calculation receipt"
        )
        detail.update(
            _log_detail(run, self._log_query)
            if self._show_log
            else _run_detail(run)
        )


def _run_option_text(run: dict[str, Any]) -> str:
    return (
        f"{str(run.get('status') or 'unknown').upper():<18} "
        f"{str(run.get('program') or 'calc').upper():<9} "
        f"{str(run.get('kind') or 'job'):<10} "
        f"{str(run.get('label') or run.get('run_id') or ''):<24} "
        f"{_elapsed(run.get('elapsed_s')):<8} "
        f"{str(run.get('stage') or '')}"
    )


def _job_option_text(job: dict[str, Any]) -> str:
    return (
        f"{str(job.get('status') or 'unknown').upper():<18} "
        f"{str(job.get('scheduler') or 'job').upper():<9} "
        f"{'scheduler':<10} {str(job.get('name') or ''):<24} "
        f"{str(job.get('runtime') or ''):<8} {str(job.get('host') or '')}"
    )


def _run_detail(run: dict[str, Any]):
    table = Table.grid(expand=True)
    table.add_column(style="dim", width=20)
    table.add_column(ratio=1)
    for label, key in (
        ("run", "run_id"),
        ("status", "status"),
        ("stage", "stage"),
        ("program / kind", "program"),
        ("project", "project"),
        ("method", "method"),
        ("basis", "basis"),
        ("input", "input_path"),
        ("output", "output_path"),
        ("PID", "pid"),
        ("elapsed", "elapsed_s"),
        ("return code", "returncode"),
        ("semantic gate", "semantic_verdict"),
        ("intent gate", "intent_verdict"),
    ):
        value = run.get(key)
        if value not in {None, ""}:
            table.add_row(
                label, _elapsed(value) if key == "elapsed_s" else str(value)
            )
    if run.get("reused_output") is True:
        table.add_row("output source", "existing completed calculation")
    energy = run.get("energy")
    if isinstance(energy, (int, float)):
        table.add_row("E(SCF)", f"{float(energy):.12f} Eh")
    chemistry_elapsed = run.get("chemistry_elapsed_s")
    if isinstance(chemistry_elapsed, (int, float)):
        table.add_row("program runtime", f"{float(chemistry_elapsed):.3f} s")
    cycles = run.get("scf_cycles")
    if isinstance(cycles, int):
        table.add_row("SCF cycles", str(cycles))
    opt_cycles = run.get("optimization_cycles")
    if isinstance(opt_cycles, int):
        table.add_row(
            "optimization",
            f"{opt_cycles} cycles · "
            + (
                "converged"
                if run.get("optimization_converged") is True
                else "not confirmed"
            ),
        )
    imag_count = run.get("imaginary_frequency_count")
    if isinstance(imag_count, int):
        values = list(run.get("imag_freqs") or [])
        table.add_row(
            "imag frequencies",
            str(imag_count)
            + (
                " · "
                + ", ".join(f"{float(value):.1f}" for value in values)
                + " cm^-1"
                if values
                else ""
            ),
        )
    scan_points = run.get("scan_points_completed")
    if isinstance(scan_points, int):
        total = run.get("scan_points_total")
        table.add_row(
            "scan points",
            (
                f"{scan_points}/{total}"
                if isinstance(total, int)
                else str(scan_points)
            ),
        )
    neb_images = run.get("neb_images")
    if isinstance(neb_images, int):
        table.add_row("NEB images", str(neb_images))
    qm_atoms = run.get("qmmm_qm_atoms")
    mm_atoms = run.get("qmmm_mm_atoms")
    if isinstance(qm_atoms, int) or isinstance(mm_atoms, int):
        table.add_row(
            "QMMM regions",
            f"QM {qm_atoms if isinstance(qm_atoms, int) else '?'} · "
            f"MM {mm_atoms if isinstance(mm_atoms, int) else '?'}",
        )
    normal = run.get("normal_termination")
    if normal is not None:
        table.add_row("normal termination", "yes" if normal else "no")
    command = str(run.get("command") or "")
    return Group(
        table, Text("\n" + command, style="dim") if command else Text()
    )


def _log_detail(run: dict[str, Any], query: str) -> Text:
    paths = [
        Path(str(run.get("output_path"))) if run.get("output_path") else None,
        Path(str(run.get("stderr_path"))) if run.get("stderr_path") else None,
        Path(str(run.get("stdout_path"))) if run.get("stdout_path") else None,
    ]
    text = ""
    source = ""
    for path in paths:
        if path is None or not path.is_file():
            continue
        source = str(path)
        text = _bounded_tail(path)
        if text:
            break
    if not text:
        return Text("No log output is available yet.", style="dim")
    lines = text.splitlines()
    if query:
        lines = [line for line in lines if query in line.lower()]
    header = f"{source}\n" + (f"filter: {query}\n" if query else "")
    return Text(header + "\n".join(lines[-200:]))


def _job_detail(job: dict[str, Any]):
    table = Table.grid(expand=True)
    table.add_column(style="dim", width=18)
    table.add_column(ratio=1)
    for key in (
        "job_id",
        "name",
        "scheduler",
        "status",
        "server_name",
        "input_path",
        "output_path",
        "runtime",
    ):
        value = job.get(key)
        if value not in {None, ""}:
            table.add_row(key, str(value))
    return table


def _bounded_tail(path: Path, limit: int = 65_536) -> str:
    try:
        with path.open("rb") as handle:
            handle.seek(0, os.SEEK_END)
            size = handle.tell()
            handle.seek(max(0, size - limit))
            payload = handle.read(limit)
    except OSError as exc:
        return f"Failed to read {path}: {exc}"
    return payload.decode("utf-8", errors="replace")


def _elapsed(value: object) -> str:
    try:
        total = max(0, int(float(value)))
    except (TypeError, ValueError):
        return ""
    hours, remainder = divmod(total, 3600)
    minutes, seconds = divmod(remainder, 60)
    if hours:
        return f"{hours}:{minutes:02d}:{seconds:02d}"
    return f"{minutes:02d}:{seconds:02d}"


__all__ = ["CalculationMonitor", "CalculationMonitorAction"]
