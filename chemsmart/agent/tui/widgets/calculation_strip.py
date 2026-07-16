"""Persistent compact status for active and recent calculations."""

from __future__ import annotations

from rich.text import Text
from textual.widgets import Static

_ACTIVE = {"validating", "starting", "running"}


class CalculationStatusStrip(Static):
    DEFAULT_CSS = """
    CalculationStatusStrip {
        display: none;
        height: 1;
        padding: 0 1;
        background: $panel;
        border-left: outer $accent;
    }
    CalculationStatusStrip.visible { display: block; }
    """

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.runs: dict[str, dict[str, object]] = {}

    def update_run(self, run: dict[str, object]) -> None:
        run_id = str(run.get("run_id") or "")
        if not run_id:
            return
        self.runs[run_id] = dict(run)
        self.add_class("visible")
        self.update(self._render_status())

    def replace_runs(self, runs: list[dict[str, object]]) -> None:
        self.runs = {
            str(run.get("run_id")): dict(run)
            for run in runs
            if run.get("run_id")
        }
        self.set_class(bool(self.runs), "visible")
        self.update(self._render_status() if self.runs else "")

    def _render_status(self) -> Text:
        ordered = sorted(
            self.runs.values(),
            key=lambda run: (
                str(run.get("status")) not in _ACTIVE,
                str(run.get("updated_at") or run.get("started_at") or ""),
            ),
        )
        current = ordered[0] if ordered else {}
        status = str(current.get("status") or "unknown").upper()
        style = (
            "accent"
            if status.lower() in _ACTIVE
            else "success"
            if status == "COMPLETED"
            else "error"
        )
        program = str(current.get("program") or "calculation").upper()
        kind = str(current.get("kind") or "job").upper()
        label = str(current.get("label") or "calculation")
        stage = str(current.get("stage") or "")
        elapsed = _format_elapsed(current.get("elapsed_s"))
        active = sum(
            str(run.get("status") or "") in _ACTIVE
            for run in self.runs.values()
        )
        failed = sum(
            str(run.get("status") or "")
            in {"chemistry_failed", "process_failed", "timeout"}
            for run in self.runs.values()
        )
        text = Text.assemble(
            (status, f"bold {style}"),
            (f" · {program} {kind} · ", "dim"),
            (label, "bold"),
        )
        if stage:
            text.append(f" · {stage}", style="dim")
        if elapsed:
            text.append(f" · {elapsed}", style="dim")
        text.append(
            f" · active {active} / failed {failed} · Ctrl+B details",
            style="dim",
        )
        return text


def _format_elapsed(value: object) -> str:
    try:
        total = max(0, int(float(value)))
    except (TypeError, ValueError):
        return ""
    hours, remainder = divmod(total, 3600)
    minutes, seconds = divmod(remainder, 60)
    if hours:
        return f"{hours}:{minutes:02d}:{seconds:02d}"
    return f"{minutes:02d}:{seconds:02d}"


__all__ = ["CalculationStatusStrip"]
