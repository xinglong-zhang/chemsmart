"""Plan/workflow transcript cell."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from rich.console import Group
from rich.text import Text
from textual.binding import Binding
from textual.events import Click

from chemsmart.agent.core import Plan, Step, _restore_json_result
from chemsmart.io.molecules.structure import Molecule

from .base import BaseCell, no_data_text


@dataclass(slots=True)
class _WorkflowRow:
    index: int
    tool: str
    args: dict[str, Any]
    detail: str
    status: str = "pending"


_STATUS_ICON = {
    "done": "✓",
    "running": "⠋",
    "failed": "✖",
    "pending": "•",
}


class PlanCell(BaseCell):
    BINDINGS = [
        Binding("enter", "toggle_expand", "Expand", show=False),
        Binding("space", "toggle_expand", "Expand", show=False),
    ]

    def __init__(self, content: Plan | str) -> None:
        self.plan = content if isinstance(content, Plan) else None
        self.plan_text = content if isinstance(content, str) else None
        self._rows = (
            [
                _row_from_step(index, step)
                for index, step in enumerate(
                    self.plan.steps,
                    start=1,
                )
            ]
            if self.plan is not None
            else []
        )
        self._last_changed: _WorkflowRow | None = None
        self.expanded = False
        super().__init__(
            self._build_renderable(),
            title="Workflow" if self.plan is not None else "Plan",
            classes="plan-cell",
        )

    def mark_started(
        self,
        step_index: int,
        tool: str,
        args: dict[str, Any] | None = None,
    ) -> None:
        row = self._find_row(step_index, tool)
        if row is None:
            return
        row.status = "running"
        row.detail = _running_detail(tool, args or row.args)
        self._last_changed = row
        self.update(self._build_renderable())

    def mark_completed(
        self,
        step_index: int,
        tool: str,
        payload: Any | None = None,
    ) -> None:
        row = self._find_row(step_index, tool)
        if row is None:
            return
        row.status = "done"
        row.detail = _completed_detail(tool, row.args, payload)
        self._last_changed = row
        self.update(self._build_renderable())

    def mark_failed(
        self,
        step_index: int,
        tool: str,
        message: str,
    ) -> None:
        row = self._find_row(step_index, tool)
        if row is None:
            return
        row.status = "failed"
        row.detail = message.strip() or "Execution was interrupted."
        self._last_changed = row
        self.expanded = True
        self.update(self._build_renderable())

    def on_click(self, _event: Click) -> None:
        self.action_toggle_expand()

    def action_toggle_expand(self) -> None:
        self.set_expanded(not self.expanded)

    def set_expanded(self, expanded: bool) -> None:
        self.expanded = bool(expanded)
        self.update(self._build_renderable())

    def _find_row(self, step_index: int, tool: str) -> _WorkflowRow | None:
        for row in self._rows:
            if row.index == step_index and row.tool == tool:
                return row
        return None

    def _build_renderable(self):
        if self.plan is None:
            text = self.plan_text or ""
            return Text(text) if text.strip() else no_data_text()
        if not self._rows:
            if self.plan.rationale:
                return Text(self.plan.rationale)
            return no_data_text()
        completed = sum(row.status == "done" for row in self._rows)
        current = next(
            (row for row in self._rows if row.status == "running"),
            next((row for row in self._rows if row.status == "pending"), None),
        )
        failed = next(
            (row for row in self._rows if row.status == "failed"), None
        )
        active = failed or current
        state = "failed" if failed is not None else (
            active.status if active is not None else "done"
        )
        header = Text.assemble(
            ("Workflow", "bold"),
            (f" · {completed}/{len(self._rows)}", "dim"),
            (
                f" · {active.tool} {state}" if active is not None else " · complete",
                _status_style(state),
            ),
            ("  [enter expand]", "dim"),
        )
        if not self.expanded:
            display = failed or (
                active if active is not None and active.status == "running" else None
            ) or self._last_changed or active
            if display is None:
                return Group(header, Text("All steps complete", style="success"))
            icon = _STATUS_ICON.get(display.status, "•")
            line = Text()
            line.append(f"{icon} ", style=_status_style(display.status))
            line.append(f"{display.index}. {display.tool}", style="bold")
            if display.detail:
                line.append("  ", style="dim")
                line.append(display.detail, style=_status_style(display.status))
            return Group(header, line)

        lines = [header]
        for row in self._rows:
            icon = _STATUS_ICON.get(row.status, "•")
            line = Text()
            style = _status_style(row.status)
            line.append(f"{icon} ", style=style)
            line.append(f"{row.index}. {row.tool}", style="bold")
            if row.detail:
                line.append("  ", style="dim")
                line.append(row.detail, style=style)
            lines.append(line)
        return Group(*lines)


def _row_from_step(index: int, step: Step) -> _WorkflowRow:
    return _WorkflowRow(
        index=index,
        tool=step.tool,
        args=dict(step.args or {}),
        detail=_pending_detail(step.tool, dict(step.args or {})),
    )


def _status_style(status: str) -> str:
    if status == "done":
        return "success"
    if status == "running":
        return "accent"
    if status == "failed":
        return "error"
    return "dim"


def _pending_detail(tool: str, args: dict[str, Any]) -> str:
    if tool == "build_molecule":
        path = str(args.get("filepath") or "")
        return f"Waiting on {path}" if path else "Waiting on structure file"
    if tool == "dry_run_input":
        return "Waiting to render input file"
    if tool == "validate_runtime":
        return "Waiting on runtime check"
    if tool == "run_local":
        return "Waiting to run"
    if tool == "submit_hpc":
        return "Waiting on submission preview"
    return "Waiting"


def _running_detail(tool: str, args: dict[str, Any]) -> str:
    if tool == "build_molecule":
        path = str(args.get("filepath") or "")
        return f"Checking {path}…" if path else "Checking structure file…"
    if tool == "build_gaussian_settings":
        return "Resolving Gaussian settings…"
    if tool == "build_job":
        return "Assembling calculation job…"
    if tool == "dry_run_input":
        return "Rendering input file…"
    if tool == "validate_runtime":
        return "Running runtime check…"
    if tool == "run_local":
        return "Starting local run…"
    if tool == "submit_hpc":
        return "Generating submission preview…"
    return "In progress…"


def _completed_detail(
    tool: str,
    args: dict[str, Any],
    payload: Any | None,
) -> str:
    restored = _restore_json_result(payload)
    if tool == "build_molecule" and isinstance(restored, Molecule):
        formula = (
            restored.formula() if hasattr(restored, "formula") else "Molecule"
        )
        count = len(restored.symbols or [])
        charge = restored.charge if restored.charge is not None else 0
        multiplicity = (
            restored.multiplicity if restored.multiplicity is not None else 1
        )
        source = Path(str(args.get("filepath") or "")).as_posix()
        details = f"{formula} · {count} atoms · q={charge} mult={multiplicity}"
        return f"{details} · {source}" if source else details
    if tool == "build_gaussian_settings" and isinstance(restored, dict):
        functional = str(restored.get("functional") or "").upper()
        basis = str(restored.get("basis") or "").upper()
        return f"{functional}/{basis}".strip("/")
    if tool == "build_gaussian_settings":
        functional = str(getattr(restored, "functional", "") or "").upper()
        basis = str(getattr(restored, "basis", "") or "").upper()
        return f"{functional}/{basis}".strip("/")
    if tool == "build_job":
        job = restored
        kind = str(args.get("kind") or "")
        label = str(getattr(job, "label", "") or args.get("label") or "")
        return " · ".join(part for part in (kind, label) if part)
    if tool == "dry_run_input" and isinstance(restored, dict):
        inputfile = str(restored.get("inputfile") or "")
        return Path(inputfile).name if inputfile else "Input file ready"
    if tool == "validate_runtime" and isinstance(restored, dict):
        return _runtime_summary(restored)
    if tool == "run_local" and isinstance(restored, dict):
        if restored.get("ok"):
            return "Local run complete"
        return f"returncode {restored.get('returncode', '?')}"
    if tool == "submit_hpc" and isinstance(restored, dict):
        if restored.get("job_id"):
            return f"job {restored['job_id']}"
        return "Submission preview ready"
    return "Done"


def _runtime_summary(validation: dict[str, Any]) -> str:
    remote_unknown = validation.get("remote_unknown") or []
    local_issues = validation.get("local_issues") or []
    if local_issues:
        return f"{len(local_issues)} local issue(s)"
    if remote_unknown:
        return f"Local OK / {len(remote_unknown)} remote item(s) needed"
    return "Local/remote run ready"
