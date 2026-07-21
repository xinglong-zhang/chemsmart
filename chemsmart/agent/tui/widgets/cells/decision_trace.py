"""Collapsible public decision trace for API-routed synthesis turns."""

from __future__ import annotations

from rich.console import Group
from rich.table import Table
from rich.text import Text
from textual.binding import Binding
from textual.events import Click

from .base import BaseCell


class DecisionTraceCell(BaseCell):
    """Render a compact, user-auditable routing trace.

    This cell intentionally does not expose hidden chain-of-thought. It shows
    the frontier provider's public routing decision, observable evidence, and
    rejected action categories so a user can audit why the agent synthesized,
    explained, critiqued, repaired, or asked for clarification.
    """

    BINDINGS = [
        Binding("enter", "toggle", "Expand", show=False),
        Binding("space", "toggle", "Expand", show=False),
    ]

    def __init__(self, trace: dict[str, object]) -> None:
        self.trace = trace
        self.expanded = False
        self.source_text = _format_trace(trace)
        super().__init__(
            _render_trace(trace, expanded=False),
            title="Decision Trace ▸",
            classes="agent-cell decision-trace-cell",
        )

    def on_click(self, _event: Click) -> None:
        self.toggle()

    def action_toggle(self) -> None:
        self.toggle()

    def toggle(self) -> None:
        self.set_expanded(not self.expanded)

    def set_expanded(self, expanded: bool) -> None:
        self.expanded = bool(expanded)
        self.border_title = (
            "Decision Trace ▾" if self.expanded else "Decision Trace ▸"
        )
        self.update(_render_trace(self.trace, expanded=self.expanded))


def _render_trace(trace: dict[str, object], *, expanded: bool) -> Group:
    action = str(trace.get("action") or "unknown")
    confidence = str(trace.get("confidence") or "unknown")
    summary = str(trace.get("decision_summary") or "No summary.")
    target = str(trace.get("target_command") or "")
    headline = Text.assemble(
        ("routing: ", "dim"),
        (action, "bold cyan"),
        (" · confidence: ", "dim"),
        (confidence, "bold"),
    )
    hint = Text(
        "click to hide details" if expanded else "click to show evidence",
        style="dim",
    )
    if not expanded:
        return Group(headline, Text(summary), hint)

    table = Table.grid(padding=(0, 2))
    table.add_column("Field", style="cyan", no_wrap=True)
    table.add_column("Value", overflow="fold")
    table.add_row("router", str(trace.get("router") or "unknown"))
    table.add_row("action", action)
    table.add_row("confidence", confidence)
    table.add_row("summary", summary)
    table.add_row("target command", target or "none")
    table.add_row(
        "default project", str(trace.get("default_project") or "none")
    )
    table.add_row(
        "memory",
        (
            "previous command available"
            if bool(trace.get("last_command_available"))
            else "no previous command"
        ),
    )
    table.add_row("request", str(trace.get("request_excerpt") or ""))

    reasoning_table = Table.grid(padding=(0, 1))
    reasoning_table.add_column("Reasoning", style="magenta", no_wrap=True)
    reasoning_table.add_column("Value", overflow="fold")
    reasoning_steps = _as_list(trace.get("reasoning"))
    for index, item in enumerate(reasoning_steps, start=1):
        reasoning_table.add_row(f"{index}.", item)

    evidence_table = Table.grid(padding=(0, 1))
    evidence_table.add_column("Evidence", style="cyan", no_wrap=True)
    evidence_table.add_column("Value", overflow="fold")
    for item in _as_list(trace.get("evidence")):
        evidence_table.add_row("-", item)

    caveat_table = Table.grid(padding=(0, 1))
    caveat_table.add_column("Caveat", style="yellow", no_wrap=True)
    caveat_table.add_column("Value", overflow="fold")
    for item in _as_list(trace.get("caveats")):
        caveat_table.add_row("!", item)

    rejected_table = Table.grid(padding=(0, 1))
    rejected_table.add_column("Rejected", style="yellow", no_wrap=True)
    rejected_table.add_column("Reason", overflow="fold")
    rejected = trace.get("rejected_actions")
    if isinstance(rejected, dict) and rejected:
        for key, value in sorted(rejected.items()):
            rejected_table.add_row(str(key), str(value))
    else:
        rejected_table.add_row("none", "No rejected action was recorded.")

    note = Text(str(trace.get("note") or ""), style="dim")
    sections: list[object] = [headline, table]
    if reasoning_steps:
        sections.append(Text("reasoning:", style="bold magenta"))
        sections.append(reasoning_table)
    sections.append(evidence_table)
    if _as_list(trace.get("caveats")):
        sections.append(Text("caveats:", style="bold yellow"))
        sections.append(caveat_table)
    sections.extend([rejected_table, note, hint])
    return Group(*sections)


def _format_trace(trace: dict[str, object]) -> str:
    lines = [
        "public decision trace:",
        f"- router: `{trace.get('router') or 'unknown'}`",
        f"- action: `{trace.get('action') or 'unknown'}`",
        f"- confidence: `{trace.get('confidence') or 'unknown'}`",
        f"- summary: {trace.get('decision_summary') or 'No summary.'}",
        f"- target command: `{trace.get('target_command') or 'none'}`",
        f"- default project: `{trace.get('default_project') or 'none'}`",
    ]
    reasoning = _as_list(trace.get("reasoning"))
    if reasoning:
        lines.append("- reasoning:")
        lines.extend(f"  - {item}" for item in reasoning)
    evidence = _as_list(trace.get("evidence"))
    if evidence:
        lines.append("- evidence:")
        lines.extend(f"  - {item}" for item in evidence)
    caveats = _as_list(trace.get("caveats"))
    if caveats:
        lines.append("- caveats:")
        lines.extend(f"  - {item}" for item in caveats)
    rejected = trace.get("rejected_actions")
    if isinstance(rejected, dict) and rejected:
        lines.append("- rejected actions:")
        lines.extend(
            f"  - {key}: {value}" for key, value in sorted(rejected.items())
        )
    if trace.get("note"):
        lines.append(f"- note: {trace['note']}")
    return "\n".join(lines)


def _as_list(value: object) -> list[str]:
    if not isinstance(value, list):
        return []
    return [str(item) for item in value if str(item).strip()]
