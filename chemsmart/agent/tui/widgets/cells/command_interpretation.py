"""Deterministic interpretation cell for model-created commands."""

from __future__ import annotations

from rich.console import Group
from rich.table import Table
from rich.text import Text
from textual.events import Click

from chemsmart.agent.model_command_parser import (
    ParsedModelCommand,
    format_parsed_model_command,
)

from .base import BaseCell


class CommandInterpretationCell(BaseCell):
    """Collapsible deterministic grounding for a synthesized command.

    Collapsed by default so the agent's composed answer stays the focus; the
    full deterministic fact table is available on demand for user audit.
    """

    def __init__(self, parsed: ParsedModelCommand) -> None:
        self.parsed = parsed
        self.expanded = False
        self.source_text = format_parsed_model_command(parsed)
        super().__init__(
            _render_collapsed_interpretation(parsed),
            title="Command Interpretation ▸",
            classes="agent-cell command-interpretation-cell",
        )

    def on_click(self, _event: Click) -> None:
        self.toggle()

    def action_toggle(self) -> None:
        self.toggle()

    def toggle(self) -> None:
        self.expanded = not self.expanded
        self.border_title = (
            "Command Interpretation ▾"
            if self.expanded
            else "Command Interpretation ▸"
        )
        self.update(
            _render_command_interpretation(self.parsed)
            if self.expanded
            else _render_collapsed_interpretation(self.parsed)
        )


def _render_collapsed_interpretation(parsed: ParsedModelCommand) -> Group:
    if parsed.parse_error:
        return Group(
            Text.assemble(
                ("grounded facts: ", "dim"),
                ("parse error", "bold red"),
            ),
            Text(parsed.parse_error, style="red"),
            Text("click to show details", style="dim"),
        )
    headline = Text.assemble(
        ("grounded facts: ", "dim"),
        (parsed.program or "unknown", "bold cyan"),
        (" · ", "dim"),
        (_job_label(parsed), "bold"),
        (" · ", "dim"),
        (parsed.filename or "runtime input", "green"),
    )
    return Group(headline, Text("click to show grounded facts", style="dim"))


def _render_command_interpretation(parsed: ParsedModelCommand) -> Group:
    if parsed.parse_error:
        return Group(
            Text("Deterministic command parser", style="bold"),
            Text(f"Parse error: {parsed.parse_error}", style="red"),
            Text(f"Workspace: {parsed.workspace}", style="dim"),
            Text(
                f"Summary: This command could not be deterministically parsed: {parsed.parse_error}.",
                style="bold yellow",
            ),
        )

    table = Table.grid(padding=(0, 2))
    table.add_column("Field", style="cyan", no_wrap=True)
    table.add_column("Value", overflow="fold")

    table.add_row("workspace", parsed.workspace)
    table.add_row("execution", _execution_label(parsed))
    table.add_row("program", parsed.program or "unknown")
    table.add_row("job", _job_label(parsed))
    table.add_row("server", _server_label(parsed))
    table.add_row("dry run", "yes" if parsed.dry_run else "no")
    table.add_row("input file", parsed.filename or "not specified")
    table.add_row("label", parsed.label or "runtime-derived")
    table.add_row(
        "charge / mult",
        f"{parsed.charge or 'runtime/default'} / {parsed.multiplicity or 'runtime/default'}",
    )
    table.add_row("project", parsed.project or "not specified")
    if parsed.project_p_flag_meaning:
        table.add_row("-p meaning", parsed.project_p_flag_meaning)
    if parsed.top_level_program:
        table.add_row(
            "top-level -p",
            f"{parsed.top_level_program} (output-file processing target, not project)",
        )
    table.add_row("method", _method_label(parsed))
    solvent = _solvent_label(parsed)
    if solvent:
        table.add_row("solvent", solvent)
    route_controls = _route_controls_label(parsed)
    if route_controls:
        table.add_row("route controls", route_controls)
    if parsed.route_parameters:
        table.add_row("route params", parsed.route_parameters)
    if parsed.opt_options:
        table.add_row("opt route opts", parsed.opt_options)
    if parsed.resources:
        table.add_row("resources", _dict_label(parsed.resources))
    if parsed.structural_options:
        table.add_row("job options", _dict_label(parsed.structural_options))

    parts: list[object] = [
        Text("Deterministic command parser", style="bold"),
        table,
    ]
    if parsed.warnings:
        warning_table = Table.grid(padding=(0, 1))
        warning_table.add_column("Marker", style="yellow", no_wrap=True)
        warning_table.add_column("Warning", overflow="fold")
        for warning in parsed.warnings:
            warning_table.add_row("warning", warning)
        parts.append(warning_table)
    parts.append(Text(_summary_sentence(parsed), style="bold green"))
    return Group(*parts)


def _execution_label(parsed: ParsedModelCommand) -> str:
    if parsed.action == "sub":
        return "sub (submit to an HPC/server queue)"
    if parsed.action == "run":
        return "run (run locally)"
    return parsed.action or "unknown"


def _job_label(parsed: ParsedModelCommand) -> str:
    if parsed.job is None:
        return "unknown"
    labels = {
        "sp": "single-point energy",
        "opt": "geometry optimization",
        "ts": "transition-state search",
        "irc": "intrinsic reaction coordinate",
        "scan": "coordinate scan",
        "modred": "constrained optimization",
        "td": "TD-DFT excited-state calculation",
        "neb": "nudged elastic band",
        "qmmm": "QM/MM calculation",
    }
    return f"{parsed.job} ({labels.get(parsed.job, parsed.job)})"


def _server_label(parsed: ParsedModelCommand) -> str:
    if parsed.server:
        return parsed.server
    if parsed.action == "sub":
        return "auto/default server"
    return "local/default"


def _method_label(parsed: ParsedModelCommand) -> str:
    parts = []
    if parsed.ab_initio:
        parts.append(f"ab initio {parsed.ab_initio}")
    if parsed.functional:
        parts.append(f"functional {parsed.functional}")
    if parsed.basis:
        parts.append(f"basis {parsed.basis}")
    if parsed.aux_basis:
        parts.append(f"auxiliary basis {parsed.aux_basis}")
    if parsed.extrapolation_basis:
        parts.append(f"extrapolation basis {parsed.extrapolation_basis}")
    return ", ".join(parts) if parts else "unresolved"


def _solvent_label(parsed: ParsedModelCommand) -> str | None:
    if not parsed.solvent_model and not parsed.solvent_id:
        return None
    return f"{parsed.solvent_model or 'model default'} / {parsed.solvent_id or 'id default'}"


def _route_controls_label(parsed: ParsedModelCommand) -> str | None:
    parts = []
    if parsed.defgrid:
        parts.append(f"defgrid={parsed.defgrid}")
    if parsed.scf_tol:
        parts.append(f"scf_tol={parsed.scf_tol}")
    if parsed.scf_algorithm:
        parts.append(f"scf_algorithm={parsed.scf_algorithm}")
    return ", ".join(parts) if parts else None


def _dict_label(values: dict[str, str]) -> str:
    return ", ".join(f"{key}={value}" for key, value in sorted(values.items()))


def _summary_sentence(parsed: ParsedModelCommand) -> str:
    text = format_parsed_model_command(parsed)
    for line in reversed(text.splitlines()):
        if line.startswith("Summary:"):
            return line
    return "Summary: This command was parsed deterministically."

