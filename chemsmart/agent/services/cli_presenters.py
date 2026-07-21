"""Rich presenters and compact projections for ``chemsmart agent`` CLI."""

from __future__ import annotations

import inspect
import json

from rich.console import Console
from rich.panel import Panel
from rich.syntax import Syntax
from rich.table import Table

from chemsmart.agent.decision_events import (
    AssistantTurnEvent,
    CriticVerdictEvent,
    DryRunInputEvent,
    ErrorEvent,
    MethodEvent,
    PlanEvent,
    RequestEvent,
    RuntimeValidationEvent,
    SubmissionPreviewEvent,
    ToolUseEvent,
    parse_decision_event,
)

_INLINE_CLI_STDERR_HINTS = (
    "agent tui exited; debug log:",
    "Debug log saved to:",
)


def format_wizard_validation_errors(errors: list[str]) -> str:
    """Render wizard validation failures for a Click exception."""
    if not errors:
        return "wizard output did not validate"
    lines = ["wizard output did not validate:"]
    lines.extend(f"- {error}" for error in errors)
    return "\n".join(lines)


def wizard_verify_table(result) -> Table:
    """Render the stable wizard verification fields."""
    table = Table(title=f"Wizard verify: {result.server_name}")
    table.add_column("Field", style="cyan", no_wrap=True)
    table.add_column("Value")
    table.add_row("server_name", result.server_name)
    table.add_row("host", result.host or "-")
    table.add_row("mode", result.mode)
    table.add_row("would_submit_via", result.would_submit_via)
    table.add_row(
        "transport_invocation",
        json.dumps(result.transport_invocation or []),
    )
    table.add_row(
        "warnings",
        "\n".join(result.warnings) if result.warnings else "-",
    )
    table.add_row(
        "errors",
        "\n".join(result.errors) if result.errors else "-",
    )
    return table


def wizard_refresh_table(result: dict[str, object]) -> Table:
    """Render the stable wizard cache refresh fields."""
    table = Table(
        title=f"Wizard refresh: {result.get('server_name') or 'server'}"
    )
    table.add_column("Field", style="cyan", no_wrap=True)
    table.add_column("Value")
    table.add_row("cache_path", str(result.get("cache_path") or "-"))
    table.add_row("status", str(result.get("status") or "-"))
    table.add_row("host", str(result.get("host") or "-"))
    table.add_row("mode", str(result.get("mode") or "-"))
    table.add_row("scheduler", str(result.get("scheduler") or "-"))
    table.add_row("probed_at", str(result.get("probed_at") or "-"))
    node_summary = result.get("node_summary")
    if not isinstance(node_summary, dict):
        node_summary = {}
    table.add_row(
        "selected_queue",
        str(node_summary.get("selected_queue") or "-"),
    )
    table.add_row(
        "resources",
        f"cpu={node_summary.get('cpu')} mem_gb={node_summary.get('mem_gb')} gpu={node_summary.get('gpu')}",
    )
    table.add_row(
        "queue_counts",
        ("total={total} enabled={enabled} started={started} gpu={gpu}").format(
            total=node_summary.get("queue_count"),
            enabled=node_summary.get("enabled_queue_count"),
            started=node_summary.get("started_queue_count"),
            gpu=node_summary.get("gpu_queue_count"),
        ),
    )
    table.add_row("project", str(node_summary.get("project") or "-"))
    table.add_row(
        "scratch",
        f"{node_summary.get('scratch_dir') or '-'} (writable={node_summary.get('scratch_writable')})",
    )
    program_candidates = result.get("program_candidates")
    if not isinstance(program_candidates, dict):
        program_candidates = {}
    programs = []
    for name in ["gaussian", "orca", "nciplot"]:
        candidate = program_candidates.get(name)
        if not isinstance(candidate, dict):
            continue
        source = candidate.get("source")
        location = candidate.get("exefolder") or ", ".join(
            str(item) for item in candidate.get("module_candidates") or []
        )
        programs.append(f"{name}: {source} {location}".strip())
    table.add_row("programs", "\n".join(programs) if programs else "-")
    table.add_row("last_error", str(result.get("last_error") or "-"))
    return table


def tool_description(tool) -> str:
    """Return the first documentation line for a registry tool."""
    doc = inspect.getdoc(tool.func) or tool.name
    return doc.splitlines()[0].strip()


def principal_args(tool) -> str:
    """Return up to three principal schema arguments for table display."""
    schema = tool.input_schema.model_json_schema()
    properties = schema.get("properties", {})
    required = schema.get("required", [])
    ordered = [name for name in required if name in properties]
    ordered.extend(name for name in properties if name not in ordered)
    top_args = ordered[:3]
    return ", ".join(top_args) if top_args else "-"


def get_gateway_url(providers_module, provider_name: str) -> str:
    """Resolve the compatibility gateway shown by the doctor command."""
    if provider_name == "anthropic":
        return providers_module._GATEWAY_URL_ANTHROPIC
    return providers_module._GATEWAY_URL_OPENAI


def first_dry_run_result(result: dict) -> dict | None:
    """Project the first dry-run result from legacy or current payloads."""
    dry_run_results = result.get("dry_run_results") or []
    return next(iter(dry_run_results), result.get("dry_run_result"))


def is_advisory_plan(plan) -> bool:
    """Return whether a plan intentionally contains no executable steps."""
    return not bool(getattr(plan, "steps", None))


def sanitize_inline_cli_output(text: str) -> str:
    """Remove stderr-only logging hints from embedded CLI output."""
    lines = [
        line
        for line in text.splitlines()
        if not any(hint in line for hint in _INLINE_CLI_STDERR_HINTS)
    ]
    return "\n".join(lines).strip()


def stream_decision_event(console: Console, entry: dict) -> None:
    """Render one public decision-log event in the plain CLI stream."""
    event = parse_decision_event(entry)
    if isinstance(event, RequestEvent):
        console.print(Panel(event.request, title="Request"))
    elif isinstance(event, PlanEvent):
        _stream_plan(console, event)
    elif isinstance(event, AssistantTurnEvent):
        if event.text.strip():
            console.print(Panel(event.text.strip(), title="Assistant"))
    elif isinstance(event, ToolUseEvent):
        body = f"{event.tool} [{event.status}]"
        if event.reason:
            body += f"\n{event.reason}"
        console.print(Panel(body, title="Tool use"))
    elif isinstance(event, MethodEvent):
        _stream_method(console, event)
    elif isinstance(event, DryRunInputEvent):
        console.print(
            Panel(
                Syntax(event.content, "text", theme="ansi_dark"),
                title=event.inputfile or "Dry run input",
            )
        )
    elif isinstance(event, RuntimeValidationEvent):
        text = json.dumps(event.validation, indent=2, sort_keys=True)
        console.print(Panel(text, title="Runtime validation"))
    elif isinstance(event, SubmissionPreviewEvent):
        preview = json.dumps(event.preview, indent=2, sort_keys=True)
        console.print(Panel(preview, title="Submission preview"))
    elif isinstance(event, CriticVerdictEvent):
        _stream_critic(console, event)
    elif isinstance(event, ErrorEvent):
        console.print(Panel(event.message, title=event.title))


def _stream_plan(console: Console, event: PlanEvent) -> None:
    if is_advisory_plan(event.plan):
        console.print(
            Panel(
                event.plan.rationale or "No tool execution required.",
                title="Advice",
            )
        )
        return
    console.print(Panel(event.text, title="Plan"))


def _stream_method(console: Console, event: MethodEvent) -> None:
    recommendation = event.recommendation
    text = (
        f"{recommendation.get('functional') or 'manual'} / "
        f"{recommendation.get('basis') or 'manual'}\n"
        f"{recommendation.get('rationale') or ''}"
    ).strip()
    console.print(Panel(text, title="Method"))


def _stream_critic(console: Console, event: CriticVerdictEvent) -> None:
    verdict = event.verdict
    issues = "\n".join(f"- {issue}" for issue in verdict.issues)
    text = f"verdict: {verdict.verdict}\nconfidence: {verdict.confidence:.2f}"
    if issues:
        text += f"\n{issues}"
    console.print(Panel(text, title="Critic"))
