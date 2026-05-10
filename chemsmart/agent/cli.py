"""Click group and subcommands for `chemsmart agent`."""

from __future__ import annotations

import inspect
import json
import logging
import os
import threading
import time
from contextlib import contextmanager
from datetime import UTC, datetime
from pathlib import Path

import click
from rich.console import Console
from rich.panel import Panel
from rich.syntax import Syntax
from rich.table import Table

from chemsmart.agent.core import (
    AgentSession,
    _default_session_root,
    render_plan,
)
from chemsmart.agent.permissions import (
    ApprovalDecision,
    PermissionMode,
    PermissionPolicy,
)
from chemsmart.agent.providers import ProviderError
from chemsmart.agent.registry import ToolRegistry
from chemsmart.agent.tui._logging import (
    _enable_console_logging,
    _flush_logging_handlers,
    _silence_console_logging,
)
from chemsmart.agent.tui.events import (
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

logger = logging.getLogger(__name__)

_INLINE_CLI_STDERR_HINTS = (
    "agent tui exited; debug log:",
    "Debug log saved to:",
)


@click.group(name="agent", invoke_without_command=True)
@click.option(
    "--plain",
    is_flag=True,
    default=False,
    help="Run the agent TUI in plain inline mode for conservative terminals.",
)
@click.option(
    "--verbose",
    "verbose",
    "-v",
    is_flag=True,
    default=False,
    help=(
        "Keep console logging at DEBUG instead of writing it only to the "
        "agent log file."
    ),
)
@click.pass_context
def agent(ctx, plain: bool, verbose: bool):
    """AI-scientist agent commands for chemsmart."""
    ctx.ensure_object(dict)
    ctx.obj["agent_verbose"] = verbose
    if ctx.invoked_subcommand is not None:
        return

    if verbose:
        _enable_console_logging()
        logger.debug("Agent verbose logging enabled")
    try:
        from chemsmart.agent.tui import launch_tui
    except ImportError as exc:
        click.echo(
            "The agent TUI requires optional dependencies. "
            'Install them with: pip install -e ".[agent-tui]"'
        )
        raise click.exceptions.Exit(1) from exc
    launch_tui(plain=plain)


@agent.command()
@click.option(
    "--no-ping",
    is_flag=True,
    default=False,
    help="Skip the provider connectivity ping.",
)
def doctor(no_ping: bool):
    """Validate api.env, AI_PROVIDER, and provider connectivity."""
    with _agent_command_logging():
        import chemsmart.agent.providers as providers

        try:
            provider = providers.get_provider()
        except ProviderError as exc:
            raise click.ClickException(str(exc)) from exc

        registry = ToolRegistry.default()
        api_key = os.environ.get("ai_api_key", "").strip()
        provider_name = os.environ.get("AI_PROVIDER", "").strip()
        gateway_url = _get_gateway_url(providers, provider_name)

        click.echo(f"AI_PROVIDER={provider_name} OK")
        click.echo(f"api.env: OK (key length={len(api_key)})")
        click.echo(f"gateway: {gateway_url}")

        if no_ping:
            click.echo("ping: skipped (--no-ping)")
        else:
            try:
                ping = provider.ping()
            except ProviderError as exc:
                click.echo(f"ping: FAILED {exc}")
                click.echo(f"tools registered: {len(registry.list_tools())}")
                if provider_name == "anthropic":
                    click.echo(
                        "WARN: this gateway tenant may not have "
                        "Anthropic access"
                    )
                raise click.exceptions.Exit(1) from exc
            click.echo(
                "ping: ok "
                f"(model={ping['resolved_model']}, "
                f"latency={ping['latency_ms']}ms)"
            )

        click.echo(f"tools registered: {len(registry.list_tools())}")
        if provider_name == "anthropic":
            click.echo(
                "WARN: this gateway tenant may not have Anthropic access"
            )


@agent.command(name="run")
@click.option(
    "--dry-submit/--execute",
    default=True,
    help="Write scripts without real remote submission, or execute submit_hpc.",
)
@click.option(
    "--mode",
    "mode_name",
    type=click.Choice(["permission", "driving"]),
    default="driving",
    show_default=True,
    help="Permission mode prompts for every tool; driving mode is autonomous.",
)
@click.option(
    "--yolo",
    is_flag=True,
    default=False,
    help="Allow run_local and submit_hpc in driving mode.",
)
@click.option(
    "--allow-remote-unknown",
    is_flag=True,
    default=False,
    help="Allow validate_runtime partial results to proceed.",
)
@click.option(
    "--allow-critic-override",
    is_flag=True,
    default=False,
    help="Allow warn verdicts that are not remote-unknown gating failures.",
)
@click.option(
    "--resume",
    "resume_id",
    default=None,
    help="Resume an existing agent session id instead of planning a new one.",
)
@click.argument("request", required=False)
def agent_run(
    request: str | None,
    dry_submit: bool,
    mode_name: str,
    yolo: bool,
    allow_remote_unknown: bool,
    allow_critic_override: bool,
    resume_id: str | None,
):
    """Plan and execute an agent workflow."""
    with _agent_command_logging():
        session = AgentSession()
        policy = _permission_policy_from_flags(mode_name, yolo)
        try:
            if resume_id:
                result = _resume_session_or_fail(
                    resume_id,
                    dry_submit=dry_submit,
                    allow_remote_unknown=allow_remote_unknown,
                    allow_critic_override=allow_critic_override,
                    policy=policy,
                    approver=_auto_deny_approver,
                )
            else:
                if not request:
                    raise click.ClickException(
                        "REQUEST is required unless --resume is used"
                    )
                result = session.run_loop(
                    request,
                    policy=policy,
                    approver=_auto_deny_approver,
                )
        except RuntimeError as exc:
            raise click.ClickException(str(exc)) from exc

        click.echo(f"session: {result['session_id']}")
        if _is_advisory_plan(result["plan"]):
            click.echo("Advice:")
            click.echo(
                result.get("assistant_output")
                or result["plan"].rationale
                or "No tool execution required."
            )
        else:
            click.echo(render_plan(result["plan"]))
            click.echo(f"approval mode: {result.get('approval_mode')}")
            click.echo(f"yolo: {bool(result.get('yolo'))}")
            dry_run_result = _first_dry_run_result(result)
            if dry_run_result:
                click.echo(f"inputfile: {dry_run_result['inputfile']}")
        if result.get("assistant_output"):
            click.echo(f"assistant: {result['assistant_output']}")
        click.echo(f"decision log: {result['session_dir']}/decision_log.jsonl")
        if result.get("limit_reason"):
            click.echo(f"limit reason: {result['limit_reason']}")


@agent.command()
@click.option(
    "--mode",
    "mode_name",
    type=click.Choice(["permission", "driving"]),
    default="driving",
    show_default=True,
    help="Permission mode prompts for every tool; driving mode is autonomous.",
)
@click.option(
    "--yolo",
    is_flag=True,
    default=False,
    help="Allow run_local and submit_hpc in driving mode.",
)
@click.argument("request")
def ask(request: str, mode_name: str, yolo: bool):
    """Run a one-shot dry-run request and stream Rich output to stdout."""
    with _agent_command_logging():
        session = AgentSession()
        console = Console()
        result_box: dict[str, dict] = {}
        error_box: dict[str, Exception] = {}
        policy = _permission_policy_from_flags(mode_name, yolo)

        def runner() -> None:
            try:
                result_box["result"] = session.run_loop(
                    request,
                    policy=policy,
                    approver=_auto_deny_approver,
                )
            except Exception as exc:  # pragma: no cover - bridged below
                error_box["error"] = exc

        thread = threading.Thread(target=runner, daemon=True)
        thread.start()

        seen = 0
        while thread.is_alive() or (
            session.session_dir
            and (session.session_dir / "decision_log.jsonl").exists()
        ):
            log_path = (
                None
                if session.session_dir is None
                else session.session_dir / "decision_log.jsonl"
            )
            if log_path and log_path.exists():
                lines = log_path.read_text(encoding="utf-8").splitlines()
                for line in lines[seen:]:
                    if line.strip():
                        _stream_event(console, json.loads(line))
                seen = len(lines)
            if not thread.is_alive():
                break
            time.sleep(0.05)
        thread.join()

        if "error" in error_box:
            raise click.ClickException(str(error_box["error"])) from error_box[
                "error"
            ]
        result = result_box["result"]
        if result.get("limit_reason"):
            click.echo(f"limit reason: {result['limit_reason']}")


@agent.command()
@click.option(
    "--dry-submit/--execute",
    default=True,
    help="Write scripts without real remote submission, or execute submit_hpc.",
)
@click.option(
    "--allow-remote-unknown",
    is_flag=True,
    default=False,
    help="Allow validate_runtime partial results to proceed.",
)
@click.option(
    "--allow-critic-override",
    is_flag=True,
    default=False,
    help="Allow warn verdicts that are not remote-unknown gating failures.",
)
@click.argument("session_id")
def resume(
    session_id: str,
    dry_submit: bool,
    allow_remote_unknown: bool,
    allow_critic_override: bool,
):
    """Resume an agent session."""
    with _agent_command_logging():
        try:
            result = _resume_session_or_fail(
                session_id,
                dry_submit=dry_submit,
                allow_remote_unknown=allow_remote_unknown,
                allow_critic_override=allow_critic_override,
            )
        except RuntimeError as exc:
            raise click.ClickException(str(exc)) from exc
        click.echo(f"session: {result['session_id']}")
        if _is_advisory_plan(result["plan"]):
            click.echo("Advice:")
            click.echo(
                result["plan"].rationale or "No tool execution required."
            )
        else:
            click.echo(render_plan(result["plan"]))
            click.echo(f"critic verdict: {result['critic_verdict'].verdict}")
            dry_run_result = _first_dry_run_result(result)
            if dry_run_result:
                click.echo(f"inputfile: {dry_run_result['inputfile']}")
        if result.get("blocked"):
            raise click.ClickException("critic gating blocked execution")


@agent.command()
@click.option(
    "--all",
    "show_all",
    is_flag=True,
    default=False,
    help="Show every saved session instead of the 10 most recent.",
)
def sessions(show_all: bool):
    """List recent agent sessions."""
    with _agent_command_logging():
        session_root = Path(_default_session_root())
        snapshots = _load_session_snapshots(session_root)
        if not show_all:
            snapshots = snapshots[:10]
        if not snapshots:
            click.echo(f"No agent sessions found in {session_root}")
            return

        table = Table(title="Agent sessions")
        table.add_column("Session ID", style="cyan", no_wrap=True)
        table.add_column("Age", no_wrap=True)
        table.add_column("Request")
        table.add_column("Status", no_wrap=True)
        for snapshot in snapshots:
            table.add_row(
                snapshot["session_id"],
                _format_age(snapshot["timestamp"]),
                _truncate_request(snapshot["request"]),
                snapshot["status"],
            )
        Console().print(table)


@agent.command()
def tools():
    """List registered agent tools."""
    with _agent_command_logging():
        registry = ToolRegistry.default()
        table = Table(title=f"Registered tools ({len(registry.list_tools())})")
        table.add_column("Tool")
        table.add_column("Description")
        table.add_column("Principal args")
        for tool in registry.list_tools():
            table.add_row(
                tool.name,
                _tool_description(tool),
                _principal_args(tool),
            )
        Console().print(table)


@contextmanager
def _agent_command_logging():
    verbose = _agent_verbose_enabled()
    if verbose:
        _enable_console_logging()
        logger.debug("Agent verbose logging enabled")
        try:
            yield
        finally:
            _flush_logging_handlers()
        return

    log_path = _silence_console_logging(_agent_log_root())
    try:
        yield
    finally:
        _flush_logging_handlers()
        if log_path.exists() and log_path.stat().st_size:
            click.echo(f"Debug log saved to: {log_path}", err=True)


def _agent_verbose_enabled() -> bool:
    ctx = click.get_current_context(silent=True)
    while ctx is not None:
        if isinstance(ctx.obj, dict) and "agent_verbose" in ctx.obj:
            return bool(ctx.obj["agent_verbose"])
        ctx = ctx.parent
    return False


def _agent_log_root() -> Path:
    return Path(_default_session_root()).parent


def _resume_session_or_fail(session_id: str, **kwargs):
    try:
        policy = kwargs.pop("policy", None)
        approver = kwargs.pop("approver", None)
        if policy is not None:
            session = AgentSession.load(session_id, **kwargs)
            request = session.state.request or "Continue."
            return session.run_loop(
                request,
                policy=policy,
                approver=approver,
            )
        return AgentSession.resume(session_id, **kwargs)
    except FileNotFoundError as exc:
        raise click.ClickException(
            f"session '{session_id}' not found. List recent sessions with: "
            "chemsmart agent sessions"
        ) from exc


def _load_session_snapshots(session_root: Path) -> list[dict[str, object]]:
    if not session_root.exists():
        return []
    snapshots = []
    for session_dir in session_root.iterdir():
        if not session_dir.is_dir():
            continue
        snapshots.append(_session_snapshot(session_dir))
    return sorted(
        snapshots,
        key=lambda snapshot: snapshot["timestamp"],
        reverse=True,
    )


def _session_snapshot(session_dir: Path) -> dict[str, object]:
    metadata = _load_json(session_dir / "session_metadata.json")
    state = _load_json(session_dir / "session.json") or _load_json(
        session_dir / "state.json"
    )
    entries = _load_decision_entries(session_dir / "decision_log.jsonl")

    request = ""
    if isinstance(metadata, dict):
        request = str(metadata.get("request") or "")
    if not request and isinstance(state, dict):
        request = str(state.get("request") or "")
    if not request:
        for entry in entries:
            if entry.get("kind") == "request":
                request = str(entry.get("payload", {}).get("request") or "")
                break

    timestamp = _coerce_timestamp(
        (metadata or {}).get("ended_at")
        or (metadata or {}).get("started_at")
        or (entries[-1].get("ts") if entries else None)
        or (state or {}).get("started_at")
    )
    if timestamp is None:
        timestamp = datetime.fromtimestamp(session_dir.stat().st_mtime, tz=UTC)

    status = _session_status(metadata, state, entries)
    return {
        "session_id": session_dir.name,
        "request": request,
        "status": status,
        "timestamp": timestamp,
    }


def _session_status(metadata, state, entries: list[dict[str, object]]) -> str:
    if isinstance(metadata, dict):
        if metadata.get("blocked"):
            return "blocked"
        if metadata.get("critic_verdict"):
            return "ok"

    summary = next(
        (
            entry.get("payload")
            for entry in reversed(entries)
            if entry.get("kind") == "session_summary"
        ),
        None,
    )
    if isinstance(summary, dict):
        if summary.get("blocked"):
            return "blocked"
        return "ok"

    if any(entry.get("kind") in {"error", "tool_error"} for entry in entries):
        return "error"

    if isinstance(state, dict):
        planned = int(state.get("total_steps_planned") or 0)
        current = int(state.get("current_step_index") or 0)
        if planned == 0 or current < planned:
            return "in-progress"
        return "ok"

    return "in-progress"


def _load_json(path: Path) -> dict | None:
    if not path.exists():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def _load_decision_entries(path: Path) -> list[dict[str, object]]:
    if not path.exists():
        return []
    return [
        json.loads(line)
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]


def _coerce_timestamp(value: object) -> datetime | None:
    if not value:
        return None
    text = str(value)
    if text.endswith("Z"):
        text = text[:-1] + "+00:00"
    try:
        return datetime.fromisoformat(text).astimezone(UTC)
    except ValueError:
        return None


def _format_age(timestamp: datetime) -> str:
    delta_seconds = max(
        0, int((datetime.now(UTC) - timestamp).total_seconds())
    )
    if delta_seconds < 60:
        return f"{delta_seconds}s ago"
    minutes = delta_seconds // 60
    if minutes < 60:
        return f"{minutes}m ago"
    hours = minutes // 60
    if hours < 24:
        return f"{hours}h ago"
    days = hours // 24
    return f"{days}d ago"


def _truncate_request(request: object, limit: int = 60) -> str:
    compact = " ".join(str(request or "").split())
    if len(compact) <= limit:
        return compact
    return f"{compact[: limit - 3]}..."


def _tool_description(tool) -> str:
    doc = inspect.getdoc(tool.func) or tool.name
    return doc.splitlines()[0].strip()


def _principal_args(tool) -> str:
    schema = tool.input_schema.model_json_schema()
    properties = schema.get("properties", {})
    required = schema.get("required", [])
    ordered = [name for name in required if name in properties]
    ordered.extend(name for name in properties if name not in ordered)
    top_args = ordered[:3]
    return ", ".join(top_args) if top_args else "-"


def _get_gateway_url(providers_module, provider_name: str) -> str:
    if provider_name == "anthropic":
        return providers_module._GATEWAY_URL_ANTHROPIC
    return providers_module._GATEWAY_URL_OPENAI


def _first_dry_run_result(result: dict) -> dict | None:
    dry_run_results = result.get("dry_run_results") or []
    return next(iter(dry_run_results), result.get("dry_run_result"))


def _is_advisory_plan(plan) -> bool:
    return not bool(getattr(plan, "steps", None))


def _permission_policy_from_flags(
    mode_name: str,
    yolo: bool,
) -> PermissionPolicy:
    return PermissionPolicy(
        mode=PermissionMode(mode_name),
        yolo=yolo,
    )


def _auto_deny_approver(_request) -> ApprovalDecision:
    return ApprovalDecision.DENY


def sanitize_inline_cli_output(text: str) -> str:
    lines = [
        line
        for line in text.splitlines()
        if not any(hint in line for hint in _INLINE_CLI_STDERR_HINTS)
    ]
    return "\n".join(lines).strip()


def _stream_event(console: Console, entry: dict) -> None:
    event = parse_decision_event(entry)
    if isinstance(event, RequestEvent):
        console.print(Panel(event.request, title="Request"))
    elif isinstance(event, PlanEvent):
        if _is_advisory_plan(event.plan):
            console.print(
                Panel(
                    event.plan.rationale or "No tool execution required.",
                    title="Advice",
                )
            )
        else:
            console.print(Panel(event.text, title="Plan"))
    elif isinstance(event, AssistantTurnEvent):
        if event.text.strip():
            console.print(Panel(event.text.strip(), title="Assistant"))
    elif isinstance(event, ToolUseEvent):
        body = f"{event.tool} [{event.status}]"
        if event.reason:
            body += f"\n{event.reason}"
        console.print(Panel(body, title="Tool use"))
    elif isinstance(event, MethodEvent):
        recommendation = event.recommendation
        text = (
            f"{recommendation.get('functional') or 'manual'} / "
            f"{recommendation.get('basis') or 'manual'}\n"
            f"{recommendation.get('rationale') or ''}"
        ).strip()
        console.print(Panel(text, title="Method"))
    elif isinstance(event, DryRunInputEvent):
        console.print(
            Panel(
                Syntax(event.content, "text", theme="ansi_dark"),
                title=event.inputfile or "Dry run input",
            )
        )
    elif isinstance(event, RuntimeValidationEvent):
        validation = event.validation
        text = json.dumps(validation, indent=2, sort_keys=True)
        console.print(Panel(text, title="Runtime validation"))
    elif isinstance(event, SubmissionPreviewEvent):
        preview = json.dumps(event.preview, indent=2, sort_keys=True)
        console.print(Panel(preview, title="Submission preview"))
    elif isinstance(event, CriticVerdictEvent):
        verdict = event.verdict
        issues = "\n".join(f"- {issue}" for issue in verdict.issues)
        text = (
            f"verdict: {verdict.verdict}\n"
            f"confidence: {verdict.confidence:.2f}"
        )
        if issues:
            text += f"\n{issues}"
        console.print(Panel(text, title="Critic"))
    elif isinstance(event, ErrorEvent):
        console.print(Panel(event.message, title=event.title))
