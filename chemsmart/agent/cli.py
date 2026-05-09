"""Click group and subcommands for `chemsmart agent`."""

from __future__ import annotations

import json
import os
import threading
import time

import click
from rich.console import Console
from rich.panel import Panel
from rich.syntax import Syntax

from chemsmart.agent.core import AgentSession, render_plan
from chemsmart.agent.providers import ProviderError
from chemsmart.agent.registry import ToolRegistry
from chemsmart.agent.tui._logging import _silence_console_logging
from chemsmart.agent.tui.events import (
    CriticVerdictEvent,
    DryRunInputEvent,
    ErrorEvent,
    MethodEvent,
    PlanEvent,
    RequestEvent,
    RuntimeValidationEvent,
    SubmissionPreviewEvent,
    parse_decision_event,
)


@click.group(name="agent", invoke_without_command=True)
@click.option(
    "--plain",
    is_flag=True,
    default=False,
    help="Run the agent TUI in plain inline mode for conservative terminals.",
)
@click.pass_context
def agent(ctx, plain: bool):
    """AI-scientist agent commands for chemsmart."""
    if ctx.invoked_subcommand is not None:
        return
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
                    "WARN: this gateway tenant may not have Anthropic access"
                )
            raise click.exceptions.Exit(1) from exc
        click.echo(
            "ping: ok "
            f"(model={ping['resolved_model']}, latency={ping['latency_ms']}ms)"
        )

    click.echo(f"tools registered: {len(registry.list_tools())}")
    if provider_name == "anthropic":
        click.echo("WARN: this gateway tenant may not have Anthropic access")


@agent.command(name="run")
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
    allow_remote_unknown: bool,
    allow_critic_override: bool,
    resume_id: str | None,
):
    """Plan and execute an agent workflow."""
    session = AgentSession()
    try:
        if resume_id:
            result = AgentSession.resume(
                resume_id,
                dry_submit=dry_submit,
                allow_remote_unknown=allow_remote_unknown,
                allow_critic_override=allow_critic_override,
            )
        else:
            if not request:
                raise click.ClickException(
                    "REQUEST is required unless --resume is used"
                )
            result = session.run(
                request,
                dry_submit=dry_submit,
                allow_remote_unknown=allow_remote_unknown,
                allow_critic_override=allow_critic_override,
            )
    except RuntimeError as exc:
        raise click.ClickException(str(exc)) from exc

    click.echo(f"session: {result['session_id']}")
    if _is_advisory_plan(result["plan"]):
        click.echo("Advice:")
        click.echo(result["plan"].rationale or "No tool execution required.")
    else:
        click.echo(render_plan(result["plan"]))
        verdict = result["critic_verdict"]
        click.echo(f"critic verdict: {verdict.verdict}")
        if verdict.issues:
            for issue in verdict.issues:
                click.echo(f"- {issue}")
        dry_run_result = _first_dry_run_result(result)
        if dry_run_result:
            click.echo(f"inputfile: {dry_run_result['inputfile']}")
    click.echo(f"decision log: {result['session_dir']}/decision_log.jsonl")
    if result.get("blocked"):
        raise click.ClickException("critic gating blocked execution")


@agent.command()
@click.argument("request")
def ask(request: str):
    """Run a one-shot dry-run request and stream Rich output to stdout."""
    session = AgentSession()
    _silence_console_logging(session.session_root.parent)
    console = Console()
    result_box: dict[str, dict] = {}
    error_box: dict[str, Exception] = {}

    def runner() -> None:
        try:
            result_box["result"] = session.run(
                request,
                dry_submit=True,
                pause_before_risky=True,
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
    if result.get("blocked"):
        raise click.ClickException("critic gating blocked execution")


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
    try:
        result = AgentSession.resume(
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
        click.echo(result["plan"].rationale or "No tool execution required.")
    else:
        click.echo(render_plan(result["plan"]))
        click.echo(f"critic verdict: {result['critic_verdict'].verdict}")
        dry_run_result = _first_dry_run_result(result)
        if dry_run_result:
            click.echo(f"inputfile: {dry_run_result['inputfile']}")
    if result.get("blocked"):
        raise click.ClickException("critic gating blocked execution")


@agent.command()
def tools():
    """List registered agent tools."""
    registry = ToolRegistry.default()
    click.echo(f"registered tools: {len(registry.list_tools())}")
    for tool in registry.list_tools():
        click.echo(f"- {tool.name}")


def _get_gateway_url(providers_module, provider_name: str) -> str:
    if provider_name == "anthropic":
        return providers_module._GATEWAY_URL_ANTHROPIC
    return providers_module._GATEWAY_URL_OPENAI


def _first_dry_run_result(result: dict) -> dict | None:
    dry_run_results = result.get("dry_run_results") or []
    return next(iter(dry_run_results), result.get("dry_run_result"))


def _is_advisory_plan(plan) -> bool:
    return not bool(getattr(plan, "steps", None))


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
            f"verdict: {verdict.verdict}\nconfidence: {verdict.confidence:.2f}"
        )
        if issues:
            text += f"\n{issues}"
        console.print(Panel(text, title="Critic"))
    elif isinstance(event, ErrorEvent):
        console.print(Panel(event.message, title=event.title))
