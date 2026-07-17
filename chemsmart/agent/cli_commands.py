"""Click group and subcommands for `chemsmart agent`."""

from __future__ import annotations

import json
import logging
import os
import warnings
from collections.abc import Callable
from contextlib import contextmanager
from pathlib import Path

import click
from rich.console import Console
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
from chemsmart.agent.services.cli_presenters import (
    first_dry_run_result,
    format_wizard_validation_errors,
    get_gateway_url,
    is_advisory_plan,
    principal_args,
)
from chemsmart.agent.services.cli_presenters import (
    sanitize_inline_cli_output as sanitize_inline_cli_output,
)
from chemsmart.agent.services.cli_presenters import (
    tool_description,
    wizard_refresh_table,
    wizard_verify_table,
)
from chemsmart.agent.services.command_logging import (
    _apply_third_party_silence,
    _enable_console_logging,
    _flush_logging_handlers,
    _silence_console_logging,
)
from chemsmart.agent.services.session_listing import (
    format_session_age,
    load_session_snapshots,
    truncate_session_request,
)
from chemsmart.agent.services.session_store import (
    SessionMigrationError,
    migrate_legacy_session,
    resolve_session_source,
)
from chemsmart.agent.wizard import (
    ProbeRunner,
    run_wizard,
    verify_server_yaml,
    write_server_yaml,
)
from chemsmart.agent.wizard.tools import wizard_refresh as run_wizard_refresh

logger = logging.getLogger(__name__)

_TuiLauncher = Callable[..., None]
_tui_launcher: _TuiLauncher | None = None


def configure_tui_launcher(launcher: _TuiLauncher | None) -> None:
    """Inject the optional TUI launcher without importing presentation code."""
    global _tui_launcher
    _tui_launcher = launcher


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
@click.option(
    "--debug",
    "debug",
    is_flag=True,
    default=False,
    help=(
        "Show DEBUG console logs including 3rd-party library chatter "
        "(numexpr/openai/anthropic/httpx). Equivalent to --verbose plus "
        "passing verbose flags through to executed chemsmart subprocesses."
    ),
)
@click.pass_context
def agent(ctx, plain: bool, verbose: bool, debug: bool):
    """AI-scientist agent commands for chemsmart."""
    ctx.ensure_object(dict)
    effective_verbose = verbose or debug
    ctx.obj["agent_verbose"] = effective_verbose
    ctx.obj["agent_debug"] = debug
    if ctx.invoked_subcommand is not None:
        return

    if effective_verbose:
        _enable_console_logging()
        logger.debug("Agent verbose logging enabled")
    if _tui_launcher is None:
        click.echo(
            "The agent TUI requires optional dependencies. "
            'Install them with: pip install -e ".[agent-tui]"'
        )
        raise click.exceptions.Exit(1)
    try:
        _tui_launcher(plain=plain)
    except ImportError as exc:
        click.echo(
            "The agent TUI requires optional dependencies. "
            'Install them with: pip install -e ".[agent-tui]"'
        )
        raise click.exceptions.Exit(1) from exc


@agent.command(name="_dump-cli-schema", hidden=True)
@click.option(
    "--out",
    "out_path",
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help="Write the CLI schema JSON to this path instead of stdout.",
)
def dump_cli_schema(out_path: Path | None) -> None:
    """Emit a JSON schema for the ChemSmart click command tree."""
    from chemsmart.agent.cli_schema import (
        build_chemsmart_cli_schema,
        dump_schema_to_json,
        schema_with_metadata,
    )

    if out_path is not None:
        dump_schema_to_json(out_path)
        return

    document = schema_with_metadata(build_chemsmart_cli_schema())
    click.echo(json.dumps(document, indent=2, sort_keys=True))


@agent.command()
@click.option(
    "--no-ping",
    is_flag=True,
    default=False,
    help="Skip the provider connectivity ping.",
)
def doctor(no_ping: bool):
    """Validate agent.yaml, provider connectivity, and registered tools."""
    with _agent_command_logging():
        import chemsmart.agent.providers as providers
        from chemsmart.agent.provider_config import (
            AgentProviderConfigError,
            load_active_provider_config,
        )

        try:
            provider_config = load_active_provider_config()
        except AgentProviderConfigError as exc:
            raise click.ClickException(str(exc)) from exc

        try:
            provider = providers.get_provider()
        except ProviderError as exc:
            raise click.ClickException(str(exc)) from exc

        registry = ToolRegistry.default()
        if provider_config is not None:
            provider_name = provider_config.type
            gateway_url = provider_config.base_url or "(default gateway)"
            click.echo(
                "agent.yaml: "
                f"{Path.home() / '.chemsmart' / 'agent' / 'agent.yaml'}"
            )
            click.echo(f"active: {provider_config.name}")
            click.echo(f"type: {provider_config.type}")
            click.echo(f"model: {provider_config.model}")
            click.echo(f"base_url: {gateway_url}")
        else:
            api_key = os.environ.get("ai_api_key", "").strip()
            provider_name = os.environ.get("AI_PROVIDER", "").strip()
            gateway_url = get_gateway_url(providers, provider_name)
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
        if is_advisory_plan(result["plan"]):
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
            dry_run_result = first_dry_run_result(result)
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
@click.option(
    "--debug",
    "debug",
    is_flag=True,
    default=False,
    help=(
        "Show DEBUG console logs and forward verbose to executed "
        "chemsmart subprocesses."
    ),
)
@click.argument("request")
@click.pass_context
def ask(ctx, request: str, mode_name: str, yolo: bool, debug: bool):
    """Synthesize one ChemSmart CLI command from a natural-language request."""
    if mode_name != "driving" or yolo:
        warnings.warn(
            "ask now uses synthesis mode; legacy flags ignored",
            DeprecationWarning,
            stacklevel=2,
        )
    if debug:
        if isinstance(ctx.obj, dict):
            ctx.obj["agent_verbose"] = True
            ctx.obj["agent_debug"] = True
    with _agent_command_logging():
        from chemsmart.agent.synthesis import SynthesisSession

        SynthesisSession(debug=debug).run_interactive(request)


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
        if is_advisory_plan(result["plan"]):
            click.echo("Advice:")
            click.echo(
                result["plan"].rationale or "No tool execution required."
            )
        else:
            click.echo(render_plan(result["plan"]))
            click.echo(f"critic verdict: {result['critic_verdict'].verdict}")
            dry_run_result = first_dry_run_result(result)
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
        snapshots = load_session_snapshots(session_root)
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
                format_session_age(snapshot["timestamp"]),
                truncate_session_request(snapshot["request"]),
                snapshot["status"],
            )
        Console().print(table)


@agent.command(name="migrate-session")
@click.argument("source")
@click.option(
    "--destination",
    type=click.Path(path_type=Path),
    default=None,
    help="Write the migrated copy to this new directory.",
)
def migrate_session(source: str, destination: Path | None) -> None:
    """Copy a legacy state.json session into the current session schema."""
    session_root = Path(_default_session_root())
    source_path = resolve_session_source(source, session_root)
    try:
        result = migrate_legacy_session(source_path, destination)
    except SessionMigrationError as exc:
        raise click.ClickException(str(exc)) from exc
    click.echo(f"source preserved: {result['source']}")
    click.echo(f"migrated session: {result['session_id']}")
    click.echo(f"destination: {result['destination']}")


@agent.command()
@click.argument("name")
@click.option(
    "--host",
    "host",
    default=None,
    help="Optional SSH host or alias to probe remotely.",
)
@click.option(
    "--write",
    is_flag=True,
    default=False,
    help="Write the rendered YAML to ~/.chemsmart/server/<name>.yaml.",
)
@click.option(
    "--overwrite",
    is_flag=True,
    default=False,
    help="Replace an existing ~/.chemsmart/server/<name>.yaml file.",
)
@click.option(
    "--yes",
    is_flag=True,
    default=False,
    help="Skip confirmation before writing the server YAML.",
)
def wizard(
    name: str,
    host: str | None,
    write: bool,
    overwrite: bool,
    yes: bool,
):
    """Probe a target server and render or write a wizard YAML config."""
    with _agent_command_logging():
        if overwrite and not write:
            raise click.ClickException("--overwrite requires --write")

        outcome = run_wizard(
            ProbeRunner,
            server_name=name,
            ssh_host_hint=host,
            write=False,
            overwrite=overwrite,
        )

        if not write:
            click.echo(outcome.plan.text, nl=False)
            if not outcome.validation.ok:
                raise click.ClickException(
                    format_wizard_validation_errors(outcome.validation.errors)
                )
            return

        if not outcome.validation.ok:
            raise click.ClickException(
                format_wizard_validation_errors(outcome.validation.errors)
            )

        target = Path.home() / ".chemsmart" / "server" / f"{name}.yaml"
        if not yes:
            click.confirm(
                f"Write wizard YAML to {target}?",
                default=False,
                abort=True,
            )

        written_path = write_server_yaml(
            name=name,
            yaml_text=outcome.plan.text,
            overwrite=overwrite,
        )
        click.echo(f"Wrote {written_path}")


@agent.command(name="wizard-verify")
@click.argument("name")
def wizard_verify(name: str):
    """Verify wizard/server transport wiring for an existing YAML."""
    with _agent_command_logging():
        result = verify_server_yaml(name)
        Console().print(wizard_verify_table(result))
        if result.errors:
            raise click.ClickException("\n".join(result.errors))


@agent.command(name="wizard-refresh")
@click.argument("name")
@click.option(
    "--force",
    is_flag=True,
    default=False,
    help="Reprobe even when the cache is still fresh.",
)
def wizard_refresh(name: str, force: bool):
    """Refresh or reuse the wizard node sidecar cache for a server."""
    with _agent_command_logging():
        result = run_wizard_refresh(name, force=force)
        Console().print(wizard_refresh_table(result))


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
                tool_description(tool),
                principal_args(tool),
            )
        Console().print(table)


@contextmanager
def _agent_command_logging():
    verbose = _agent_verbose_enabled()
    debug = _agent_debug_enabled()
    if verbose:
        _enable_console_logging()
        _apply_third_party_silence(logging.DEBUG if debug else logging.WARNING)
        logger.debug("Agent verbose logging enabled")
        try:
            yield
        finally:
            _flush_logging_handlers()
        return

    log_path = _silence_console_logging(_agent_log_root())
    _apply_third_party_silence(logging.WARNING)
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


def _agent_debug_enabled() -> bool:
    ctx = click.get_current_context(silent=True)
    while ctx is not None:
        if isinstance(ctx.obj, dict) and "agent_debug" in ctx.obj:
            return bool(ctx.obj["agent_debug"])
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
