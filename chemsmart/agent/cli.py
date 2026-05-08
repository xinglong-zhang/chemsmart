"""Click group and subcommands for `chemsmart agent`."""

from __future__ import annotations

import os

import click

from chemsmart.agent.core import AgentSession, render_plan
from chemsmart.agent.providers import ProviderError
from chemsmart.agent.registry import ToolRegistry


@click.group(name="agent")
def agent():
    """AI-scientist agent commands for chemsmart."""
    pass


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
    if dry_run_results:
        return dry_run_results[0]
    return result.get("dry_run_result")
