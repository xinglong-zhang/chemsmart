"""
Click group and subcommands for `chemsmart agent`.

Wave 1 ships only the `doctor` subcommand.
run / resume / tools are stubs for subsequent waves.
"""

from __future__ import annotations

import os
import sys

import click

_REGISTERED_TOOL_NAMES = ("build_molecule", "recommend_method")


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
    import chemsmart.agent.providers as _p

    try:
        provider = _p.get_provider()
    except _p.ProviderError as exc:
        raise click.ClickException(str(exc))

    api_key = os.environ.get("ai_api_key", "").strip()
    provider_name = os.environ.get("AI_PROVIDER", "").strip()
    gateway_url = _get_gateway_url(_p, provider_name)

    click.echo(f"AI_PROVIDER={provider_name} OK")
    click.echo(f"api.env: OK (key length={len(api_key)})")
    click.echo(f"gateway: {gateway_url}")

    if no_ping:
        click.echo("ping: skipped (--no-ping)")
    else:
        try:
            ping = provider.ping()
        except _p.ProviderError as exc:
            click.echo(f"ping: FAILED {exc}")
            click.echo(f"tools registered: {len(_REGISTERED_TOOL_NAMES)}")
            if provider_name == "anthropic":
                click.echo(
                    "WARN: this gateway tenant may not have "
                    "Anthropic access"
                )
            raise click.exceptions.Exit(1)

        click.echo(
            "ping: ok "
            f"(model={ping['resolved_model']}, "
            f"latency={ping['latency_ms']}ms)"
        )

    click.echo(f"tools registered: {len(_REGISTERED_TOOL_NAMES)}")
    if provider_name == "anthropic":
        click.echo("WARN: this gateway tenant may not have Anthropic access")


@agent.command(name="run")
@click.argument("request", required=False)
def agent_run(request):
    """Plan and execute an agent workflow (not yet implemented)."""
    click.echo("not yet implemented (see bin/plan.md)")
    sys.exit(1)


@agent.command()
@click.argument("session_id", required=False)
def resume(session_id):
    """Resume an agent session (not yet implemented)."""
    click.echo("not yet implemented (see bin/plan.md)")
    sys.exit(1)


@agent.command()
def tools():
    """List registered agent tools (not yet implemented)."""
    click.echo("not yet implemented (see bin/plan.md)")
    sys.exit(1)


def _get_gateway_url(providers_module, provider_name: str) -> str:
    if provider_name == "anthropic":
        return providers_module._GATEWAY_URL_ANTHROPIC
    return providers_module._GATEWAY_URL_OPENAI
