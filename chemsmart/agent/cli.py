"""
Click group and subcommands for `chemsmart agent`.

Wave 1 ships only the `doctor` subcommand.
run / resume / tools are stubs for subsequent waves.
"""

import os
import sys

import click


@click.group(name="agent")
def agent():
    """AI-scientist agent commands for chemsmart."""
    pass


@agent.command()
def doctor():
    """Validate api.env, AI_PROVIDER, and provider connectivity."""
    import chemsmart.agent.providers as _p

    try:
        _p.get_provider()
    except _p.ProviderError as exc:
        raise click.ClickException(str(exc))

    api_key = os.environ.get("ai_api_key", "").strip()
    provider_name = os.environ.get("AI_PROVIDER", "").strip()

    click.echo(f"AI_PROVIDER={provider_name} OK")
    click.echo(f"api.env: OK (key length={len(api_key)})")
    click.echo("tools registered: 0")


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
