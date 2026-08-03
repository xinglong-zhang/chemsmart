"""Provider-neutral natural-language entrypoint for ChemSmart Agent."""

from __future__ import annotations

from pathlib import Path

import click

from chemsmart.utils.cli import MyGroup


def _task_options(function):
    options = (
        click.option("--task", type=str, default=None, help="Natural-language task."),
        click.option(
            "--task-file",
            type=click.Path(exists=True, dir_okay=False, path_type=Path),
            default=None,
            help="UTF-8 file containing the natural-language task.",
        ),
        click.option(
            "--provider",
            type=click.Choice(("deepseek",), case_sensitive=False),
            required=True,
        ),
        click.option(
            "--secret-file",
            type=click.Path(exists=True, dir_okay=False, path_type=Path),
            required=True,
            help="Secret assignment file parsed as data; it is never sourced.",
        ),
        click.option(
            "--workspace",
            type=click.Path(exists=True, file_okay=False, path_type=Path),
            required=True,
            help="Disposable task workspace containing user-approved artifacts.",
        ),
    )
    for option in reversed(options):
        function = option(function)
    return function


def _read_task(task: str | None, task_file: Path | None) -> str:
    if (task is None) == (task_file is None):
        raise click.UsageError("provide exactly one of --task or --task-file")
    value = task if task is not None else task_file.read_text(encoding="utf-8")
    value = str(value).strip()
    if not value:
        raise click.UsageError("the task must not be empty")
    return value


@click.group(name="agent", cls=MyGroup)
def agent():
    """Plan or run computational chemistry through typed ChemSmart tools."""


@agent.command("plan")
@_task_options
def plan(task, task_file, provider, secret_file, workspace):
    """Create and safely preview a command-compiled research workflow."""

    from chemsmart.agent.live_session import run_live_agent_session

    result = run_live_agent_session(
        task=_read_task(task, task_file),
        provider=provider.lower(),
        secret_file=secret_file,
        workspace=workspace,
        execution_enabled=False,
        approval_file=None,
    )
    click.echo(result.public_summary_json())


@agent.command("run")
@_task_options
@click.option(
    "--approval-file",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Host workflow approval. Omission keeps the session preview-only.",
)
def run(task, task_file, provider, secret_file, workspace, approval_file):
    """Plan and, when approved, execute host-compiled workflow nodes."""

    from chemsmart.agent.live_session import run_live_agent_session

    result = run_live_agent_session(
        task=_read_task(task, task_file),
        provider=provider.lower(),
        secret_file=secret_file,
        workspace=workspace,
        execution_enabled=approval_file is not None,
        approval_file=approval_file,
    )
    click.echo(result.public_summary_json())


__all__ = ["agent"]
