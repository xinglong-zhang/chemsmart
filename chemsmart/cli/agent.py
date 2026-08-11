"""Provider-neutral natural-language entrypoint for ChemSmart Agent."""

from __future__ import annotations

import json
from pathlib import Path

import click

from chemsmart.utils.cli import MyGroup


def _task_options(function):
    options = (
        click.option(
            "--task", type=str, default=None, help="Natural-language task."
        ),
        click.option(
            "--task-file",
            type=click.Path(exists=True, dir_okay=False, path_type=Path),
            default=None,
            help="UTF-8 file containing the natural-language task.",
        ),
        click.option(
            "--provider",
            type=str,
            default=None,
            help="Named agent.yaml profile; defaults to its active profile.",
        ),
        click.option(
            "--provider-config",
            type=click.Path(exists=True, dir_okay=False, path_type=Path),
            default=None,
            help="Provider YAML; defaults to ~/.chemsmart/agent/agent.yaml.",
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
        click.option(
            "--analysis-completion-file",
            type=click.Path(exists=True, dir_okay=False, path_type=Path),
            default=None,
            help=(
                "Host-owned JSON policy for receipt-based numerical-analysis "
                "completion."
            ),
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
def plan(
    task,
    task_file,
    provider,
    provider_config,
    secret_file,
    workspace,
    analysis_completion_file,
):
    """Create and safely preview a command-compiled research workflow."""

    from chemsmart.agent.live_session import run_live_agent_session

    result = run_live_agent_session(
        task=_read_task(task, task_file),
        provider=provider.lower() if provider else None,
        provider_config_file=provider_config,
        secret_file=secret_file,
        workspace=workspace,
        execution_enabled=False,
        approval_file=None,
        analysis_completion_file=analysis_completion_file,
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
def run(
    task,
    task_file,
    provider,
    provider_config,
    secret_file,
    workspace,
    analysis_completion_file,
    approval_file,
):
    """Plan and, when approved, execute host-compiled workflow nodes."""

    from chemsmart.agent.live_session import run_live_agent_session

    result = run_live_agent_session(
        task=_read_task(task, task_file),
        provider=provider.lower() if provider else None,
        provider_config_file=provider_config,
        secret_file=secret_file,
        workspace=workspace,
        execution_enabled=approval_file is not None,
        approval_file=approval_file,
        analysis_completion_file=analysis_completion_file,
    )
    click.echo(result.public_summary_json())


@agent.command("execute")
@click.option(
    "--approval-file",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Host workflow approval to execute. Required: this command runs an "
    "already-approved plan and never produces one.",
)
@click.option(
    "--workspace",
    required=True,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    help="Approved workspace holding the geometry and project artifacts.",
)
@click.option(
    "--run-directory",
    required=True,
    type=click.Path(file_okay=False, path_type=Path),
    help="Where the append-only event stream and receipts are written.",
)
@click.option(
    "--task-spec-sha256",
    default="",
    help="Task specification digest the approval was frozen against.",
)
def execute(approval_file, workspace, run_directory, task_spec_sha256):
    """Execute an approved workflow deterministically, with no model.

    Every scientific choice -- program, project YAML, method, node graph,
    charge and multiplicity -- was made when the plan was written and is
    frozen in the approval. What remains is bookkeeping with one right
    answer, so this command takes no provider and no credential: it cannot
    reach a model, and its options say so.
    """

    from chemsmart.agent.executor import execute_approved_workflow

    result = execute_approved_workflow(
        approval_file=approval_file,
        workspace=workspace,
        run_directory=run_directory,
        task_spec_sha256=task_spec_sha256,
    )
    click.echo(
        json.dumps(
            {
                "workflow_id": result.workflow_id,
                "plan_sha256": result.plan_sha256,
                "approval_sha256": result.approval_sha256,
                "status": result.status,
                "provider_calls": result.provider_calls,
                "run_directory": result.run_directory,
                "nodes": [
                    {
                        "node_id": node.node_id,
                        "program": node.program,
                        "jobtype": node.jobtype,
                        "state": node.state,
                        "invocation_identity_sha256": (
                            node.invocation_identity_sha256
                        ),
                        "execution_receipt_sha256": (
                            node.execution_receipt_sha256
                        ),
                        "failure": node.failure,
                    }
                    for node in result.nodes
                ],
            },
            indent=2,
            sort_keys=True,
        )
    )


@agent.command("experiment")
@click.option(
    "--case-id",
    required=True,
    help="Preregistered Qwen/PySCF development or transfer case ID.",
)
@click.option("--repeat-index", type=click.IntRange(min=0), default=0)
@click.option(
    "--decomposition/--no-decomposition",
    default=False,
    help="Enable the bounded specialist factor when its dispatcher is available.",
)
@click.option(
    "--feedback-projection",
    type=click.Choice(("full-v1", "causal-v1"), case_sensitive=False),
    default="full-v1",
    show_default=True,
)
@click.option(
    "--critic/--no-critic",
    default=False,
    help="Enable the fresh read-only critic when its dispatcher is available.",
)
@click.option(
    "--max-concurrency",
    type=click.IntRange(min=1, max=4),
    default=1,
    show_default=True,
)
@click.option(
    "--provider-config",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Provider YAML; defaults to ~/.chemsmart/agent/agent.yaml.",
)
@click.option(
    "--secret-file",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    required=True,
    help="Secret assignment file parsed as data; it is never sourced.",
)
@click.option(
    "--workspace",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    required=True,
    help="Disposable task workspace containing an approved XYZ artifact.",
)
def experiment(
    case_id,
    repeat_index,
    decomposition,
    feedback_projection,
    critic,
    max_concurrency,
    provider_config,
    secret_file,
    workspace,
):
    """Run one preview-only Qwen/PySCF D/F/C campaign episode."""

    from chemsmart.agent.experiments.qwen_pyscf_dfc import (
        build_qwen_dfc_arm,
    )
    from chemsmart.agent.experiments.qwen_pyscf_fixtures import (
        qwen_pyscf_cases_v1,
    )
    from chemsmart.agent.live_session import run_live_agent_session

    cases = {item.case_id: item for item in qwen_pyscf_cases_v1()}
    case = cases.get(str(case_id).strip().upper())
    if case is None:
        raise click.UsageError("unknown Qwen/PySCF campaign case ID")
    arm = build_qwen_dfc_arm(
        decomposition=decomposition,
        feedback_projection=feedback_projection,
        critic=critic,
        max_concurrency=max_concurrency,
    )
    result = run_live_agent_session(
        task=case.task,
        provider="alibaba-token-plan",
        provider_config_file=provider_config,
        secret_file=secret_file,
        workspace=workspace,
        execution_enabled=False,
        approval_file=None,
        experiment_arm=arm,
        experiment_case=case,
        experiment_repeat_index=repeat_index,
    )
    click.echo(result.public_summary_json())


__all__ = ["agent"]
