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
        click.option(
            "--identity-manifest",
            type=click.Path(exists=True, dir_okay=False, path_type=Path),
            default=None,
            help=(
                "User-approved names, roles, states, and workspace-relative "
                "geometry bindings for one or more molecular inputs."
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
@click.option(
    "--execution-envelope",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help=(
        "Resource/program bounds used only to prepare an exact review; "
        "this command never executes an engine."
    ),
)
@click.option(
    "--review-file",
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help="Write the inert exact execution review packet to this JSON file.",
)
def plan(
    task,
    task_file,
    provider,
    provider_config,
    secret_file,
    workspace,
    analysis_completion_file,
    identity_manifest,
    execution_envelope,
    review_file,
):
    """Create and safely preview a command-compiled research workflow."""

    from chemsmart.agent.live_session import run_live_agent_session
    from chemsmart.agent.identity import load_approved_molecular_input_manifest

    approved_inputs = (
        load_approved_molecular_input_manifest(
            identity_manifest, workspace=workspace
        )
        if identity_manifest is not None
        else ()
    )

    result = run_live_agent_session(
        task=_read_task(task, task_file),
        provider=provider.lower() if provider else None,
        provider_config_file=provider_config,
        secret_file=secret_file,
        workspace=workspace,
        execution_enabled=False,
        approval_file=None,
        execution_envelope_file=execution_envelope,
        analysis_completion_file=analysis_completion_file,
        approved_molecular_inputs=approved_inputs,
        review_file=review_file.resolve() if review_file is not None else None,
    )
    click.echo(result.public_summary_json())


@agent.command("run")
@_task_options
@click.option(
    "--approval-file",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help=(
        "Legacy option retained only to fail closed; approvals are consumed "
        "by 'chemsmart agent execute'."
    ),
)
@click.option(
    "--execution-envelope",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help=(
        "Science-free resource/program bounds for planning, safe preview, "
        "and an inert exact review packet. It grants no execution authority."
    ),
)
@click.option(
    "--review-file",
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help="Write the inert exact execution review packet to this JSON file.",
)
def run(
    task,
    task_file,
    provider,
    provider_config,
    secret_file,
    workspace,
    analysis_completion_file,
    identity_manifest,
    approval_file,
    execution_envelope,
    review_file,
):
    """Plan and safely preview; real execution is a separate approved step."""

    from chemsmart.agent.live_session import run_live_agent_session
    from chemsmart.agent.identity import load_approved_molecular_input_manifest

    if approval_file is not None:
        raise click.UsageError(
            "--approval-file cannot grant authority to a provider session; "
            "use 'chemsmart agent execute --approval-file ...'"
        )

    approved_inputs = (
        load_approved_molecular_input_manifest(
            identity_manifest, workspace=workspace
        )
        if identity_manifest is not None
        else ()
    )

    result = run_live_agent_session(
        task=_read_task(task, task_file),
        provider=provider.lower() if provider else None,
        provider_config_file=provider_config,
        secret_file=secret_file,
        workspace=workspace,
        execution_enabled=False,
        approval_file=None,
        execution_envelope_file=execution_envelope,
        analysis_completion_file=analysis_completion_file,
        approved_molecular_inputs=approved_inputs,
        review_file=review_file.resolve() if review_file is not None else None,
    )
    click.echo(result.public_summary_json())


@agent.command("approve", hidden=True)
@click.option(
    "--review-file",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Inert workflow execution review produced by agent plan or run.",
)
@click.option(
    "--reviewed-sha256",
    required=True,
    help="The full review digest the human inspected and is resolving.",
)
@click.option(
    "--decision",
    required=True,
    type=click.Choice(
        ("approve", "deny", "revise", "quit"), case_sensitive=False
    ),
    help="One exact human resolution; there are no session or prefix grants.",
)
@click.option(
    "--actor", required=True, help="Human actor recorded with the decision."
)
@click.option(
    "--approval-id",
    default=None,
    help="Optional public approval ID; otherwise it is derived from the digest.",
)
@click.option(
    "--output",
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help="Required for approve: exact one-shot approval bundle JSON.",
)
@click.option(
    "--decision-log",
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help="Append-only decision event stream; defaults beside the review file.",
)
def approve(
    review_file,
    reviewed_sha256,
    decision,
    actor,
    approval_id,
    output,
    decision_log,
):
    """Legacy file-based approval compatibility command."""

    from chemsmart.agent._contracts import ContractError
    from chemsmart.agent.live_session import resolve_workflow_execution_review

    normalized = decision.lower()
    if normalized == "approve" and output is None:
        raise click.UsageError("--output is required when --decision=approve")
    if normalized != "approve" and output is not None:
        raise click.UsageError(
            "--output is valid only when --decision=approve"
        )
    try:
        resolution, bundle = resolve_workflow_execution_review(
            review_file=review_file,
            reviewed_sha256=reviewed_sha256,
            decision=normalized,
            actor=actor,
            output_file=output.resolve() if output is not None else None,
            decision_log=(
                decision_log.resolve() if decision_log is not None else None
            ),
            approval_id=str(approval_id or ""),
        )
    except ContractError as exc:
        raise click.ClickException(str(exc)) from exc
    click.echo(
        json.dumps(
            {
                "decision": resolution.decision,
                "review_sha256": resolution.review_sha256,
                "resolution_sha256": resolution.resolution_sha256,
                "approval_id": resolution.approval_id,
                "bundle_sha256": bundle.bundle_sha256 if bundle else "",
                "one_shot": bool(bundle),
            },
            indent=2,
            sort_keys=True,
        )
    )


@agent.command("tui")
@click.option(
    "--provider",
    type=str,
    default=None,
    help="Named agent.yaml profile; defaults to its active profile.",
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
    "--execution-envelope",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Optional local resource/program bounds for an execution review.",
)
@click.option(
    "--review-file",
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help=(
        "Optional audit export of the prepared workflow; it is not required "
        "for TUI approval or execution."
    ),
)
@click.option(
    "--workspace",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    required=True,
    help="Task workspace containing user-approved artifacts.",
)
@click.option(
    "--analysis-completion-file",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Host-owned JSON policy for receipt-based numerical analysis.",
)
@click.option(
    "--identity-manifest",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="User-approved molecular identity and geometry bindings.",
)
@click.option(
    "--plain",
    is_flag=True,
    help="Use an inline, no-mouse terminal mode for SSH and simple terminals.",
)
def tui(
    provider,
    provider_config,
    secret_file,
    execution_envelope,
    review_file,
    workspace,
    analysis_completion_file,
    identity_manifest,
    plain,
):
    """Open interactive planning, analysis, review, and approved execution."""

    from chemsmart.agent.tui import launch_tui

    try:
        launch_tui(
            workspace=workspace,
            secret_file=secret_file,
            execution_envelope_file=execution_envelope,
            review_file=review_file,
            provider=provider.lower() if provider else None,
            provider_config_file=provider_config,
            identity_manifest=identity_manifest,
            analysis_completion_file=analysis_completion_file,
            plain=plain,
        )
    except RuntimeError as exc:
        raise click.ClickException(str(exc)) from exc


@agent.command("execute", hidden=True)
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
    """Legacy file-based execution compatibility command.

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
                "non_executable_node_ids": result.non_executable_node_ids,
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


__all__ = ["agent"]
