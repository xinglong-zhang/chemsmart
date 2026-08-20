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
    """Operate the ChemSmart Agent pipeline: plan, review, run, or tui.

    ``plan`` creates and safely previews a workflow, ``review`` re-presents a
    stored review for the one human decision, ``run`` executes an approved
    bundle provider-free, and ``tui`` is the interactive terminal over the
    same chain.
    """


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


@agent.command("review")
@click.option(
    "--review-file",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="A stored workflow execution review to re-present unchanged.",
)
@click.option(
    "--workspace",
    required=True,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    help="The exact workspace the review was produced in.",
)
@click.option(
    "--decision",
    default=None,
    type=click.Choice(("approve", "deny"), case_sensitive=False),
    help="Omit to preview only; supply one to record a fresh human decision.",
)
@click.option(
    "--actor",
    default=None,
    help="Human actor recorded with the decision; required with --decision.",
)
@click.option(
    "--approval-id",
    default=None,
    help="Decision identity; a fresh replay id is minted when omitted.",
)
@click.option(
    "--task-spec-sha256",
    default=None,
    help="Optional: refuse unless the review targets this exact task.",
)
def review(
    review_file, workspace, decision, actor, approval_id, task_spec_sha256
):
    """Re-present a stored workflow execution review for one human decision.

    Without ``--decision`` this only displays the stored review. With one, it
    records the fresh human resolution under its own approval id and decision
    log, and an approval writes a one-shot bundle for ``chemsmart agent run``.
    A spent approval is never reused; deciding again is a new decision.
    """

    from chemsmart.agent._contracts import ContractError
    from chemsmart.agent.live_session import (
        inspect_workflow_execution_replay,
        replay_approval_id,
        resolve_workflow_execution_review,
    )

    normalized = (decision or "").lower()
    if normalized and not actor:
        raise click.UsageError("--actor is required with --decision")
    try:
        report = inspect_workflow_execution_replay(
            review_file=review_file,
            workspace=workspace,
            task_spec_sha256=str(task_spec_sha256 or ""),
        )
    except ContractError as exc:
        raise click.ClickException(str(exc)) from exc

    summary = {
        key: value
        for key, value in report.items()
        if key != "canonical_review"
    }
    if not normalized:
        summary["decision"] = "not taken; pass --decision to decide"
        click.echo(json.dumps(summary, indent=2, sort_keys=True))
        return

    if report["missing_approved_artifacts"]:
        # Resolving inputs would fail before anything is dispatched, so a
        # decision here would be spent on a run that cannot start.
        raise click.ClickException(
            "these approved bytes are no longer under the workspace, so this "
            "workflow cannot execute as approved: "
            f"{report['missing_approved_artifacts']}. Plan it again rather "
            "than deciding on a run that must fail."
        )
    chosen = str(approval_id or "").strip() or replay_approval_id()
    scope = (
        Path(workspace).resolve() / ".chemsmart-agent" / "replays" / chosen
    )
    scope.mkdir(parents=True, exist_ok=True, mode=0o700)
    try:
        resolution, bundle = resolve_workflow_execution_review(
            review_file=review_file,
            reviewed_sha256=report["review_sha256"],
            decision=normalized,
            actor=actor,
            # Scope both per decision: the default decision log is shared by
            # every approval of one review, and a second append under the same
            # resolution id is refused before a bundle can be written.
            output_file=(
                scope / "bundle.json" if normalized == "approve" else None
            ),
            decision_log=scope / "decisions.jsonl",
            approval_id=chosen,
        )
    except ContractError as exc:
        raise click.ClickException(str(exc)) from exc
    summary.update(
        {
            "decision": resolution.decision,
            "approval_id": resolution.approval_id,
            "resolution_sha256": resolution.resolution_sha256,
            "bundle_sha256": bundle.bundle_sha256 if bundle else "",
            "bundle_file": str(scope / "bundle.json") if bundle else "",
            "next": (
                f"chemsmart agent run --approval-file {scope}/bundle.json"
                f" --workspace {Path(workspace).resolve()}"
                f" --run-directory {scope}/run"
                if bundle
                else ""
            ),
        }
    )
    click.echo(json.dumps(summary, indent=2, sort_keys=True))


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


@agent.command("run")
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
def run(approval_file, workspace, run_directory, task_spec_sha256):
    """Execute an approved workflow bundle provider-free.

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
                "analysis_status": result.analysis_status,
                "analysis_report_path": result.analysis_report_path,
                "analysis_completion_receipt_sha256s": (
                    result.analysis_completion_receipt_sha256s
                ),
                "analysis_nodes": [
                    {
                        "node_id": node.node_id,
                        "analysis_kind": node.analysis_kind,
                        "state": node.state,
                        "reason": node.reason,
                        "receipt_sha256s": node.receipt_sha256s,
                    }
                    for node in result.analysis_nodes
                ],
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
